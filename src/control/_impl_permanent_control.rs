use biodivine_lib_param_bn::symbolic_async_graph::{SymbolicAsyncGraph, GraphColoredVertices, GraphColors};
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::{BooleanNetwork, VariableId, FnUpdate};
use std::collections::HashMap;
use std::convert::TryFrom;
use biodivine_lib_param_bn::biodivine_std::bitvector::{ArrayBitVector, BitVector};
use crate::control::_algo_utils::{reach_fwd_with_saturation, reach_bwd_with_saturation, strong_basin, add_auto_regulations};

pub struct PermanentControl {
    controlled_network: BooleanNetwork,
    uncontrolled_network: BooleanNetwork,
    controlled_graph: SymbolicAsyncGraph,
    uncontrolled_graph: SymbolicAsyncGraph,
    source_set: GraphColoredVertices,
    target_set: GraphColoredVertices,
    controlled_basin: GraphColoredVertices,
    var_is_controlled: HashMap<VariableId, GraphColors>,
    var_is_not_controlled: HashMap<VariableId, GraphColors>,
}

impl PermanentControl {
    pub fn new(network: BooleanNetwork, source: &ArrayBitVector, target: &ArrayBitVector) -> PermanentControl {
        let normalized_network = add_auto_regulations(network);

        let mut controlled_network = BooleanNetwork::new(normalized_network.as_graph().clone());
        let mut uncontrolled_network = BooleanNetwork::new(normalized_network.as_graph().clone());

        // First, copy existing parameters from the original network:
        for p in normalized_network.parameters() {
            let parameter = normalized_network.get_parameter(p);
            let controlled_parameter = controlled_network.add_parameter(parameter.get_name(), parameter.get_arity()).unwrap();
            let uncontrolled_parameter = uncontrolled_network.add_parameter(parameter.get_name(), parameter.get_arity()).unwrap();
            if controlled_parameter != p || uncontrolled_parameter != p {
                panic!("This should not happen. Encodings should be the same.");
            }
        }

        // Then add control parameters and modify update functions:
        for v in normalized_network.variables().rev() {
            let v_controlled = format!("{}_controlled", normalized_network.get_variable_name(v));
            let controlled_parameter = controlled_network.add_parameter(v_controlled.as_str(), 0).unwrap();
            let uncontrolled_parameter = uncontrolled_network.add_parameter(v_controlled.as_str(), 0).unwrap();
            if controlled_parameter != uncontrolled_parameter {
                panic!("This should not happen. Encodings should be the same.");
            }
            let parameter_id = controlled_parameter;
            if let Some(function) = normalized_network.get_update_function(v) {
                // A little trick to avoid always cloning the value of fn_parameter...
                let fn_parameter = || FnUpdate::mk_param(parameter_id, &[]);

                // Set controlled function to (v_controlled => v) && (!v_controlled => f(...))
                let controlled_implies_v = fn_parameter().implies(FnUpdate::mk_var(v));
                let not_controlled_implies_f = FnUpdate::mk_not(fn_parameter()).implies(function.clone());
                let controlled_function = controlled_implies_v.and(not_controlled_implies_f);
                controlled_network.add_update_function(v, controlled_function).unwrap();

                // Set uncontrolled function to (v_controlled || !v_controlled) && f(...)
                let control_tautology = fn_parameter().or(FnUpdate::mk_not(fn_parameter()));
                let uncontrolled_function = control_tautology.and(function.clone());
                uncontrolled_network.add_update_function(v, uncontrolled_function).unwrap();
            } else {
                panic!(r"
                We assume there are no implicit unknown functions in the network.
                If there are, replace them with explicit parameters.
            ");
            }
        }

        let controlled_graph = SymbolicAsyncGraph::new(controlled_network.clone()).unwrap();
        let uncontrolled_graph = SymbolicAsyncGraph::new(uncontrolled_network.clone()).unwrap();

        // Now, we can also make some helper sets that represent when variables are controlled:
        let mut var_is_controlled = HashMap::new();
        let mut var_is_not_controlled = HashMap::new();
        for v in normalized_network.variables() {
            let v_name = controlled_graph.as_network().get_variable_name(v);
            let parameter = controlled_graph.as_network().find_parameter(format!("{}_controlled", v_name).as_str()).unwrap();
            let bdd_v_is_controlled = controlled_graph.symbolic_context().mk_uninterpreted_function_is_true(parameter, &[]);
            let bdd_v_is_not_controlled = bdd_v_is_controlled.not();
            // This "copy" method is a workaround that "internally" creates a symbolic set with an arbitrary Bdd.
            // However, we are responsible for the Bdd being a valid symbolic representation of what we want.
            // It takes the metadata from the `unit_colors` symbolic set, but swaps the Bdd with our value.
            var_is_controlled.insert(v, controlled_graph.unit_colors().copy(bdd_v_is_controlled));
            var_is_not_controlled.insert(v, controlled_graph.unit_colors().copy(bdd_v_is_not_controlled));
        }

        let target_set = controlled_graph.vertex(target);
        // let target_set = target.iter().fold(uncontrolled_graph.mk_unit_colored_vertices(), |a, (v, value)| {
        //     a.intersect(&uncontrolled_graph.fix_network_variable(*v, *value))
        // });

        let source_set = controlled_graph.vertex(source);
        // let source_set = source.iter().fold(uncontrolled_graph.mk_unit_colored_vertices(), |a, (v, value)| {
        //     a.intersect(&uncontrolled_graph.fix_network_variable(*v, *value))
        // });

        let target_fwd = reach_fwd_with_saturation(&controlled_graph, target_set.clone());
        let target_bwd = reach_bwd_with_saturation(&controlled_graph, target_set.clone());
        let target_not_attractor = target_fwd.minus(&target_bwd).colors();

        // let source_fwd = reach_fwd_with_saturation(&uncontrolled_graph, source_set.clone());
        // let source_bwd = reach_bwd_with_saturation(&uncontrolled_graph, source_set.clone());
        // let source_not_attractor = source_fwd.minus(&source_bwd).colors();

        // We have to also remove these colors, because these appear in situations where source and
        // target are in the same attractor.
        // let target_can_reach_source = target_fwd.intersect(&source_fwd).colors();
        let valid_colours = uncontrolled_graph.unit_colors()
            // .minus(&source_not_attractor)
            // .minus(&target_can_reach_source)
            .minus(&target_not_attractor);

        // And make a new target set but only with these valid colors:
        let target_set = target_set.clone().intersect_colors(&valid_colours);

        // Consequently, we can now compute a strong basin of the set in the CONTROLLED graph,
        // which will give us the state that will eventually always reach this basin of "good" states.
        let mut controlled_basin = strong_basin(&controlled_graph, target_set.clone());

        // This gives us triples `(state, control, parametrisation)` such that on all paths starting
        // in `state`, there is a point after which if `control` is released, the network will converge
        // to `target`.

        // Now what we need to do is to only keep states that are equivalent to source under control.
        // This boils down to removing any state where uncontrolled variables do not match the value
        // in `source` (since control application will fix all controlled variables).

        for (v, source_value) in normalized_network.clone().variables().zip(source.values()) {
            let v_different_than_source = uncontrolled_graph.fix_network_variable(v, !source_value);
            let v_not_controlled = var_is_not_controlled.get(&v).unwrap();

            let v_not_controlled_and_different_than_source = v_different_than_source.intersect_colors(v_not_controlled);
            controlled_basin = controlled_basin.minus(&v_not_controlled_and_different_than_source);
        }

        println!("{}", controlled_basin.symbolic_size());

        PermanentControl {
            controlled_network: controlled_network.clone(),
            uncontrolled_network,
            controlled_graph,
            uncontrolled_graph,
            source_set,
            target_set: target_set.clone(),
            controlled_basin,
            var_is_controlled,
            var_is_not_controlled
        }
    }

    pub fn one_size_controls(&self) {
        for v in self.controlled_network.variables() {
                let variables_other_then_v_are_not_controlled = self.var_is_not_controlled.iter()
                    .filter(|(id, _)| **id != v)
                    .fold(self.uncontrolled_graph.mk_unit_colors(), |a, (_, b)| a.intersect(b));
                let only_v_is_controlled = self.var_is_controlled.get(&v).unwrap().intersect(&variables_other_then_v_are_not_controlled);

                // In this symbolic set, all values of `v_controlled` parameters are fixed. We could do
                // this more easily if we saved `BddVariables` that the `v_controlled` parameters correspond
                // to, but for now this is ok.
                let controls_where_only_v_is_controlled = self.controlled_basin.intersect_colors(&only_v_is_controlled);

                let v_true = self.uncontrolled_graph.fix_network_variable(v, true);
                let v_false = self.controlled_graph.fix_network_variable(v, false);

                println!("Control [{}=true] works for {} colors.", self.controlled_network.get_variable_name(v), v_true.intersect(&controls_where_only_v_is_controlled).colors().approx_cardinality());
                println!("Control [{}=false] works for {} colors.", self.controlled_network.get_variable_name(v), v_false.intersect(&controls_where_only_v_is_controlled).colors().approx_cardinality());
            }
    }
}

// fn powerset(s: &VariableIdIterator) -> Vec<Vec<VariableId>> {
//     (0..2usize.pow(s.len() as u32)).map(|i| {
//         s.enumerate().filter(|&(t, _)| (i >> t) % 2 == 1)
//             .map(|(_, element)| element.clone())
//             .collect()
//     }).collect()
// }

    //test
    // let witness = uncontrolled_graph.pick_witness(uncontrolled_graph.unit_colors());
    // let witness_graph = SymbolicAsyncGraph::new(witness).unwrap();
    //
    // // Variables in this set will be controlled to true
    // let true_controls = powerset(&witness.clone().variables());
    // for true_control in true_controls {
    //
    //     // Variables in this set will be controlled to false
    //     let false_controls = powerset(&witness.clone().variables());
    //     for false_control in false_controls {
    //         let intersection = true_control.clone().iter().any(|x| false_control.contains(x));
    //         if intersection {
    //             // One variable can not be set to same value
    //             continue
    //         }
    //
    //         let mut immediate_successor = HashMap::new();
    //         for v in witness.clone().variables() {
    //             if true_control.clone().contains(&v) {
    //                 immediate_successor.insert(v, true);
    //             } else if false_control.contains(&v) {
    //                 immediate_successor.insert(v, false);
    //             } else {
    //                 immediate_successor.insert(v, *source.get(&v).unwrap());
    //             }
    //         }
    //         let immediate_successor_set = immediate_successor.iter().fold(uncontrolled_graph.mk_unit_colored_vertices(), |a, (v, value)| {
    //             a.intersect(&uncontrolled_graph.fix_network_variable(*v, *value))
    //         });
    //
    //         let control = true_control.clone().into_iter().chain(false_control.into_iter())
    //                                              .map(|name| var_is_controlled.get(&name).unwrap())
    //                                              .fold(uncontrolled_graph.mk_unit_colors(),
    //                                                    |a, b| a.intersect(b));
    //
    //         let final_color = control.intersect(&witness_graph.unit_colors());
    //         let controlled_witness = controlled_graph.pick_witness(&final_color);
    //         let controlled_witness_graph = SymbolicAsyncGraph::new(controlled_witness).unwrap();
    //         let stabilized_after_control = find_reachable_attractors(&controlled_witness_graph, immediate_successor_set);
    //
    //         let relevant_basin = basin_of_target.intersect_colors(witness_graph.unit_colors());
    //         let not_in_basin_vertices = stabilized_after_control.vertices().minus(&relevant_basin.vertices());
    //         assert!(not_in_basin_vertices.is_empty());
    //     }
    // }
//}


// #[cfg(test)]
// mod tests {
//     use biodivine_lib_param_bn::BooleanNetwork;
//     use std::convert::TryFrom;
//     use std::collections::HashMap;
//     use crate::control::_algo_temporary_control::temporary_source_target_control;
//
//     #[test]
//     #[allow(non_snake_case)]
//     fn test_g2a_control() {
//         // This is the G2A network which is often used as an example, but here, we have added
//         // an a self-loop to all variables that otherwise don't have one (because in a controlled
//         // graph, technically every variable has a self-loop).
//
//         // We have also given each variable an explicit update function.
//         let g2a_network = r"
//             $CtrA:f_5(CtrA, GcrA, CcrM, SciP)
//             CtrA -> CtrA
//             GcrA -> CtrA
//             CcrM -| CtrA
//             SciP -| CtrA
//             $GcrA:f_4(CtrA, DnaA)
//             GcrA -?? GcrA
//             CtrA -| GcrA
//             DnaA -> GcrA
//             $CcrM:f_3(CtrA, CcrM, SciP)
//             CtrA -> CcrM
//             CcrM -| CcrM
//             SciP -| CcrM
//             $SciP:f_1(CtrA, DnaA)
//             CtrA -> SciP
//             DnaA -| SciP
//             SciP -?? SciP
//             $DnaA:f_2(CtrA, GcrA, DnaA, CcrM)
//             CtrA -> DnaA
//             GcrA -| DnaA
//             DnaA -| DnaA
//             CcrM -> DnaA
//         ";
//
//         let network = BooleanNetwork::try_from(g2a_network).unwrap();
//         let CcrM = network.as_graph().find_variable("CcrM").unwrap();
//         let CtrA = network.as_graph().find_variable("CtrA").unwrap();
//         let DnaA = network.as_graph().find_variable("DnaA").unwrap();
//         let GcrA = network.as_graph().find_variable("GcrA").unwrap();
//         let SciP = network.as_graph().find_variable("SciP").unwrap();
//
//         let mut source = HashMap::new();
//         source.insert(CcrM, true);
//         source.insert(CtrA, true);
//         source.insert(DnaA, true);
//         source.insert(GcrA, false);
//         source.insert(SciP, false);
//
//         let mut target = HashMap::new();
//         target.insert(CcrM, true);
//         target.insert(CtrA, false);
//         target.insert(DnaA, false);
//         target.insert(GcrA, false);
//         target.insert(SciP, false);
//
//         let result = temporary_source_target_control(network, &source, &target);
//     }

//}