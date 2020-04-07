use biodivine_lib_param_bn::async_graph::{AsyncGraph};
use std::collections::HashMap;
use biodivine_lib_param_bn::bdd_params::BddParams;
use biodivine_lib_std::param_graph::{Graph, EvolutionOperator, Params};
use biodivine_lib_std::{IdState};

fn find_strong_basin(graph: &AsyncGraph, attractor: IdState, params: &BddParams) -> HashMap<IdState, BddParams> {
    let bwd = graph.bwd();
    let fwd = graph.fwd();
    let mut basin = HashMap::new();
    basin.insert(attractor, params.clone());

    // Just a quick sanity check to verify that the given `attractor` state is really a sink for
    // all parameters requested in `params`. If a successor has non-empty intersection with
    // `params`, then we have a problem.
    let successor_count = fwd.step(attractor).filter(|(_, p)| {
         !params.intersect(p).is_empty()
    }).count();
    if successor_count != 0 {
        panic!("Given state ({:?}) is not an attractor. It has {} successor(s).", attractor, successor_count);
    }


    let empty_params = &graph.empty_params();
    let mut found_something = true;

    while found_something {
        found_something = false;
        let current_basin: Vec<IdState> = basin.keys().cloned().collect();

        for node in current_basin {
            for (parent, parent_edge_params) in bwd.step(node) {
                // Reference to parameters for which parent is currently in basin.
                let parent_in_basin = basin.get(&parent).unwrap_or(empty_params);
                // Parameters which can be potentially still added to the basin.
                let can_become_basin = parent_edge_params.minus(parent_in_basin);

                if can_become_basin.is_empty() {    // skip useless states
                    continue;
                }

                // We start with all potential parameters and remove any parametrisation that
                // has a successor which is not in the basin.
                let mut add_to_basin = can_become_basin;
                for (successor, successor_edge_params) in fwd.step(parent) {
                    if add_to_basin.is_empty() {
                        break;  // We might as well end this if we know there is nothing to add.
                    }
                    // Reference to parameters for which successor is currently in basin.
                    let successor_in_basin = basin.get(&successor).unwrap_or(empty_params);
                    // Parameters for which successor is reachable, but not in basin.
                    let reachable_but_not_in_basin = successor_edge_params.minus(successor_in_basin);
                    add_to_basin = add_to_basin.minus(&reachable_but_not_in_basin);
                }

                // If something still remains in add_to_basin, actually add it to basin!
                if !add_to_basin.is_empty() {
                    found_something = true;
                    match basin.get_mut(&parent) {
                        None => {
                            basin.insert(parent, add_to_basin);
                        },
                        Some(value) => {
                            *value = value.union(&add_to_basin);
                        }
                    }
                }
            }
        }
    }
    return basin;
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use biodivine_lib_param_bn::BooleanNetwork;
    use biodivine_lib_bdd::{Bdd, BddVariableSet};
    use biodivine_lib_param_bn::bdd_params::BddParameterEncoder;
    use std::env::var;
    use std::path::Path;
    use std::fmt::Debug;
    use std::convert::TryFrom;

    const TEST_MODEL: &str = "\
        # Since the model is fully instantiated, we don't really care about the exact
        # types of regulations. But whatever...

        CtrA -> CtrA
        GcrA -> CtrA
        CcrM -| CtrA
        SciP -| CtrA
        $CtrA: !CcrM | CtrA | (GcrA & !SciP)

        CtrA -| GcrA
        DnaA -> GcrA
        $GcrA: !CtrA & DnaA

        CtrA -> SciP
        DnaA -| SciP
        $SciP: CtrA & !DnaA

        CtrA -> CcrM
        CcrM -| CcrM
        SciP -| CcrM
        $CcrM: !CcrM | CtrA | !SciP

        CcrM -> DnaA
        CtrA -> DnaA
        GcrA -| DnaA
        DnaA -| DnaA
        $DnaA: (CcrM & CtrA) | (!CcrM & CtrA & (!DnaA | !GcrA))
    ";

    #[test]
    fn test_reach() {
        //let sbml_str = fs::read_to_string("src/g2b.sbml").unwrap();
        //let (model, _) = BooleanNetwork::from_sbml(&sbml_str).unwrap();

        let model = BooleanNetwork::try_from(TEST_MODEL).unwrap();

        //println!("{:?}", model.graph());

        let graph = &AsyncGraph::new(model).unwrap();
        // We want to start from state
        // CcrM | CtrA | DnaA | GcrA | SciP
        //  1   |  1   |  1   |  0   |  0
        // Note that bits are reversed because first variable corresponds to least significant bit.
        let state = IdState::from(0b00111 as usize);

        let basin = find_strong_basin(graph, state, graph.unit_params());
        println!("Strong basin has {} states.", basin.len());
        for (id, param) in basin {
            println!("{} {:?}", id, param.cardinality());
        }
    }
}
