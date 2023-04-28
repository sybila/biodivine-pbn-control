use std::collections::HashMap;
use std::time::SystemTime;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::GraphVertices;
use biodivine_lib_param_bn::symbolic_async_graph::projected_iteration::{RawProjection, RawSymbolicIterator};
use crate::perturbation::PerturbationGraph;
use crate::phenotype_control::_symbolic_utils::mk_bdd_of_bound;

pub fn bounded_phenotype_control(
    graph: &PerturbationGraph,
    phenotype: &GraphVertices,
    size_bound: usize,
) {

    let not_phenotype = graph.unit_colored_vertices()
        .vertices().minus(&phenotype);

    // A map which gives us the symbolic variable of the perturbation parameter.
    let perturbation_var_map = graph.variables()
        .filter_map(|var| graph.get_perturbation_parameter(var).map(|it| (var, it)))
        .map(|(var, param)| {
            (var, graph.as_symbolic_context()
                .get_explicit_function_table(param)
                .symbolic_variables()[0])
        })
        .collect::<HashMap<_, _>>();

    let bdd_vars = graph.as_symbolic_context().bdd_variable_set();
    // The list of symbolic variables of perturbation parameters.
    let perturbation_vars = {
        let mut values = perturbation_var_map
            .values()
            .cloned()
            .collect::<Vec<_>>();
        values.sort();
        values
    };

    for perturbation_size in 0..(size_bound + 1) {
        let start = SystemTime::now();
        println!("Perturbation size: {}", perturbation_size);
        let admissible_perturbations = mk_bdd_of_bound(bdd_vars, &perturbation_vars, perturbation_size);
        {
            let factor = 2.0f64.powi(bdd_vars.num_vars() as i32 - perturbation_vars.len() as i32);
            println!("[{}] >> Admissible fixed(Q) sets: {}", perturbation_size, admissible_perturbations.cardinality() / factor);
        }
        let admissible_perturbations = graph.empty_colors().copy(admissible_perturbations);

        // This is a trap set of all state-color-perturbation combinations that are
        // guaranteed to stay in the phenotype.
        let mut trap = graph.unit_colored_vertices()
            .intersect_colors(&admissible_perturbations)
            .intersect_vertices(&phenotype);

        'trap: loop {
            for var in graph.variables().rev() {
                let can_leave = graph.as_perturbed().var_can_post_out(var, &trap);
                if !can_leave.is_empty() {
                    trap = trap.minus(&can_leave);
                    if trap.symbolic_size() > 100_000 {
                        println!("[{}] >> Trap phenotype progress: {} / {}", perturbation_size, trap.symbolic_size(), trap.approx_cardinality());
                    }
                    continue 'trap;
                }
            }
            break;
        }

        let mut trap = graph.unit_colored_vertices()
            .intersect_colors(&admissible_perturbations)
            .minus(&trap);

        'trap: loop {
            for var in graph.variables().rev() {
                let can_leave = graph.as_perturbed().var_can_post_out(var, &trap);
                if !can_leave.is_empty() {
                    trap = trap.minus(&can_leave);
                    if trap.symbolic_size() > 100_000 {
                        println!("[{}] >> Trap non-phenotype progress: {} / {}", perturbation_size, trap.symbolic_size(), trap.approx_cardinality());
                    }
                    continue 'trap;
                }
            }
            break;
        }

        let mut inverse_control = trap.into_bdd();
        for var in graph.variables() {
            let state_var = graph.as_symbolic_context().get_state_variable(var);
            if let Some(perturbation_var) = perturbation_var_map.get(&var) {
                // If the variable can be perturbed, we split into two cases and eliminate
                // it in the unperturbed cases.

                let is_perturbed = inverse_control
                    .var_select(*perturbation_var, true);
                let is_not_perturbed = inverse_control
                    .var_select(*perturbation_var, false)
                    .var_project(state_var);
                inverse_control = is_perturbed.or(&is_not_perturbed);
            } else {
                // If the variable cannot be perturbed, we can eliminate it everywhere.
                inverse_control = inverse_control.var_project(state_var);
            }
        }

        let inverse_control_map = graph.empty_colored_vertices().copy(inverse_control);

        // Control map consists of admissible state-color pairs that are also admissible
        // for our perturbation size and are not in the inverse map.
        let control_map = graph.unit_colored_vertices()
            .intersect_colors(&admissible_perturbations)
            .minus(&inverse_control_map);

        {
            // Compute the number of valuations of the perturbation parameters.
            let factor = 2.0f64.powi(bdd_vars.num_vars() as i32 - perturbation_vars.len() as i32);
            let mut only_perturbation_parameters = control_map.clone().into_bdd();
            for var in bdd_vars.variables() {
                if !perturbation_vars.contains(&var) {
                    only_perturbation_parameters = only_perturbation_parameters.var_project(var);
                }
            }
            println!("[{}] >> fixed(Q) sets in control map: {}", perturbation_size, only_perturbation_parameters.cardinality() / factor);
        }

        let control_map_bdd = control_map.clone().into_bdd();
        let perturbation_vars_projection = RawProjection::new(perturbation_vars.clone(), &control_map_bdd);

        let all_colors_size = graph.unit_colors().approx_cardinality() / 2.0f64.powi(perturbation_vars.len() as i32);
        let mut best_robustness = 0.0;
        let mut with_best_robustness = 0;

        for is_perturbed_vector in perturbation_vars_projection.iter() {
            let perturbed_variables = perturbation_var_map
                .iter()
                .filter_map(|(var, p_var)| {
                    if is_perturbed_vector.get_value(*p_var).unwrap() {
                        Some(*var)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();

            // This should remove all perturbation symbolic variables from the set.
            let control_subset = control_map_bdd.restrict(&is_perturbed_vector.to_values());

            let mut perturbed_state_vars = perturbed_variables
                .iter()
                .map(|var| graph.as_symbolic_context().get_state_variable(*var))
                .collect::<Vec<_>>();

            let state_vars_projection = RawProjection::new(perturbed_state_vars, &control_subset);
            for state_vector in state_vars_projection.iter() {
                // This should remove all remaining perturbed state variables. Unperturbed
                // variables should not appear in the control map add all.
                let control_colors = control_subset.restrict(&state_vector.to_values());
                let mut map = HashMap::new();
                for var in &perturbed_variables {
                    let state_var = graph.as_symbolic_context().get_state_variable(*var);
                    map.insert(
                        graph.as_original().as_network().get_variable_name(*var).clone(),
                        state_vector.get_value(state_var).unwrap(),
                    );
                }
                let factor = 2.0f64.powi(
                    (graph.as_symbolic_context().state_variables().len() + graph.num_perturbation_parameters()) as i32
                );

                let working_colors = control_colors.cardinality() / factor;
                let robustness = working_colors / all_colors_size;
                if robustness >= best_robustness {
                    if robustness != best_robustness {
                        with_best_robustness = 0;
                    }
                    best_robustness = robustness;
                    with_best_robustness += 1;
                    println!("[{}] >>>> {:?}: {}; rho = {:.2}", perturbation_size, map, working_colors, robustness);
                }
            }
        }

        println!("[{}] Elapsed: {}ms", perturbation_size, start.elapsed().unwrap().as_millis());

        println!("[{}] Best robustness {} for {} perturbations.", perturbation_size, best_robustness, with_best_robustness);

        // if best_robustness == 1.0 {
        //     println!("Sufficient robustness achieved for perturbation size {}.", perturbation_size);
        //     return;
        // }
    }

    println!("Sufficient robustness not achieved with perturbation size {}.", size_bound);
}