use crate::perturbation::PerturbationGraph;
use crate::phenotype_control::PhenotypeControlMap;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColors, GraphVertices};
use biodivine_lib_param_bn::VariableId;
use itertools::Itertools;
use std::collections::HashMap;
use std::time::SystemTime;
use biodivine_lib_param_bn::symbolic_async_graph::projected_iteration::RawProjection;
use crate::aeon::reachability::backward_within;
use crate::phenotype_control::_symbolic_utils::mk_bdd_of_bound;

impl PerturbationGraph {
    pub fn phenotype_permanent_control(
        &self,
        phenotype: GraphVertices,
        admissible_perturbations: GraphColors,
        allow_oscillation: bool
    ) -> PhenotypeControlMap {
        let perturbation_bbd_vars_mapping = self.variables()
            .filter_map(|var| self.get_perturbation_parameter(var).map(|it| (var, it)))
            .map(|(var, param)| {
                (var, self.as_symbolic_context()
                    .get_explicit_function_table(param)
                    .symbolic_variables()[0])
            })
            .collect::<HashMap<_, _>>();

        let mut trap;
        let control_universe = self.unit_colored_vertices().intersect_colors(&admissible_perturbations);
        let phenotype_coloured_vertices = control_universe.intersect_vertices(&phenotype);
        if allow_oscillation {
            let bwd_reach = backward_within(self.as_perturbed(), &phenotype_coloured_vertices, &control_universe);
            trap = bwd_reach
        } else {
            trap = phenotype_coloured_vertices;
        }

        'trap: loop {
            for var in self.variables().rev() {
                let can_leave = self.as_perturbed().var_can_post_out(var, &trap);
                if !can_leave.is_empty() {
                    trap = trap.minus(&can_leave);
                    if trap.symbolic_size() > 100_000 {
                        println!(" Trap phenotype progress: {} / {}",  trap.symbolic_size(), trap.approx_cardinality());
                    }
                    continue 'trap;
                }
            }
            break;
        }

        let mut trap = self.unit_colored_vertices()
            .intersect_colors(&admissible_perturbations)
            .minus(&trap);

        'trap: loop {
            for var in self.variables().rev() {
                let can_leave = self.as_perturbed().var_can_post_out(var, &trap);
                if !can_leave.is_empty() {
                    trap = trap.minus(&can_leave);
                    if trap.symbolic_size() > 100_000 {
                        println!(" Trap non-phenotype progress: {} / {}", trap.symbolic_size(), trap.approx_cardinality());
                    }
                    continue 'trap;
                }
            }
            break;
        }

        let mut inverse_control = trap.into_bdd();
        for var in self.variables() {
            let state_var = self.as_symbolic_context().get_state_variable(var);
            if let Some(perturbation_var) = perturbation_bbd_vars_mapping.get(&var) {
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

        let inverse_control_map = self.empty_colored_vertices().copy(inverse_control);

        // Control map consists of admissible state-color pairs that are also admissible
        // for our perturbation size and are not in the inverse map.
        let control_map = self.unit_colored_vertices()
            .intersect_colors(&admissible_perturbations)
            .minus(&inverse_control_map);

        PhenotypeControlMap {
            perturbation_set: control_map,
            context: self.clone()
        }
    }

    pub fn ceiled_phenotype_permanent_control(
        &self,
        phenotype: GraphVertices,
        size_bound: usize,
        perturbation_variables: Vec<VariableId>,
        stop_early: bool,
        allow_oscillation: bool
    ) -> PhenotypeControlMap {
        // A map which gives us the symbolic variable of the perturbation parameter.
        let perturbation_bbd_vars_mapping = perturbation_variables.iter()
            .filter_map(|var| self.get_perturbation_parameter(var.clone()).map(|it| (var.clone(), it)))
            .map(|(var, param)| {
                (var, self.as_symbolic_context()
                    .get_explicit_function_table(param)
                    .symbolic_variables()[0])
            })
            .collect::<HashMap<_, _>>();

        let bdd_vars = self.as_symbolic_context().bdd_variable_set();
        // The list of symbolic variables of perturbation parameters.
        let perturbation_bdd_vars = {
            let mut values = perturbation_bbd_vars_mapping
                .values()
                .cloned()
                .collect::<Vec<_>>();
            values.sort();
            values
        };

        let mut control_map_all = self.mk_empty_colored_vertices();

        for perturbation_size in 0..(size_bound + 1) {
            let start = SystemTime::now();
            println!("Perturbation size: {}", perturbation_size);
            let admissible_perturbations = mk_bdd_of_bound(bdd_vars, &perturbation_bdd_vars, perturbation_size);
            {
                let factor = 2.0f64.powi(bdd_vars.num_vars() as i32 - perturbation_bdd_vars.len() as i32);
                println!("[{}] >> Admissible fixed(Q) sets: {}", perturbation_size, admissible_perturbations.cardinality() / factor);
            }
            let admissible_perturbations = self.empty_colors().copy(admissible_perturbations);
            let control_map = self.phenotype_permanent_control(phenotype.clone(), admissible_perturbations, allow_oscillation).perturbation_set;

            // Compute the number of valuations of the perturbation parameters.
            let factor = 2.0f64.powi(bdd_vars.num_vars() as i32 - perturbation_bdd_vars.len() as i32);
            let mut only_perturbation_parameters = control_map.clone().into_bdd();
            for var in bdd_vars.variables() {
                if !perturbation_bdd_vars.contains(&var) {
                    only_perturbation_parameters = only_perturbation_parameters.var_project(var);
                }
            }
            println!("[{}] >> fixed(Q) sets in control map: {}", perturbation_size, only_perturbation_parameters.cardinality() / factor);


            let control_map_bdd = control_map.clone().into_bdd();
            let perturbation_vars_projection = RawProjection::new(perturbation_bdd_vars.clone(), &control_map_bdd);

            let all_colors_size = self.unit_colors().approx_cardinality() / 2.0f64.powi(perturbation_bdd_vars.len() as i32);
            let mut best_robustness = 0.0;
            let mut with_best_robustness = 0;

            for is_perturbed_vector in perturbation_vars_projection.iter() {
                let perturbed_variables = perturbation_bbd_vars_mapping
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
                    .map(|var| self.as_symbolic_context().get_state_variable(*var))
                    .collect::<Vec<_>>();

                let state_vars_projection = RawProjection::new(perturbed_state_vars, &control_subset);
                for state_vector in state_vars_projection.iter() {
                    // This should remove all remaining perturbed state variables. Unperturbed
                    // variables should not appear in the control map add all.
                    let control_colors = control_subset.restrict(&state_vector.to_values());
                    let mut map = HashMap::new();
                    for var in &perturbed_variables {
                        let state_var = self.as_symbolic_context().get_state_variable(*var);
                        map.insert(
                            self.as_original().as_network().get_variable_name(*var).clone(),
                            state_vector.get_value(state_var).unwrap(),
                        );
                    }
                    let factor = 2.0f64.powi(
                        (self.as_symbolic_context().state_variables().len() + self.num_perturbation_parameters()) as i32
                    );

                    let working_colors = control_colors.cardinality() / factor;
                    let robustness = working_colors / all_colors_size;
                    if robustness >= best_robustness {
                        if robustness != best_robustness {
                            with_best_robustness = 0;
                        }
                        best_robustness = robustness;
                        with_best_robustness += 1;
                    }
                    println!("[{}] >>>> {:?}: {}; rho = {:.3}", perturbation_size, map, working_colors, robustness);
                }
            }

            println!("[{}] Elapsed: {}ms", perturbation_size, start.elapsed().unwrap().as_millis());

            println!("[{}] Best robustness {} for {} perturbations.", perturbation_size, best_robustness, with_best_robustness);

            control_map_all = control_map_all.union(&control_map);
            if best_robustness == 1.0 && stop_early {
                println!("Sufficient robustness achieved for perturbation size {}.", perturbation_size);
                return PhenotypeControlMap {
                    perturbation_set: control_map_all,
                    context: self.clone(),
                }
            }
        }

        println!("Sufficient robustness not achieved with perturbation size {}.", size_bound);

        return PhenotypeControlMap {
            perturbation_set: control_map_all,
            context: self.clone(),
        }
    }



    pub fn ceiled_phenotype_permanent_control_with_true_oscillation(
        &self,
        phenotype: GraphVertices,
        size_bound: usize,
        perturbation_variables: Vec<VariableId>,
        stop_early: bool
    ) -> PhenotypeControlMap {
        // A map which gives us the symbolic variable of the perturbation parameter.
        let perturbation_bbd_vars_mapping = perturbation_variables.iter()
            .filter_map(|var| self.get_perturbation_parameter(var.clone()).map(|it| (var.clone(), it)))
            .map(|(var, param)| {
                (var, self.as_symbolic_context()
                    .get_explicit_function_table(param)
                    .symbolic_variables()[0])
            })
            .collect::<HashMap<_, _>>();

        let bdd_vars = self.as_symbolic_context().bdd_variable_set();
        // The list of symbolic variables of perturbation parameters.
        let perturbation_bdd_vars = {
            let mut values = perturbation_bbd_vars_mapping
                .values()
                .cloned()
                .collect::<Vec<_>>();
            values.sort();
            values
        };

        let mut control_map_all = self.mk_empty_colored_vertices();

        for perturbation_size in 0..(size_bound + 1) {
            let start = SystemTime::now();
            println!("Perturbation size: {}", perturbation_size);
            let admissible_perturbations = mk_bdd_of_bound(bdd_vars, &perturbation_bdd_vars, perturbation_size);
            {
                let factor = 2.0f64.powi(bdd_vars.num_vars() as i32 - perturbation_bdd_vars.len() as i32);
                println!("[{}] >> Admissible fixed(Q) sets: {}", perturbation_size, admissible_perturbations.cardinality() / factor);
            }
            let admissible_perturbations = self.empty_colors().copy(admissible_perturbations);


            // oscillation with phenotype set to true
            let control_map_in_phenotype = self.phenotype_permanent_control(phenotype.clone(), admissible_perturbations.clone(), true).perturbation_set;

            // oscillation with outside of phenotype
            let outside_phenotype = self.mk_unit_colored_vertices().minus_vertices(&phenotype).vertices();
            let control_map_outside_phenotype = self.phenotype_permanent_control(outside_phenotype, admissible_perturbations, true).perturbation_set;

            let control_map = control_map_in_phenotype.intersect(&control_map_outside_phenotype);


            // Compute the number of valuations of the perturbation parameters.
            let factor = 2.0f64.powi(bdd_vars.num_vars() as i32 - perturbation_bdd_vars.len() as i32);
            let mut only_perturbation_parameters = control_map.clone().into_bdd();
            for var in bdd_vars.variables() {
                if !perturbation_bdd_vars.contains(&var) {
                    only_perturbation_parameters = only_perturbation_parameters.var_project(var);
                }
            }
            println!("[{}] >> fixed(Q) sets in control map: {}", perturbation_size, only_perturbation_parameters.cardinality() / factor);


            let control_map_bdd = control_map.clone().into_bdd();
            let perturbation_vars_projection = RawProjection::new(perturbation_bdd_vars.clone(), &control_map_bdd);

            let all_colors_size = self.unit_colors().approx_cardinality() / 2.0f64.powi(perturbation_bdd_vars.len() as i32);
            let mut best_robustness = 0.0;
            let mut with_best_robustness = 0;

            for is_perturbed_vector in perturbation_vars_projection.iter() {
                let perturbed_variables = perturbation_bbd_vars_mapping
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
                    .map(|var| self.as_symbolic_context().get_state_variable(*var))
                    .collect::<Vec<_>>();

                let state_vars_projection = RawProjection::new(perturbed_state_vars, &control_subset);
                for state_vector in state_vars_projection.iter() {
                    // This should remove all remaining perturbed state variables. Unperturbed
                    // variables should not appear in the control map add all.
                    let control_colors = control_subset.restrict(&state_vector.to_values());
                    let mut map = HashMap::new();
                    for var in &perturbed_variables {
                        let state_var = self.as_symbolic_context().get_state_variable(*var);
                        map.insert(
                            self.as_original().as_network().get_variable_name(*var).clone(),
                            state_vector.get_value(state_var).unwrap(),
                        );
                    }
                    let factor = 2.0f64.powi(
                        (self.as_symbolic_context().state_variables().len() + self.num_perturbation_parameters()) as i32
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

            control_map_all = control_map_all.union(&control_map);
            if best_robustness == 1.0 && stop_early {
                println!("Sufficient robustness achieved for perturbation size {}.", perturbation_size);
                return PhenotypeControlMap {
                    perturbation_set: control_map_all,
                    context: self.clone(),
                }
            }
        }

        println!("Sufficient robustness not achieved with perturbation size {}.", size_bound);

        return PhenotypeControlMap {
            perturbation_set: control_map_all,
            context: self.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::aeon::phentoype::build_phenotype;
    use crate::perturbation::PerturbationGraph;
    use biodivine_lib_param_bn::BooleanNetwork;
    use std::collections::HashMap;
    use std::convert::TryFrom;

    #[test]
    pub fn test_trivial_permanent_myeloid() {
        let model_string = &std::fs::read_to_string("models/myeloid_witness.aeon").unwrap();
        let model = BooleanNetwork::try_from(model_string.as_str()).unwrap();
        println!(
            "========= {}({}) =========",
            "models/myeloid_witness.aeon",
            model.num_vars()
        );
        let perturbations = PerturbationGraph::new(&model);
        let erythrocyte_phenotype = build_phenotype(
            perturbations.as_perturbed(),
            HashMap::from([("EKLF", true)]),
        );
        let control = perturbations.phenotype_permanent_control(
            erythrocyte_phenotype,
            perturbations.mk_unit_colors(),
            false
        );

        // Trivial working control
        let working_colors =
            control.perturbation_working_colors(&HashMap::from([(String::from("EKLF"), true)]));

        assert_eq!(1.0, working_colors.approx_cardinality());

        // Trivial not-working control
        let not_working_colors =
            control.perturbation_working_colors(&HashMap::from([(String::from("EKLF"), false)]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());

        // Non-trivial working control
        let working_colors = control.perturbation_working_colors(&HashMap::from([
            (String::from("Fli1"), false),
            (String::from("GATA1"), true),
            (String::from("GATA2"), true),
        ]));

        assert_eq!(1.0, working_colors.approx_cardinality());

        let not_working_colors = control.perturbation_working_colors(&HashMap::from([
            (String::from("CEBPa"), true),
            (String::from("PU1"), false),
        ]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());

        let not_working_colors = control.perturbation_working_colors(&HashMap::from([]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());
    }

    #[test]
    pub fn test_ceiled_permanent_myeloid() {
        let model_string = &std::fs::read_to_string("models/myeloid_witness.aeon").unwrap();
        let model = BooleanNetwork::try_from(model_string.as_str()).unwrap();
        println!(
            "========= {}({}) =========",
            "models/myeloid_witness.aeon",
            model.num_vars()
        );
        let mut all_vars = Vec::new();
        for v in model.variables() {
            all_vars.push(v.clone());
        }

        let perturbations = PerturbationGraph::new(&model);
        let erythrocyte_phenotype = build_phenotype(
            perturbations.as_perturbed(),
            HashMap::from([("EKLF", true)]),
        );
        let control =
            perturbations.ceiled_phenotype_permanent_control(erythrocyte_phenotype, 3, model.variables().collect(), false, false);

        // Trivial working control
        let working_colors =
            control.perturbation_working_colors(&HashMap::from([(String::from("EKLF"), true)]));

        assert_eq!(1.0, working_colors.approx_cardinality());

        // Trivial not-working control
        let not_working_colors =
            control.perturbation_working_colors(&HashMap::from([(String::from("EKLF"), false)]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());

        // Non-trivial working control
        let working_colors = control.perturbation_working_colors(&HashMap::from([
            (String::from("Fli1"), false),
            (String::from("GATA1"), true),
            (String::from("GATA2"), true),
        ]));

        assert_eq!(1.0, working_colors.approx_cardinality());

        let not_working_colors = control.perturbation_working_colors(&HashMap::from([
            (String::from("CEBPa"), true),
            (String::from("PU1"), false),
        ]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());

        let not_working_colors = control.perturbation_working_colors(&HashMap::from([]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());
    }
}
