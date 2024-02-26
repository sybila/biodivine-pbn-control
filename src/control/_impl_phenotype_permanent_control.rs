use crate::aeon::reachability::backward_within;
use crate::perturbation::PerturbationGraph;

use crate::control::{ControlMap, PhenotypeControlMap, PhenotypeOscillationType};
use biodivine_lib_bdd::BddVariable;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColors, GraphVertices};
use biodivine_lib_param_bn::VariableId;
use std::collections::HashMap;
use std::time::SystemTime;

impl PerturbationGraph {
    pub fn phenotype_permanent_control(
        &self,
        phenotype: GraphVertices,
        oscillation: PhenotypeOscillationType,
        verbose: bool,
    ) -> PhenotypeControlMap {
        self.phenotype_permanent_control_internal(
            phenotype,
            self.as_original().mk_unit_colors(),
            oscillation,
            verbose,
        )
    }

    fn phenotype_permanent_control_internal(
        &self,
        phenotype: GraphVertices,
        admissible_colors_perturbations: GraphColors,
        oscillation: PhenotypeOscillationType,
        verbose: bool,
    ) -> PhenotypeControlMap {
        let start = SystemTime::now();

        // println!("PERTURBED {}",self.as_perturbed().unit_colors().as_bdd().num_vars());
        //
        // for v in self.as_perturbed().symbolic_context().parameter_variables() {
        //     println!("{}", self.as_perturbed().symbolic_context().bdd_variable_set().name_of(*v));
        // }
        //
        // println!("UNPERTURBED {}", self.as_non_perturbable().unit_colors().as_bdd().num_vars());
        // for v in self.as_non_perturbable().symbolic_context().parameter_variables() {
        //     println!("{}", self.as_non_perturbable().symbolic_context().bdd_variable_set().name_of(*v));
        // }
        //

        let allowed_colors = self
            .as_perturbed()
            .transfer_colors_from(
                self.as_non_perturbable().unit_colors(),
                &self.as_non_perturbable(),
            )
            .unwrap()
            .intersect(&admissible_colors_perturbations);

        if verbose {
            println!(">>>>>>>>>>>>> Computing phenotype control with allowed perturbations of cardinality {}", allowed_colors.approx_cardinality())
        }
        let perturbation_bbd_vars_mapping =
            self.get_perturbation_bdd_mapping(self.perturbable_variables());
        let bdd_vars = self.as_symbolic_context().bdd_variable_set();

        // The list of symbolic variables of perturbation parameters.
        let perturbation_bdd_vars = Self::get_perturbation_bdd_vars(&perturbation_bbd_vars_mapping);
        let mut trap;
        let control_universe = self
            .unit_colored_vertices()
            .intersect_colors(&allowed_colors);
        let phenotype_coloured_vertices = control_universe.intersect_vertices(&phenotype);
        match oscillation {
            PhenotypeOscillationType::Required => {
                return self.phenotype_permanent_control_with_true_oscillation_internal(
                    phenotype,
                    allowed_colors,
                    verbose,
                )
            }
            PhenotypeOscillationType::Allowed => {
                let bwd_reach = backward_within(
                    self.as_perturbed(),
                    &phenotype_coloured_vertices,
                    &control_universe,
                    verbose
                );
                trap = bwd_reach
            }
            PhenotypeOscillationType::Forbidden => {
                trap = phenotype_coloured_vertices;
            }
        }

        'trap: loop {
            for var in self.variables().rev() {
                let can_leave = self.as_perturbed().var_can_post_out(var, &trap);
                if !can_leave.is_empty() {
                    trap = trap.minus(&can_leave);
                    if trap.symbolic_size() > 100_000 {
                        println!(
                            " Trap phenotype progress: {} / {}",
                            trap.symbolic_size(),
                            trap.approx_cardinality()
                        );
                    }
                    continue 'trap;
                }
            }
            break;
        }

        if verbose {
            println!("Trap cardinality {}", trap.approx_cardinality());
        }

        if verbose {
            println!(
                "Cardinality of phenotype trap set: {}",
                trap.approx_cardinality()
            )
        }

        let mut trap = self
            .unit_colored_vertices()
            .intersect_colors(&admissible_colors_perturbations)
            .minus(&trap);

        'trap: loop {
            for var in self.variables().rev() {
                let can_leave = self.as_perturbed().var_can_post_out(var, &trap);
                if !can_leave.is_empty() {
                    trap = trap.minus(&can_leave);
                    if trap.symbolic_size() > 100_000 {
                        println!(
                            " Trap non-phenotype progress: {} / {}",
                            trap.symbolic_size(),
                            trap.approx_cardinality()
                        );
                    }
                    continue 'trap;
                }
            }
            break;
        }

        if verbose {
            println!("Inverse trap cardinality {}", trap.approx_cardinality());
        }

        if verbose {
            println!(
                "Cardinality of inversed trap set: {}",
                trap.approx_cardinality()
            )
        }

        let mut inverse_control = trap.into_bdd();
        for var in self.variables() {
            let state_var = self.as_symbolic_context().get_state_variable(var);
            if let Some(perturbation_var) = perturbation_bbd_vars_mapping.get(&var) {
                // If the variable can be perturbed, we split into two cases and eliminate
                // it in the unperturbed cases.

                let is_perturbed = inverse_control.var_select(*perturbation_var, true);
                let is_not_perturbed = inverse_control
                    .var_select(*perturbation_var, false)
                    .var_exists(state_var);
                inverse_control = is_perturbed.or(&is_not_perturbed);
            } else {
                // If the variable cannot be perturbed, we can eliminate it everywhere.
                inverse_control = inverse_control.var_exists(state_var);
            }
        }

        let inverse_control_map = self.empty_colored_vertices().copy(inverse_control);
        if verbose {
            println!(
                "Inverse control cardinality {}",
                inverse_control_map.approx_cardinality()
            );
        }

        if verbose {
            println!(
                "Cardinality of inversed control map: {}",
                inverse_control_map.approx_cardinality()
            )
        }

        // Control map consists of admissible state-color pairs that are also admissible
        // for our perturbation size and are not in the inverse map.
        let control_map = self
            .unit_colored_vertices()
            .intersect_colors(&admissible_colors_perturbations)
            .minus(&inverse_control_map);

        if verbose {
            println!(
                "Cardinality of control map: {}",
                control_map.approx_cardinality()
            )
        }

        let result = PhenotypeControlMap {
            perturbation_variables: self.perturbable_variables().clone(),
            perturbation_set: control_map.clone(),
            context: self.clone(),
        };

        if verbose {
            // Compute the number of valuations of the perturbation parameters.
            let _factor =
                2.0f64.powi(bdd_vars.num_vars() as i32 - perturbation_bdd_vars.len() as i32);
            let mut only_perturbation_parameters = control_map.into_bdd();
            for var in bdd_vars.variables() {
                if !perturbation_bdd_vars.contains(&var) {
                    only_perturbation_parameters = only_perturbation_parameters.var_exists(var);
                }
            }

            result.working_perturbations(0.0, true, false);
            println!(">>>>>> Elapsed: {}ms", start.elapsed().unwrap().as_millis());
        }

        result
    }

    pub fn phenotype_permanent_control_of_specific_size(
        &self,
        phenotype: GraphVertices,
        perturbation_size: usize,
        allow_oscillation: PhenotypeOscillationType,
        verbose: bool,
    ) -> PhenotypeControlMap {
        let admissible_perturbations = self.create_perturbation_colors(perturbation_size, verbose);
        let control_map = self
            .phenotype_permanent_control_internal(
                phenotype.clone(),
                admissible_perturbations,
                allow_oscillation,
                verbose,
            )
            .perturbation_set;

        let map = PhenotypeControlMap {
            perturbation_variables: self.perturbable_variables().clone(),
            perturbation_set: control_map,
            context: self.clone(),
        };

        return map;
    }

    pub fn ceiled_phenotype_permanent_control(
        &self,
        phenotype: GraphVertices,
        size_bound: usize,
        allow_oscillation: PhenotypeOscillationType,
        stop_early: bool,
        verbose: bool,
    ) -> PhenotypeControlMap {
        // A map which gives us the symbolic variable of the perturbation parameter.
        let mut control_map_all = self.mk_empty_colored_vertices();

        for perturbation_size in 0..(size_bound + 1) {
            let _start = SystemTime::now();
            let admissible_perturbations =
                self.create_perturbation_colors(perturbation_size, verbose);

            // perturbation_graph.as_perturbed().transfer_colors_from(normal_graph.unit_colors(), &normal_graph).unwrap();

            let control_map = self
                .phenotype_permanent_control_internal(
                    phenotype.clone(),
                    admissible_perturbations,
                    allow_oscillation.clone(),
                    verbose,
                )
                .perturbation_set;

            control_map_all = control_map_all.union(&control_map);

            let working_perturbations = PhenotypeControlMap {
                perturbation_variables: self.perturbable_variables().clone(),
                perturbation_set: control_map,
                context: self.clone(),
            }.working_perturbations(1.0, verbose, false);

            if verbose || stop_early {
                let mut best_robustness = 0.0;
                for (_perturbation, working_colors) in working_perturbations {
                    let robustness = working_colors.approx_cardinality()
                        / self.as_non_perturbable().unit_colors().approx_cardinality();
                    if robustness > best_robustness {
                        best_robustness = robustness
                    }
                }

                if verbose {
                    println!(
                        "Robustness {} achieved for perturbation size {}.",
                        best_robustness, perturbation_size
                    );
                }

                if best_robustness >= 1.0 && stop_early {
                    if verbose {
                        println!(
                            "Robustness {} achieved for perturbation size {}.",
                            best_robustness, perturbation_size
                        );
                    }
                    return PhenotypeControlMap {
                        perturbation_variables: self.perturbable_variables().clone(),
                        perturbation_set: control_map_all,
                        context: self.clone(),
                    };
                }
            }
        }

        return PhenotypeControlMap {
            perturbation_variables: self.perturbable_variables().clone(),
            perturbation_set: control_map_all,
            context: self.clone(),
        };
    }

    fn phenotype_permanent_control_with_true_oscillation_internal(
        &self,
        phenotype: GraphVertices,
        admissible_colors: GraphColors,
        verbose: bool,
    ) -> PhenotypeControlMap {
        // oscillation with phenotype set to true
        let control_map_in_phenotype = self
            .phenotype_permanent_control_internal(
                phenotype.clone(),
                admissible_colors.clone(),
                PhenotypeOscillationType::Allowed,
                verbose,
            )
            .perturbation_set;

        // oscillation with outside of phenotype
        let outside_phenotype = self
            .mk_unit_colored_vertices()
            .minus_vertices(&phenotype)
            .vertices();
        let control_map_outside_phenotype = self
            .phenotype_permanent_control_internal(
                outside_phenotype,
                admissible_colors,
                PhenotypeOscillationType::Allowed,
                verbose,
            )
            .perturbation_set;

        let control_map = control_map_in_phenotype.intersect(&control_map_outside_phenotype);

        return PhenotypeControlMap {
            perturbation_variables: self.perturbable_variables().clone(),
            perturbation_set: control_map,
            context: self.clone(),
        };
    }

    pub fn get_perturbation_bdd_vars(
        perturbation_bbd_vars_mapping: &HashMap<VariableId, BddVariable>,
    ) -> Vec<BddVariable> {
        let perturbation_bdd_vars = {
            let mut values = perturbation_bbd_vars_mapping
                .values()
                .cloned()
                .collect::<Vec<_>>();
            values.sort();
            values
        };
        perturbation_bdd_vars
    }

    pub fn get_perturbation_bdd_mapping(
        &self,
        perturbation_variables: &Vec<VariableId>,
    ) -> HashMap<VariableId, BddVariable> {
        let perturbation_bbd_vars_mapping = perturbation_variables
            .iter()
            .filter_map(|var| {
                self.get_perturbation_parameter(var.clone())
                    .map(|it| (var.clone(), it))
            })
            .map(|(var, param)| {
                (
                    var,
                    self.as_symbolic_context()
                        .get_explicit_function_table(param)
                        .symbolic_variables()[0],
                )
            })
            .collect::<HashMap<_, _>>();
        perturbation_bbd_vars_mapping
    }
}

#[cfg(test)]
mod tests {
    use crate::aeon::phentoype::build_phenotype;
    use crate::perturbation::PerturbationGraph;
    use crate::control::_impl_phenotype_permanent_control::PhenotypeOscillationType;
    use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
    use biodivine_lib_param_bn::BooleanNetwork;
    use std::collections::HashMap;
    use std::convert::TryFrom;
    use crate::control::ControlMap;

    #[test]
    pub fn test_standard_permanent_myeloid() {
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
            PhenotypeOscillationType::Forbidden,
            false,
        );

        let working_perturbations = control.working_perturbations(1.0, false, false);

        // println!("{:?}", control.working_perturbations(1.0, true));

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

        let perturbations = PerturbationGraph::new(&model);
        let erythrocyte_phenotype = build_phenotype(
            perturbations.as_perturbed(),
            HashMap::from([("EKLF", true)]),
        );
        let control = perturbations.ceiled_phenotype_permanent_control(
            erythrocyte_phenotype,
            3,
            PhenotypeOscillationType::Forbidden,
            false,
            true,
        );

        let working_perturbations = control.working_perturbations(1.0, false, false);
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
    pub fn test_restricted_vars_permanent_myeloid() {
        let model_string = &std::fs::read_to_string("models/myeloid_witness.aeon").unwrap();
        let model = BooleanNetwork::try_from(model_string.as_str()).unwrap();
        println!(
            "========= {}({}) =========",
            "models/myeloid_witness.aeon",
            model.num_vars()
        );

        let mut perturbable_vars = vec![];
        for v in model.variables() {
            if model.get_variable_name(v) != "EKLF" {
                perturbable_vars.push(v);
            }
        }

        let perturbations = PerturbationGraph::with_restricted_variables(&model, perturbable_vars);
        let erythrocyte_phenotype = build_phenotype(
            perturbations.as_perturbed(),
            HashMap::from([("EKLF", true)]),
        );
        let control = perturbations.ceiled_phenotype_permanent_control(
            erythrocyte_phenotype,
            1,
            PhenotypeOscillationType::Forbidden,
            false,
            true,
        );

        println!("{:?}", control.working_perturbations(1.0, false, false).len());
        // Trivial working control
        let working_colors =
            control.perturbation_working_colors(&HashMap::from([(String::from("EKLF"), true)]));

        assert_eq!(0.0, working_colors.approx_cardinality());
        assert_eq!(0, control.working_perturbations(1.0, false, false).len());
    }

    #[test]
    pub fn test_implicit_parameter_compatibility() {
        let bn = BooleanNetwork::try_from(
            r"
            v_1 -> v_2
            v_3 -| v_2
            v_2 -> v_3
            $v_3: v_2
        ",
        )
        .unwrap();

        let stg = SymbolicAsyncGraph::new(&bn).unwrap();
        let perturbed = PerturbationGraph::new(&bn);

        assert!(perturbed
            .as_perturbed()
            .transfer_colors_from(stg.unit_colors(), &stg)
            .is_some());
    }
}
