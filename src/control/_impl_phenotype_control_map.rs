use crate::perturbation::PerturbationGraph;
use crate::control::{ControlMap, PhenotypeControlMap};
use biodivine_lib_bdd::Bdd;
use biodivine_lib_param_bn::symbolic_async_graph::projected_iteration::RawProjection;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, GraphColors};
use std::collections::HashMap;

impl ControlMap for PhenotypeControlMap {
    fn new(
        context: PerturbationGraph,
        perturbation_set: GraphColoredVertices,
    ) -> PhenotypeControlMap {
        return PhenotypeControlMap {
            perturbation_variables: context.perturbable_variables().clone(),
            context,
            perturbation_set,
        };
    }

    fn as_bdd(&self) -> &Bdd {
        self.perturbation_set.as_bdd()
    }

    fn as_colored_vertices(&self) -> &GraphColoredVertices {
        &self.perturbation_set
    }

    fn working_perturbations(
        &self,
        min_robustness: f64,
        verbose: bool,
    ) -> Vec<(HashMap<String, bool>, GraphColors)> {
        if min_robustness < 0.0 || min_robustness > 1.0 {
            panic!("Min robustness must be in range between 0.0 and 1.0")
        }

        let perturbation_bbd_vars_mapping = self
            .context
            .get_perturbation_bdd_mapping(&self.perturbation_variables);
        let perturbation_bdd_vars =
            PerturbationGraph::get_perturbation_bdd_vars(&perturbation_bbd_vars_mapping);

        let control_map_bdd = self.as_bdd();
        let perturbation_vars_projection =
            RawProjection::new(perturbation_bdd_vars.clone(), &control_map_bdd);

        let all_colors_size = self
            .context
            .as_non_perturbable()
            .unit_colors()
            .approx_cardinality();
        let mut best_robustness = 0.0;
        let mut with_best_robustness = 0;

        let mut result = vec![];

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

            let perturbed_state_vars = perturbed_variables
                .iter()
                .map(|var| self.context.as_symbolic_context().get_state_variable(*var))
                .collect::<Vec<_>>();

            let state_vars_projection = RawProjection::new(perturbed_state_vars, &control_subset);
            for state_vector in state_vars_projection.iter() {
                // This should remove all remaining perturbed state variables. Unperturbed
                // variables should not appear in the control map add all.
                let control_colors = control_subset.restrict(&state_vector.to_values());
                let converted_colors = self
                    .context
                    .as_non_perturbable()
                    .transfer_colors_from(
                        &self
                            .context
                            .as_original()
                            .empty_colors()
                            .copy(control_colors),
                        self.context.as_original(),
                    )
                    .unwrap();
                let mut map = HashMap::new();
                for var in &perturbed_variables {
                    let state_var = self.context.as_symbolic_context().get_state_variable(*var);
                    map.insert(
                        self.context
                            .as_original()
                            .as_network()
                            .unwrap()
                            .get_variable_name(*var)
                            .clone(),
                        state_vector.get_value(state_var).unwrap(),
                    );
                }
                let working_colors = converted_colors.approx_cardinality();
                let robustness = working_colors / all_colors_size;
                if robustness >= best_robustness {
                    if robustness != best_robustness {
                        with_best_robustness = 0;
                    }
                    best_robustness = robustness;
                    with_best_robustness += 1;
                }
                if robustness >= min_robustness {
                    result.push((map.clone(), converted_colors))
                }

                if verbose {
                    println!(">>> {:?}: {}; rho = {:.3}", map, working_colors, robustness);
                }
            }
        }
        if verbose {
            println!(
                "Best robustness {} observed in {} perturbations.",
                best_robustness, with_best_robustness
            );
        }
        result
    }

    fn perturbation_working_colors(&self, perturbation: &HashMap<String, bool>) -> GraphColors {
        let mut perturbation_bdd = self.perturbation_set.as_bdd().clone();
        // Obtain BDD having given variables perturbed to the specified value and remaining variables having unperturbed
        for v in self.context.as_perturbed().variables() {
            // println!("{:?}", GraphColoredVertices::new(perturbation_bdd.clone(), &self.context.as_perturbed().symbolic_context()).colors().approx_cardinality());
            let var_name = self
                .context
                .as_perturbed()
                .as_network()
                .unwrap()
                .get_variable_name(v);
            let bdd_var = self
                .context
                .as_perturbed()
                .symbolic_context()
                .get_state_variable(v);
            if perturbation.contains_key(var_name) {
                // println!("{:?}", var_name);
                let perturbation_value = perturbation.get(var_name).unwrap();
                // Fix states & params to the perturbation value
                perturbation_bdd = perturbation_bdd.and(
                    &self
                        .context
                        .fix_perturbation(v, Some(&perturbation_value))
                        .into_bdd(),
                );
                perturbation_bdd = perturbation_bdd.var_exists(bdd_var);
            } else {
                // Discard non-perturbed var
                perturbation_bdd = perturbation_bdd.and(&self.context.not_perturbed(v).into_bdd());
                perturbation_bdd = perturbation_bdd.var_for_all(bdd_var);
            }
        }

        // TODO switch to following, when projection works
        // let colors = self.context.as_non_perturbable().transfer_colors_from(&self.context.as_original().empty_colors().copy(perturbation_bdd), self.context.as_original()).unwrap();
        let colors = self
            .context
            .as_original()
            .empty_colors()
            .copy(perturbation_bdd);
        colors
    }
}
