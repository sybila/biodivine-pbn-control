use std::collections::HashMap;
use std::iter::zip;
use crate::control::{AttractorControlMap, ControlMap, PhenotypeControlMap};
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, GraphColors};
use biodivine_lib_param_bn::VariableId;
use itertools::Itertools;
use crate::perturbation::PerturbationGraph;

impl ControlMap for AttractorControlMap {
    fn new(context: PerturbationGraph, perturbation_set: GraphColoredVertices, ) -> Self {
        return AttractorControlMap {
            perturbation_variables: context.perturbable_variables().clone(),
            context,
            perturbation_set,
        };
    }

    fn as_bdd(&self) -> &biodivine_lib_bdd::Bdd {
        self.perturbation_set.as_bdd()
    }

    fn as_colored_vertices(&self) -> &GraphColoredVertices {
        &self.perturbation_set
    }

    fn working_perturbations(&self, min_robustness: f64, verbose: bool) -> Vec<(HashMap<String, bool>, GraphColors)> {
        let mut perturbations = vec![];
        let mut max_achieved_rob = 0.0;
        for i in 0..self.perturbation_variables.len() {
            for combination in self.clone().perturbation_variables.into_iter().combinations(i) {
                for variation in (0..i).map(|_| vec![true, false]).multi_cartesian_product() {
                    let mut perturbation = HashMap::new();
                    for (var, val) in zip(combination.clone(), variation) {
                        let var_name = self.context.as_perturbed().get_variable_name(var);
                        perturbation.insert(var_name, val);
                    }
                    let colors = self.perturbation_working_colors(&perturbation);
                    if colors.clone().approx_cardinality() > min_robustness {
                        perturbations.push((perturbation, colors.clone()));
                    }
                    if colors.clone().approx_cardinality() > max_achieved_rob {
                        max_achieved_rob = colors.approx_cardinality();
                    }
                }
            }
            if max_achieved_rob == 1.0 {
                break;
            }
        }
        return perturbations
    }

    fn perturbation_working_colors(&self, perturbation: &HashMap<String, bool>) -> GraphColors {
        let mut controlled_cols = self.perturbation_set.clone();

        for v in self.context.as_perturbed().variables() {
            // println!("{:?}", GraphColoredVertices::new(perturbation_bdd.clone(), &self.context.as_perturbed().symbolic_context()).colors().approx_cardinality());
            let var_name = self
                .context
                .as_perturbed()
                .as_network()
                .unwrap()
                .get_variable_name(v);

            if perturbation.contains_key(var_name) {
                let require = self.context.fix_perturbation(v, perturbation.get(var_name));
                controlled_cols = controlled_cols.intersect(&require);
            }
        }

        return controlled_cols.colors();
    }

}
