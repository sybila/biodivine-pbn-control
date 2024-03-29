use std::collections::HashMap;
use std::iter::zip;
use biodivine_lib_bdd::Bdd;
use crate::control::{AttractorControlMap, ControlMap};
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, GraphColors};
use biodivine_lib_param_bn::VariableId;

use itertools::Itertools;
use crate::perturbation::PerturbationGraph;

impl ControlMap for AttractorControlMap {
    fn new(context: PerturbationGraph, perturbation_set: GraphColoredVertices) -> Self {
        return AttractorControlMap {
            perturbation_variables: context.perturbable_variables().clone(),
            context,
            perturbation_set
        };
    }

    fn as_bdd(&self) -> &Bdd {
        self.perturbation_set.as_bdd()
    }

    fn as_colored_vertices(&self) -> &GraphColoredVertices {
        &self.perturbation_set
    }

    fn working_perturbations(&self, min_robustness: f64, _verbose: bool, return_all: bool) -> Vec<(HashMap<String, bool>, GraphColors)> {
        if min_robustness < 0.0 || min_robustness > 1.0 {
            panic!("Min robustness must be in range between 0.0 and 1.0")
        }
        let mut max_robustness = 0.0;
        let mut best_control_size = 0;
        let all_colors_size = self
            .context
            .as_non_perturbable()
            .unit_colors()
            .approx_cardinality();
        let mut perturbations = vec![] ;
        let var_names = self.clone().perturbation_variables.into_iter().map(|v| { self.context.as_original().as_network().unwrap().get_variable_name(v.clone()).clone() });
        for combination in var_names.powerset() {
            if combination.contains(&"Fli1".to_string()) && combination.contains(&"PU1".to_string()) && combination.contains(&"GATA1".to_string()) {
                println!("OK");
            }
            if max_robustness >= min_robustness && combination.len() > best_control_size && !return_all {
                return perturbations
            }
            for true_values in combination.clone().iter().powerset()
            {
                let perturbation = HashMap::from_iter(combination.clone().into_iter().map(|var| {(var.clone(), true_values.contains(&&var))}));
                let control_colors = self.perturbation_working_colors(&perturbation);

                let robustness = control_colors.approx_cardinality() / all_colors_size;
                if max_robustness < robustness {
                    best_control_size = perturbation.len();
                    max_robustness = robustness;
                }
                if robustness > min_robustness {
                    perturbations.push((perturbation.clone(), control_colors.clone()));
                }
            }
        }
        return perturbations

    }

    fn perturbation_working_colors(&self, perturbation: &HashMap<String, bool>) -> GraphColors {
        let mut perturbed = self.clone();

        for v in self.context.variables() {
            let var_name = self
                .context
                .as_perturbed()
                .as_network()
                .unwrap()
                .get_variable_name(v.clone());

            if perturbation.contains_key(var_name) {
                perturbed.require_perturbation(v.clone(), perturbation.get(var_name));
            } else {
                perturbed.exclude_perturbation(v.clone(), None);
            }
        }
        perturbed.controllable_colors()
    }

}

impl AttractorControlMap {
    /// Remove from this control map any results that *do not* perturb `variable`.
    /// If `value` is given, only keep perturbations which result in this value.
    pub fn require_perturbation(&mut self, variable: VariableId, value: Option<&bool>) {
        let require = self.context.fix_perturbation(variable, value);
        self.perturbation_set = self.perturbation_set.intersect(&require);
    }

    /// Remove from this control map any results that perturb `variable`. If `value` is given,
    /// only remove perturbations which result in this value.
    pub fn exclude_perturbation(&mut self, variable: VariableId, value: Option<&bool>) {
        let exclude = self.context.fix_perturbation(variable, value);
        self.perturbation_set = self.perturbation_set.minus(&exclude);
    }

    pub fn as_bdd(&self) -> &Bdd {
        self.perturbation_set.as_bdd()
    }

    pub fn as_colored_vertices(&self) -> &GraphColoredVertices {
        &self.perturbation_set
    }

     fn controllable_colors(&self) -> GraphColors {
         let bdd_context = self.context.as_symbolic_context();
         let mut bdd = self.perturbation_set.colors().into_bdd();
         for v in self.context.variables() {
             let parameter = self.context.get_perturbation_parameter(v);
             if let Some(parameter) = parameter {
                 for (_, bdd_var) in bdd_context.get_explicit_function_table(parameter) {
                     // There should be only one because arity is zero
                     bdd = bdd.var_exists(bdd_var);
                 }
             }
         }

         // Now just fix the control parameters and variables to true so that the final cardinality is correct.
         for s_var in bdd_context.state_variables() {
             bdd = bdd.var_for_all(*s_var);
         }
         for v in self.context.variables() {
             let parameter = self.context.get_perturbation_parameter(v);
             if let Some(parameter) = parameter {
                 for (_, bdd_var) in bdd_context.get_explicit_function_table(parameter) {
                     // There should be only one because arity is zero
                     bdd = bdd.var_for_all(bdd_var);
                 }
             }
         }

        self.context.as_non_perturbable().transfer_colors_from(&self.context.as_perturbed().empty_colors().copy(bdd), self.context.as_perturbed()).unwrap()
    }

    /// Compute the set of original colors (without perturbations) that can be controlled
    /// using this control map.
    pub fn controllable_colors_cardinality(&self) -> f64 {
        self.controllable_colors().approx_cardinality()
    }

    /// Compute the number of vertices the source can jump to due to different perturbations.
    pub fn jump_vertices(&self) -> f64 {
        self.perturbation_set.vertices().approx_cardinality()
    }
}
