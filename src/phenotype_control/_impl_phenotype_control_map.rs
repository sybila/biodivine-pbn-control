use std::collections::HashMap;
use std::hash::Hash;
use std::time::Instant;
use biodivine_lib_bdd::Bdd;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, GraphColors};
use biodivine_lib_param_bn::VariableId;
use itertools::min;
use crate::phenotype_control::PhenotypeControlMap;

impl PhenotypeControlMap {
    pub fn as_bdd(&self) -> &Bdd {
        self.perturbation_set.as_bdd()
    }

    pub fn as_colored_vertices(&self) -> &GraphColoredVertices {
        &self.perturbation_set
    }

    pub fn perturbation_working_colors(&self, perturbation: &HashMap<String, bool>) -> GraphColors {
        let mut perturbation_bdd = self.perturbation_set.as_bdd().clone();
        // Obtain BDD having given variables perturbed to the specified value and remaining variables having unperturbed
        for v in self.context.as_perturbed().as_network().variables() {
            let var_name = self.context.as_perturbed().as_network().get_variable_name(v);
            if perturbation.contains_key(var_name) {
                let perturbation_value = perturbation.get(var_name).unwrap();
                let bdd_var = self.context.as_perturbed().symbolic_context().get_state_variable(v);
                // Fix states & params to the perturbation value
                perturbation_bdd = perturbation_bdd.and(&self.context.fix_perturbation(v, Some(perturbation_value.clone())).into_bdd());
                perturbation_bdd = perturbation_bdd.var_project(bdd_var);
            } else {
                // Require the variable to NOT be perturbed
                perturbation_bdd = perturbation_bdd.and(&self.context.not_perturbed(v).into_bdd());
            }
        }

        // Do universal projection across all non-perturbed variables
        for v in self.context.variables() {
            let var_name = self.context.as_perturbed().as_network().get_variable_name(v);
            if !perturbation.contains_key(var_name) {
                let var = self.context.as_symbolic_context().get_state_variable(v);
                perturbation_bdd = Bdd::fused_binary_flip_op(
                    (&perturbation_bdd, None),
                    (&perturbation_bdd, Some(var)),
                    None,
                    biodivine_lib_bdd::op_function::and
                )
            }
        }

        let colors = GraphColoredVertices::new(perturbation_bdd, &self.context.as_perturbed().symbolic_context()).colors();
        colors
    }

    // Returns all perturbations of working for at least min_cardinality colours with size up to max_size
    pub fn ceiled_size_perturbation_working_colors(&self, max_size: i32, min_cardinality: f64, controllable_vars: &Vec<VariableId>) -> Vec<HashMap<String, bool>> {
        let mut perturbations = Vec::new();
        for i in 1..(max_size+1) {
            let now = Instant::now();
            println!("Exploring perturbations of size {:?}", i);
            let mut controls = self.rec_ceiled_size_perturbation_working_colors(i, min_cardinality, controllable_vars, self.perturbation_set.as_bdd().clone(), HashMap::new());
            println!("Perturbations working for at least {:?} colors : {:?}", min_cardinality, controls.len());
            let duration = now.elapsed();
            println!("Exploring perturbations of size {:?} took {:?}", i, duration);
            perturbations.append(&mut controls);
        }
        return perturbations
    }

    fn rec_ceiled_size_perturbation_working_colors(&self, remaining_size: i32, min_cardinality: f64, controllable_vars: &Vec<VariableId>, current_perturbation_bdd: Bdd, current_perturbation: HashMap<String, bool>) -> Vec<HashMap<String, bool>> {
        if remaining_size == 0 {
            let mut result_bdd = current_perturbation_bdd.clone();
            for v in controllable_vars {
                let var_name = self.context.as_perturbed().as_network().get_variable_name(v.clone());
                if !current_perturbation.contains_key(var_name) {
                    // Require the variable to NOT be perturbed
                    result_bdd = result_bdd.and(&self.context.not_perturbed(v.clone()).into_bdd());
                }
            }

            for v in self.context.variables() {
                let var_name = self.context.as_perturbed().as_network().get_variable_name(v);
                if !current_perturbation.contains_key(var_name) {
                    let var = self.context.as_symbolic_context().get_state_variable(v);
                    result_bdd = Bdd::fused_binary_flip_op(
                        (&result_bdd, None),
                        (&result_bdd, Some(var)),
                        None,
                        biodivine_lib_bdd::op_function::and
                    )
                }
            }

            let gc = GraphColoredVertices::new(result_bdd, &self.context.as_perturbed().symbolic_context()).colors();
            // println!("{:?} {:?}", current_perturbation, gc.approx_cardinality());
            if gc.approx_cardinality() >= min_cardinality {
                return vec![current_perturbation]
            } else {
                return Vec::new()
            }
        }

        let mut values = Vec::new();
        for var in controllable_vars {
            for value in [true, false] {
                let mut new_perturbation = current_perturbation.clone();
                new_perturbation.insert(self.context.as_perturbed().as_network().get_variable_name(var.clone()).clone(), value);
                let bdd_var = self.context.as_perturbed().symbolic_context().get_state_variable(var.clone());
                let mut result_bdd = current_perturbation_bdd.clone();
                result_bdd = result_bdd.and(&self.context.fix_perturbation(var.clone(), Some(value.clone())).into_bdd());
                result_bdd = result_bdd.var_project(bdd_var);
                let mut other = self.rec_ceiled_size_perturbation_working_colors(remaining_size-1, min_cardinality, controllable_vars, result_bdd, new_perturbation);
                values.append(&mut other);
            }
        }

        return values;
    }

}
