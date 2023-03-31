use std::collections::HashMap;
use biodivine_lib_bdd::Bdd;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, GraphColors};
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

}
