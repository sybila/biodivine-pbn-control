use crate::phenotype_control::PhenotypeControlMap;
use biodivine_lib_bdd::Bdd;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, GraphColors};
use biodivine_lib_param_bn::VariableId;
use std::collections::HashMap;
use std::time::Instant;

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
        for v in self.context.as_perturbed().variables() {
            // println!("{:?}", GraphColoredVertices::new(perturbation_bdd.clone(), &self.context.as_perturbed().symbolic_context()).colors().approx_cardinality());
            let var_name = self.context.as_perturbed().get_variable_name(v);
            if perturbation.contains_key(var_name.as_str()) {
                // println!("{:?}", var_name);
                let perturbation_value = perturbation.get(var_name.as_str()).unwrap();
                let bdd_var = self
                    .context
                    .as_perturbed()
                    .symbolic_context()
                    .get_state_variable(v);
                // Fix states & params to the perturbation value
                perturbation_bdd = perturbation_bdd.and(
                    &self
                        .context
                        .fix_perturbation(v, Some(*perturbation_value))
                        .into_bdd(),
                );
                perturbation_bdd = perturbation_bdd.var_exists(bdd_var);
            } else {
                // Require the variable to NOT be perturbed
                perturbation_bdd = perturbation_bdd.and(&self.context.not_perturbed(v).into_bdd());
            }
        }

        // let mut file = File::create("./dot_result.txt").unwrap();
        // file.write_all(perturbation_bdd.to_dot_string(self.context.as_symbolic_context().bdd_variable_set(), true).as_bytes()).unwrap();

        // Do universal projection across all non-perturbed variables
        for v in self.context.variables() {
            // println!("{:?}", GraphColoredVertices::new(perturbation_bdd.clone(), &self.context.as_perturbed().symbolic_context()).colors().approx_cardinality());
            let var_name = self.context.as_perturbed().get_variable_name(v);
            if !perturbation.contains_key(var_name.as_str()) {
                let var = self.context.as_symbolic_context().get_state_variable(v);
                perturbation_bdd = Bdd::fused_binary_flip_op(
                    (&perturbation_bdd, None),
                    (&perturbation_bdd, Some(var)),
                    None,
                    biodivine_lib_bdd::op_function::and,
                )
            }
        }

        let colors = GraphColoredVertices::new(
            perturbation_bdd,
            self.context.as_perturbed().symbolic_context(),
        )
        .colors();
        colors
    }

    // Returns all perturbations of working for at least min_cardinality colours with size up to max_size
    pub fn ceiled_size_perturbation_working_colors(
        &self,
        max_size: usize,
        min_cardinality: f64,
        controllable_vars: &Vec<VariableId>,
        stop_early: bool,
        verbose: bool,
    ) -> Vec<HashMap<String, bool>> {
        let mut perturbations = Vec::new();
        let mut minimal_found = false;
        for i in 1..(max_size + 1) {
            let now = Instant::now();
            println!("Exploring perturbations of size {:?}", i);
            let mut controls = self.rec_ceiled_size_perturbation_working_colors(
                i,
                min_cardinality,
                controllable_vars,
                self.perturbation_set.as_bdd().clone(),
                HashMap::new(),
                verbose,
            );
            println!(
                "Perturbations working for at least {:?} colors : {:?}",
                min_cardinality,
                controls.len()
            );
            let duration = now.elapsed();
            println!(
                "Exploring perturbations of size {:?} took {:?}",
                i, duration
            );

            perturbations.append(&mut controls);

            if !perturbations.is_empty() && !minimal_found {
                println!("Minaml perturbations: {:?}", perturbations);
                minimal_found = true;
                if stop_early {
                    return perturbations;
                }
            }
        }

        perturbations
    }

    fn rec_ceiled_size_perturbation_working_colors(
        &self,
        max_size: usize,
        min_cardinality: f64,
        controllable_vars: &Vec<VariableId>,
        current_perturbation_bdd: Bdd,
        current_perturbation: HashMap<String, bool>,
        verbose: bool,
    ) -> Vec<HashMap<String, bool>> {
        if current_perturbation.len() == max_size {
            let mut result_bdd = current_perturbation_bdd.clone();
            for v in controllable_vars {
                let var_name = self.context.as_perturbed().get_variable_name(*v);
                if !current_perturbation.contains_key(var_name.as_str()) {
                    // Require the variable to NOT be perturbed
                    result_bdd = result_bdd.and(&self.context.not_perturbed(*v).into_bdd());
                }
            }

            for v in self.context.variables() {
                let var_name = self.context.as_perturbed().get_variable_name(v);
                if !current_perturbation.contains_key(var_name.as_str()) {
                    let var = self.context.as_symbolic_context().get_state_variable(v);
                    result_bdd = Bdd::fused_binary_flip_op(
                        (&result_bdd, None),
                        (&result_bdd, Some(var)),
                        None,
                        biodivine_lib_bdd::op_function::and,
                    )
                }
            }

            let gc = GraphColoredVertices::new(
                result_bdd,
                self.context.as_perturbed().symbolic_context(),
            )
            .colors();

            return if gc.approx_cardinality() >= min_cardinality {
                if verbose {
                    // println!("{:?}: {:?}", current_perturbation, gc.to_dot_string(self.context.as_symbolic_context()));
                    println!("{:?}: ", current_perturbation);
                }
                vec![current_perturbation]
            } else {
                Vec::new()
            };
        }

        let mut values = Vec::new();
        let mut found_vars = Vec::new();
        let mut pointer_ok = current_perturbation.is_empty();
        for var in controllable_vars {
            let var_name = self.context.as_perturbed().get_variable_name(*var);
            if pointer_ok {
                for value in [true, false] {
                    let mut new_perturbation = current_perturbation.clone();
                    new_perturbation.insert(var_name.to_string(), value);
                    let bdd_var = self
                        .context
                        .as_perturbed()
                        .symbolic_context()
                        .get_state_variable(*var);
                    let mut result_bdd = current_perturbation_bdd.clone();
                    result_bdd = result_bdd
                        .and(&self.context.fix_perturbation(*var, Some(value)).into_bdd());
                    result_bdd = result_bdd.var_exists(bdd_var);
                    let mut other = self.rec_ceiled_size_perturbation_working_colors(
                        max_size,
                        min_cardinality,
                        controllable_vars,
                        result_bdd,
                        new_perturbation,
                        verbose,
                    );
                    values.append(&mut other);
                }
            } else if current_perturbation.contains_key(var_name.as_str()) {
                found_vars.push(var_name);
                if current_perturbation.len() == found_vars.len() {
                    pointer_ok = true;
                }
            }
        }

        values
    }
}
