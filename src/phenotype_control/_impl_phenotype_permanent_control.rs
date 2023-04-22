use crate::aeon::attractors;
use crate::perturbation::PerturbationGraph;
use crate::phenotype_control::PhenotypeControlMap;
use biodivine_lib_bdd::{Bdd, BddPartialValuation, BddVariable};
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::fixed_points::FixedPoints;
use biodivine_lib_param_bn::symbolic_async_graph::reachability::Reachability;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColors, GraphVertices};
use biodivine_lib_param_bn::{ParameterId, VariableId};
use chrono::Local;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

impl PerturbationGraph {
    pub fn ceiled_phenotype_permanent_control(
        &self,
        phenotype: GraphVertices,
        max_size: usize,
        perturbation_variables: Vec<VariableId>,
        attractor_search_method: &str
    ) -> PhenotypeControlMap {
        assert!(!perturbation_variables.is_empty());

        let now = Instant::now();
        println!("Starting phenotype permanent control ceiled to size {} controlling {} different variables, started at: {}", max_size, perturbation_variables.len(), Local::now());
        // A color set which will (eventually) hold the perturbations over
        // `perturbation_variables` up to a certain size.
        let mut admissible_perturbations = self.as_perturbed().mk_empty_colors();

        // Then, make a map which gives us a BDD variable of the "perturbation parameter"
        // for each network variable.
        let symbolic_context = self.as_symbolic_context();
        let mut perturbation_bdd_vars = HashMap::new();
        for var in self.variables() {
            if let Some(p) = self.get_perturbation_parameter(var) {
                let table = symbolic_context.get_explicit_function_table(p);
                assert_eq!(0, table.arity);
                perturbation_bdd_vars.insert(var, table.symbolic_variables()[0]);
            }
        }

        // A partial valuation that is false exactly for all perturbation parameters
        // (i.e. this is the "non-perturbed" color).
        let all_false_valuation = perturbation_bdd_vars
            .values()
            .map(|v| (*v, false))
            .collect::<Vec<_>>();
        let all_false_valuation = BddPartialValuation::from_values(&all_false_valuation);

        // Resolve perturbation_variables into their BDD counterparts.
        let bdd_perturbation_variables = perturbation_variables
            .iter()
            .map(|var| {
                *perturbation_bdd_vars
                    .get(var)
                    .expect("Variable does not have a perturbation parameter.")
            })
            .collect::<Vec<_>>();

        let mut admissible_bdd = symbolic_context.bdd_variable_set().mk_conjunctive_clause(&all_false_valuation);
        for _i in 0..max_size {
            let mut new_bdd = admissible_bdd.clone();
            for var in &bdd_perturbation_variables {
                // First, filter perturbations in `admissible_bdd` to those which do not perturb
                // `var`. Then flip the perturbation parameter such that the perturbations in the
                // resulting BDD actually do perturb `var`.
                // In other words, take perturbations of size `k` where `var` is not perturbed,
                // and then perturb `var` in all of them, resulting in perturbations of size `k+1`.
                let var_not_perturbed = symbolic_context.bdd_variable_set().mk_literal(*var, false);
                let also_perturbs_var = Bdd::fused_binary_flip_op(
                    (&admissible_bdd, None),
                    (&var_not_perturbed, None),
                    Some(*var),
                    biodivine_lib_bdd::op_function::and
                );
                new_bdd = new_bdd.or(&also_perturbs_var);
            }
            admissible_bdd = new_bdd;
        }

        let colors = self.empty_colors().copy(admissible_bdd);
        admissible_perturbations = admissible_perturbations.union(&colors);

        let result = self.phenotype_permanent_control(phenotype, admissible_perturbations, attractor_search_method);
        let duration = now.elapsed();
        println!("Control map computation finished at {:?} ", Local::now());
        println!("Time elapsed for computing control map: {:?}", duration);
        result
    }

    pub fn phenotype_permanent_control(
        &self,
        phenotype: GraphVertices,
        admissible_perturbations: GraphColors,
        attractor_search_method: &str
    ) -> PhenotypeControlMap {
        println!(
            "all space {}",
            self.as_perturbed()
                .unit_colored_vertices()
                .approx_cardinality()
        );
        println!(
            "all vertices {}",
            self.as_perturbed()
                .unit_colored_vertices()
                .vertices()
                .approx_cardinality()
        );
        println!("phenotype vertices {}", phenotype.approx_cardinality());


        let selected_attractor_search_method;
        if attractor_search_method  == "heuristic" {
            selected_attractor_search_method= self.get_attractor_type_in_unperturbed_network()
        } else {
            selected_attractor_search_method = attractor_search_method;
        }

        let mut phenotype_violating_attractors = self.mk_empty_colored_vertices();
        if selected_attractor_search_method == "sinks" {
            println!("------- Using fixed points implementation");
            let phenotype_violating_space = self
                .as_perturbed()
                .unit_colored_vertices()
                .minus_vertices(&phenotype)
                .intersect_colors(&admissible_perturbations);
            println!(
                "space to explore attractors {}",
                phenotype_violating_space.approx_cardinality()
            );
            phenotype_violating_attractors =
                FixedPoints::symbolic(self.as_perturbed(), &phenotype_violating_space);
        } else if selected_attractor_search_method == "complex" {
            println!("------- Using all attractors implementation");
            let complex_attractors = attractors::compute_restricted(self.as_original(), self.mk_unit_colored_vertices().intersect_colors(&admissible_perturbations));
            for ca in complex_attractors {
                let states_in_ca_but_not_phenotype = ca.minus_vertices(&phenotype);
                let colors_with_states_outside_phenotype = states_in_ca_but_not_phenotype.colors();
                let relevant_atts = ca.intersect_colors(&colors_with_states_outside_phenotype);
                let violating_attractors = ca.intersect(&relevant_atts);
                phenotype_violating_attractors =
                    phenotype_violating_attractors.union(&violating_attractors);
            }
        } else {
            panic!("Unknown attractor search method {:?}", attractor_search_method);
        }

        println!(
            "violating atts cardinality {}",
            phenotype_violating_attractors.approx_cardinality()
        );

        let now = Instant::now();

        let phenotype_violating_space =
            Reachability::reach_bwd(self.as_perturbed(), &phenotype_violating_attractors);
        println!(
            "violating space {}",
            phenotype_violating_space.approx_cardinality()
        );

        let phenotype_respecting_space = self
            .as_perturbed()
            .unit_colored_vertices()
            .minus(&phenotype_violating_space);
        println!(
            "ok space {}",
            phenotype_respecting_space.approx_cardinality()
        );

        print!("space computation took: {:?}", now.elapsed());

        PhenotypeControlMap {
            perturbation_set: phenotype_respecting_space,
            context: self.clone(),
        }
    }

    fn get_attractor_type_in_unperturbed_network(&self) -> &str{
        let selected_attractor_search_method: &str;

        let unperturbed_attractors = attractors::compute(self.as_original());
        let mut unperturbed_attractors_all = self.as_original().mk_empty_vertices();
        for ua in unperturbed_attractors {
            unperturbed_attractors_all = unperturbed_attractors_all.union(&ua);
        }
        let unperturbed_attractors_fps = FixedPoints::symbolic(
            self.as_original(),
            self.as_original().unit_colored_vertices(),
        );

        println!(
            "FPs {:?}",
            unperturbed_attractors_fps.vertices().approx_cardinality()
        );
        println!(
            "ALL {:?}",
            unperturbed_attractors_all.vertices().approx_cardinality()
        );
        if unperturbed_attractors_all
            .vertices()
            .is_subset(&unperturbed_attractors_fps.vertices()) {
            "sinks"
        } else {
            "complex"
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
            perturbations.as_perturbed().mk_unit_colors(),
            "sinks"
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
            perturbations.ceiled_phenotype_permanent_control(erythrocyte_phenotype, 3, all_vars, "sinks");

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
