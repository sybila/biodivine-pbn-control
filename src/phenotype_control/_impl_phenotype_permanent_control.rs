use std::collections::HashSet;
use std::hash::Hash;
use std::time::Instant;
use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::fixed_points::FixedPoints;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices,GraphVertices,GraphColors};
use biodivine_lib_param_bn::symbolic_async_graph::reachability::Reachability;
use biodivine_lib_param_bn::{VariableId, VariableIdIterator};
use chrono::{DateTime, Local};
use crate::phenotype_control::PhenotypeControlMap;
use itertools::Itertools;

impl PerturbationGraph {
    pub fn ceiled_phenotype_permanent_control(
        &self,
        phenotype: GraphVertices,
        max_size: i32,
        pert_variables: Vec<VariableId>
    ) -> PhenotypeControlMap {
        let now = Instant::now();
        println!("Starting phenotype permanent control ceiled to size {} controlling {} different variables, started at: {}", max_size, pert_variables.len(), Local::now());
        let mut admissible_perturbations = self.as_perturbed().mk_empty_colors();
        for i in 1..(max_size+1) {
            for c in pert_variables.iter().combinations(i as usize) {
                let mut perturbed = HashSet::new();
                let mut color_combination = self.as_perturbed().mk_unit_colors();
                for perturbation in c {
                    perturbed.insert(perturbation.clone());
                }
                for v in self.as_perturbed().as_network().variables() {
                    if perturbed.contains(&v) {
                        color_combination = color_combination.intersect(
                            &self.fix_perturbation(v, None).colors());
                    } else {
                        color_combination = color_combination.intersect(
                            &self.not_perturbed(v)
                        )
                    }
                }
                admissible_perturbations = admissible_perturbations.union(&color_combination);
            }
        }

        let result = self.phenotype_permanent_control(phenotype, admissible_perturbations);
        let duration = now.elapsed();
        println!("Control map computation finished at {:?} ", Local::now());
        println!("Time elapsed for computing control map: {:?}", duration);
        result
    }


    pub fn phenotype_permanent_control(
        &self,
        phenotype: GraphVertices,
        admissible_perturbations: GraphColors
    ) -> PhenotypeControlMap {
        println!("all space {}", self.as_perturbed().unit_colored_vertices().approx_cardinality());
        println!("all vertices {}", self.as_perturbed().unit_colored_vertices().vertices().approx_cardinality());
        println!("phenotype vertices {}", phenotype.approx_cardinality());

        let phenotype_violating_space = self.as_perturbed()
                                            .unit_colored_vertices()
                                            .minus_vertices(&phenotype)
                                            .intersect_colors(&admissible_perturbations);
        println!("space to explore attractors {}", phenotype_violating_space.approx_cardinality());

        let phenotype_violating_attractors = FixedPoints::symbolic(self.as_perturbed(), &phenotype_violating_space);
        println!("violating atts {}", phenotype_violating_attractors.approx_cardinality());

        let phenotype_violating_space = Reachability::reach_bwd(self.as_perturbed(), &phenotype_violating_attractors);
        println!("violating space {}", phenotype_violating_space.approx_cardinality());

        let phenotype_respecting_space = self.as_perturbed().unit_colored_vertices().minus(&phenotype_violating_space);
        println!("ok space {}", phenotype_respecting_space.approx_cardinality());

        PhenotypeControlMap {
            perturbation_set: phenotype_respecting_space,
            context: self.clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use crate::perturbation::PerturbationGraph;
    use biodivine_lib_param_bn::BooleanNetwork;
    use std::convert::TryFrom;
    use crate::aeon::phentoype::build_phenotype;


    #[test]
    pub fn test_trivial_permanent_myeloid() {
        let model_string = &std::fs::read_to_string("models/myeloid_witness.aeon").unwrap();
        let model = BooleanNetwork::try_from(model_string.as_str()).unwrap();
        println!("========= {}({}) =========", "models/myeloid_witness.aeon", model.num_vars());
        let perturbations = PerturbationGraph::new(&model);
        let erythrocyte_phenotype = build_phenotype(
                perturbations.as_perturbed(),
            HashMap::from([
                ("EKLF", true),
            ]));
        let control = perturbations.phenotype_permanent_control(
            erythrocyte_phenotype,
            perturbations.as_perturbed().mk_unit_colors()
        );


        // Trivial working control
        let working_colors = control.perturbation_working_colors(
            &HashMap::from([
                (String::from("EKLF"), true)
            ]));

        assert_eq!(1.0, working_colors.approx_cardinality());

        // Trivial not-working control
        let not_working_colors = control.perturbation_working_colors(
            &HashMap::from([
                (String::from("EKLF"), false)
            ]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());


        // Non-trivial working control
        let working_colors = control.perturbation_working_colors(
            &HashMap::from([
                (String::from("Fli1"), false),
                (String::from("GATA1"), true),
                (String::from("GATA2"), true)
            ]));

        assert_eq!(1.0, working_colors.approx_cardinality());


        let not_working_colors = control.perturbation_working_colors(
            &HashMap::from([
                (String::from("CEBPa"), true),
                (String::from("PU1"), false),
            ]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());
    }
}
