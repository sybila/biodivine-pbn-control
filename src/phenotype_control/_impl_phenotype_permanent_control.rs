use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::fixed_points::FixedPoints;
use biodivine_lib_param_bn::symbolic_async_graph::GraphColoredVertices;
use crate::aeon;
use crate::phenotype_control::PhenotypeControlMap;

impl PerturbationGraph {
    pub fn phenotype_permanent_control(
        &self,
        phenotype: GraphColoredVertices,
    ) -> PhenotypeControlMap {
        println!("all space {}", self.as_perturbed().unit_colored_vertices().approx_cardinality());
        println!("phenotype space {}", phenotype.approx_cardinality());

        let attractors = FixedPoints::symbolic(self.as_perturbed(), self.as_perturbed().unit_colored_vertices());
        println!("all atts {}", attractors.approx_cardinality());
        let phenotype_violating_attractors = attractors.minus(&phenotype);
        println!("violating atts {}", phenotype_violating_attractors.approx_cardinality());
        let phenotype_violating_space = aeon::reachability::backward(self.as_perturbed(), &phenotype_violating_attractors);
        println!("violating space {}", phenotype_violating_space.approx_cardinality());
        let phenotype_respecting_space = self.as_perturbed().unit_colored_vertices().minus(&phenotype_violating_space);
        println!("ok space {}", phenotype_respecting_space.approx_cardinality());

        PhenotypeControlMap {
            perturbation_set: phenotype_respecting_space,
            context: self.clone(),
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
            erythrocyte_phenotype
        );


        // Trivial working control
        let working_colors = control.perturbation_working_colors(
            HashMap::from([
                (String::from("EKLF"), true)
            ]));

        assert_eq!(1.0, working_colors.approx_cardinality());

        // Trivial not-working control
        let not_working_colors = control.perturbation_working_colors(
            HashMap::from([
                (String::from("EKLF"), false)
            ]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());


        // Non-trivial working control
        let working_colors = control.perturbation_working_colors(
            HashMap::from([
                (String::from("Fli1"), false),
                (String::from("GATA1"), true),
                (String::from("GATA2"), true)
            ]));

        assert_eq!(1.0, working_colors.approx_cardinality());


        let not_working_colors = control.perturbation_working_colors(
            HashMap::from([
                (String::from("CEBPa"), true),
                (String::from("PU1"), false),
            ]));

        assert_eq!(0.0, not_working_colors.approx_cardinality());
    }
}
