use crate::control::ControlMap;
use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;

impl PerturbationGraph {
    /// Compute permanent control map. That is, controls which work when a perturbation is
    /// applied and then never lifted.
    pub fn permanent_control(
        &self,
        source: &ArrayBitVector,
        target: &ArrayBitVector,
    ) -> ControlMap {
        /*
           Permanent control works exactly as one-step, but in the perturbed graph instead of original.
        */
        let target_set = self.vertex(target);
        let weak_basin = crate::aeon::reachability::backward(self.as_perturbed(), &target_set);
        let strong_basin =
            crate::aeon::reachability::forward_closed(self.as_perturbed(), &weak_basin);
        let can_jump_to = self.post_perturbation(source, &strong_basin);
        ControlMap {
            perturbation_set: can_jump_to,
            context: &self,
        }
    }


    pub fn permanent_control_basin(
        &self,
        target: &ArrayBitVector,
    ) -> ControlMap {
        /*
           Permanent control works exactly as one-step, but in the perturbed graph instead of original.
        */
        let target_set = self.vertex(target);
        let weak_basin = crate::aeon::reachability::backward(self.as_perturbed(), &target_set);
        let strong_basin =
            crate::aeon::reachability::forward_closed(self.as_perturbed(), &weak_basin);
        ControlMap {
            perturbation_set: strong_basin,
            context: &self,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::perturbation::PerturbationGraph;
    use biodivine_lib_param_bn::biodivine_std::bitvector::{ArrayBitVector, BitVector};
    use biodivine_lib_param_bn::BooleanNetwork;
    use std::convert::TryFrom;

    // Test that in non-parametrised models, trivial permanent control always leads to target,
    // and that there are also other controls we can use except trivial.
    //
    // Note that because the attractors in test models are single-state, we don't have to worry
    // about attractor disappearing due to control.
    //
    // Also, due to the nature of permanent control, the trivial-ish control we considered in
    // one-step does not work here, because any trivial control must lead directly to an attractor,
    // since it always creates a new attractor in the jump state.
    fn test_trivial_permanent_control(model_file: &str) {
        let model_string =
            &std::fs::read_to_string(format!("models/{}_witness.aeon", model_file)).unwrap();
        let model = BooleanNetwork::try_from(model_string.as_str()).unwrap();
        println!("========= {}({}) =========", model_file, model.num_vars());
        let perturbations = PerturbationGraph::new(&model);

        let attractors = crate::aeon::attractors::compute(perturbations.as_original());
        // We are using first attractor as source and remaining attractors as targets.
        let source_state: ArrayBitVector = attractors[0]
            .vertices()
            .materialize()
            .iter()
            .next()
            .unwrap();
        for target in attractors.iter().skip(1) {
            let target_state = target.vertices().materialize().iter().next().unwrap();

            let control = perturbations.permanent_control(&source_state, &target_state);
            println!(
                "Control from {:?} to {:?} cardinality: {}",
                source_state,
                target_state,
                control.as_bdd().cardinality()
            );
            println!(
                "Control jump targets: {}",
                control
                    .as_colored_vertices()
                    .vertices()
                    .approx_cardinality()
            );

            // Most trivial control that goes straight into target state. There should be exactly one.
            let mut trivial_control = control.clone();
            for v in perturbations.variables() {
                trivial_control.require_perturbation(v, Some(target_state.get(usize::from(v))));
            }
            assert_eq!(1.0, trivial_control.as_bdd().cardinality());

            // A slightly less restrictive control that *requires* perturbation of all
            // variables, but they don't have to exactly match target, just something in
            // the strong basin.
            //
            // This should still be just the target vertex, because trivial permanent control
            // will always lock the state to a sink.
            let mut all_perturbed = control.clone();
            for v in perturbations.variables() {
                all_perturbed.require_perturbation(v, None);
            }
            assert_eq!(1.0, all_perturbed.as_bdd().cardinality());

            // However, there should still be some other vertices we can jump to aside from target.
            assert!(
                control
                    .as_colored_vertices()
                    .vertices()
                    .approx_cardinality()
                    > 1.0
            );
        }
    }

    #[test]
    pub fn test_trivial_permanent_myeloid() {
        test_trivial_permanent_control("myeloid");
    }

    #[test]
    pub fn test_trivial_permanent_cardiac() {
        test_trivial_permanent_control("cardiac");
    }

    #[test]
    pub fn test_trivial_permanent_erbb() {
        test_trivial_permanent_control("erbb");
    }

    /*
        #[test]
        pub fn test_trivial_permanent_all() {
            // Larger models are too slow for debug builds, but should be ok if running with --release
            for model_file in ["myeloid", "cardiac", "erbb"/*, "tumour", "mapk", "hgf"*/].iter() {
                test_trivial_permanent_control(model_file);
            }
        }
    */
}
