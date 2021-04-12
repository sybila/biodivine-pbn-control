use crate::aeon::reachability::{backward, forward_closed};
use crate::control::ControlMap;
use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;

impl PerturbationGraph {
    /// Compute temporary control map. That is, controls which work when a perturbation is applied,
    /// then held "as long as necessary" and eventually released.
    pub fn temporary_control(
        &self,
        source: &ArrayBitVector,
        target: &ArrayBitVector,
    ) -> ControlMap {
        /*
           Temporary control is the most challenging, because the control jump needs to be into
           the perturbed basin of a normal basin of target.
        */
        let target_set = self.vertex(target);
        let original_weak_basin = backward(self.as_original(), &target_set);
        let original_strong_basin = forward_closed(self.as_original(), &original_weak_basin);
        let perturbed_weak_basin = backward(self.as_perturbed(), &original_strong_basin);
        let perturbed_strong_basin = forward_closed(self.as_perturbed(), &perturbed_weak_basin);
        let can_jump_and_hold = self.post_perturbation(source, &perturbed_strong_basin);
        ControlMap {
            perturbation_set: can_jump_and_hold,
            context: &self,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::aeon::reachability::{backward, forward_closed};
    use crate::perturbation::PerturbationGraph;
    use biodivine_lib_param_bn::biodivine_std::bitvector::{ArrayBitVector, BitVector};
    use biodivine_lib_param_bn::biodivine_std::traits::Set;
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
    fn test_trivial_temporary_control(model_file: &str) {
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

            let control = perturbations.temporary_control(&source_state, &target_state);
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
            // This should be the vertices of the strong basin.
            let mut all_perturbed = control.clone();
            for v in perturbations.variables() {
                all_perturbed.require_perturbation(v, None);
            }
            assert!(all_perturbed.as_bdd().cardinality() > 1.0);

            let target_set = perturbations.vertex(&target_state);
            let weak_basin = backward(perturbations.as_original(), &target_set);
            let strong_basin = forward_closed(perturbations.as_original(), &weak_basin);
            assert_eq!(
                all_perturbed.as_bdd().cardinality(),
                strong_basin.vertices().approx_cardinality()
            );

            // However, there should still be some cases we can jump to aside from target.
            // Note that this does not necessarily have to be new vertices we jump to, it may just
            // be simpler perturbations that can be used instead of the "trivial-ish".
            assert!(control.as_bdd().cardinality() > all_perturbed.as_bdd().cardinality());

            // Finally, we know temporary control should usually work in more instances than
            // one-step control. Sadly, our test models contain one instance where temporary and
            // one step are equal :/ However, we can still test that one is a subset of the other.
            let one_step_control = perturbations.one_step_control(&source_state, &target_state);
            let extra_controls = control
                .as_colored_vertices()
                .minus(one_step_control.as_colored_vertices());
            println!("Extra controls: {}", extra_controls.approx_cardinality());
            assert!(one_step_control
                .as_colored_vertices()
                .is_subset(control.as_colored_vertices()));
        }
    }

    #[test]
    pub fn test_trivial_temporary_myeloid() {
        test_trivial_temporary_control("myeloid");
    }

    #[test]
    pub fn test_trivial_temporary_cardiac() {
        test_trivial_temporary_control("cardiac");
    }

    #[test]
    pub fn test_trivial_temporary_erbb() {
        test_trivial_temporary_control("erbb");
    }

    /*
        #[test]
        pub fn test_trivial_temporary_all() {
            // Larger models are too slow for debug builds, but should be ok if running with --release
            for model_file in ["myeloid", "cardiac", "erbb"/*, "tumour", "mapk", "hgf"*/].iter() {
                test_trivial_temporary_control(model_file);
            }
        }
    */
}
