use crate::control::ControlMap;
use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;
use biodivine_lib_param_bn::symbolic_async_graph::GraphColors;

impl PerturbationGraph {
    /// Compute one-step control map. That is, controls which work by applying the perturbation
    /// immediately for a single time step.
    pub fn one_step_control(
        &self,
        source: &ArrayBitVector,
        target: &ArrayBitVector,
        compute_params: &GraphColors,
    ) -> ControlMap {
        /*
           To eventually stabilize in target, we have to reach its strong basin using a
           perturbation. We thus first compute the basin and then compute which perturbations
           can jump into the basin from source.

           Note that colors where target is not in an attractor will be eliminated using the
           strong basin procedure.
        */
        let target_set = self.vertex(target).intersect_colors(compute_params);
        let weak_basin = crate::aeon::reachability::backward(self.as_original(), &target_set);
        let strong_basin =
            crate::aeon::reachability::forward_closed(self.as_original(), &weak_basin);
        let can_jump_to = self.post_perturbation(source, &strong_basin);
        ControlMap {
            perturbation_set: can_jump_to,
            context: self.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::perturbation::PerturbationGraph;
    use biodivine_lib_param_bn::biodivine_std::bitvector::{ArrayBitVector, BitVector};
    use biodivine_lib_param_bn::BooleanNetwork;
    use std::convert::TryFrom;

    // Test that in non-parametrised models, trivial one-step control always leads to target,
    // and that we can also reach the whole strong basin using "trivial-ish" control.
    fn test_trivial_one_step_control(model_file: &str) {
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

            let control = perturbations.one_step_control(
                &source_state,
                &target_state,
                perturbations.unit_colors(),
            );
            println!(
                "Control from {:?} to {:?} cardinality: {}",
                source_state,
                target_state,
                control.as_bdd().cardinality()
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
            // This should actually match the size of the strong basin, because every state
            // there will have its own "trivial" control.
            let mut all_perturbed = control.clone();
            for v in perturbations.variables() {
                all_perturbed.require_perturbation(v, None);
            }
            println!(
                "All perturbed control: {}",
                all_perturbed.as_bdd().cardinality()
            );
            assert!(all_perturbed.as_bdd().cardinality() > 1.0);

            let target_set = perturbations.vertex(&target_state);
            let weak_basin =
                crate::aeon::reachability::backward(perturbations.as_original(), &target_set);
            let strong_basin =
                crate::aeon::reachability::forward_closed(perturbations.as_original(), &weak_basin);
            assert_eq!(
                all_perturbed.as_bdd().cardinality(),
                strong_basin.vertices().approx_cardinality()
            );
        }
    }

    #[test]
    pub fn test_trivial_one_step_myeloid() {
        test_trivial_one_step_control("myeloid");
    }

    #[test]
    pub fn test_trivial_one_step_cardiac() {
        test_trivial_one_step_control("cardiac");
    }

    #[test]
    pub fn test_trivial_one_step_erbb() {
        test_trivial_one_step_control("erbb");
    }

    /*
    #[test]
    pub fn test_trival_one_step_all() {
        // Larger models are too slow for debug builds, but should be ok if running with --release
        for model_file in ["myeloid", "cardiac", "erbb"/*, "tumour", "mapk", "hgf"*/].iter() {
            test_trivial_one_step_control(model_file);
        }
    }
    */
}
