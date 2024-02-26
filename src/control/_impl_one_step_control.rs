use crate::control::AttractorControlMap;
use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;
use biodivine_lib_param_bn::symbolic_async_graph::GraphColors;
use itertools::Itertools;

impl PerturbationGraph {
    /// Compute one-step control map. That is, controls which work by applying the perturbation
    /// immediately for a single time step.
    pub fn one_step_control(
        &self,
        source: &ArrayBitVector,
        target: &ArrayBitVector,
        compute_params: &GraphColors,
        verbose: bool
    ) -> AttractorControlMap {
        /*
           To eventually stabilize in target, we have to reach its strong basin using a
           perturbation. We thus first compute the basin and then compute which perturbations
           can jump into the basin from source.

           Note that colors where target is not in an attractor will be eliminated using the
           strong basin procedure.
        */
        let allowed_colors = self
            .as_perturbed()
            .transfer_colors_from(
                self.as_non_perturbable().unit_colors(),
                &self.as_non_perturbable(),
            )
            .unwrap();

        let target_set = self.vertex(target).intersect_colors(compute_params).intersect_colors(&allowed_colors);
        let weak_basin = crate::aeon::reachability::backward(self.as_original(), &target_set, verbose);
        let strong_basin =
            crate::aeon::reachability::forward_closed(self.as_original(), &weak_basin, verbose);
        let can_jump_to = self.post_perturbation(source, &strong_basin);
        AttractorControlMap {
            perturbation_set: can_jump_to,
            context: self.clone(),
            perturbation_variables: self.variables().collect_vec()
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use crate::perturbation::PerturbationGraph;
    use biodivine_lib_param_bn::biodivine_std::bitvector::{ArrayBitVector, BitVector};
    use biodivine_lib_param_bn::BooleanNetwork;
    use std::convert::TryFrom;
    use std::iter::zip;
    use crate::control::ControlMap;

    // Test that in non-parametrised models, trivial one-step control always leads to target,
    // and that we can also reach the whole strong basin using "trivial-ish" control.
    fn test_trivial_one_step_control(model_file: &str) {
        let model_string =
            &std::fs::read_to_string(format!("models/{}_witness.aeon", model_file)).unwrap();
        let model = BooleanNetwork::try_from(model_string.as_str()).unwrap();
        println!("========= {}({}) =========", model_file, model.num_vars());
        let perturbations = PerturbationGraph::new(&model);

        let attractors = crate::aeon::attractors::compute(perturbations.as_original(), false);
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
                false
            );
            println!(
                "Control from {:?} to {:?} cardinality: {}",
                source_state,
                target_state,
                control.as_bdd().cardinality()
            );

            // Most trivial control that goes straight into target state.
            let mut pert = HashMap::new();
            for (v, val) in zip(perturbations.variables(), target_state.values()) {
                let v_name = perturbations.as_perturbed().get_variable_name(v);
                pert.insert(v_name, val);
            }
            assert_eq!(1.0, control.perturbation_working_colors(&pert).approx_cardinality());

            // // Control number should be the same as the strong basin size
            // let mut all_working = control.working_perturbations(1.0, false, true);
            //
            // let target_set = perturbations.vertex(&target_state);
            // let weak_basin =
            //     crate::aeon::reachability::backward(perturbations.as_original(), &target_set, false);
            // let strong_basin =
            //     crate::aeon::reachability::forward_closed(perturbations.as_original(), &weak_basin, false);
            // assert_eq!(
            //     all_working.len(),
            //     strong_basin.vertices().approx_cardinality() as usize
            // );
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
