#[cfg(test)]
mod tests {
    use biodivine_lib_param_bn::BooleanNetwork;
    use std::convert::TryFrom;
    use std::fs;
    use biodivine_lib_std::{IdState};
    use biodivine_aeon_server::scc::StateSet;
    use crate::controlled_async_graph::ControlledAsyncGraph;
    use crate::async_graph_with_control::AsyncGraphWithControl;
    use crate::strong_basin_control_experimental::_algo_utils::get_all_params_with_attractor;

    #[test]
    fn test_witness_permanent() {
        let aeon_str: &str = &fs::read_to_string("models/4vars_test.aeon").unwrap();
        let model = BooleanNetwork::try_from(aeon_str).unwrap();

        let mut graph = AsyncGraphWithControl::new(model);
        // We want to start from state
        // CcrM | CtrA | DnaA | GcrA | SciP
        //  0   |  0   |  0   |  0   |  1
        // And find control into
        // CcrM | CtrA | DnaA | GcrA | SciP
        //  1   |  0   |  1   |  1   |  0
        // Note that bits are reversed because first variable corresponds to least significant bit.
        let source = IdState::from(0b0000 as usize);
        let target = IdState::from(0b1111 as usize);

        let relevant_params = get_all_params_with_attractor(&graph, target);
        let target_attractor = &StateSet::new_with_fun(graph.num_states(), |s| if s.eq(&target) { Some(relevant_params.clone()) } else { None });
        let controls = graph.find_permanent_control(source, target_attractor);

        assert_eq!(controls.len(), 5);
    }


    #[test]
    fn test_witness_temporary() {
        let aeon_str: &str = &fs::read_to_string("models/4vars_test.aeon").unwrap();
        let model = BooleanNetwork::try_from(aeon_str).unwrap();

        let mut graph = AsyncGraphWithControl::new(model);
        // We want to start from state
        // CcrM | CtrA | DnaA | GcrA | SciP
        //  0   |  0   |  0   |  0   |  1
        // And find control into
        // CcrM | CtrA | DnaA | GcrA | SciP
        //  1   |  0   |  1   |  1   |  0
        // Note that bits are reversed because first variable corresponds to least significant bit.
        let source = IdState::from(0b0000 as usize);
        let target = IdState::from(0b1100 as usize);

        let relevant_params = get_all_params_with_attractor(&graph, target);
        let target_attractor = &StateSet::new_with_fun(graph.num_states(), |s| if s.eq(&target) { Some(relevant_params.clone()) } else { None });
        let controls = graph.find_temporary_control(source, target_attractor);

        assert_eq!(controls.len(), 6);
    }
}
