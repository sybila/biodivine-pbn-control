#[cfg(test)]
mod tests {
    use biodivine_lib_param_bn::BooleanNetwork;
    use std::convert::TryFrom;
    use std::fs;
    use biodivine_lib_std::{IdState};
    use crate::strong_basin::_algo_utils::get_all_params_with_attractor;
    use crate::strong_basin::_algo_sb_parallel_fixed_point::find_strong_basin;
    use crate::async_graph_with_control::AsyncGraphWithControl;
    use biodivine_aeon_server::scc::StateSet;

    #[test]
    fn test_witness() {
        let aeon_str: &str = &fs::read_to_string("models/g2b_2stable_attractors.aeon").unwrap();
        let model = BooleanNetwork::try_from(aeon_str).unwrap();

        let graph = &AsyncGraphWithControl::new(model);
        // We want to start from state
        // CcrM | CtrA | DnaA | GcrA | SciP
        //  1   |  1   |  1   |  0   |  0
        // Note that bits are reversed because first variable corresponds to least significant bit.
        let state = IdState::from(0b00111 as usize);

        let relevant_params = get_all_params_with_attractor(graph, state);
        let seed = &StateSet::new_with_fun(graph.num_states(), |s| if s.eq(&state) { Some(relevant_params.clone()) } else { None });
        let basin = find_strong_basin(graph, seed);

        assert_eq!(basin.len(), 16);
    }

    #[test]
    fn test_full_model() {
        let aeon_str: &str = &fs::read_to_string("models/g2b.aeon").unwrap();
        let model = BooleanNetwork::try_from(aeon_str).unwrap();

        let graph = &AsyncGraphWithControl::new(model);
        // We want to start from state
        // CcrM | CtrA | DnaA | GcrA | SciP
        //  1   |  1   |  1   |  0   |  0
        // Note that bits are reversed because first variable corresponds to least significant bit.
        let state = IdState::from(0b00111 as usize);

        let relevant_params = get_all_params_with_attractor(graph, state);
        let seed = &StateSet::new_with_fun(graph.num_states(), |s| if s.eq(&state) { Some(relevant_params.clone()) } else { None });
        let basin = find_strong_basin(graph, seed);

        assert_eq!(basin.len(), 32);
    }
}