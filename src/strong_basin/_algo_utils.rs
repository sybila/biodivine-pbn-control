use crate::controlled_async_graph::ControlledAsyncGraph;
use biodivine_lib_param_bn::async_graph::{AsyncGraph, DefaultEdgeParams};
use biodivine_lib_param_bn::bdd_params::BddParams;
use biodivine_lib_std::param_graph::{EvolutionOperator, Graph, Params};
use biodivine_lib_std::IdState;

pub fn get_all_params_with_attractor(graph: &ControlledAsyncGraph, state: IdState) -> BddParams {
    let fwd = graph.fwd();
    let successors = fwd.step(state);

    let mut bad_params = graph.empty_params().clone();
    for (_, par) in successors {
        //if !succ.eq(&state) {
        bad_params = bad_params.union(&par);
        //}
    }
    bad_params = bad_params.intersect(graph.unit_params());

    return graph.unit_params().minus(&bad_params);
}

pub fn find_attractors(graph: &AsyncGraph<DefaultEdgeParams>) -> Vec<IdState> {
    let fwd = graph.fwd();
    let nodes = graph.num_states();
    let mut attractors = Vec::new();

    for n in 0..nodes {
        let state = IdState::from(n);
        let has_successor = fwd
            .step(state)
            .fold(graph.empty_params().clone(), |a, (_, b)| a.union(&b));
        if has_successor.is_empty() {
            attractors.push(state);
        }
    }

    return attractors;
}
