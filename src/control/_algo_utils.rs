use biodivine_lib_param_bn::symbolic_async_graph::{SymbolicAsyncGraph, GraphColoredVertices};
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::BooleanNetwork;
use crate::control::_algo_tgr_reduction::tgr_reduction;

pub fn reach_bwd_with_saturation(graph: &SymbolicAsyncGraph, initial: GraphColoredVertices) -> GraphColoredVertices {
    let mut result = initial;

    loop {
        let mut stop = true;
        for var in graph.as_network().variables().rev() { // The order is important to update Bdd based on the "easiest" variables first.
            let step = graph.var_pre(var, &result).minus(&result);

            if !step.is_empty() {
                result = result.union(&step);
                stop = false;
                break;
            }
        }
        if stop {
            return result;
        }
    }
}

pub fn reach_fwd_with_saturation(graph: &SymbolicAsyncGraph, initial: GraphColoredVertices) -> GraphColoredVertices {
    let mut result = initial;

    loop {
        let mut stop = true;
        for var in graph.as_network().variables().rev() { // The order is important to update Bdd based on the "easiest" variables first.
            let step = graph.var_post(var, &result).minus(&result);

            if !step.is_empty() {
                result = result.union(&step);
                stop = false;
                break;
            }
        }
        if stop {
            return result;
        }
    }
}


pub fn strong_basin(graph: &SymbolicAsyncGraph, initial: GraphColoredVertices) -> GraphColoredVertices {
    let mut basin = reach_bwd_with_saturation(graph, initial);
    loop {
        let mut stop = true;
        for var in graph.as_network().variables().rev() {
            let outside_successors = graph.var_post(var, &basin).minus(&basin);
            let can_go_out = graph.var_pre(var, &outside_successors).intersect(&basin);

            if !can_go_out.is_empty() {
                basin = basin.minus(&can_go_out);
                stop = false;
                break;
            }
        }
        if stop {
            return basin;
        }
    }
}

pub fn add_auto_regulations(model: BooleanNetwork) -> BooleanNetwork {
    let mut result = model.as_graph().clone();
    for v in result.variables() {
        if result.find_regulation(v, v).is_none() {
            let name = result.get_variable_name(v).clone();
            result.add_regulation(name.as_str(), name.as_str(), false, None);
        }
    }

    let mut result = BooleanNetwork::new(result);
    for v in result.variables() {
        // Technically, the models should have equivalent variable ids!
        if let Some(function) = model.get_update_function(v) {
            result.add_update_function(v, function.clone());
        }
    }

    result
}

pub fn find_attractors(graph: &SymbolicAsyncGraph) -> Vec<GraphColoredVertices> {
    let mut result = Vec::new();
    let mut universe = tgr_reduction(graph, graph.mk_unit_colored_vertices());
    while !universe.is_empty(){
        let pivot = universe.pick_vertex();
        let fwd = reach_fwd_with_saturation(graph, pivot.clone());
        let bwd = reach_bwd_with_saturation(graph, pivot.clone());
        let scc = fwd.intersect(&bwd);
        let not_attractor_colors = fwd.minus(&scc).colors();
        let attractor = scc.minus_colors(&not_attractor_colors);
        if !attractor.is_empty() {
            result.push(attractor);
        }
        universe = universe.minus(&bwd);
    }

    return result;
}


pub fn find_reachable_attractors(graph: &SymbolicAsyncGraph, pivot: GraphColoredVertices) -> GraphColoredVertices{
    let pivot = graph.unit_colored_vertices().pick_vertex();
    let fwd = reach_fwd_with_saturation(graph, pivot.clone());
    let bwd = reach_bwd_with_saturation(graph, pivot.clone());
    let scc = fwd.intersect(&bwd);
    let not_attractor_colors = fwd.minus(&scc).colors();
    let attractor = scc.minus_colors(&not_attractor_colors);
    return attractor;
}
