use biodivine_lib_param_bn::symbolic_async_graph::{SymbolicAsyncGraph, GraphColoredVertices};
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::{BooleanNetwork, FnUpdate};
use crate::control::_algo_tgr_reduction::tgr_reduction;
use std::convert::TryFrom;

/// Compute the coloured set of all backward reachable states from the `initial` set.
///
/// The algorithm applies individual variable transitions in the order which should affect
/// the structure of the BDD the least, i.e. the application of individual transitions should
/// be much more efficient as if they were just applied in a row.
pub fn reach_bwd_with_saturation(graph: &SymbolicAsyncGraph, initial: GraphColoredVertices) -> GraphColoredVertices {
    let mut result = initial;

    loop {
        let mut stop = true;
        // The order is important to update Bdd based on the "easiest" variables first.
        for var in graph.as_network().variables().rev() {
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

/// Compute the coloured set of all forward reachable states from the `initial` set.
///
/// The algorithm applies individual variable transitions in the order which should affect
/// the structure of the BDD the least, i.e. the application of individual transitions should
/// be much more efficient as if they were just applied in a row.
pub fn reach_fwd_with_saturation(graph: &SymbolicAsyncGraph, initial: GraphColoredVertices) -> GraphColoredVertices {
    let mut result = initial;

    loop {
        let mut stop = true;
        // The order is important to update Bdd based on the "easiest" variables first.
        for var in graph.as_network().variables().rev() {
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

/// Compute the so called *strong basin* of the given `initial` set. Strong basin is the largest
/// forward-closed subset of the vertices that can reach the `initial` set.
///
/// A forward-closed set must contain all forward reachable vertices of every vertex it contains.
/// In particular, attractor is forward closed. Any vertex together with all paths to its
/// reachable attractors is also forward-closed. Essentially, a forward-closed set must always
/// contain some attractors and states that can only reach these attractors.
///
/// In particular, if `initial` is some attractor vertex, the result contains the attractor plus
/// the vertices that can only reach this particular attractor.
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

/// Create a copy of the given `model`, but add an unspecified auto-regulation to each variable
/// (that does not have an autoregulation already). Additionally, convert every implicit update
/// function to and explicit parameter.
pub fn add_auto_regulations(model: &BooleanNetwork) -> BooleanNetwork {
    let mut result = model.as_graph().clone();
    // Copy variables and regulations.
    for v in result.variables() {
        if result.find_regulation(v, v).is_none() {
            let name = result.get_variable_name(v).clone();
            result.add_regulation(name.as_str(), name.as_str(), false, None).unwrap();
        }
    }

    let mut result = BooleanNetwork::new(result);

    // Copy parameters.
    for p in model.parameters() {
        let parameter = &model[p];
        result.add_parameter(parameter.get_name(), parameter.get_arity()).unwrap();
    }

    // Copy update functions.
    for v in result.variables() {
        // Technically, the models should have equivalent variable ids!
        if let Some(function) = model.get_update_function(v) {
            result.add_update_function(v, function.clone()).unwrap();
        } else {    // Create an explicit parameter to replace the implicit function.
            let regulators = result.regulators(v);
            let parameter = result.add_parameter(
                format!("update_{}", result.get_variable_name(v)).as_str(),
                    u32::try_from(regulators.len()).unwrap(),
            ).unwrap();
            result.add_update_function(v, FnUpdate::Param(parameter, regulators)).unwrap();
        }
    }

    result
}

/// Use transition guided reduction and Xie-Beerel algorithm to uncover all attractors
/// of the given graph.
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