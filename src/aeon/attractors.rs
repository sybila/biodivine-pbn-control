use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, SymbolicAsyncGraph};

/// Use transition guided reduction and Xie-Beerel algorithm to uncover all attractors
/// of the given graph.
pub fn compute(graph: &SymbolicAsyncGraph) -> Vec<GraphColoredVertices> {
    let mut result = Vec::new();
    let mut universe = super::tgr::reduction(graph, graph.mk_unit_colored_vertices());
    while !universe.is_empty() {
        let pivot = universe.pick_vertex();
        let fwd = super::reachability::forward(graph, &pivot);
        let bwd = super::reachability::backward(graph, &pivot);
        let scc = fwd.intersect(&bwd);
        let not_attractor_colors = fwd.minus(&scc).colors();
        let attractor = scc.minus_colors(&not_attractor_colors);
        if !attractor.is_empty() {
            result.push(attractor);
        }
        universe = universe.minus(&bwd);
    }

    result
}