use biodivine_lib_param_bn::symbolic_async_graph::{SymbolicAsyncGraph, GraphColoredVertices};
use biodivine_lib_param_bn::biodivine_std::traits::Set;

/// Returns a subset of the state space that is guaranteed to contain all attractors of the graph.
///
/// This is based on the interleaved-transition-guided-reduction, but it has no interleaving and
/// complete elimination of variables is removed.
pub fn tgr_reduction(
    graph: &SymbolicAsyncGraph,
    mut universe: GraphColoredVertices,
) -> GraphColoredVertices {
    for var in graph.as_network().variables() {
        let var_can_post = graph.var_can_post(var, &universe);
        let reach_from_post = crate::control::_algo_utils::reach_fwd_with_saturation(graph, var_can_post.clone()).intersect(&universe);

        // Remove basin of the reachable area.
        if reach_from_post != universe {
            let reach_basin = crate::control::_algo_utils::reach_bwd_with_saturation(graph, reach_from_post.clone()).intersect(&universe)
                .minus(&reach_from_post);
            if !reach_basin.is_empty() {
                universe = universe.minus(&reach_basin);
            }
        }

        let post_extended_component =
            crate::control::_algo_utils::reach_bwd_with_saturation(graph, var_can_post.clone()).intersect(&reach_from_post);
        let bottom_region = reach_from_post.minus(&post_extended_component);

        // Remove basin of the bottom area.
        if !bottom_region.is_empty() {
            let bottom_basin = crate::control::_algo_utils::reach_bwd_with_saturation(graph, bottom_region.clone()).intersect(&universe)
                .minus(&bottom_region);
            if !bottom_basin.is_empty() {
                universe = universe.minus(&bottom_basin);
            }
        }
    }
    universe
}