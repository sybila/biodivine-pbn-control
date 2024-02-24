use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, SymbolicAsyncGraph};

/// Returns a subset of the state space that is guaranteed to contain all attractors of the graph.
///
/// This is based on the interleaved-transition-guided-reduction, but it has no interleaving and
/// complete elimination of variables is removed.
pub fn reduction(
    graph: &SymbolicAsyncGraph,
    mut universe: GraphColoredVertices,
    verbose: bool
) -> GraphColoredVertices {
    for var in graph.variables() {
        let var_can_post = graph.var_can_post(var, &universe);
        let reach_from_post =
            super::reachability::forward(graph, &var_can_post, verbose).intersect(&universe);

        // Remove basin of the reachable area.
        if reach_from_post != universe {
            let reach_basin = super::reachability::backward(graph, &reach_from_post, false)
                .intersect(&universe)
                .minus(&reach_from_post);
            if !reach_basin.is_empty() {
                universe = universe.minus(&reach_basin);
            }
        }

        let post_extended_component =
            super::reachability::backward(graph, &var_can_post, false).intersect(&reach_from_post);
        let bottom_region = reach_from_post.minus(&post_extended_component);

        // Remove basin of the bottom area.
        if !bottom_region.is_empty() {
            let bottom_basin = super::reachability::backward(graph, &bottom_region, false)
                .intersect(&universe)
                .minus(&bottom_region);
            if !bottom_basin.is_empty() {
                universe = universe.minus(&bottom_basin);
            }
        }
    }
    universe
}
