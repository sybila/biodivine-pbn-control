use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, SymbolicAsyncGraph};

/// Compute the coloured set of all backward reachable states from the `initial` set.
pub fn backward(
    graph: &SymbolicAsyncGraph,
    initial: &GraphColoredVertices,
) -> GraphColoredVertices {
    let mut result = initial.clone();

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
        if cfg!(feature = "print_progress") && result.as_bdd().size() > 100_000 {
            println!("Backward progress: {}", result.as_bdd().size())
        }
        if stop {
            return result;
        }
    }
}

/// Compute the coloured set of all forward reachable states from the `initial` set.
pub fn forward(graph: &SymbolicAsyncGraph, initial: &GraphColoredVertices) -> GraphColoredVertices {
    let mut result = initial.clone();

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
        if cfg!(feature = "print_progress") && result.as_bdd().size() > 100_000 {
            println!("Forward progress: {}", result.as_bdd().size())
        }
        if stop {
            return result;
        }
    }
}

/// Compute the largest forward closed subset of the `initial` set. That is, the states that can
/// only reach states inside `initial`.
///
/// In particular, if the initial set is a weak basin, the result is a strong basin.
pub fn forward_closed(
    graph: &SymbolicAsyncGraph,
    initial: &GraphColoredVertices,
) -> GraphColoredVertices {
    let mut basin = initial.clone();
    loop {
        let mut stop = true;
        for var in graph.as_network().variables().rev() {
            let can_go_out = graph.var_can_post_out(var, &basin);
            if !can_go_out.is_empty() {
                basin = basin.minus(&can_go_out);
                stop = false;
                break;
            }
        }
        if cfg!(feature = "print_progress") && basin.as_bdd().size() > 100_000 {
            println!("Forward closed progress: {}", basin.as_bdd().size())
        }
        if stop {
            return basin;
        }
    }
}
