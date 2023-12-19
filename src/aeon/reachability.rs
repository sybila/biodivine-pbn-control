use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, SymbolicAsyncGraph};

/// Compute the coloured set of all backward reachable states from the `initial` set that
/// all paths are within the given `bounds` set.
pub fn backward_within(
    graph: &SymbolicAsyncGraph,
    initial: &GraphColoredVertices,
    bounds: &GraphColoredVertices,
) -> GraphColoredVertices {
    assert!(initial.is_subset(bounds));
    let mut result = initial.clone();

    loop {
        let mut stop = true;
        for var in graph.variables().rev() {
            let step = graph.var_pre(var, &result).intersect(bounds).minus(&result);

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

/// Compute the coloured set of all backward reachable states from the `initial` set.
pub fn backward(
    graph: &SymbolicAsyncGraph,
    initial: &GraphColoredVertices,
) -> GraphColoredVertices {
    let mut result = initial.clone();

    loop {
        let mut stop = true;
        // The order is important to update Bdd based on the "easiest" variables first.
        for var in graph.variables().rev() {
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
        for var in graph.variables().rev() {
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

/// Compute the coloured set of all forward reachable states from the `initial` set such that
/// the paths are all within the given `bounds` set.
pub fn forward_within(
    graph: &SymbolicAsyncGraph,
    initial: &GraphColoredVertices,
    bounds: &GraphColoredVertices,
) -> GraphColoredVertices {
    let mut result = initial.clone();

    loop {
        let mut stop = true;
        for var in graph.variables().rev() {
            let step = graph
                .var_post(var, &result)
                .intersect(bounds)
                .minus(&result);

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

/// Compute the colored set of all forward reachable states from `initial` within `bounds` where
/// the full set is also forward closed. I.e. colors where the forward reachable set is not
/// closed within `bounds` will be removed completely.
pub fn forward_closed_within(
    graph: &SymbolicAsyncGraph,
    initial: &GraphColoredVertices,
    bounds: &GraphColoredVertices,
) -> GraphColoredVertices {
    let mut result = initial.clone();

    loop {
        let mut stop = true;
        for var in graph.variables().rev() {
            let step = graph.var_post(var, &result).minus(&result);

            if !step.is_empty() {
                let goes_outside = step.minus(bounds).colors();
                result = result.union(&step).minus_colors(&goes_outside);
                stop = false;
                break;
            }
        }
        if cfg!(feature = "print_progress") && result.as_bdd().size() > 100_000 {
            println!("Forward-closed-within progress: {}", result.as_bdd().size())
        }
        if stop {
            return result;
        }
    }
}

/// Remove all colors for which the `set` isn't backward closed.
pub fn backward_closed_subset(
    graph: &SymbolicAsyncGraph,
    set: &GraphColoredVertices,
) -> GraphColoredVertices {
    let mut result = set.clone();

    for var in graph.variables().rev() {
        let bwd_step = graph.var_pre(var, &result).minus(&result);
        result = result.minus_colors(&bwd_step.colors());
    }

    result
}

/// Remove all colors for which the `set` isn't forward closed.
pub fn forward_closed_subset(
    graph: &SymbolicAsyncGraph,
    set: &GraphColoredVertices,
) -> GraphColoredVertices {
    let mut result = set.clone();

    for var in graph.variables().rev() {
        let fwd_step = graph.var_post(var, &result).minus(&result);
        result = result.minus_colors(&fwd_step.colors());
    }

    result
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
        for var in graph.variables().rev() {
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
