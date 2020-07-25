use biodivine_lib_param_bn::async_graph::{AsyncGraph};
use std::collections::{HashMap, HashSet};
use biodivine_lib_param_bn::bdd_params::{BddParams};
use biodivine_lib_std::param_graph::{Graph, EvolutionOperator, Params, InvertibleEvolutionOperator};
use biodivine_lib_std::{IdState};
use biodivine_aeon_server::scc::{StateSet};
use std::clone::Clone;
use std::io;
use std::io::Write;
use biodivine_aeon_server::scc::algo_reach::guarded_reach;
use std::borrow::Borrow;
use rayon::prelude::*;

fn all_possible_predecessors<F>(bwd: &F, set: &HashSet<IdState>) -> HashSet<IdState>
    where
        F: EvolutionOperator<State = IdState, Params = BddParams> + Send + Sync,
{
    return set
        .par_iter()
        .flat_map(|s| bwd.step(*s).map(|(t, _)| t).collect::<Vec<_>>())
        .collect();
}

trait FoldUnion {
    fn fold_union(self) -> Option<BddParams>;
}

impl<I> FoldUnion for I
    where
        I: Iterator<Item = Option<BddParams>>,
{
    fn fold_union(self) -> Option<BddParams> {
        return self.fold(None, |a, b| match (a, b) {
            (Some(a), Some(b)) => Some(a.union(&b)),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        });
    }
}

pub fn find_strong_basin(graph: &AsyncGraph, attractor: IdState, params: BddParams) -> HashMap<IdState, BddParams>
{
    let fwd = graph.fwd();
    let bwd = graph.bwd();
    // Just a quick sanity check to verify that the given `attractor` state is really a sink for
    // all parameters requested in `params`. If a successor has non-empty intersection with
    // `params`, then we have a problem.
    let successor_count = fwd.step(attractor).filter(|(_, p)| {
        !params.intersect(p).is_empty()
    }).count();
    if successor_count != 0 {
        panic!("Given state ({:?}) is not an attractor. It has {} successor(s).", attractor, successor_count);
    }

    let empty_params = graph.empty_params();

    let state_count = graph.states().count();
    let seed = StateSet::new_with_fun(state_count, |s| if s.eq(&attractor) { Some(params.clone()) } else { None });
    let no_guard = StateSet::new_with_initial(state_count, graph.unit_params());
    let backward_reach = guarded_reach(&bwd, &seed, &no_guard);
    let mut basin = HashMap::new();
    for (n, p) in backward_reach.iter() {
        basin.insert(n, p.clone());
    }

    println!("Weak basin has {} states.", basin.len());

    let mut changed: HashSet<IdState> = basin.keys().cloned().collect();
    loop {
        let recompute: HashSet<IdState> = all_possible_predecessors(&bwd, &changed);

        let to_remove: Vec<(&IdState, BddParams)> = recompute
            .par_iter()
            .filter_map(|t| {
                let params_t = basin.get(t);
                if let Some(params_t) = params_t {
                    let remove_from_t: Option<BddParams> = fwd
                        .step(*t)
                        .map(|(s, edge)| {
                            // Under what parameters is successor a part of the basin
                            let params_s = basin.get(&s).unwrap_or(&empty_params);

                            // The parameters which lead outside of basin. These are the params which node
                            // currently has assigned, for which there is an edge outside (and suc_params)
                            // and these are NOT in the basin (and !basin_params).
                            let difference = params_t.intersect(&edge).minus(params_s);
                            if !difference.is_empty() {
                                // There are params which lead outside the basin
                                Some(difference)
                            } else {
                                // All parameters lead to basin
                                None
                            }
                        })
                        .fold_union();
                    remove_from_t.map(|p| (t, p))
                } else {
                    // Can't remove node which is not in the basin
                    None
                }
            })
            .collect::<Vec<_>>();

        if to_remove.is_empty() {
            break;
        }

        changed.clear();
        for (node, params) in to_remove {
            *basin.get_mut(&node).unwrap() = basin[&node].minus(&params);
            changed.insert(*node);

            if basin.get(&node).unwrap().is_empty() {
                // Clean up node which does not belong to basin under any parameter
                basin.remove(&node);
            }
        }
    }

    return basin;
}
