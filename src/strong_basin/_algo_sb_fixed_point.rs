use biodivine_lib_param_bn::async_graph::{AsyncGraph};
use std::collections::HashMap;
use biodivine_lib_param_bn::bdd_params::{BddParams};
use biodivine_lib_std::param_graph::{Graph, EvolutionOperator, Params};
use biodivine_lib_std::{IdState};
use biodivine_aeon_server::scc::{StateSet};
use std::clone::Clone;
use std::io;
use std::io::Write;
use biodivine_aeon_server::scc::algo_reach::guarded_reach;
use std::borrow::Borrow;


pub fn find_strong_basin(graph: &AsyncGraph, attractor: IdState, params: BddParams) -> HashMap<IdState, BddParams> {
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
    loop {
        let mut to_remove = HashMap::new();
        let current_basin: Vec<IdState> = basin.keys().cloned().collect();

        for node in current_basin {
            // Find all nodes+params which have successors outside the basin
            // therefore they are not part of the strong basin
            let node_in_basin = basin.get(&node).unwrap(); // nodes are keys in basin
            let successors = fwd.step(node);
            for (suc, suc_params) in successors {
                // Under what parameters is successor a part of the basin
                let basin_params = basin.get(&suc).unwrap_or(&empty_params);

                // The parameters which lead outside of basin. These are the params which node
                // currently has assigned, for which there is an edge outside (and suc_params)
                // and these are NOT in the basin (and !basin_params).
                let difference = node_in_basin.intersect(&suc_params).minus(basin_params);
                if !difference.is_empty() {
                    // There are params which lead outside the basin
                    match to_remove.get_mut(&node) {
                        None => {
                            to_remove.insert(node, difference);
                        },
                        Some(value) => {
                            *value = value.union(&difference);
                        }
                    }
                }
            }
        }

        let mut total = 0.0;
        for (s, p) in basin.iter() {
            total += p.cardinality();
        }
        total = 0.0;
        for (s, p) in to_remove.iter() {
            total += p.cardinality();
        }
        io::stdout().flush();
        if to_remove.is_empty() {
            break;
        }

        for (node, params) in to_remove {
            *basin.get_mut(&node).unwrap() = basin[&node].minus(&params);

            if basin.get(&node).unwrap().is_empty() {
                // Clean up node which does not belong to basin under some parameter
                basin.remove(&node);
            }
        }
    }

    return basin;
}


fn find_weak_basin(graph: &AsyncGraph, attractor: IdState, params: &BddParams) -> HashMap<IdState, BddParams> {
    let bwd = graph.bwd();
    let mut basin = HashMap::new();
    basin.insert(attractor, params.clone());

    let empty_params = &graph.empty_params();
    let mut found_something = true;

    while found_something {
        found_something = false;
        let current_basin: Vec<IdState> = basin.keys().cloned().collect();

        for node in current_basin {
            for (parent, parent_edge_params) in bwd.step(node) {
                let node_in_basin = basin.get(&node).unwrap();  // unwrap ok because nodes are keys from basin
                // Reference to parameters for which parent is currently in basin.
                let parent_in_basin = basin.get(&parent).unwrap_or(empty_params);
                // Parameters which can be potentially still added to the basin.
                let can_become_basin = node_in_basin.intersect(&parent_edge_params).minus(parent_in_basin);

                if can_become_basin.is_empty() {    // skip useless states
                    continue;
                }

                // If something still remains in add_to_basin, actually add it to basin!
                found_something = true;
                match basin.get_mut(&parent) {
                    None => {
                        basin.insert(parent, can_become_basin);
                    },
                    Some(value) => {
                        *value = value.union(&can_become_basin);
                    }
                }

                println!("Weak basin size: {}", basin.len())
            }
        }
    }
    return basin;
}