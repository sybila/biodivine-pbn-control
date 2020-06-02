use biodivine_lib_param_bn::async_graph::{AsyncGraph};
use std::collections::{HashMap, HashSet, VecDeque};
use biodivine_lib_param_bn::bdd_params::{BddParams};
use biodivine_lib_std::param_graph::{Graph, EvolutionOperator, Params};
use biodivine_lib_std::{IdState};
use biodivine_aeon_server::scc::{StateSet};
use std::clone::Clone;
use std::io;
use std::io::Write;
use biodivine_aeon_server::scc::algo_reach::guarded_reach;
use std::borrow::Borrow;

pub fn paremeterless_find_strong_basin(graph: &AsyncGraph, attractor: IdState) -> HashSet<IdState> {
    let fwd = graph.fwd();

    let mut basin= find_weak_basin(graph, attractor);
    println!("Weak basin has {} states.", basin.len());

    loop {
        let mut removed_something = false;

        for node in basin.clone() {
            // Find all nodes which have successors outside the basin
            // therefore they are not part of the strong basin
            for (suc, p) in fwd.step(node) {
                if !basin.contains(&suc) && !p.is_empty(){
                    // Node has successor outside of the basin, the node is not a part of the SB
                    removed_something = true;
                    basin.remove(&node);
                    break;
                }
            }
        }

        if !removed_something {
            break;
        }
    }

    return basin;
}


fn find_weak_basin(graph: &AsyncGraph, attractor: IdState) -> HashSet<IdState> {
    let bwd = graph.bwd();
    let mut basin = HashSet::new();
    let mut queue = VecDeque::new();

    queue.push_back(attractor);

    while !queue.is_empty() {
        let current = queue.pop_front().unwrap();
        if basin.contains(&current) {
            continue;
        }

        basin.insert(current);
        for (next, p) in bwd.step(current) {
            if !p.is_empty() {
                queue.push_back(next);
            }
        }
    }
    return basin;
}