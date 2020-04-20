use biodivine_lib_param_bn::async_graph::{AsyncGraph};
use std::collections::HashMap;
use biodivine_lib_param_bn::bdd_params::{BddParams};
use biodivine_lib_std::param_graph::{Graph, EvolutionOperator, Params};
use biodivine_lib_std::{IdState};


fn find_strong_basin(graph: &AsyncGraph, attractor: IdState, params: &BddParams) -> HashMap<IdState, BddParams> {
    let fwd = graph.fwd();
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
    let mut basin = find_weak_basin(graph, attractor, params);

    loop {
        let mut to_remove = HashMap::new();
        let current_basin: Vec<IdState> = basin.keys().cloned().collect();

        for node in current_basin {
            // Find all nodes+params which have successors outside the basin
            // therefore they are not part of the strong basin
            let successors = fwd.step(node);
            for (suc, suc_params) in successors {
                // Under what parameters is successor a part of the basin
                let basin_params = basin.get(&suc).unwrap_or(&empty_params);

                let difference = suc_params.minus(basin_params);
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
                // Reference to parameters for which parent is currently in basin.
                let parent_in_basin = basin.get(&parent).unwrap_or(empty_params);
                // Parameters which can be potentially still added to the basin.
                let can_become_basin = parent_edge_params.minus(parent_in_basin);

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
            }
        }
    }
    return basin;
}


#[cfg(test)]
mod tests {
    use super::*;
    use biodivine_lib_param_bn::BooleanNetwork;
    use std::convert::TryFrom;

    const TEST_MODEL: &str = "\
        CtrA -> CtrA
        GcrA -> CtrA
        CcrM -| CtrA
        SciP -| CtrA

        CtrA -| GcrA
        DnaA -> GcrA

        CtrA -> SciP
        DnaA -| SciP

        CtrA -> CcrM
        CcrM -| CcrM
        SciP -| CcrM

        CcrM -> DnaA
        CtrA -> DnaA
        GcrA -| DnaA
        DnaA -| DnaA
    ";

    fn get_all_params_with_attractor(graph: &AsyncGraph, state: IdState) -> BddParams {
        let fwd = graph.fwd();
        let successors = fwd.step(state);

        let mut bad_params = graph.empty_params();
        for (succ, par) in successors {
            //if !succ.eq(&state) {
                bad_params = bad_params.union(&par);
            //}
        }
        bad_params = bad_params.intersect(graph.unit_params());

        return graph.unit_params().minus(&bad_params);
    }

    #[test]
    fn test_reach() {
        //let sbml_str = fs::read_to_string("src/g2b.sbml").unwrap();
        //let (model, _) = BooleanNetwork::from_sbml(&sbml_str).unwrap();

        let model = BooleanNetwork::try_from(TEST_MODEL).unwrap();

        //println!("{:?}", model.graph());

        let graph = &AsyncGraph::new(model).unwrap();
        // We want to start from state
        // CcrM | CtrA | DnaA | GcrA | SciP
        //  1   |  1   |  1   |  0   |  0
        // Note that bits are reversed because first variable corresponds to least significant bit.
        let state = IdState::from(0b00111 as usize);

        let relevant_params = get_all_params_with_attractor(graph, state);
        println!("Attractor parameter set cardinality: {}", relevant_params.cardinality());
        let basin = find_strong_basin(graph, state, &relevant_params);
        println!("Strong basin has {} states.", basin.len());
        for (id, param) in basin {
            println!("{} {:?}", id, param.cardinality());
        }
    }
}
