use biodivine_lib_param_bn::async_graph::{AsyncGraph};
use std::collections::HashMap;
use biodivine_lib_param_bn::bdd_params::BddParams;
use biodivine_lib_std::param_graph::{Graph, EvolutionOperator, Params};
use biodivine_lib_std::{IdState};

fn find_strong_basin(graph: &AsyncGraph, attractor: IdState, params: &BddParams) -> HashMap<IdState, BddParams> {
    let bwd = graph.bwd();
    let fwd = graph.fwd();
    let mut basin = HashMap::new();
    basin.insert(attractor, params.clone());

    let mut found_something = true;

    while found_something {
        found_something = false;
        let mut current_basin = Vec::new();
        for (node, _) in &basin {
            current_basin.push(node.clone());
        }

        for node in current_basin {
            let node_parents = bwd.step(node);
            for (parent, param) in node_parents {
                if basin.contains_key(&parent) {
                    if !basin[&parent].intersect(&param).is_empty() {
                        continue
                    }
                }

                let mut is_in_basin = true;
                for (parent_successor, parent_param) in fwd.step(parent) {
                    if parent_param != param {
                        continue
                    }
                    match basin.get(&parent_successor) {
                        Some(basin_params) => {
                            if !basin_params.intersect(&parent_param).is_empty() {
                                is_in_basin = false;
                                break;
                            }
                        },
                        None => {
                            is_in_basin = false;
                            break;
                        }
                    }
                }

                if is_in_basin {
                    found_something = true;
                    match basin.get_mut(&parent) {
                        Some(basin_params) => {
                            let params = basin_params.union(&param);
                            basin.entry(parent).or_insert(params);
                            break;
                        },
                        None => {
                            basin.insert(parent, param);
                        }
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
    use std::fs;
    use biodivine_lib_param_bn::BooleanNetwork;
    use biodivine_lib_bdd::{Bdd, BddVariableSet};
    use biodivine_lib_param_bn::bdd_params::BddParameterEncoder;
    use std::env::var;
    use std::path::Path;
    use std::fmt::Debug;

    #[test]
    fn test_reach() {
        let sbml_str = fs::read_to_string("src/g2b.sbml").unwrap();
        let (model, _) = BooleanNetwork::from_sbml(&sbml_str).unwrap();
        let graph = &AsyncGraph::new(model).unwrap();
        let state = IdState::from(0b10101 as usize);

        let basin = find_strong_basin(graph, state, graph.unit_params());
        for (id, param) in basin {
            println!("{} {:?}", id, param);
        }
    }
}
