use biodivine_lib_param_bn::async_graph::{AsyncGraph};
use std::collections::HashMap;
use biodivine_lib_param_bn::bdd_params::BddParams;
use biodivine_lib_std::param_graph::{Graph, EvolutionOperator, Params};
use biodivine_lib_std::{IdState};

fn find_strong_basin(graph: &AsyncGraph, attractor: IdState, params: BddParams) -> HashMap<IdState, BddParams> {
    let bwd = graph.bwd();
    let fwd = graph.fwd();
    let mut basin = HashMap::new();
    basin.insert(attractor, params);

    let mut foundSomething = true;

    while foundSomething {
        foundSomething = false;
        let mut current_basin = Vec::new();
        for (node, _) in &basin {
            current_basin.push(node.clone());
        }

        for node in current_basin {
            let nodeParents = bwd.step(node);
            for (parent, param) in nodeParents {
                if basin.contains_key(&parent) {
                    if !basin[&parent].intersect(&param).is_empty() {
                        continue
                    }
                }

                let mut isInBasin = true;
                for (parentSuccessor, parentParam) in fwd.step(parent) {
                    if parentParam != param {
                        continue
                    }
                    match basin.get(&parentSuccessor) {
                        Some(basin_params) => {
                            if !basin_params.intersect(&parentParam).is_empty() {
                                isInBasin = false;
                                break;
                            }
                        },
                        None => {
                            isInBasin = false;
                            break;
                        }
                    }
                }

                if isInBasin {
                    foundSomething = true;
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
        let sbmlStr = fs::read_to_string("src/g2b.sbml").expect("Something went wrong");
        let (model, _) = BooleanNetwork::from_sbml(&sbmlStr).expect("Something went wrong");
        let graph = &AsyncGraph::new(model).expect("tadida");
        let state = IdState::from(0b10101 as usize);

        let params = graph.unit_params().clone();
        let basin = find_strong_basin(graph, state, params);
        for (id, param) in basin {
            println!("{} {:?}", id, param);
        }
    }
}
