use biodivine_lib_param_bn::BooleanNetwork;
use std::convert::TryFrom;
use std::fs;
use biodivine_lib_param_bn::async_graph::AsyncGraph;
use biodivine_lib_std::{IdState};
use biodivine_pbn_control::strong_basin::_algo_utils::{get_all_params_with_attractor, find_attractors};
use biodivine_pbn_control::strong_basin::_algo_sb_fixed_point::find_strong_basin;
use std::time::{SystemTime, Instant};
use std::ops::Add;

fn main() {
    let aeon_str: &str = &fs::read_to_string("models/myeloid_witness.aeon").unwrap();
    let model = BooleanNetwork::try_from(aeon_str).unwrap();

    let graph = &AsyncGraph::new(model).unwrap();
    let state = IdState::from(1284);

    // let attractors = find_attractors(graph);
    // println!("{:?}", attractors);
    let relevant_params = get_all_params_with_attractor(graph, state);
    println!("Attractor parameter set cardinality: {}", relevant_params.cardinality());

    let begin = Instant::now();

    let basin = find_strong_basin(graph, state, &relevant_params);
    println!("Strong basin has {} states.", basin.len());
    for (id, param) in basin {
        println!("{} {:?}", id, param.cardinality());
    }
    println!("Time elapsed (s): {:?}", begin.elapsed().as_secs());
}