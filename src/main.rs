use biodivine_lib_param_bn::BooleanNetwork;
use std::convert::TryFrom;
use std::fs;
use biodivine_lib_param_bn::async_graph::AsyncGraph;
use biodivine_lib_std::{IdState};
use biodivine_pbn_control::strong_basin::_algo_utils::{get_all_params_with_attractor, find_attractors};
use biodivine_pbn_control::parameterless_experimental::_algo_sb_fixed_point::paremeterless_find_strong_basin;
use std::time::{SystemTime, Instant};
use std::ops::Add;
use biodivine_lib_param_bn::bdd_params::BddParameterEncoder;
use std::fs::File;
use std::io::{LineWriter, Write};
use biodivine_pbn_control::strong_basin::_algo_sb_fixed_point::find_strong_basin;
use biodivine_pbn_control::control::_algo_control_to_basin::{find_smallest_control_to_basin, find_robust_control_to_basin, control_dist};

fn main() {
    let cell_fate_witness: &str = &fs::read_to_string("models/cell_fate_7stable_attractors.aeon").unwrap();
    let cf_model = BooleanNetwork::try_from(cell_fate_witness).unwrap();
    let cf_graph = &AsyncGraph::new(cf_model).unwrap();
    let cf_attractors = &find_attractors(cf_graph);
    analyse_model("models/cell_fate_7stable_attractors.aeon", cf_attractors);
    analyse_model("models/cell_fate_7stable_attractors_2params.aeon", cf_attractors);
    analyse_model("models/cell_fate_7stable_attractors_4params.aeon", cf_attractors);

    let myeloid_witness: &str = &fs::read_to_string("models/myeloid_witness.aeon").unwrap();
    let m_model = BooleanNetwork::try_from(myeloid_witness).unwrap();
    let m_graph = &AsyncGraph::new(m_model).unwrap();
    let m_attractors = &find_attractors(m_graph);
    analyse_model("models/myeloid_witness.aeon", m_attractors);
    analyse_model("models/myeloid_4params.aeon", m_attractors);
    analyse_model("models/myeloid_8params.aeon", m_attractors);
    analyse_model("models/myeloid_11params.aeon", m_attractors);

    source_target_controls("models/myeloid_11params.aeon", IdState::from(553), IdState::from(1285))
}


fn analyse_model(model_file_name: &str, attractors: &Vec<IdState>) {
    println!("Model : {}", model_file_name);

    let aeon_str: &str = &fs::read_to_string(model_file_name).unwrap();
    let model = BooleanNetwork::try_from(aeon_str).unwrap();

    let graph = &AsyncGraph::new(model).unwrap();

    for a in attractors.clone() {
        println!("Attractor : {}", a);
        let relevant_params = get_all_params_with_attractor(graph, a);
        let relevant_params_cardinality =  relevant_params.cardinality();
        println!("Attractor parameter set cardinality: {}", relevant_params_cardinality);

        let begin = Instant::now();

        let basin = find_strong_basin(graph, a, relevant_params);
        println!("Strong basin has {} states.", basin.len());
        println!("Strong basin computation time (ms): {:?}", begin.elapsed().as_millis());

        println!();
    }
    println!();
}

fn source_target_controls(model_file_name: &str, source: IdState, target: IdState) {
    println!("Analysing source-target from {:?} to {:?} of model {}", source, target, model_file_name);

    let aeon_str: &str = &fs::read_to_string(model_file_name).unwrap();
    let model = BooleanNetwork::try_from(aeon_str).unwrap();

    let graph = &AsyncGraph::new(model).unwrap();

    let relevant_params = get_all_params_with_attractor(graph, target);
    let relevant_params_cardinality =  relevant_params.cardinality();
    println!("Attractor parameter set cardinality: {}", relevant_params_cardinality);

    let begin = Instant::now();

    let basin = find_strong_basin(graph, target, relevant_params);
    println!("Strong basin has {} states.", basin.len());
    println!("Strong basin computation time (ms): {:?}", begin.elapsed().as_millis());

    let controls = find_smallest_control_to_basin(&source, &basin);
    for (s, p) in controls {
        println!("Control into {:?}, size {:?}, robustness {:?}.", s, control_dist(&s, &target), p.cardinality() / relevant_params_cardinality);
    }
}
