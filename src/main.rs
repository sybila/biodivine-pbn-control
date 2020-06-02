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
use biodivine_pbn_control::control::_algo_control_to_basin::{find_best_control_to_basin, find_smallest_control_to_basin, find_robust_control_to_basin, control_dist};

fn main() {
    analyse_model_parameterless("models/cell_fate_7stable_attractors.aeon");
    analyse_model_parameterless("models/myeloid_witness.aeon");

    let aeon_str: &str = &fs::read_to_string("models/cell_fate_7stable_attractors.aeon").unwrap();
    let model = BooleanNetwork::try_from(aeon_str).unwrap();
    let graph = &AsyncGraph::new(model).unwrap();
    let attractors = &find_attractors(graph);
    analyse_model_with_pars("models/cell_fate_7stable_attractors.aeon", attractors);
    analyse_model_with_pars("models/cell_fate_7stable_attractors_2params.aeon", attractors);
    analyse_model_with_pars("models/cell_fate_7stable_attractors_4params.aeon", attractors);

    let aeon_strr: &str = &fs::read_to_string("models/myeloid_witness.aeon").unwrap();
    let modell = BooleanNetwork::try_from(aeon_strr).unwrap();
    let graphh = &AsyncGraph::new(modell).unwrap();
    let attractorss = &find_attractors(graphh);
    analyse_model_with_pars("models/myeloid_witness.aeon", attractorss);
    analyse_model_with_pars("models/myeloid_4params.aeon", attractorss);
    analyse_model_with_pars("models/myeloid_8params.aeon", attractorss);
    analyse_model_with_pars("models/myeloid_11params.aeon", attractorss);
}


fn analyse_model_with_pars(model_file_name: &str, attractors: &Vec<IdState>) {
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

        for b in attractors.clone() {
            let begin = Instant::now();
            let control = find_smallest_control_to_basin(&b, &basin);

            println!("Smallest control from {:?} to {:?}: {:?}, distance: {:?}, cardinality {:?}/{:?}", b, a, control[0].0, control_dist(&b, &a), control[0].1.cardinality(), relevant_params_cardinality);
            println!("Time elapsed to compute the control from basin (ms): {:?}", begin.elapsed().as_millis());

            let begin2 = Instant::now();
            let control = find_robust_control_to_basin(&basin);

            println!("Robust control from {:?} to {:?}: {:?}, distance: {:?}, cardinality {:?}/{:?}", b, a, control[0].0, control_dist(&b, &a), control[0].1.cardinality(), relevant_params_cardinality);
            println!("Time elapsed to compute the control from basin (ms): {:?}", begin2.elapsed().as_millis());
        }

        println!();
    }
    println!();
}

fn analyse_model_parameterless(model_file_name: &str) {
    println!("Model : {}", model_file_name);

    let aeon_str: &str = &fs::read_to_string(model_file_name).unwrap();
    let model = BooleanNetwork::try_from(aeon_str).unwrap();

    let graph = &AsyncGraph::new(model).unwrap();

    let attractors = find_attractors(graph);
    for a in attractors.clone() {
        println!("Attractor : {}", a);

        let begin = Instant::now();

        let basin = paremeterless_find_strong_basin(graph, a);
        println!("Strong basin has {} states.", basin.len());
        println!("Strong basin computation time (ms): {:?}", begin.elapsed().as_millis());
        println!();

        for b in attractors.clone() {
            let begin = Instant::now();
            let control = find_best_control_to_basin(&b, &basin);

            println!("Control from {:?} to {:?}: {:?}", b, a, control);
            println!("Time elapsed to compute the control from basin (ms): {:?}", begin.elapsed().as_millis())
        }
    }
    println!();
}



