use biodivine_lib_param_bn::BooleanNetwork;
use std::convert::TryFrom;
use std::fs;
use biodivine_lib_param_bn::async_graph::AsyncGraph;
use biodivine_lib_std::{IdState};
use biodivine_pbn_control::strong_basin::_algo_utils::{get_all_params_with_attractor, find_attractors};
use std::time::{Instant};
use biodivine_lib_param_bn::bdd_params::BddParameterEncoder;
use std::fs::File;
use std::io::{LineWriter, Write};
use biodivine_pbn_control::strong_basin::_algo_sb_parallel_fixed_point::find_strong_basin;
use biodivine_pbn_control::control::_algo_control_to_basin::{find_smallest_control_to_basin, find_robust_control_to_basin, control_dist};
use biodivine_aeon_server::scc::StateSet;
use biodivine_pbn_control::async_graph_with_control::AsyncGraphWithControl;
use biodivine_pbn_control::controlled_async_graph::ControlledAsyncGraph;

fn main() {
    // let cell_fate_witness: &str = &fs::read_to_string("models/cell_fate_7stable_attractors.aeon").unwrap();
    // let cf_model = BooleanNetwork::try_from(cell_fate_witness).unwrap();
    // let cf_graph = &AsyncGraph::new(cf_model).unwrap();
    // let cf_attractors = &find_attractors(cf_graph);
    // analyse_model("models/cell_fate_7stable_attractors.aeon", cf_attractors);
    // analyse_model("models/cell_fate_7stable_attractors_2params.aeon", cf_attractors);
    // analyse_model("models/cell_fate_7stable_attractors_4params.aeon", cf_attractors);

    let myeloid_witness: &str = &fs::read_to_string("models/myeloid_witness.aeon").unwrap();
    let m_model = BooleanNetwork::try_from(myeloid_witness).unwrap();
    let m_graph = &AsyncGraph::new(m_model).unwrap();
    let m_attractors = &find_attractors(m_graph);
    // analyse_model("models/myeloid_witness.aeon", m_attractors);
    // analyse_model("models/myeloid_4params.aeon", m_attractors);
    analyse_model("models/myeloid_8params.aeon", m_attractors);
    // analyse_model("models/myeloid_11params.aeon", m_attractors);

    //source_target_controls("models/myeloid_11params.aeon", IdState::from(553), IdState::from(1285))
}


fn analyse_model(model_file_name: &str, attractors: &Vec<IdState>) {
    println!("Model : {}", model_file_name);

    let aeon_str: &str = &fs::read_to_string(model_file_name).unwrap();
    let model = BooleanNetwork::try_from(aeon_str).unwrap();

    let mut graph = ControlledAsyncGraph::new(model);

    for source in attractors.clone() {
        println!("Source attractor : {}", source);
        for target in attractors.clone() {
        println!("Target Attractor : {}", target);
            let target_seed = &StateSet::new_with_fun(graph.num_states(), |s| if s.eq(&target) { Some(graph.unit_params().clone()) } else { None });

            let mut begin = Instant::now();
            let temporary = graph.find_temporary_control(source, target_seed);
            println!("Temporary control can be done in {} ways.", temporary.len());
            println!("Temporary control computation time (ms): {:?}", begin.elapsed().as_millis());

            begin = Instant::now();
            let permanent = graph.find_permanent_control(source, target_seed);
            println!("Permanent control can be done in {} ways.", permanent.len());
            println!("Permanent control computation time (ms): {:?}", begin.elapsed().as_millis());

            println!();
        }
        println!();
    }

}

// fn source_target_controls(model_file_name: &str, source: IdState, target: IdState) {
//     println!("Analysing source-target from {:?} to {:?} of model {}", source, target, model_file_name);
//
//     let aeon_str: &str = &fs::read_to_string(model_file_name).unwrap();
//     let model = BooleanNetwork::try_from(aeon_str).unwrap();
//
//     let graph = &AsyncGraphWithControl::new(model);
//
//     let relevant_params = get_all_params_with_attractor(graph, target);
//     let relevant_params_cardinality =  relevant_params.cardinality();
//     println!("Attractor parameter set cardinality: {}", relevant_params_cardinality);
//
//     let begin = Instant::now();
//
//     let seed = &StateSet::new_with_fun(graph.num_states(), |s| if s.eq(&target) { Some(relevant_params.clone()) } else { None });
//     let basin = find_strong_basin(graph, seed, graph.unit_params());
//     println!("Strong basin has {} states.", basin.len());
//     println!("Strong basin computation time (ms): {:?}", begin.elapsed().as_millis());
//
//     let controls = find_smallest_control_to_basin(&source, &basin);
//     for (s, p) in controls {
//         println!("Control into {:?}, size {:?}, robustness {:?}.", s, control_dist(&s, &target), p.cardinality() / relevant_params_cardinality);
//     }
// }
