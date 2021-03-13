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
    for m in ["cell_fate", "myeloid", "drosophila", "proinflamatory_tumor", "apoptosis",].iter() {
        let cell_fate_witness: &str = &fs::read_to_string(format!("models/{}_witness.aeon", m)).unwrap();
        let cf_model = BooleanNetwork::try_from(cell_fate_witness).unwrap();
        let cf_graph = &AsyncGraph::new(cf_model).unwrap();

        println!("{}", m);
        println!("{}", cf_graph.num_states());
        println!("{}", cf_graph.unit_params().cardinality());

        let cf_attractors = &find_attractors(cf_graph);

        println!("{:?}", cf_attractors);
    }

}
