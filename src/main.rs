use std::fs;
use biodivine_lib_param_bn::async_graph::AsyncGraph;
use biodivine_pbn_control::control::_algo_temporary_control::find_attractors;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_lib_param_bn::BooleanNetwork;
use std::convert::TryFrom;

fn main() {
    for m in ["myeloid", "cell_fate", "proinflamatory_tumor", "apoptosis", "cholesterol"].iter() {
        let witness_str: &str = &fs::read_to_string(format!("models/{}_witness.aeon", m)).unwrap();
        let model = BooleanNetwork::try_from(witness_str).unwrap();
        let graph = &SymbolicAsyncGraph::new(model.clone()).unwrap();


        println!("Model: {}", m);
        println!("# variables: {}", model.variables().len());
        //println!("{}", graph.unit_colors().approx_cardinality());

        let attractors = &find_attractors(graph);

        println!("# attractors: {}", attractors.len());
        println!("{:?}", attractors);
    }

}
