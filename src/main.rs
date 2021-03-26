use std::fs;
use biodivine_pbn_control::control::_algo_temporary_control::find_attractors;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_lib_param_bn::BooleanNetwork;
use std::convert::TryFrom;
use biodivine_lib_param_bn::biodivine_std::traits::Set;

fn main() {
    for m in ["myeloid", "cell_fate", "proinflamatory_tumor", "apoptosis", "cholesterol"].iter() {
        let witness_str: &str = &fs::read_to_string(format!("models/{}_witness.aeon", m)).unwrap();
        let model = BooleanNetwork::try_from(witness_str).unwrap();
        let graph = &SymbolicAsyncGraph::new(model.clone()).unwrap();


        println!("Model: {}", m);
        println!("# variables: {}", model.variables().len());
        //println!("{}", graph.unit_colors().approx_cardinality());

        let attractors = &find_attractors(graph);

        let all = attractors.iter().fold(graph.mk_empty_vertices(), |a, b| a.union(b));
        println!("# attractors: {}", attractors.len());
        println!("# attractor vertices: {}", all.vertices().approx_cardinality());

        for vertex in all.vertices().materialize().iter() {
            println!("State: {:?}", vertex);
            let singleton_set = graph.vertex(&vertex);
            println!("Set vertices: {}", singleton_set.vertices().materialize().iter().count());
        }
        //println!("{:?}", attractors);
    }

}
