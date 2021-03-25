use std::fs;
use biodivine_lib_param_bn::async_graph::AsyncGraph;

fn main() {
    for m in ["cell_fate", "myeloid", "drosophila", "proinflamatory_tumor", "apoptosis",].iter() {
        let cell_fate_witness: &str = &fs::read_to_string(format!("models/{}_witness.aeon", m)).unwrap();
        // let cf_model = BooleanNetwork::try_from(cell_fate_witness).unwrap();
        // let cf_graph = &AsyncGraph::new(cf_model).unwrap();
        //
        // println!("{}", m);
        // println!("{}", cf_graph.num_states());
        // println!("{}", cf_graph.unit_params().cardinality());

        // let cf_attractors = &find_attractors(cf_graph);

        // println!("{:?}", cf_attractors);
    }

}
