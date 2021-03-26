use biodivine_lib_param_bn::async_graph::AsyncGraph;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_pbn_control::strong_basin::_algo_utils::find_attractors;
use std::convert::TryFrom;
use std::fs;

fn main() {
    for m in [
        "cell_fate",
        "myeloid",
        "drosophila",
        "proinflamatory_tumor",
        "apoptosis",
    ]
    .iter()
    {
        let cell_fate_witness: &str =
            &fs::read_to_string(format!("models/{}_witness.aeon", m)).unwrap();
        let cf_model = BooleanNetwork::try_from(cell_fate_witness).unwrap();
        let cf_graph = &AsyncGraph::new(cf_model).unwrap();

        println!("{}", m);
        println!("{}", cf_graph.num_states());
        println!("{}", cf_graph.unit_params().cardinality());

        let cf_attractors = &find_attractors(cf_graph);

        println!("{:?}", cf_attractors);
    }
}
