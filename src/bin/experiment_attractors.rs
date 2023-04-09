use biodivine_pbn_control::experiment_utils::{parse_experiment, run_control_experiment};
use biodivine_pbn_control::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::fixed_points::FixedPoints;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use std::time::Instant;


fn main() {
    let models = [
        "[id-010]__[var-13]__[in-2]__[CARDIAC-DEVELOPMENT].aeon",
        "[id-089]__[var-13]__[in-4]__[MAPK-REDUCED-1].aeon",
        "[id-096]__[var-19]__[in-1]__[ERBB-REGULATED-G1-S-TRANSITION].aeon",
        "[id-065]__[var-30]__[in-2]__[TUMOUR-CELL-INVASION-AND-MIGRATION].aeon",
        "[id-150]__[var-31]__[in-2]__[CELL-FATE-DECISION-MULTISCALE].aeon",
        "[id-179]__[var-46]__[in-10]__[MICROENVIRONMENT-CONTROL].aeon",
        "[id-070]__[var-49]__[in-4]__[MAPK-CANCER-CELL-FATE].aeon",
        "[id-014]__[var-54]__[in-7]__[T-LGL-SURVIVAL-NETWORK-2008].aeon"
    ];

    for m in models {
        let model_string = std::fs::read_to_string(format!("./models_phenotype/{}",m)).unwrap();
        let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
        let graph_tmp = SymbolicAsyncGraph::new(bn.clone()).unwrap();
        let witness = graph_tmp.pick_witness(graph_tmp.unit_colors());
        let graph = SymbolicAsyncGraph::new(witness.clone()).unwrap();

        let now = Instant::now();

        let _attractors = FixedPoints::symbolic(&graph, graph.unit_colored_vertices());

        let duration = now.elapsed();
        println!("{:?}: Time elapsed for attractor search: {:?}", m, duration);


        let now = Instant::now();

        let _attractors2 = biodivine_pbn_control::aeon::attractors::compute(&graph);

        let duration = now.elapsed();
        println!("{:?}: Time elapsed for attractor search: {:?}", m, duration);
    }
}
