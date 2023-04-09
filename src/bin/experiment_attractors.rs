use biodivine_pbn_control::experiment_utils::{parse_experiment, run_control_experiment};
use biodivine_pbn_control::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::fixed_points::FixedPoints;
use biodivine_lib_param_bn::{BooleanNetwork, VariableId};
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use std::time::Instant;


fn main() {
    let models = [
        "cardiac",
        "reduced_mapk",
        "erbb",
        "tumour",
        "cell_fate",
        "full_mapk",
        "t_lgl"
    ];
    let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();

    for m in models {
        let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();
        let model_file = config[m]["file"].as_str().unwrap();
        let controllable_vars = get_controllable_vars(m, model_file);

        let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", model_file)).unwrap();
        let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
        let graph = PerturbationGraph::with_restricted_variables(&bn, &controllable_vars);

        let now = Instant::now();

        let _attractors = FixedPoints::symbolic(&graph.as_perturbed(), graph.unit_colored_vertices());

        let duration = now.elapsed();
        println!("{:?}: Time elapsed for attractor search: {:?}", m, duration);


        let now = Instant::now();

        let _attractors2 = biodivine_pbn_control::aeon::attractors::compute(&graph.as_perturbed());

        let duration = now.elapsed();
        println!("{:?}: Time elapsed for attractor search: {:?}", m, duration);
    }
}


fn get_controllable_vars(model_name: &str, model_file: &str) -> Vec<VariableId> {
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", model_file)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();

    let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
    let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();

    let mut controllable_vars = Vec::new();
    let uncontrollable = config[model_name]["uncontrollable"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
    let inputs = config[model_name]["inputs"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
    for v in bn.variables() {
        if !uncontrollable.contains(&bn.get_variable_name(v).as_str()) && !inputs.contains(&bn.get_variable_name(v).as_str()) && !extra_forbidden.contains(&bn.get_variable_name(v).as_str()) {
            controllable_vars.push(v);
        }
    }

    return controllable_vars
}
