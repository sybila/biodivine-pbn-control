use biodivine_lib_param_bn::fixed_points::FixedPoints;
use biodivine_lib_param_bn::symbolic_async_graph::reachability::Reachability;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_pbn_control::aeon::phentoype::build_phenotype;
use std::collections::HashMap;
use std::time::Instant;

fn main() {
    let phenotypes = [
        "epithelial",
        "hybrid_1",
        "hybrid_2",
        "hybrid_3",
        "mensenchymal_1",
        "mensenchymal_2",
        "mensenchymal_3",
        "undefined",
    ];
    let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
    let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();

    for phenotype in phenotypes {
        println!("{}", phenotype);
        let model_file = config["emt"]["file"].as_str().unwrap();
        // let controllable_vars = get_controllable_vars("emt", model_file);

        let model_string =
            std::fs::read_to_string(format!("./models_phenotype/{}", model_file)).unwrap();
        let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
        let graph = SymbolicAsyncGraph::new(&bn).unwrap();

        let phenotype_map = config["emt"]["targets"][phenotype].as_object().unwrap();
        let mut phenotype_vals = HashMap::new();
        for (k, v) in phenotype_map {
            phenotype_vals.insert(k.as_str(), v.as_bool().unwrap());
        }
        let p_space = build_phenotype(&graph, phenotype_vals);
        let space = graph.unit_colored_vertices().intersect_vertices(&p_space);

        let now = Instant::now();

        let attractors = FixedPoints::symbolic(&graph, &space);
        let _bwd = Reachability::reach_bwd(&graph, &attractors);

        let duration = now.elapsed();
        println!(
            "{:?}: Time elapsed for attractor search: {:?}",
            phenotype, duration
        );

        // let now = Instant::now();
        //
        // let _attractors2 = biodivine_pbn_control::aeon::attractors::compute(&graph.as_perturbed());
        //
        // let duration = now.elapsed();
        // println!("{:?}: Time elapsed for attractor search: {:?}", m, duration);
    }
}

/*
fn get_controllable_vars(model_name: &str, model_file: &str) -> Vec<VariableId> {
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", model_file)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();

    let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
    let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();

    let mut controllable_vars = Vec::new();
    let uncontrollable = config[model_name]["uncontrollable"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
    let inputs = config[model_name]["inputs"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
    for v in bn.variables() {
        if !uncontrollable.contains(&bn.get_variable_name(v).as_str()) && !inputs.contains(&bn.get_variable_name(v).as_str()){
            controllable_vars.push(v);
        }
    }

    return controllable_vars
}*/
