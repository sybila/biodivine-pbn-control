use std::collections::HashMap;
use biodivine_lib_param_bn::{BooleanNetwork, VariableId};
use biodivine_lib_param_bn::symbolic_async_graph::GraphVertices;
use crate::aeon::phentoype::build_phenotype;
use crate::perturbation::PerturbationGraph;


pub(crate) fn get_controllable_vars(model_file: &str, model_name: &str, extra_forbidden: Vec<&str>) -> Vec<VariableId> {
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", model_file)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();

    let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
    let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();

    let mut controllable_vars = Vec::new();
    let uncontrollable = config[model_name]["uncontrollable"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
    // let inputs = config[model_name]["inputs"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
    for v in bn.variables() {
        if !uncontrollable.contains(&bn.get_variable_name(v).as_str()) && !extra_forbidden.contains(&bn.get_variable_name(v).as_str()) {
            controllable_vars.push(v);
        }
    }

    return controllable_vars
}

pub(crate) fn get_trivial_phenotype(model_name: &str, phenotype_name: &str, stg: &PerturbationGraph) -> GraphVertices {
    let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
    let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();
    let phenotype_map  = config[model_name]["targets"][phenotype_name].as_object().unwrap();
    let mut phenotype_vals = HashMap::new();
    for (k,v) in phenotype_map {
        phenotype_vals.insert(k.as_str(), v.as_bool().unwrap());
    }
    build_phenotype(stg.as_perturbed(),phenotype_vals)
}
