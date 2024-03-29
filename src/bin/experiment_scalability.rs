use std::collections::HashMap;

use biodivine_lib_param_bn::BooleanNetwork;
use std::time::Instant;

use biodivine_pbn_control::aeon::phentoype::build_phenotype;
use biodivine_pbn_control::perturbation::PerturbationGraph;
use chrono::Local;

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let model = args[1].as_str();
    let phenotype = args[2].as_str();
    let max_control_size: usize = args[3].parse::<usize>().unwrap();
    let max_control_vars: usize = args[4].parse::<usize>().unwrap();
    let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
    let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();
    let model_name = config[model]["file"].as_str().unwrap();
    let model_string =
        std::fs::read_to_string(format!("./models_phenotype/{}", model_name)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();

    let mut controllable_vars = Vec::new();
    let uncontrollable = config[model]["uncontrollable"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap())
        .collect::<Vec<&str>>();
    let inputs = config[model]["inputs"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap())
        .collect::<Vec<&str>>();
    for v in bn.variables() {
        if !uncontrollable.contains(&bn.get_variable_name(v).as_str())
            && !inputs.contains(&bn.get_variable_name(v).as_str())
        {
            controllable_vars.push(v);
        }
    }

    let mut p_vars = Vec::new();
    let mut i = 0;
    for v in controllable_vars.clone() {
        p_vars.push(v);
        i += 1;
        if i >= max_control_vars {
            break;
        }
    }
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &p_vars);

    let phenotype_map = config[model]["targets"][phenotype].as_object().unwrap();
    let mut phenotype_vals = HashMap::new();
    for (k, v) in phenotype_map {
        phenotype_vals.insert(k.as_str(), v.as_bool().unwrap());
    }
    let phenotype = build_phenotype(perturbation_graph.as_perturbed(), phenotype_vals);

    // Print model statistics:
    let model_variables = bn.num_vars();
    assert!(i32::try_from(model_variables).is_ok());

    // All colors considered by the perturbation graph
    let all_colors = perturbation_graph.unit_colors().approx_cardinality();
    // The (combinatorial) portion of colours that appear due to perturbation parameters.
    let perturbation_colors = 2.0f64.powi(p_vars.len() as i32);
    // The (combinatorial) portion of colours that are carried over from the original model.
    let model_colors = all_colors / perturbation_colors;

    println!("Variables: {}", model_variables);
    println!("Inputs: {}", inputs.len());
    println!("Controllable variables: {}", p_vars.len());
    println!("Uncertainty colors: {}", model_colors);
    println!("Perturbation colors: {}", perturbation_colors);
    println!("All colors: {}", all_colors);
    println!(
        "Unknown update functions: {}",
        bn.variables()
            .filter(|it| { bn.get_update_function(*it).is_none() })
            .count()
    );
    println!(
        "Lowest cardinality: {}",
        bn.variables()
            .filter(|it| { bn.get_update_function(*it).is_none() })
            .map(|it| { bn.as_graph().regulators(it).len() })
            .min()
            .unwrap()
    );
    println!(
        "Highest cardinality: {}",
        bn.variables()
            .filter(|it| { bn.get_update_function(*it).is_none() })
            .map(|it| { bn.as_graph().regulators(it).len() })
            .max()
            .unwrap()
    );

    let result = PerturbationGraph::ceiled_phenotype_permanent_control(
        &perturbation_graph,
        phenotype,
        max_control_size,
        p_vars.clone(),
        "heuristic",
    );

    let zero_perturbation_working_colors = result.perturbation_working_colors(&HashMap::from([]));
    println!(
        "No perturbation working for {:?}",
        zero_perturbation_working_colors.approx_cardinality()
    );

    let now = Instant::now();
    println!("Starting control enumeration at: {}", Local::now());

    result.ceiled_size_perturbation_working_colors(
        max_control_size,
        model_colors,
        &p_vars,
        false,
        false,
    );

    let duration = now.elapsed();
    println!("Control enumeration finished at {:?} ", Local::now());
    println!("Time elapsed for control enumeration: {:?}", duration);
}

/*
fn powerset(s: &[VariableId]) -> Vec<Vec<VariableId>> {
    let mut subsets: Vec<Vec<VariableId>> = vec![];
    let empty: Vec<VariableId> = vec![];
    subsets.push(empty);

    let mut updated: Vec<Vec<VariableId>> = vec![];

    for ele in s {
        for mut sub in subsets.clone() {
            sub.push(*ele);
            updated.push(sub);
        }
        subsets.append(&mut updated);
    }
    subsets
}
*/
