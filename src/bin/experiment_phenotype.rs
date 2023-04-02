use std::borrow::Borrow;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::time::Instant;
use std::vec;
use biodivine_lib_param_bn::{BooleanNetwork, VariableId};
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use chrono::Local;
use itertools::Itertools;
use serde_json::Value;
use biodivine_pbn_control::aeon::phentoype::build_phenotype;
use biodivine_pbn_control::experiment_utils::{parse_experiment, run_control_experiment};
use biodivine_pbn_control::perturbation::PerturbationGraph;

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let model = args[1].as_str();
    let phenotype = args[2].as_str();
    let max_control_size: i32 = args[3].parse::<i32>().unwrap();
    // let max_control_vars: usize = args[4].parse::<usize>().unwrap();
    let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
    let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();
    let model_name = config[model]["file"].as_str().unwrap();
    let model_string = std::fs::read_to_string( format!("./models_phenotype/{}", model_name)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();

    // let mut p_vars = Vec::new();
    // let mut i = 0;
    // for v in bn.variables() {
    //     if i == 0 {
    //         i += 1;
    //         continue;
    //     }
    //     if i > max_control_vars {
    //         break;
    //     }
    //     p_vars.push(v);
    //     i += 1;
    // }
    // let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &p_vars);

    let mut controllable_vars = Vec::new();
    let uncontrollable = config[model]["uncontrollable"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
    let inputs = config[model]["inputs"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
    for v in bn.variables() {
        if !uncontrollable.contains(&bn.get_variable_name(v).as_str()) && !inputs.contains(&bn.get_variable_name(v).as_str()) {
            controllable_vars.push(v);
        }
    }

    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &controllable_vars);


    let phenotype_map  = config[model]["targets"][phenotype].as_object().unwrap();
    let mut phenotype_vals = HashMap::new();
    for (k,v) in phenotype_map {
        phenotype_vals.insert(k.as_str(), v.as_bool().unwrap());
    }
    let phenotype = build_phenotype(perturbation_graph.as_perturbed(),
                                    phenotype_vals);

    // Print model statistics:
    let model_variables = bn.num_vars();
    assert!(i32::try_from(model_variables).is_ok());

    // All colors considered by the perturbation graph
    let all_colors = perturbation_graph.unit_colors().approx_cardinality();
    // The (combinatorial) portion of colours that appear due to perturbation parameters.
    let perturbation_colors = 2.0f64.powi(controllable_vars.len() as i32);
    // The (combinatorial) portion of colours that are carried over from the original model.
    let model_colors = all_colors / perturbation_colors;

    println!("Variables: {}", model_variables);
    println!("Inputs: {}", inputs.len());
    println!("Controllable variables: {}", controllable_vars.len());
    println!("Uncertainty colors: {}", model_colors);
    println!("Perturbation colors: {}", perturbation_colors);
    println!("All colors: {}", all_colors);
    println!(
        "Unknown update functions: {}",
        bn
            .variables()
            .filter(|it| { bn.get_update_function(*it).is_none() })
            .count()
    );
    println!(
        "Lowest cardinality: {}",
        bn
            .variables()
            .filter(|it| { bn.get_update_function(*it).is_none() })
            .map(|it| { bn.as_graph().regulators(it).len() })
            .min()
            .unwrap()
    );
    println!(
        "Highest cardinality: {}",
        bn
            .variables()
            .filter(|it| { bn.get_update_function(*it).is_none() })
            .map(|it| { bn.as_graph().regulators(it).len() })
            .max()
            .unwrap()
    );


    let result = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype, max_control_size, controllable_vars.clone());


    let now = Instant::now();
    println!("Starting control enumeration at at: {}", Local::now());

    for i in 1..(max_control_size+1) {
        let mut max_cardinality = 0.0;
        let mut max_cardinality_control = HashMap::new();
        let mut non_zero_working = 1;
        // let mut union_working_colors = perturbation_graph.as_original().mk_empty_colors();
        for controlled in controllable_vars.iter().combinations(i as usize) {
            let mut controlled_cpy = Vec::new();
            for c in controlled {
                controlled_cpy.push(c.clone());
            }
            for over_expressed in powerset(&(controlled_cpy.clone())) {
                let mut control = HashMap::new();
                for v in controlled_cpy.clone() {
                    if over_expressed.contains(&v) {
                        control.insert(bn.get_variable_name(v).clone(), true);
                    } else {
                        control.insert(bn.get_variable_name(v).clone(), false);
                    }
                }

                let working_colors = result.perturbation_working_colors(&control);
                if working_colors.approx_cardinality() > 0.0 {
                    non_zero_working += 1;
                    if working_colors.approx_cardinality() > max_cardinality {
                        max_cardinality = working_colors.approx_cardinality();
                        max_cardinality_control = control;
                    }
                    // union_working_colors = union_working_colors.union(&working_colors);
                }
            }
        }

        println!("Max cardinality control for size {:?}: {:?} working for {:?} colors", i, max_cardinality_control, max_cardinality);
        println!("Controls of size {:?} working for some colors: {:?}", i, non_zero_working);
        // SOME VAR PROJECT IS MISSING :(
        // println!("Union working colors cardinality of size {:?}: {:?}", i, union_working_colors.approx_cardinality());
    }

    let duration = now.elapsed();
    println!("Control enumeration finished at {:?} ", Local::now());
    println!("Time elapsed for control enumeration: {:?}", duration);
}

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
