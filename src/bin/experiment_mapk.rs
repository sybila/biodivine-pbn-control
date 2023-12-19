use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_pbn_control::aeon::phentoype::build_phenotype;
use biodivine_pbn_control::perturbation::PerturbationGraph;
use biodivine_pbn_control::phenotype_control::_simplified_algorithm::bounded_phenotype_control;
use std::collections::HashMap;

fn main() {
    for model_name in [
        "[id-070]__[var-49]__[in-4]__[MAPK-CANCER-CELL-FATE]_uncertain_DUSP1_reg.aeon",
        "[id-070]__[var-49]__[in-4]__[MAPK-CANCER-CELL-FATE]_uncertain_FRS2.aeon",
        "[id-070]__[var-49]__[in-4]__[MAPK-CANCER-CELL-FATE]_uncertain_MEK1_2.aeon",
    ] {
        for phenotype in ["apoptosis", "proliferation", "no_decision", "growth_arrest"] {
            println!(">>>>>>>>>>>>>> {:?} {:?}", model_name, phenotype);
            let max_control_size = 3;
            let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
            let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();
            let model_string =
                std::fs::read_to_string(format!("./models_phenotype/{}", model_name)).unwrap();
            let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();

            let mut controllable_vars = Vec::new();
            let uncontrollable = config["full_mapk"]["uncontrollable"]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_str().unwrap())
                .collect::<Vec<&str>>();
            let inputs = config["full_mapk"]["inputs"]
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

            let perturbation_graph =
                PerturbationGraph::with_restricted_variables(&bn, &controllable_vars);

            let phenotype_map = config["full_mapk"]["targets"][phenotype]
                .as_object()
                .unwrap();
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

            bounded_phenotype_control(&perturbation_graph, &phenotype, max_control_size);
        }
    }

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

    // let ag = SymbolicAsyncGraph::new(bn.clone()).unwrap();
    // let w = ag.pick_witness(ag.unit_colors());
    // let perturbation_graph = PerturbationGraph::with_restricted_variables(&w, &controllable_vars);

    /*
    let result = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype, max_control_size, controllable_vars.clone(), "complex");

    let zero_perturbation_working_colors = result.perturbation_working_colors(&HashMap::from([]));
    println!("No perturbation working for {:?}", zero_perturbation_working_colors.approx_cardinality());

    let now = Instant::now();
    println!("Starting control enumeration at: {}", Local::now());

    result.ceiled_size_perturbation_working_colors(max_control_size, model_colors, &controllable_vars, false, false);

    let duration = now.elapsed();
    println!("Control enumeration finished at {:?} ", Local::now());
    println!("Time elapsed for control enumeration: {:?}", duration);
     */
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
