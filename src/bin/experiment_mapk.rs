use std::collections::HashMap;
use biodivine_pbn_control::experiment_utils::{parse_experiment, run_control_experiment};
use biodivine_pbn_control::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::fixed_points::FixedPoints;
use biodivine_lib_param_bn::{BooleanNetwork, VariableId};
use biodivine_lib_param_bn::symbolic_async_graph::{GraphVertices, SymbolicAsyncGraph};
use std::time::Instant;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_pbn_control::aeon::phentoype::build_phenotype;


fn main() {
    let mapk_full = "[id-070]__[var-49]__[in-4]__[MAPK-CANCER-CELL-FATE].aeon";
    let mapk_reduced = "[id-089]__[var-13]__[in-4]__[MAPK-REDUCED-1]__egfr_oe_witness_0000.aeon";

    let mapk_controllable = get_controllable_vars(mapk_full, "full_mapk", vec![]);
    let mapk_reduced_controllable = get_controllable_vars(mapk_reduced, "reduced_mapk", vec![]);
    // let mapk_lof_p53_controllable = get_controllable_vars(mapk_lof_p53, "reduced_mapk", vec!["v_p53"]);
    // let mapk_lof_pten_p14 = get_controllable_vars(mapk_lof_pten_p14, "reduced_mapk", vec!["v_PTEN, v_p14"]);

    let mut apoptosis = HashMap::new();
    apoptosis.insert("v_Apoptosis", true);
    apoptosis.insert("v_Growth_Arrest", true);
    apoptosis.insert("v_Proliferation", false );

    let mut no_decision = HashMap::new();
    no_decision.insert("v_Apoptosis", false);
    no_decision.insert("v_Growth_Arrest", false);
    no_decision.insert("v_Proliferation", false );

    let mut proliferation = HashMap::new();
    proliferation.insert("v_Apoptosis", false);
    proliferation.insert("v_Growth_Arrest", false);
    proliferation.insert("v_Proliferation", true );

    let phenotypes = vec![apoptosis.clone(), no_decision.clone(), proliferation.clone()];
    // let phenotypes = vec![apoptosis];
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", mapk_reduced)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &mapk_reduced_controllable);

    // bn.find_parameter()


    for p in phenotypes {
        println!("{:?}", p);
        // Find controls in normal networks
        let phenotype = build_phenotype(perturbation_graph.as_perturbed(), p);
        control_ceiled(mapk_reduced, &mapk_reduced_controllable, phenotype);
        // control_ceiled_to_3(mapk_full, &mapk_controllable, p.clone());
    }


    println!("{:?}", "Not apoptosis");
    // Find controls in normal networks
    let apoptosis_space = build_phenotype(perturbation_graph.as_perturbed(), apoptosis);
    let not_apoptosis = perturbation_graph.mk_unit_colored_vertices().vertices().minus(&apoptosis_space);
    control_ceiled(mapk_reduced, &mapk_reduced_controllable, not_apoptosis);


    println!("{:?}", "Not no_decision");
    // Find controls in normal networks
    let no_decision_space = build_phenotype(perturbation_graph.as_perturbed(), no_decision);
    let not_no_decision = perturbation_graph.mk_unit_colored_vertices().vertices().minus(&no_decision_space);
    control_ceiled(mapk_reduced, &mapk_reduced_controllable, not_no_decision);


    println!("{:?}", "Not proliferation");
    // Find controls in normal networks
    let proliferation_space = build_phenotype(perturbation_graph.as_perturbed(), proliferation);
    let no_proliferation_space = perturbation_graph.mk_unit_colored_vertices().vertices().minus(&proliferation_space);
    control_ceiled(mapk_reduced, &mapk_reduced_controllable, no_proliferation_space);
}

const MAX_CONTROL: usize = 1;

fn control_ceiled(model_file: &str, controllable_vars: &Vec<VariableId>, phenotype: GraphVertices) {
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", model_file)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &controllable_vars);
    let all_colors = perturbation_graph.unit_colors().approx_cardinality();
    let perturbation_colors = 2.0f64.powi(controllable_vars.len() as i32);
    let model_colors = all_colors / perturbation_colors;

    println!("{:?}", controllable_vars);
    let result = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype, MAX_CONTROL, controllable_vars.clone(), "heuristic");

    println!("{:?}", result.as_bdd().to_dot_string(&perturbation_graph.as_symbolic_context().bdd_variable_set(), true));

    let zero_perturbation_working_colors = result.perturbation_working_colors(&HashMap::from([]));
    println!("No perturbation working for {:?}: {:?}", zero_perturbation_working_colors.approx_cardinality(),  zero_perturbation_working_colors.to_dot_string(perturbation_graph.as_symbolic_context()));

    result.ceiled_size_perturbation_working_colors(MAX_CONTROL, model_colors, &controllable_vars.clone(), false, true);
}

fn get_controllable_vars(model_file: &str, model_name: &str, extra_forbidden: Vec<&str>) -> Vec<VariableId> {
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
