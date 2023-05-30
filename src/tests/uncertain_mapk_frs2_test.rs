// use rstest::rstest;
// use std::collections::HashMap;
// use biodivine_lib_param_bn::{BooleanNetwork, VariableId};
// use biodivine_lib_param_bn::biodivine_std::traits::Set;
// use biodivine_lib_param_bn::FnUpdate::Var;
// use crate::aeon::config::{get_controllable_vars, get_trivial_phenotype};
// use biodivine_lib_param_bn::symbolic_async_graph::{GraphVertices, SymbolicAsyncGraph};
// use serde_json::{Map, Value};
// use crate::aeon::phentoype::build_phenotype;
// use crate::perturbation::PerturbationGraph;
// use crate::phenotype_control::_simplified_algorithm::bounded_phenotype_control;
//
// const MAX_CONTROL: usize = 3;
//
// // Figure 5,6
// #[rstest]
// #[case(false, "growth_arrest", "full_mapk_uncertain_FRS2")]
// #[case(false, "no_decision", "full_mapk_uncertain_FRS2")]
// #[case(true, "apoptosis", "full_mapk_uncertain_FRS2")]
// #[case(true, "proliferation", "full_mapk_uncertain_FRS2")]
// #[case(true, "no_decision", "full_mapk_uncertain_FRS2")]
// fn find_all_controls_frs2(#[case] stay: bool, #[case] phenotype: &str, #[case] model:&str) {
//     let max_control_size = 3;
//
//     println!(">>>>>>> TEST: {:?} {:?} {:?}", stay, phenotype, model);
//     let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
//     let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();
//     let model_name = config[model]["file"].as_str().unwrap();
//     let model_string = std::fs::read_to_string( format!("./models_phenotype/{}", model_name)).unwrap();
//     let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
//     let mut controllable_vars = Vec::new();
//     let uncontrollable = config[model]["uncontrollable"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
//     let inputs = config[model]["inputs"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
//     for v in bn.variables() {
//         if !uncontrollable.contains(&bn.get_variable_name(v).as_str()) && !inputs.contains(&bn.get_variable_name(v).as_str()) {
//             controllable_vars.push(v);
//         }
//     }
//
//     let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &(controllable_vars.clone()));
//
//     let mut phenotype_space = get_trivial_phenotype(model, phenotype, &perturbation_graph);
//     if !stay {
//         // Don't stay -> avoid
//         phenotype_space = perturbation_graph.mk_unit_colored_vertices().minus_vertices(&phenotype_space).vertices();
//     }
//
//     bounded_phenotype_control(&perturbation_graph, &phenotype_space, max_control_size);
//
//     assert_eq!(0.0, 1.0);
// }
//
//
// #[rstest]
// #[case(false, "apoptosis", "full_mapk_uncertain_MEK1_2")]
// #[case(true, "apoptosis", "full_mapk_uncertain_MEK1_2")]
// #[case(true, "proliferation", "full_mapk_uncertain_MEK1_2")]
// #[case(true, "no_decision", "full_mapk_uncertain_DUSP1")]
// #[case(true, "growth_arrest", "full_mapk_uncertain_DUSP1")]
// #[case(true, "apoptosis", "full_mapk_uncertain_DUSP1")]
// #[case(true, "proliferation", "full_mapk_uncertain_DUSP1")]
// #[case(true, "no_decision", "full_mapk_uncertain_DUSP1")]
// fn find_all_controls_mek1_2(#[case] stay: bool, #[case] phenotype: &str, #[case] model:&str) {
//     let max_control_size = 3;
//
//     println!(">>>>>>> TEST: {:?} {:?} {:?}", stay, phenotype, model);
//     let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
//     let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();
//     let model_name = config[model]["file"].as_str().unwrap();
//     let model_string = std::fs::read_to_string( format!("./models_phenotype/{}", model_name)).unwrap();
//     let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
//     let mut controllable_vars = Vec::new();
//     let uncontrollable = config[model]["uncontrollable"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
//     let inputs = config[model]["inputs"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
//     for v in bn.variables() {
//         if !uncontrollable.contains(&bn.get_variable_name(v).as_str()) && !inputs.contains(&bn.get_variable_name(v).as_str()) {
//             controllable_vars.push(v);
//         }
//     }
//
//     let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &(controllable_vars.clone()));
//
//     let mut phenotype_space = get_trivial_phenotype(model, phenotype, &perturbation_graph);
//     if !stay {
//         // Don't stay -> avoid
//         phenotype_space = perturbation_graph.mk_unit_colored_vertices().minus_vertices(&phenotype_space).vertices();
//     }
//
//     bounded_phenotype_control(&perturbation_graph, &phenotype_space, max_control_size);
//
//     assert_eq!(0.0, 1.0);
// }
//
// #[rstest]
// #[case(true, "no_decision", "full_mapk_uncertain_DUSP1")]
// #[case(true, "growth_arrest", "full_mapk_uncertain_DUSP1")]
// #[case(true, "apoptosis", "full_mapk_uncertain_DUSP1")]
// #[case(true, "proliferation", "full_mapk_uncertain_DUSP1")]
// #[case(true, "no_decision", "full_mapk_uncertain_DUSP1")]
// fn find_all_controls_mek_dusp1(#[case] stay: bool, #[case] phenotype: &str, #[case] model:&str) {
//     let max_control_size = 3;
//
//     println!(">>>>>>> TEST: {:?} {:?} {:?}", stay, phenotype, model);
//     let config_str = std::fs::read_to_string("./models_phenotype/benchmark.json").unwrap();
//     let config: serde_json::Value = serde_json::from_str(config_str.as_str()).unwrap();
//     let model_name = config[model]["file"].as_str().unwrap();
//     let model_string = std::fs::read_to_string( format!("./models_phenotype/{}", model_name)).unwrap();
//     let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
//     let mut controllable_vars = Vec::new();
//     let uncontrollable = config[model]["uncontrollable"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
//     let inputs = config[model]["inputs"].as_array().unwrap().into_iter().map(|x| x.as_str().unwrap()).collect::<Vec<&str>>();
//     for v in bn.variables() {
//         if !uncontrollable.contains(&bn.get_variable_name(v).as_str()) && !inputs.contains(&bn.get_variable_name(v).as_str()) {
//             controllable_vars.push(v);
//         }
//     }
//
//     let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &(controllable_vars.clone()));
//
//     let mut phenotype_space = get_trivial_phenotype(model, phenotype, &perturbation_graph);
//     if !stay {
//         // Don't stay -> avoid
//         phenotype_space = perturbation_graph.mk_unit_colored_vertices().minus_vertices(&phenotype_space).vertices();
//     }
//
//     bounded_phenotype_control(&perturbation_graph, &phenotype_space, max_control_size);
//
//     assert_eq!(0.0, 1.0);
// }
//
//
//
//
// fn parse_simple_json_dict(perturbations: &str) -> HashMap<String, bool> {
//     let parsed: HashMap<String,serde_json::Value> = serde_json::from_str(perturbations).unwrap();
//     let mut perturbation_vals = HashMap::new();
//     for (k,v) in parsed.iter() {
//         perturbation_vals.insert(k.clone(), v.as_bool().unwrap());
//     }
//     perturbation_vals
// }
