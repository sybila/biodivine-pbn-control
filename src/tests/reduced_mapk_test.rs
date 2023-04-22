use rstest::rstest;
use std::collections::HashMap;
use biodivine_lib_param_bn::{BooleanNetwork, VariableId};
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use crate::aeon::config::{get_controllable_vars, get_trivial_phenotype};
use biodivine_lib_param_bn::symbolic_async_graph::GraphVertices;
use serde_json::{Map, Value};
use crate::aeon::phentoype::build_phenotype;
use crate::perturbation::PerturbationGraph;

static MAPK_REDUCED_KEY: &str = "reduced_mapk";
const MAX_CONTROL: usize = 3;

// Figure 4C, EGFR over-expression - hardcoded witness model
#[rstest]
#[case(false, "no_decision", "{}", "sinks")] // r3
#[case(true, "proliferation", "{\"v_p53\": false}", "sinks")] // r5
#[case(true, "apoptosis", "{\"v_DNA_damage\": true}", "sinks")] // r7
#[case(true, "apoptosis", "{\"v_TGFBR_stimulus\": true}", "sinks")] // r13
#[case(false, "apoptosis", "{\"v_AKT\": true, \"v_PI3K\": true}", "sinks")] // r11
#[case(true, "proliferation", "{\"v_p14\": false}", "sinks")] // r9
#[case(false, "apoptosis", "{\"v_PTEN\": false}", "sinks")] // r15
#[case(false, "no_decision", "{}", "complex")] // r3
#[case(true, "proliferation", "{\"v_p53\": false}", "complex")] // r5
#[case(true, "apoptosis", "{\"v_DNA_damage\": true}", "complex")] // r7
#[case(true, "apoptosis", "{\"v_TGFBR_stimulus\": true}", "complex")] // r13
#[case(false, "apoptosis", "{\"v_AKT\": true, \"v_PI3K\": true}", "complex")] // r11
#[case(true, "proliferation", "{\"v_p14\": false}", "complex")] // r9
#[case(false, "apoptosis", "{\"v_PTEN\": false}", "complex")] // r15
fn mapk_egfr_oe_witness_phenotype(#[case] stay: bool,
                                  #[case] phenotype: &str,
                                  #[case] perturbed_vals: &str,
                                  #[case] attractor_serach: &str) {
    const  MAPK_REDUCED_FILE : &str = "[id-089]__[var-13]__[in-4]__[MAPK-REDUCED-1]__egfr_oe_witness_0000.aeon";
    println!(">>>>>>> TEST: {:?} {:?} {:?} {:?}", stay, phenotype, perturbed_vals, attractor_serach);
    let mapk_reduced_controllable = get_controllable_vars(MAPK_REDUCED_FILE, MAPK_REDUCED_KEY, vec![]);
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", MAPK_REDUCED_FILE)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &mapk_reduced_controllable);

    let mut phenotype_space = get_trivial_phenotype(MAPK_REDUCED_KEY, phenotype, &perturbation_graph);
    if !stay {
        // Don't stay -> avoid
        phenotype_space = perturbation_graph.mk_unit_colored_vertices().minus_vertices(&phenotype_space).vertices();
    }

    let control_map = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype_space, MAX_CONTROL, mapk_reduced_controllable, attractor_serach);

    let perturbation = parse_simple_json_dict(perturbed_vals);
    let working_colors = control_map.perturbation_working_colors(&perturbation).approx_cardinality();
    assert_eq!(working_colors, 1.0);
}


// Figure 4C, FGFR3 gain-of-function - hardcoded witness model
#[rstest]
#[case(false, "apoptosis", "{\"v_p53\": false}", "sinks")] // r6
#[case(true, "apoptosis", "{\"v_DNA_damage\": true}", "sinks")] // r8
#[case(true, "apoptosis", "{\"v_TGFBR_stimulus\": true}", "sinks")] // r14
#[case(false, "apoptosis", "{\"v_AKT\": true, \"v_PI3K\": true}", "sinks")] // r12
#[case(false, "apoptosis", "{\"v_p14\": false}", "sinks")] // r10
#[case(false, "apoptosis", "{\"v_PTEN\": false}", "sinks")] // r16
#[case(false, "apoptosis", "{\"v_p53\": false}", "complex")] // r6
#[case(true, "apoptosis", "{\"v_DNA_damage\": true}", "complex")] // r8
#[case(true, "apoptosis", "{\"v_TGFBR_stimulus\": true}", "complex")] // r14
#[case(false, "apoptosis", "{\"v_AKT\": true, \"v_PI3K\": true}", "complex")] // r12
#[case(false, "apoptosis", "{\"v_p14\": false}", "complex")] // r10
#[case(false, "apoptosis", "{\"v_PTEN\": false}", "complex")] // r16
fn mapk_fgfr3_oe_witness_phenotype(#[case] stay: bool,
                                   #[case] phenotype: &str,
                                   #[case] perturbed_vals: &str,
                                   #[case] attractor_search: &str) {
    const  MAPK_REDUCED_FILE: &str = "[id-089]__[var-13]__[in-4]__[MAPK-REDUCED-1]__egfr_oe_witness_0000.aeon";
    println!(">>>>>>> TEST: {:?} {:?} {:?} {:?}", stay, phenotype, perturbed_vals, attractor_search);
    let mapk_reduced_controllable = get_controllable_vars(MAPK_REDUCED_FILE, MAPK_REDUCED_KEY, vec![]);
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", MAPK_REDUCED_FILE)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &mapk_reduced_controllable);

    let mut phenotype_space = get_trivial_phenotype(MAPK_REDUCED_KEY, phenotype, &perturbation_graph);
    if !stay {
        // Don't stay -> avoid
        phenotype_space = perturbation_graph.mk_unit_colored_vertices().minus_vertices(&phenotype_space).vertices();
    }

    let control_map = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype_space, MAX_CONTROL, mapk_reduced_controllable, attractor_search);

    let perturbation = parse_simple_json_dict(perturbed_vals);
    let working_colors = control_map.perturbation_working_colors(&perturbation).approx_cardinality();
    assert_eq!(working_colors, 1.0);
}


// Figure 4C, EGFR over-expression - witness model (all inputs 0), EGFR normal and needs to be perturbed
#[rstest]
#[case(false, "no_decision", "{\"v_EGFR\": true}")] // r3
#[case(true, "proliferation", "{\"v_p53\": false, \"v_EGFR\": true}")] // r5
#[case(true, "apoptosis", "{\"v_DNA_damage\": true, \"v_EGFR\": true}")] // r7
#[case(true, "apoptosis", "{\"v_TGFBR_stimulus\": true, \"v_EGFR\": true}")] // r13
#[case(false, "apoptosis", "{\"v_AKT\": true, \"v_PI3K\": true, \"v_EGFR\": true}")] // r11
#[case(true, "proliferation", "{\"v_p14\": false, \"v_EGFR\": true}")] // r9
#[case(false, "apoptosis", "{\"v_PTEN\": false, \"v_EGFR\": true}")] // r15
fn mapk_witness_phenotype_with_egfr_oe_perturbation(#[case] stay: bool,
                                  #[case] phenotype: &str,
                                  #[case] perturbed_vals: &str) {
    const  MAPK_REDUCED_FILE: &str = "[id-089]__[var-13]__[in-4]__[MAPK-REDUCED-1]__witness_0000.aeon";
    println!(">>>>>>> TEST: {:?} {:?} {:?}", stay, phenotype, perturbed_vals);
    let mapk_reduced_controllable = get_controllable_vars(MAPK_REDUCED_FILE, MAPK_REDUCED_KEY, vec![]);
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", MAPK_REDUCED_FILE)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &mapk_reduced_controllable);

    let mut phenotype_space = get_trivial_phenotype(MAPK_REDUCED_KEY, phenotype, &perturbation_graph);
    if !stay {
        // Don't stay -> avoid
        phenotype_space = perturbation_graph.mk_unit_colored_vertices().minus_vertices(&phenotype_space).vertices();
    }

    let control_map = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype_space, MAX_CONTROL, mapk_reduced_controllable, "complex");

    let perturbation = parse_simple_json_dict(perturbed_vals);
    let working_colors = control_map.perturbation_working_colors(&perturbation).approx_cardinality();
    assert_eq!(working_colors, 1.0);
}


// Figure 4C, FGFR3 gain of function - witness model (all inputs 0), FGFR3 normal and needs to be perturbed
#[rstest]
#[case(false, "apoptosis", "{\"v_p53\": false, \"v_FGFR3\": true}")] // r6
#[case(true, "apoptosis", "{\"v_DNA_damage\": true, \"v_FGFR3\": true}")] // r8
#[case(true, "apoptosis", "{\"v_TGFBR_stimulus\": true, \"v_FGFR3\": true}")] // r14
#[case(false, "apoptosis", "{\"v_AKT\": true, \"v_PI3K\": true, \"v_FGFR3\": true}")] // r12
#[case(false, "apoptosis", "{\"v_p14\": false, \"v_FGFR3\": true}")] // r10
#[case(false, "apoptosis", "{\"v_PTEN\": false, \"v_FGFR3\": true}")] // r16
fn mapk_witness_phenotype_require_fgfr3_oe_perturbation(#[case] stay: bool,
                                   #[case] phenotype: &str,
                                   #[case] perturbed_vals: &str) {
    const  MAPK_REDUCED_FILE: &str = "[id-089]__[var-13]__[in-4]__[MAPK-REDUCED-1]__witness_0000.aeon";
    println!(">>>>>>> TEST: {:?} {:?} {:?}", stay, phenotype, perturbed_vals);
    let mapk_reduced_controllable = get_controllable_vars(MAPK_REDUCED_FILE, MAPK_REDUCED_KEY, vec![]);
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", MAPK_REDUCED_FILE)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &mapk_reduced_controllable);

    let mut phenotype_space = get_trivial_phenotype(MAPK_REDUCED_KEY, phenotype, &perturbation_graph);
    if !stay {
        // Don't stay -> avoid
        phenotype_space = perturbation_graph.mk_unit_colored_vertices().minus_vertices(&phenotype_space).vertices();
    }

    let control_map = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype_space, MAX_CONTROL, mapk_reduced_controllable, "complex");

    let perturbation = parse_simple_json_dict(perturbed_vals);
    let working_colors = control_map.perturbation_working_colors(&perturbation).approx_cardinality();
    assert_eq!(working_colors, 1.0);
}


// Figure 4C, EGFR over-expression - model with inputs, EGFR normal and needs to be perturbed, working just for colors with correct inputs
#[rstest]
#[case(false, "no_decision", "{}", "{\"v_EGFR\": true}")] // r3
#[case(true, "proliferation", "{\"v_p53\": false, }", "{\"v_EGFR\": true}")] // r5
#[case(true, "apoptosis", "{}", "{\"v_EGFR\": true, \"v_DNA_damage\": true}")] // r7
#[case(true, "apoptosis", "{}", "{\"v_EGFR\": true, \"v_TGFBR_stimulus\": true}")] // r13
#[case(false, "apoptosis", "{\"v_AKT\": true, \"v_PI3K\": true, }", "{\"v_EGFR\": true}")] // r11
#[case(true, "proliferation", "{\"v_p14\": false, }", "{\"v_EGFR\": true}")] // r9
#[case(false, "apoptosis", "{\"v_PTEN\": false, }", "{\"v_EGFR\": true}")] // r15
fn mapk_phenotype_require_egfr_oe_perturbation(#[case] stay: bool,
                                               #[case] phenotype: &str,
                                               #[case] perturbed_vals: &str,
                                               #[case] color_vals: &str) {
    const  MAPK_REDUCED_FILE: &str = "[id-089]__[var-13]__[in-4]__[MAPK-REDUCED-1].aeon";
    println!(">>>>>>> TEST: {:?} {:?} {:?}", stay, phenotype, perturbed_vals);
    let mapk_reduced_controllable = get_controllable_vars(MAPK_REDUCED_FILE, MAPK_REDUCED_KEY, vec![]);
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", MAPK_REDUCED_FILE)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &mapk_reduced_controllable);

    let mut phenotype_space = get_trivial_phenotype(MAPK_REDUCED_KEY, phenotype, &perturbation_graph);
    if !stay {
        // Don't stay -> avoid
        phenotype_space = perturbation_graph.mk_unit_colored_vertices().minus_vertices(&phenotype_space).vertices();
    }

    let control_map = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype_space, MAX_CONTROL, mapk_reduced_controllable, "complex");

    let perturbation = parse_simple_json_dict(perturbed_vals);
    let working_colors = control_map.perturbation_working_colors(&perturbation);

    let parsed_color = parse_simple_json_dict(color_vals);
    let tested_color = perturbation_graph.build_colors_with_values(&bn, parsed_color);
    assert_eq!(working_colors.intersect(&tested_color).approx_cardinality(), 1.0);
}


// Figure 4C, FGFR3 gain of function - model with inputs, FGFR3 normal and needs to be perturbed, working just for colors with correct input values
#[rstest]
#[case(false, "apoptosis", "{\"v_p53\": false, }", "{\"v_FGFR3\": true}")] // r6
#[case(true, "apoptosis", "{}", "{\"v_FGFR3\": true, \"v_DNA_damage\": true}")] // r8
#[case(true, "apoptosis", "{}", "{\"v_FGFR3\": true, \"v_TGFBR_stimulus\": true}")] // r14
#[case(false, "apoptosis", "{\"v_AKT\": true, \"v_PI3K\": true, }", "{\"v_FGFR3\": true}")] // r12
#[case(false, "apoptosis", "{\"v_p14\": false, }", "{\"v_FGFR3\": true}")] // r10
#[case(false, "apoptosis", "{\"v_PTEN\": false, }", "{\"v_FGFR3\": true}")] // r16
fn mapk_phenotype_require_fgfr3_oe_perturbation(#[case] stay: bool,
                                                #[case] phenotype: &str,
                                                #[case] perturbed_vals: &str,
                                                #[case] color_vals: &str) {
    const  MAPK_REDUCED_FILE: &str = "[id-089]__[var-13]__[in-4]__[MAPK-REDUCED-1]__witness_0000.aeon";
    println!(">>>>>>> TEST: {:?} {:?} {:?}", stay, phenotype, perturbed_vals);
    let mapk_reduced_controllable = get_controllable_vars(MAPK_REDUCED_FILE, MAPK_REDUCED_KEY, vec![]);
    let model_string = std::fs::read_to_string(format!("./models_phenotype/{}", MAPK_REDUCED_FILE)).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &mapk_reduced_controllable);

    let mut phenotype_space = get_trivial_phenotype(MAPK_REDUCED_KEY, phenotype, &perturbation_graph);
    if !stay {
        // Don't stay -> avoid
        phenotype_space = perturbation_graph.mk_unit_colored_vertices().minus_vertices(&phenotype_space).vertices();
    }

    let control_map = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype_space, MAX_CONTROL, mapk_reduced_controllable, "complex");

    let perturbation = parse_simple_json_dict(perturbed_vals);
    let working_colors = control_map.perturbation_working_colors(&perturbation);

    let parsed_color = parse_simple_json_dict(color_vals);
    let tested_color = perturbation_graph.build_colors_with_values(&bn, parsed_color);
    assert_eq!(working_colors.intersect(&tested_color).approx_cardinality(), 1.0);
}

fn parse_simple_json_dict(perturbations: &str) -> HashMap<String, bool> {
    let parsed: HashMap<String,serde_json::Value> = serde_json::from_str(perturbations).unwrap();
    let mut perturbation_vals = HashMap::new();
    for (k,v) in parsed.iter() {
        perturbation_vals.insert(k.clone(), v.as_bool().unwrap());
    }
    perturbation_vals
}
