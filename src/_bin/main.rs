#![allow(dead_code)] // In this file, we want to allow unused functions.

use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColors, SymbolicAsyncGraph};
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_pbn_control::aeon::reachability::backward;
use biodivine_pbn_control::control::ControlMap;
use biodivine_pbn_control::perturbation::PerturbationGraph;
use chrono::Utc;
use itertools::Itertools;
use std::convert::TryFrom;
use std::time::{Duration, Instant};

const MODELS: [&str; 5] = ["myeloid", "cardiac", "erbb", "tumour", "mapk"];
// const models1: [&str; 2] = ["myeloid", "cardiac"];
const SUFFIXES: [&str; 7] = [
    "witness", "1unknown", "2unknown", "3unknown", "4unknown", "5unknown", "6unknown",
];
const MAX_CONTROL_SIZE: usize = 5;

fn main() {
    main_models_experiment();

    // main_scalability_experiment();
    //
    // main_robustness_experiment()
}

fn main_models_experiment() {
    //main_control_template(&MODELS, &["4unknown"], PerturbationGraph::one_step_control, "one-step");
    //main_control_template(&MODELS, &["4unknown"], PerturbationGraph::temporary_control, "permanent");
    main_control_template(
        &["tumour"],
        &["4unknown"],
        PerturbationGraph::permanent_control,
        "temporary",
    );
}

fn main_scalability_experiment() {
    main_control_template(
        &["tumour"],
        &SUFFIXES,
        PerturbationGraph::one_step_control,
        "one-step",
    );
    main_control_template(
        &["tumour"],
        &SUFFIXES,
        PerturbationGraph::temporary_control,
        "temporary",
    );
    main_control_template(
        &["tumour"],
        &SUFFIXES,
        PerturbationGraph::permanent_control,
        "permanent",
    );
}

fn main_robustness_experiment() {
    for s in 0..3 {
        for t in 0..3 {
            if s == t {
                continue;
            }
            main_all_control_types_robustness("cardiac", s, t);
        }
    }
    main_all_control_types_robustness("erbb", 0, 1);
    main_all_control_types_robustness("erbb", 1, 0);
    for s in 0..3 {
        for t in 0..3 {
            if s == t {
                continue;
            }
            main_all_control_types_robustness("tumour", s, t);
        }
    }
}

fn main_all_control_types_robustness(m: &str, source_ix: usize, target_ix: usize) {
    main_control_robustness_template(
        m,
        source_ix,
        target_ix,
        PerturbationGraph::one_step_control,
        "one-step",
    );
    main_control_robustness_template(
        m,
        source_ix,
        target_ix,
        PerturbationGraph::temporary_control,
        "temporary",
    );
    main_control_robustness_template(
        m,
        source_ix,
        target_ix,
        PerturbationGraph::permanent_control,
        "permanent",
    );
}

fn main_control_template<F>(
    models: &[&str],
    suffixes: &[&str],
    control_function: F,
    control_type: &str,
) where
    F: for<'a> Fn(
        &'a PerturbationGraph,
        &'a ArrayBitVector,
        &'a ArrayBitVector,
        &'a GraphColors,
    ) -> ControlMap,
{
    println!(
        ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> {} CONTROL",
        control_type
    );
    for model_name in models {
        for model_suffix in suffixes {
            let model_path = format!("models/{}_{}.aeon", model_name, model_suffix);
            let model_string = std::fs::read_to_string(model_path).unwrap();
            let model = BooleanNetwork::try_from(model_string.as_str()).unwrap();
            let perturbation_graph = PerturbationGraph::new(&model);

            {
                // Print model statistics:
                let model_variables = model.num_vars();
                assert!(i32::try_from(model_variables).is_ok());

                // All colors considered by the perturbation graph
                let all_colors = perturbation_graph.unit_colors().approx_cardinality();
                // The (combinatorial) portion of colours that appear due to perturbation parameters.
                let perturbation_colors = 2.0f64.powi(model_variables as i32);
                // The (combinatorial) portion of colours that are carried over from the original model.
                let model_colors = all_colors / perturbation_colors;

                println!("========= {} (variables {})(uncertainty colours {})(perturbation colours {})(all colours {}) =========",
                         model_name,
                         model_variables,
                         model_colors,
                         perturbation_colors,
                         all_colors);
            }

            // Compute attractor states in the witness model.
            let start_attractors = Instant::now();
            let attractors = find_witness_attractors(model_name);
            let all_attractors_colors = attractors
                .iter()
                .map(|it| get_all_params_with_attractor(&perturbation_graph, it))
                .collect::<Vec<_>>();
            println!(
                "Attractors ready in {}ms, attractors count: {}, starting control...",
                start_attractors.elapsed().as_millis(),
                attractors.len()
            );

            let mut control_times: Vec<Duration> = Vec::new();
            for (t_i, target) in attractors.iter().enumerate() {
                for (s_i, source) in attractors.iter().enumerate() {
                    if s_i == t_i {
                        continue;
                    }

                    let attractor_colors = all_attractors_colors[t_i].clone();
                    let start = Instant::now();
                    let control =
                        control_function(&perturbation_graph, source, target, &attractor_colors);
                    println!(
                        "Control from attr. #{:?} (source) to attr. #{:?} (target) exists for {} color(s), jumping through {} vertices.",
                        s_i,
                        t_i,
                        control.controllable_colors_cardinality(),
                        control.jump_vertices()
                    );
                    let elapsed = start.elapsed();
                    println!("Elapsed: {} ms", elapsed.as_millis());
                    control_times.push(elapsed);
                }
            }

            let total_time_ms = control_times.iter().sum::<Duration>().as_millis();
            let average_time_ms = (total_time_ms as f32) / (control_times.len() as f32);
            println!(
                "<<<<<<<<<<<< Average {} control time for {} model with suffix `{}` is {:.2} ms",
                control_type, model_name, model_suffix, average_time_ms,
            );
        }
    }
}

fn main_control_robustness_template<F>(
    m: &str,
    source_ix: usize,
    target_ix: usize,
    control_function: F,
    control_type: &str,
) where
    F: for<'a> Fn(
        &'a PerturbationGraph,
        &'a ArrayBitVector,
        &'a ArrayBitVector,
        &'a GraphColors,
    ) -> ControlMap,
{
    println!(
        "Robustness of {} control in model {}, source: {}, target: {}",
        control_type, m, source_ix, target_ix
    );
    assert_ne!(source_ix, target_ix);

    let model_string: &str =
        &std::fs::read_to_string(format!("models/{}_4unknown.aeon", m)).unwrap();
    let model = BooleanNetwork::try_from(model_string).unwrap();
    let perturbations = PerturbationGraph::new(&model);
    println!(
        "========= {}(v{})(p{}) =========",
        m,
        model.num_vars(),
        perturbations
            .as_original()
            .unit_colors()
            .approx_cardinality()
    );
    let attractors = find_witness_attractors(m);
    println!("Attractors count: {}", attractors.len());
    let source = attractors.get(source_ix).unwrap();
    let target = attractors.get(target_ix).unwrap();
    let att_colors = get_all_params_with_attractor(&perturbations, target);
    println!(
        "Attractor params cardinality: {:?}",
        att_colors.approx_cardinality()
    );
    let start = Instant::now();
    let control = control_function(&perturbations, source, target, &att_colors);
    println!(
        "Control from Attractor {:?} (source) to Attractor {:?} (target) works for {} color(s), jumping through {} vertices.",
        source_ix,
        target_ix,
        control.controllable_colors_cardinality(),
        control.jump_vertices()
    );
    println!("Elapsed: {}ms", start.elapsed().as_millis());

    let robustness_all = control.controllable_colors_cardinality();

    let mut current_iter = 0;
    let mut max_robustness = 0.0;
    let mut union_robustness = perturbations
        .empty_colors()
        .as_bdd()
        .or(perturbations.empty_colors().as_bdd());
    for combination in model.variables().powerset() {
        let control_size = combination.len();
        if control_size > current_iter {
            println!(
                "[{}] Control of size {:?} max robustnesss: {:?}, union robustness: {:?}",
                Utc::now().format("%H:%M:%S %d.%m.%Y"),
                current_iter,
                max_robustness / robustness_all,
                union_robustness.cardinality() / robustness_all
            );
            current_iter = control_size;
            max_robustness = 0.0;
            union_robustness = perturbations
                .empty_colors()
                .as_bdd()
                .or(perturbations.empty_colors().as_bdd());
        }
        if control_size > MAX_CONTROL_SIZE {
            break;
        }
        let mut local_control = control.clone();
        for var_id in model.variables() {
            if combination.contains(&var_id) {
                local_control.require_perturbation(var_id, None);
            } else {
                local_control.exclude_perturbation(var_id, None);
            }
        }
        let card = local_control.controllable_colors_cardinality();
        if max_robustness < card {
            max_robustness = card;
        }

        union_robustness = union_robustness.or(&local_control.controllable_colors())
    }
}

/// Compute possible source-target attractor pairs for a network.
///
/// At the moment, we don't consider every pair of attractor states, we just pick one state
/// from the first attractor as target, and then different source states from the remaining attractors.
fn compute_attractor_pairs(network: &BooleanNetwork) -> Vec<(ArrayBitVector, ArrayBitVector)> {
    let graph = SymbolicAsyncGraph::new(network).unwrap();
    let attractors = biodivine_pbn_control::aeon::attractors::compute(&graph);
    let target: ArrayBitVector = attractors[0]
        .vertices()
        .materialize()
        .iter()
        .next()
        .unwrap();
    let mut result = Vec::new();
    for source in attractors.iter().skip(1) {
        let source = source.vertices().materialize().iter().next().unwrap();
        result.push((source, target.clone()));
    }

    result
}

fn find_witness_attractors(m: &str) -> Vec<ArrayBitVector> {
    let model_string: &str =
        &std::fs::read_to_string(format!("models/{}_witness.aeon", m)).unwrap();
    let model = BooleanNetwork::try_from(model_string).unwrap();
    let graph = SymbolicAsyncGraph::new(&model).unwrap();
    let attractors = biodivine_pbn_control::aeon::attractors::compute(&graph);
    let mut vertices = Vec::new();
    for a in attractors {
        vertices.push(a.pick_vertex());
    }

    vertices
        .into_iter()
        .map(|x| x.vertices().materialize().iter().next().unwrap())
        .collect()
}

pub fn get_all_params_with_attractor(
    graph: &PerturbationGraph,
    state: &ArrayBitVector,
) -> GraphColors {
    let seed = graph.vertex(state);
    let bwd = backward(graph.as_original(), &seed);
    let mut attractor = seed;
    'forward: loop {
        if attractor.as_bdd().size() > 10_000 {
            println!("FWD: {}", attractor.as_bdd().size());
        }
        for var in graph.as_original().variables().rev() {
            let step = graph
                .as_original()
                .var_post(var, &attractor)
                .minus(&attractor);

            if !step.is_empty() {
                // Any color that is found for a state not in BWD is not attractor color.
                let not_attractor = step.minus(&bwd).colors();
                attractor = attractor.union(&step).minus_colors(&not_attractor);
                continue 'forward;
            }
        }
        // No new steps found. Attractor is now truly the set of attractor states:
        return attractor.colors();
    }
}
