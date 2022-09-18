#![allow(dead_code)]    // In this file, we want to allow unused functions.

use std::borrow::{Borrow, BorrowMut};
use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColors, SymbolicAsyncGraph};
use biodivine_lib_param_bn::{BooleanNetwork, VariableId};
use biodivine_pbn_control::perturbation::PerturbationGraph;
use std::convert::TryFrom;
use std::time::Instant;
use biodivine_lib_param_bn::biodivine_std::traits::{Graph, Set};
use itertools::Itertools;
use biodivine_pbn_control::control::ControlMap;
use chrono::{DateTime, Utc};
use biodivine_pbn_control::aeon::reachability::{backward, forward};

// const models1: [&str; 5] = ["myeloid", "cardiac", "erbb", "tumour", "mapk"];
const models1: [&str; 2] = ["myeloid", "cardiac"];
const suffixes2: [&str; 7] = ["witness", "1unknown", "2unknown", "3unknown", "4unknown", "5unknown", "6unknown"];
const MAX_CONTROL_SIZE: usize = 5;


fn main() {
    main_models_experiment();

    // main_scalability_experiment();
    //
    // main_robustness_experiment()
}

fn main_models_experiment() {
    main_control_template(models1.to_vec(), vec!["4unknown"], PerturbationGraph::one_step_control, "one-step");
    main_control_template(models1.to_vec(), vec!["4unknown"], PerturbationGraph::temporary_control, "temporary");
    main_control_template(models1.to_vec(), vec!["4unknown"], PerturbationGraph::permanent_control, "permanent");
}

fn main_scalability_experiment() {
    main_control_template(vec!["tumour"], suffixes2.to_vec(), PerturbationGraph::one_step_control, "one-step");
    main_control_template(vec!["tumour"], suffixes2.to_vec(), PerturbationGraph::temporary_control, "temporary");
    main_control_template(vec!["tumour"], suffixes2.to_vec(), PerturbationGraph::permanent_control, "permanent");
}

fn main_robustness_experiment() {
    for s in 0..3 {
        for t in 0..3 {
            if s == t {
                continue
            }
            main_all_control_types_robustness("cardiac", s, t);
        }
    }
    main_all_control_types_robustness("erbb", 0, 1);
    main_all_control_types_robustness("erbb", 1, 0);
    for s in 0..3 {
        for t in 0..3 {
            if s == t {
                continue
            }
            main_all_control_types_robustness("tumour", s, t);
        }
    }
}

fn main_all_control_types_robustness(m: &str, source_ix: usize, target_ix: usize) {
    main_control_robustness_template(m, source_ix, target_ix, PerturbationGraph::one_step_control, "one-step");
    main_control_robustness_template(m, source_ix, target_ix, PerturbationGraph::temporary_control, "temporary");
    main_control_robustness_template(m, source_ix, target_ix, PerturbationGraph::permanent_control, "permanent");
}

fn main_control_template<F>(models: Vec<&str>, suffixes: Vec<&str>, control_function: F, control_type: &str)
    where F: for <'a> Fn(&'a PerturbationGraph, &'a ArrayBitVector, &'a ArrayBitVector, &'a GraphColors) -> ControlMap
{
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> {} CONTROL", control_type);
    for m in models.clone().iter() {
        let mut results = Vec::new();
        for s in suffixes.clone().iter() {
            let model_string: &str =
                &std::fs::read_to_string(format!("models/{}_{}.aeon", m, s)).unwrap();
            let model = BooleanNetwork::try_from(model_string).unwrap();
            let perturbations = PerturbationGraph::new(&model);
            let all_colours_car = perturbations.as_original().unit_colors().approx_cardinality();
            println!("========= {} (variables {})(uncertainty colours {})(perturbation colours {})(all colours {}) =========",
                     m,
                     model.num_vars(),
                     all_colours_car / f64::from(i32::pow(i32::try_from(model.num_vars()).unwrap(), 2)),
                     i32::pow(i32::try_from(model.num_vars()).unwrap(), 2),
                     all_colours_car);
            let start = Instant::now();
            let attractors = find_witness_attractors(m);
            println!(
                "Attractors ready in {}ms, attractors count: {}, starting control...",
                start.elapsed().as_millis(),
                attractors.len()
            );

            for (t_i, target) in attractors.clone().iter().enumerate() {
                for (s_i, source) in attractors.clone().iter().enumerate() {
                    if s_i == t_i {
                        continue
                    }
                    let att_colors = get_all_params_with_attractor(perturbations.borrow(), target);
                    let start = Instant::now();
                    let control = control_function(&perturbations, source, target, &att_colors);
                    println!(
                        "Control from Attractor {:?} (source) to Attractor {:?} (target) works for {} color(s), jumping through {} vertices.",
                        s_i,
                        t_i,
                        control.controllable_colors_cardinality(),
                        control.jump_vertices()
                    );
                    let el_time = start.elapsed().as_millis();
                    println!("Elapsed: {} ms", el_time);
                    results.push(el_time);
                }
            }
        }
        let avg_time = results.iter().sum::<u128>() as f32 / results.len() as f32;
        println!("<<<<<<<<<<<< Average {} control time for {} model is {} ms", control_type, m, avg_time);
    }
}

fn main_control_robustness_template<F>(m: &str, source_ix: usize, target_ix: usize, control_function: F, control_type: &str)
    where F: for <'a> Fn(&'a PerturbationGraph, &'a ArrayBitVector, &'a ArrayBitVector, &'a GraphColors) -> ControlMap
{
    println!("Robustness of {} control in model {}, source: {}, target: {}", control_type, m, source_ix, target_ix);
    assert_ne!(source_ix, target_ix);

    let model_string: &str = &std::fs::read_to_string(format!("models/{}_4unknown.aeon", m)).unwrap();
    let model = BooleanNetwork::try_from(model_string).unwrap();
    let perturbations = PerturbationGraph::new(&model);
    println!("========= {}(v{})(p{}) =========", m, model.num_vars(), perturbations.as_original().unit_colors().approx_cardinality());
    let attractors = find_witness_attractors(m);
    println!("Attractors count: {}", attractors.len());
    let source = attractors.get(source_ix).unwrap();
    let target = attractors.get(target_ix). unwrap();
    let att_colors = get_all_params_with_attractor(perturbations.borrow(), target);
    println!("Attractor params cardinality: {:?}", att_colors.approx_cardinality());
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
    let mut union_robustness = perturbations.empty_colors().as_bdd().or(perturbations.empty_colors().as_bdd());
    for combination in model.variables().into_iter().powerset() {
        let control_size = combination.len();
        if control_size > current_iter {
            println!("[{}] Control of size {:?} max robustnesss: {:?}, union robustness: {:?}",
                     Utc::now().format("%H:%M:%S %d.%m.%Y"),
                     current_iter,
                     max_robustness/robustness_all,
                     union_robustness.cardinality()/robustness_all);
            current_iter = control_size;
            max_robustness = 0.0;
            union_robustness = perturbations.empty_colors().as_bdd().or(perturbations.empty_colors().as_bdd());
        }
        if control_size > MAX_CONTROL_SIZE {
            break
        }
        let mut local_control = control.clone();
        for varId in model.variables() {
            if combination.contains(&varId) {
                local_control.require_perturbation(varId, None);
            } else {
                local_control.exclude_perturbation(varId, None);
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
    let graph = SymbolicAsyncGraph::new(network.clone()).unwrap();
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
    let graph = SymbolicAsyncGraph::new(model.clone()).unwrap();
    let attractors = biodivine_pbn_control::aeon::attractors::compute(&graph);
    let mut vertices = Vec::new();
    for a in attractors {
        vertices.push(a.pick_vertex());
    }

    vertices.into_iter().map(|x| x.vertices().materialize().iter().next().unwrap()).collect()
}


pub fn get_all_params_with_attractor(graph: &PerturbationGraph, state: &ArrayBitVector) -> GraphColors {
    let seed = graph.vertex(state);
    let fwd = forward(graph.as_original(), seed.borrow());
    let bwd = backward(graph.as_original(), seed.borrow());
    let scc = fwd.intersect(&bwd);
    let not_attractor_colors = fwd.minus(&scc).colors();
    let attractor = scc.minus_colors(&not_attractor_colors);
    return attractor.colors();
}