#![allow(dead_code)]    // In this file, we want to allow unused functions.

use std::borrow::{Borrow, BorrowMut};
use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_lib_param_bn::{BooleanNetwork, VariableId};
use biodivine_pbn_control::perturbation::PerturbationGraph;
use std::convert::TryFrom;
use std::time::Instant;
use itertools::Itertools;
use biodivine_pbn_control::control::ControlMap;
use chrono::{DateTime, Utc};

// const models1: [&str; 5] = ["myeloid", "cardiac", "erbb", "tumour", "mapk"];
// const models1: [&str; 3] = ["myeloid", "cardiac", "erbb"];
// const suffixes2: [&str; 7] = ["witness", "1unknown", "2unknown", "3unknown", "4unknown", "5unknown", "6unknown"];
const MAX_CONTROL_SIZE: usize = 5;


fn main() {
    // main_one_step_old_benchmark("myeloid", ["witness", "4params", "8params"].to_vec());
    // main_one_step_old_benchmark("cell_fate", ["7stable_attractors", "7stable_attractors_2params", "7stable_attractors_4params"].to_vec());

    // main_one_step(models1.to_vec(), vec!["4unknown"]);
    // main_permanent(models1.to_vec(), vec!["4unknown"]);
    // main_temporary(models1.to_vec(), vec!["4unknown"]);
    //
    // main_one_step(vec!["tumour"],suffixes2.to_vec());
    // main_permanent(vec!["tumour"],suffixes2.to_vec());
    // main_temporary(vec!["tumour"],suffixes2.to_vec());

    for s in 0..3 {
        for t in 0..3 {
            if s == t {
                continue
            }
            main_all_robustness("cardiac", s, t);
        }
    }
    main_all_robustness("erbb", 0, 1);
    main_all_robustness("erbb", 1, 0);
    for s in 0..3 {
        for t in 0..3 {
            if s == t {
                continue
            }
            main_all_robustness("tumour", s, t);
        }
    }
}

fn main_all_robustness(m: &str, source_ix: usize, target_ix: usize) {
    main_robustness(m, source_ix, target_ix, PerturbationGraph::one_step_control, "one-step");
    main_robustness(m, source_ix, target_ix, PerturbationGraph::temporary_control, "temporary");
    main_robustness(m, source_ix, target_ix, PerturbationGraph::permanent_control, "permanent");
}

fn main_robustness<F>(m: &str, source_ix: usize, target_ix: usize, control_function: F, control_type: &str)
    where F: for <'a> Fn(&'a PerturbationGraph, &'a ArrayBitVector, &'a ArrayBitVector) -> ControlMap<'a>
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
    let start = Instant::now();
    let control = control_function(&perturbations, source, target);
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
        if control_size > MAX_CONTROL_SIZE {
            break
        }
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


fn main_one_step_old_benchmark(m: &str, suffixes: Vec<&str>) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ONE STEP CONTROL - OLD BENCHMARK");
    
    for s in suffixes.clone().iter() {
        println!("models/old_benchmark/{}_{}.aeon", m, s);
        let model_string: &str =
            &std::fs::read_to_string(format!("models/old_benchmark/{}_{}.aeon", m, s)).unwrap();
        let model = BooleanNetwork::try_from(model_string).unwrap();
        let perturbations = PerturbationGraph::new(&model);
        println!("========= {}(v{})(p{}) =========", m, model.num_vars(), perturbations.as_original().unit_colors().approx_cardinality());
        let start = Instant::now();
        let attractors = find_witness_attractors(m);
        println!(
            "Attractors ready in {}ms, starting control...",
            start.elapsed().as_millis()
        );
        for (t_i, target) in attractors.clone().iter().enumerate() {
            for (s_i, source) in attractors.clone().iter().enumerate() {
                if s_i == t_i {
                    continue
                }
                let start = Instant::now();
                let basin = perturbations.strong_basin(&target);
               
                println!("Elapsed: {}ms", start.elapsed().as_millis());
            }
        }
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

fn main_one_step(models: Vec<&str>, suffixes: Vec<&str>) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ONE STEP CONTROL");
    for m in models.clone().iter() {
        for s in suffixes.clone().iter() {
            let model_string: &str =
                &std::fs::read_to_string(format!("models/{}_{}.aeon", m, s)).unwrap();
            let model = BooleanNetwork::try_from(model_string).unwrap();
            let perturbations = PerturbationGraph::new(&model);
            println!("========= {}(v{})(p{}) =========", m, model.num_vars(), perturbations.as_original().unit_colors().approx_cardinality());
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
                    let start = Instant::now();
                    let control = perturbations.one_step_control(&source, &target);
                    println!(
                        "Control from Attractor {:?} (source) to Attractor {:?} (target) works for {} color(s), jumping through {} vertices.",
                        s_i,
                        t_i,
                        control.controllable_colors_cardinality(),
                        control.jump_vertices()
                    );
                    println!("Elapsed: {}ms", start.elapsed().as_millis());
                }
            }
        }
    }
}

fn main_permanent(models: Vec<&str>, suffixes: Vec<&str>) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  PERMANENT CONTROL");
    for m in models.clone().iter() {
        for s in suffixes.clone().iter() {
            let model_string: &str =
                &std::fs::read_to_string(format!("models/{}_{}.aeon", m, s)).unwrap();
            let model = BooleanNetwork::try_from(model_string).unwrap();
            let perturbations = PerturbationGraph::new(&model);
            println!("========= {}(v{})(p{}) =========", m, model.num_vars(), perturbations.as_original().unit_colors().approx_cardinality());
            let start = Instant::now();
            let attractors = find_witness_attractors(m);
            println!(
                "Attractors ready in {}ms, starting control...",
                start.elapsed().as_millis()
            );
            for (t_i, target) in attractors.clone().iter().enumerate() {
                for (s_i, source) in attractors.clone().iter().enumerate() {
                    if s_i == t_i {
                        continue
                    }
                    let start = Instant::now();
                    let control = perturbations.permanent_control(&source, &target);
                    println!(
                        "Control from Attractor {:?} (source) to Attractor {:?} (target) works for {} color(s), jumping through {} vertices.",
                        s_i,
                        t_i,
                        control.controllable_colors_cardinality(),
                        control.jump_vertices()
                    );
                    println!("Elapsed: {}ms", start.elapsed().as_millis());
                }
            }
        }
    }
}

fn main_temporary(models: Vec<&str>, suffixes: Vec<&str>) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  TEMPORARY CONTROL");
    for m in models.clone().iter() {
        for s in suffixes.clone().iter() {
            let model_string: &str =
                &std::fs::read_to_string(format!("models/{}_{}.aeon", m, s)).unwrap();
            let model = BooleanNetwork::try_from(model_string).unwrap();
            let perturbations = PerturbationGraph::new(&model);
            println!("========= {}(v{})(p{}) =========", m, model.num_vars(), perturbations.as_original().unit_colors().approx_cardinality());
            let start = Instant::now();
            let attractors = find_witness_attractors(m);
            println!(
                "Attractors ready in {}ms, starting control...",
                start.elapsed().as_millis()
            );
            for (t_i, target) in attractors.clone().iter().enumerate() {
                for (s_i, source) in attractors.clone().iter().enumerate() {
                    if s_i == t_i {
                        continue
                    }
                    let start = Instant::now();
                    let control = perturbations.temporary_control(&source, &target);
                    println!(
                        "Control from Attractor {:?} (source) to Attractor {:?} (target) works for {} color(s), jumping through {} vertices.",
                        s_i,
                        t_i,
                        control.controllable_colors_cardinality(),
                        control.jump_vertices()
                    );
                    println!("Elapsed: {}ms", start.elapsed().as_millis());
                }
            }
        }
    }
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