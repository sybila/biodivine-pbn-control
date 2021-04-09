#![allow(dead_code)]    // In this file, we want to allow unused functions.

use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_pbn_control::perturbation::PerturbationGraph;
use std::convert::TryFrom;
use std::time::Instant;

//const models: [&str; 6] = ["myeloid", "cardiac", "erbb", "tumour", "mapk", "hgf"];
const models: [&str; 5] = ["myeloid", "cardiac", "erbb", "tumour", "mapk"];

fn main() {
    main_one_step("4unknown");
    main_permanent("4unknown");
    main_temporary("4unknown");
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

fn main_one_step(suffix: &str) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ONE STEP CONTROL ({})", suffix);
    for m in models.clone().iter() {
        let model_string: &str =
            &std::fs::read_to_string(format!("models/{}_{}.aeon", m, suffix)).unwrap();
        let model = BooleanNetwork::try_from(model_string).unwrap();
        let perturbations = PerturbationGraph::new(&model);
        println!("========= {}(v{})(p{}) =========", m, model.num_vars(), perturbations.as_original().unit_colors().approx_cardinality());
        let start = Instant::now();
        let attractors = find_withness_attractors(&m);
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

fn main_permanent(suffix: &str) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PERMANENT CONTROL ({})", suffix);
    for m in models.clone().iter() {
        let model_string: &str =
            &std::fs::read_to_string(format!("models/{}_{}.aeon", m, suffix)).unwrap();
        let model = BooleanNetwork::try_from(model_string).unwrap();
        let perturbations = PerturbationGraph::new(&model);
        println!("========= {}(v{})(p{}) =========", m, model.num_vars(), perturbations.as_original().unit_colors().approx_cardinality());
        let start = Instant::now();
        let attractors = find_withness_attractors(&m);
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

fn main_temporary(suffix: &str) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TEMPORARY CONTROL ({})", suffix);
    for m in models.clone().iter() {
        let model_string: &str =
            &std::fs::read_to_string(format!("models/{}_{}.aeon", m, suffix)).unwrap();
        let model = BooleanNetwork::try_from(model_string).unwrap();
        let perturbations = PerturbationGraph::new(&model);
        println!("========= {}(v{})(p{}) =========", m, model.num_vars(), perturbations.as_original().unit_colors().approx_cardinality());
        let start = Instant::now();
        let attractors = find_withness_attractors(&m);
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

fn find_withness_attractors(m: &str) -> Vec<ArrayBitVector> {
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