#![allow(dead_code)]    // In this file, we want to allow unused functions.

use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_pbn_control::perturbation::PerturbationGraph;
use std::convert::TryFrom;
use std::time::Instant;

const models1: [&str; 5] = ["myeloid", "cardiac", "erbb", "tumour", "mapk"];
const suffixes2: [&str; 7] = ["witness", "1unknown", "2unknown", "3unknown", "4unknown", "5unknown", "6unknown"];


fn main() {
    main_one_step_old_benchmark();

    main_one_step(models1.to_vec(), vec!["4unknown"]);
    main_permanent(models1.to_vec(), vec!["4unknown"]);
    main_temporary(models1.to_vec(), vec!["4unknown"]);

    main_one_step(vec!["tumour"],suffixes2.to_vec());
    main_permanent(vec!["tumour"],suffixes2.to_vec());
    main_temporary(vec!["tumour"],suffixes2.to_vec());
}

fn main_one_step_old_benchmark() {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ONE STEP CONTROL - OLD BENCHMARK");

    for m in ["myeloid, cell_fate"].iter() {
        let suffixes: [&str; 4];
        let suffixes = if m == "cell_fate" {["7_stable_attractors", "7_stable_attractors_2params", "7_stable_attractors_4params"].to_vec()} else {["witness", "4_params", "8_params", "11params"].to_vec()};
        
        for s in suffixes.clone().iter() {
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