#![allow(dead_code)]    // In this file, we want to allow unused functions.

use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_pbn_control::perturbation::PerturbationGraph;
use std::convert::TryFrom;
use std::time::Instant;

fn main() {
    main_one_step("4unknown");
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
    for m in ["myeloid", "cardiac", "erbb", "tumour", "mapk", "hgf"].iter() {
        let model_string: &str =
            &std::fs::read_to_string(format!("models/{}_{}.aeon", m, suffix)).unwrap();
        let model = BooleanNetwork::try_from(model_string).unwrap();
        println!("========= {}({}) =========", m, model.num_vars());
        let perturbations = PerturbationGraph::new(&model);
        let start = Instant::now();
        let attractors = compute_attractor_pairs(&model);
        println!(
            "Attractors ready in {}ms, starting control...",
            start.elapsed().as_millis()
        );
        for (source, target) in attractors {
            let start = Instant::now();
            let control = perturbations.one_step_control(&source, &target);
            println!(
                "Control from {:?} to {:?} works for {} color(s), jumping through {} vertices.",
                source,
                target,
                control.controllable_colors_cardinality(),
                control.jump_vertices()
            );
            println!("Elapsed: {}ms", start.elapsed().as_millis());
        }
    }
}

fn main_permanent(suffix: &str) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PERMANENT CONTROL ({})", suffix);
    for m in ["myeloid", "cardiac", "erbb", "tumour", "mapk", "hgf"].iter() {
        let model_string: &str =
            &std::fs::read_to_string(format!("models/{}_{}.aeon", m, suffix)).unwrap();
        let model = BooleanNetwork::try_from(model_string).unwrap();
        println!("========= {}({}) =========", m, model.num_vars());
        let perturbations = PerturbationGraph::new(&model);
        let start = Instant::now();
        let attractors = compute_attractor_pairs(&model);
        println!(
            "Attractors ready in {}ms, starting control...",
            start.elapsed().as_millis()
        );
        for (source, target) in attractors {
            let start = Instant::now();
            let control = perturbations.permanent_control(&source, &target);
            println!(
                "Control from {:?} to {:?} works for {} color(s), jumping through {} vertices.",
                source,
                target,
                control.controllable_colors_cardinality(),
                control.jump_vertices()
            );
            println!("Elapsed: {}ms", start.elapsed().as_millis());
        }
    }
}

fn main_temporary(suffix: &str) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TEMPORARY CONTROL ({})", suffix);
    for m in ["myeloid", "cardiac", "erbb", "tumour", "mapk", "hgf"].iter() {
        let model_string: &str =
            &std::fs::read_to_string(format!("models/{}_{}.aeon", m, suffix)).unwrap();
        let model = BooleanNetwork::try_from(model_string).unwrap();
        println!("========= {}({}) =========", m, model.num_vars());
        let perturbations = PerturbationGraph::new(&model);
        let start = Instant::now();
        let attractors = compute_attractor_pairs(&model);
        println!(
            "Attractors ready in {}ms, starting control...",
            start.elapsed().as_millis()
        );
        for (source, target) in attractors {
            let start = Instant::now();
            let control = perturbations.temporary_control(&source, &target);
            println!(
                "Control from {:?} to {:?} works for {} color(s), jumping through {} vertices.",
                source,
                target,
                control.controllable_colors_cardinality(),
                control.jump_vertices()
            );
            println!("Elapsed: {}ms", start.elapsed().as_millis());
        }
    }
}

/*
fn main() {
    for m in ["myeloid", "cardiac", "erbb", "tumour", "mapk", "hgf"].iter() {
        let witness_str: &str = &fs::read_to_string(format!("models/{}_witness.aeon", m)).unwrap();
        let witness_model = BooleanNetwork::try_from(witness_str).unwrap();
        let witness_graph = &SymbolicAsyncGraph::new(witness_model.clone()).unwrap();

        println!("Model: {}", m);
        println!("# variables: {}", witness_model.clone().variables().len());

        let attractors = &find_attractors(witness_graph);
        println!("# attractors: {}", attractors.clone().len());
        let mut vertices = Vec::new();
        for a in attractors {
            vertices.push(a.pick_vertex());
        }

        let param_str: &str = &fs::read_to_string(format!("models/{}_4unknown.aeon", m)).unwrap();
        let param_model = BooleanNetwork::try_from(param_str).unwrap();
        let param_graph = &SymbolicAsyncGraph::new(param_model.clone()).unwrap();
        println!("# params: {}", param_graph.unit_colors().approx_cardinality());

        println!("WITNESS TEMPORARY");

        for (t_i, t) in vertices.clone().iter().enumerate() {
            if t_i > 3 {
                break;
            }
            let begin = Instant::now();
            for (s_i, s) in vertices.clone().iter().enumerate() {
                if s_i > 3 {
                    break;
                }
                if t_i == s_i {
                    continue;
                }
                for source in s.vertices().materialize().iter() {
                    for target in t.vertices().materialize().iter() {
                        let _ = TemporaryControl::new(BooleanNetwork::try_from(witness_str).unwrap(), &source, &target);
                    }
                }
            }
            let elapsed = begin.elapsed().as_millis();
            println!("To attractor {} elapsed {} ms", t_i, elapsed);
        }

        println!("WITNESS PERMANENT");

        for (t_i, t) in vertices.clone().iter().enumerate() {
            if t_i > 3 {
                break;
            }
            let begin = Instant::now();
            for (s_i, s) in vertices.clone().iter().enumerate() {
                if s_i > 3 {
                    break;
                }
                if t_i == s_i {
                    continue;
                }
                for source in s.vertices().materialize().iter() {
                    for target in t.vertices().materialize().iter() {
                        let _ = PermanentControl::new(BooleanNetwork::try_from(witness_str).unwrap(), &source, &target);
                    }
                }
            }
            let elapsed = begin.elapsed().as_millis();
            println!("To attractor {} elapsed {} ms", t_i, elapsed);
        }

        println!("4 PARAMS TEMPORARY");
        for (t_i, t) in vertices.clone().iter().enumerate() {
            if t_i > 3 {
                break;
            }
            let begin = Instant::now();
            for (s_i, s) in vertices.clone().iter().enumerate() {
                if s_i > 3 {
                    break;
                }
                if t_i == s_i {
                    continue;
                }
                for source in s.vertices().materialize().iter() {
                    for target in t.vertices().materialize().iter() {
                        let _ = TemporaryControl::new(BooleanNetwork::try_from(param_str).unwrap(), &source, &target);
                    }
                }
            }
            let elapsed = begin.elapsed().as_millis();
            println!("To attractor {} elapsed {} ms", t_i, elapsed);
        }

        println!("4 PARAMS PERMANENT");
        for (t_i, t) in vertices.clone().iter().enumerate() {
            if t_i > 3 {
                break;
            }
            let begin = Instant::now();
            for (s_i, s) in vertices.clone().iter().enumerate() {
                if s_i > 3 {
                    break;
                }
                if t_i == s_i {
                    continue;
                }
                for source in s.vertices().materialize().iter() {
                    for target in t.vertices().materialize().iter() {
                        let _ = PermanentControl::new(BooleanNetwork::try_from(param_str).unwrap(), &source, &target);
                    }
                }
            }
            let elapsed = begin.elapsed().as_millis();
            println!("To attractor {} elapsed {} ms", t_i, elapsed);
        }
    }
}
*/
