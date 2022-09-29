use std::convert::TryFrom;
use std::time::Instant;
use biodivine_lib_param_bn::async_graph::{AsyncGraph, DefaultEdgeParams};
use biodivine_lib_param_bn::bdd_params::BddParams;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_lib_std::IdState;
use biodivine_lib_std::param_graph::{EvolutionOperator, Graph, Params};
use biodivine_pbn_control::control::_algo_control_to_basin::{find_smallest_control_to_basin};
use biodivine_pbn_control::old_aeon_server::algo_reach::reach;
use biodivine_pbn_control::old_aeon_server::StateSet;
use biodivine_pbn_control::strong_basin::_algo_sb_fixed_point::find_strong_basin;

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let model_path = &args[1];
    let model_string = std::fs::read_to_string(model_path).unwrap();
    let (source, target, model) = parse_experiment(model_string.as_str());
    run_control_experiment(source, target, model);
}

pub fn run_control_experiment(source: IdState, target: IdState, model: BooleanNetwork) {
    println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ONE-STEP SEMI-SYMBOLIC CONTROL");
    let graph = &AsyncGraph::new(model).unwrap();

    // We are not printing model stats... it would be a pain to port this to the old code
    // reliably... you can use newer experiment code for that instead...

    let start_attractor = Instant::now();
    let attractor_colors = get_all_params_with_attractor(graph, target);
    println!("Attractor colors computed in {} ms.", start_attractor.elapsed().as_millis());

    let start = Instant::now();
    let basin = find_strong_basin(graph, target, attractor_colors.clone());
    println!("Strong basin has {} states.", basin.len());
    println!("Strong basin computation time (ms): {:?}", start.elapsed().as_millis());

    let controls = find_smallest_control_to_basin(&source, &basin);
    println!("Controls: {}", controls.len());
    /*for (s, p) in controls {
        println!("Control into {:?}, size {:?}, robustness {:?}.", s, control_dist(&s, &target), p.cardinality() / attractor_colors.cardinality());
    }*/

    let elapsed = start.elapsed();
    println!("Elapsed: {} ms", elapsed.as_millis());
    let elapsed = elapsed.as_millis() as f64 / 1000.0;
    println!("WARNING: This is not a real output of UNIX `time` utility. This is just a similar format for compatibility reasons.");
    print!("real {}\nuser ??\nsys ??\n", elapsed);
}

pub fn parse_experiment(file: &str) -> (IdState, IdState, BooleanNetwork) {
    let lines = file.lines().collect::<Vec<_>>();
    let source_line = &lines[0];
    let target_line = &lines[1];
    assert!(source_line.starts_with("#source:"));
    assert!(target_line.starts_with("#target:"));
    let source = string_to_state(&source_line[8..]);
    let target = string_to_state(&target_line[8..]);
    (source, target, BooleanNetwork::try_from(file).unwrap())
}

/// Convert a bit-vector string to an actual bit-vector.
pub fn string_to_state(string_vector: &str) -> IdState {
    assert!(string_vector.starts_with("["));
    assert!(string_vector.ends_with("]"));
    let string_vector = &string_vector[1..(string_vector.len() - 1)];
    let mut state = IdState::from(0);
    for (i, segment) in string_vector.split(",").enumerate() {
        //print!("{}", segment);
        if segment.trim() == "True" {
           state = state.flip_bit(i);
        } else if segment.trim() == "False" {
            // Do nothing. Already false.
        } else {
            panic!("Unexpected: {}", segment);
        }
    }
    //println!("\n{}", state);
    state
}

/// Compute all graph colors for which the given state is an attractor state.
pub fn get_all_params_with_attractor(graph: &AsyncGraph<DefaultEdgeParams>, state: IdState) -> BddParams {
    //println!("All colors: {}", graph.unit_params().cardinality());
    let mut pivot = StateSet::new(graph.num_states());
    pivot.union_key(state, graph.unit_params());
    let fwd = reach(&graph.fwd(), &pivot);
    let bwd = reach(&graph.bwd(), &pivot);
    let non_terminal = fwd.minus(&bwd);
    let non_terminal_colors = non_terminal.fold_union().unwrap_or(graph.empty_params().clone());
    graph.unit_params().minus(&non_terminal_colors)
}
