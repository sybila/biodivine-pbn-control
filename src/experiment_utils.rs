use crate::aeon::reachability::backward;
use crate::control::AttractorControlMap;
use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::biodivine_std::bitvector::ArrayBitVector;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::GraphColors;
use biodivine_lib_param_bn::BooleanNetwork;
use std::convert::TryFrom;
use std::time::Instant;

pub fn run_control_experiment<F>(
    source: ArrayBitVector,
    target: ArrayBitVector,
    model: BooleanNetwork,
    control_function: F,
    control_type: &str,
) where
    F: Fn(
        &PerturbationGraph,
        &ArrayBitVector,
        &ArrayBitVector,
        &GraphColors,
        bool
    ) -> AttractorControlMap,
{
    println!(
        ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> {} CONTROL",
        control_type
    );
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

        println!("Variables: {}", model_variables);
        println!("Uncertainty colors: {}", model_colors);
        println!("Perturbation colors: {}", perturbation_colors);
        println!("All colors: {}", all_colors);
        println!(
            "Unknown update functions: {}",
            model
                .variables()
                .filter(|it| { model.get_update_function(*it).is_none() })
                .count()
        );
        println!(
            "Lowest cardinality: {}",
            model
                .variables()
                .filter(|it| { model.get_update_function(*it).is_none() })
                .map(|it| { model.as_graph().regulators(it).len() })
                .min()
                .unwrap()
        );
        println!(
            "Highest cardinality: {}",
            model
                .variables()
                .filter(|it| { model.get_update_function(*it).is_none() })
                .map(|it| { model.as_graph().regulators(it).len() })
                .max()
                .unwrap()
        );
    }

    // Compute attractor states in the witness model.

    let start_attractor = Instant::now();
    let attractor_colors = get_all_params_with_attractor(&perturbation_graph, &target);
    println!(
        "Attractor colors: {}",
        attractor_colors.approx_cardinality()
    );
    println!(
        "Attractor colors computed in {} ms.",
        start_attractor.elapsed().as_millis()
    );
    let start = Instant::now();
    let control = control_function(&perturbation_graph, &source, &target, &attractor_colors, false);
    println!(
        "Control exists jumping through {} vertices.",
        // control.controllable_colors_cardinality(),
        control.jump_vertices()
    );
    let elapsed = start.elapsed();
    // println!("Elapsed: {} ms", elapsed.as_millis());
    let _elapsed = elapsed.as_millis() as f64 / 1000.0;
    // println!("WARNING: This is not a real output of UNIX `time` utility. This is just a similar format for compatibility reasons.");
    // print!("real {}\nuser ??\nsys ??\n", elapsed);
}

pub fn parse_experiment(file: &str) -> (ArrayBitVector, ArrayBitVector, BooleanNetwork) {
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
pub fn string_to_state(string_vector: &str) -> ArrayBitVector {
    assert!(string_vector.starts_with('['));
    assert!(string_vector.ends_with(']'));
    let string_vector = &string_vector[1..(string_vector.len() - 1)];
    let mut vector = Vec::new();
    for segment in string_vector.split(',') {
        if segment.trim() == "True" {
            vector.push(true);
        } else if segment.trim() == "False" {
            vector.push(false);
        } else {
            panic!("Unexpected: {}", segment);
        }
    }
    ArrayBitVector::from(vector)
}

/// Compute all graph colors for which the given state is an attractor state.
pub fn get_all_params_with_attractor(
    graph: &PerturbationGraph,
    state: &ArrayBitVector,
) -> GraphColors {
    let seed = graph.vertex(state);
    let bwd = backward(graph.as_original(), &seed, false);
    let mut attractor = seed;
    'forward: loop {
        if cfg!(feature = "print_progress") && attractor.as_bdd().size() > 100_000 {
            println!("FWD-attractor: {}", attractor.as_bdd().size());
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
