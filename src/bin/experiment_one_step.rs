use biodivine_pbn_control::experiment_utils::{parse_experiment, run_control_experiment};
use biodivine_pbn_control::perturbation::PerturbationGraph;

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let model_path = &args[1];
    let model_string = std::fs::read_to_string(model_path).unwrap();
    let (source, target, model) = parse_experiment(model_string.as_str());
    run_control_experiment(
        source,
        target,
        model,
        PerturbationGraph::one_step_control,
        "one-step",
    );
}
