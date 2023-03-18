use std::borrow::Borrow;
use std::collections::HashMap;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_pbn_control::aeon::phentoype::build_phenotype;
use biodivine_pbn_control::experiment_utils::{parse_experiment, run_control_experiment};
use biodivine_pbn_control::perturbation::PerturbationGraph;

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let model_path = &args[1];
    let max_control_size: i32 = args[2].to_string().parse::<i32>().unwrap();
    let model_string = std::fs::read_to_string(model_path).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::new(&bn);
    let phenotype_var = bn.get_variable_name(bn.variables().into_iter().next().unwrap());
    let phenotype = build_phenotype(perturbation_graph.as_perturbed(),
                                    HashMap::from([(phenotype_var.as_str(), true)]));
    let _result = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype, max_control_size);
}
