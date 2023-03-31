use std::borrow::Borrow;
use std::collections::HashMap;
use biodivine_lib_param_bn::{BooleanNetwork, VariableId};
use biodivine_pbn_control::aeon::phentoype::build_phenotype;
use biodivine_pbn_control::experiment_utils::{parse_experiment, run_control_experiment};
use biodivine_pbn_control::perturbation::PerturbationGraph;

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let model_path = &args[1];
    let max_control_size: i32 = args[2].parse::<i32>().unwrap();
    let max_control_vars: usize = args[3].parse::<usize>().unwrap();
    let model_string = std::fs::read_to_string(model_path).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let mut p_vars = Vec::new();
    let mut i = 0;
    for v in bn.variables() {
        if i == 0 {
            i += 1;
            continue;
        }
        if i > max_control_vars {
            break;
        }
        p_vars.push(v);
        i += 1;
    }
    let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, &p_vars);
    let phenotype_var = bn.get_variable_name(bn.variables().into_iter().next().unwrap());
    let phenotype = build_phenotype(perturbation_graph.as_perturbed(),
                                    HashMap::from([(phenotype_var.as_str(), true)]));
    let _result = PerturbationGraph::ceiled_phenotype_permanent_control(&perturbation_graph, phenotype, max_control_size, &p_vars);
}
