use std::collections::HashMap;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_pbn_control::aeon::phentoype::resolve_var_id;
use biodivine_pbn_control::control::ControlMap;
use biodivine_pbn_control::perturbation::PerturbationGraph;

fn main() {
    let model_path = "./models/myeloid_witness.aeon";
    let model_string = std::fs::read_to_string(model_path).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let graph = SymbolicAsyncGraph::new(&bn).unwrap();
    let perturbation_graph = PerturbationGraph::new(&bn);

    let attractors = biodivine_pbn_control::aeon::attractors::compute(&graph, false);
    let mut vertices = Vec::new();
    for a in attractors {
        vertices.push(a.pick_vertex().clone().vertices());
    }

    let eklf = resolve_var_id(&graph, "EKLF").unwrap();
    let erythrocyte = graph.mk_unit_vertices().fix_network_variable(eklf, true);
    let erythrocyte_att = vertices.iter().filter(|a| a.intersect(&erythrocyte).approx_cardinality() > 0.0).next().unwrap();

    let fli1 = resolve_var_id(&graph, "Fli1").unwrap();
    let megakaryocyte = graph.mk_unit_vertices().fix_network_variable(fli1, true);
    let megakaryocyte_att = vertices.iter().filter(|a| a.intersect(&megakaryocyte).approx_cardinality() > 0.0).next().unwrap();

    let source = erythrocyte_att.clone().into_iter().next().unwrap();
    let target = megakaryocyte_att.clone().into_iter().next().unwrap();
    let res = perturbation_graph.one_step_control(&source, &target, perturbation_graph.unit_colors(), true);

    println!("All control {:?}", res.working_perturbations(0.99, true, false));

    println!("Empty control {:?}", res.perturbation_working_colors(&HashMap::new()).approx_cardinality());
    println!("Single control {:?}", res.perturbation_working_colors(&HashMap::from_iter(vec![("Fli1".to_string(), true), ("EKLF".to_string(), false)].into_iter())).approx_cardinality());

}