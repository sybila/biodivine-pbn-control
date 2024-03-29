use std::collections::HashMap;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_pbn_control::aeon::phentoype::resolve_var_id;
use biodivine_pbn_control::control::ControlMap;
use biodivine_pbn_control::control::PhenotypeOscillationType::Forbidden;
use biodivine_pbn_control::perturbation::PerturbationGraph;

fn main() {
    let model_path = "./models/myeloid_witness.aeon";
    let model_string = std::fs::read_to_string(model_path).unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let graph = SymbolicAsyncGraph::new(&bn).unwrap();
    let perturbation_graph = PerturbationGraph::new(&bn);

    let attractors = biodivine_pbn_control::aeon::attractors::compute(&perturbation_graph.as_original(), false);
    let mut vertices = Vec::new();
    for a in attractors {
        vertices.push(a.pick_vertex().clone().vertices());
    }

    let cjun = resolve_var_id(&graph, "cJun").unwrap();
    let monocyte = perturbation_graph.mk_unit_colored_vertices().vertices().fix_network_variable(cjun, true);
    let monocyte_att = vertices.iter().filter(|a| a.intersect(&monocyte).approx_cardinality() > 0.0).next().unwrap();


    let eklf = resolve_var_id(&graph, "EKLF").unwrap();
    let erythrocyte = perturbation_graph.mk_unit_colored_vertices().vertices().fix_network_variable(eklf, true);
    let erythrocyte_att = vertices.iter().filter(|a| a.intersect(&erythrocyte).approx_cardinality() > 0.0).next().unwrap();


    let source = monocyte_att.clone().into_iter().next().unwrap();
    let target = erythrocyte_att.clone().into_iter().next().unwrap();
    let res = perturbation_graph.permanent_control(&source, &target, perturbation_graph.unit_colors(), true);

    println!("source: {:?}", source);
    println!("target: {:?}", target);

    println!("All control {:?}", res.working_perturbations(0.99, true, false));

    // {'Fli1': False, 'PU1': False, 'GATA1': True}

    println!("Empty control {:?}", res.perturbation_working_colors(&HashMap::new()).approx_cardinality());
    println!("Single control {:?}", res.perturbation_working_colors(&HashMap::from_iter(vec![("Fli1".to_string(), false), ("PU1".to_string(), false), ("GATA1".to_string(), true)].into_iter())).approx_cardinality());

    // perturbation_graph.ceiled_phenotype_permanent_control(megakaryocyte_att.clone(), 5, Forbidden, true, true);
}