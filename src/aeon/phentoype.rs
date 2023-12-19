use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphVertices, SymbolicAsyncGraph};
use biodivine_lib_param_bn::VariableId;
use std::collections::HashMap;

// Obtain subspace from human-readable variable names and their expected values (as conjunction)
pub fn build_phenotype(
    graph: &SymbolicAsyncGraph,
    phenotype: HashMap<&str, bool>,
) -> GraphVertices {
    let mut result = graph.unit_colored_vertices().clone();
    for (var, value) in phenotype {
        let var_id = resolve_var_id(graph, var).unwrap();
        let subspace = graph.fix_network_variable(var_id, value);
        result = result.intersect(&subspace);
    }

    result.vertices().clone()
}

pub fn resolve_var_id(graph: &SymbolicAsyncGraph, var: &str) -> Option<VariableId> {
    let mut v_name = String::new();
    for v in graph.variables() {
        // Resolve variable name
        v_name = graph.get_variable_name(v);
        if var == v_name {
            return Some(v);
        }
    }
    assert_eq!(var, v_name, "Unknown variable");
    None
}
