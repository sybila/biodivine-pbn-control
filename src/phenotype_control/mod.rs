use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::symbolic_async_graph::GraphColoredVertices;
use biodivine_lib_param_bn::VariableId;

pub mod _impl_phenotype_control_map;
pub mod _impl_phenotype_permanent_control;
mod _symbolic_utils;

/// A mapping between admissible perturbations and colors for which the perturbation controls
/// the network.
///
/// The map holds a reference to the `PerturbationGraph` from which it was created and thus
/// cannot outlive the graph.
#[derive(Clone)]
pub struct PhenotypeControlMap {
    perturbation_variables: Vec<VariableId>,
    context: PerturbationGraph,
    perturbation_set: GraphColoredVertices,
}

#[derive(Clone)]
pub enum PhenotypeOscillationType {
    Forbidden,
    Allowed,
    Required,
}
