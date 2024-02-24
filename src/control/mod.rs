use std::collections::HashMap;
use biodivine_lib_bdd::Bdd;
use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::symbolic_async_graph::{GraphColoredVertices, GraphColors};
use biodivine_lib_param_bn::VariableId;

pub mod _impl_one_step_control;
pub mod _impl_permanent_control;
pub mod _impl_temporary_control;

pub mod _impl_phenotype_control_map;
pub mod _impl_phenotype_permanent_control;
pub mod _impl_attractor_control_map;

pub mod _symbolic_utils;

/// A mapping between admissible perturbations and colors for which the perturbation controls
/// the network.
///
/// The map holds a reference to the `PerturbationGraph` from which it was created and thus
/// cannot outlive the graph.
#[derive(Clone)]
pub struct AttractorControlMap {
    perturbation_variables: Vec<VariableId>,
    context: PerturbationGraph,
    perturbation_set: GraphColoredVertices,
}

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

pub trait ControlMap {
    fn new(context: PerturbationGraph, perturbation_set: GraphColoredVertices) -> Self;
    fn as_bdd(&self) -> &Bdd;
    fn as_colored_vertices(&self) -> &GraphColoredVertices;
    fn working_perturbations(
        &self,
        min_robustness: f64,
        verbose: bool,
    ) -> Vec<(HashMap<String, bool>, GraphColors)>;
    fn perturbation_working_colors(&self, perturbation: &HashMap<String, bool>) -> GraphColors;
}
