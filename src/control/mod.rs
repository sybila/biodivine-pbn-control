use crate::perturbation::PerturbationGraph;
use biodivine_lib_param_bn::symbolic_async_graph::GraphColoredVertices;

pub mod _impl_one_step_control;
pub mod _impl_permanent_control;
pub mod _impl_temporary_control;

mod _impl_control_map;

/// A mapping between admissible perturbations and colors for which the perturbation controls
/// the network.
///
/// The map holds a reference to the `PerturbationGraph` from which it was created and thus
/// cannot outlive the graph.
#[derive(Clone)]
pub struct ControlMap<'a> {
    context: &'a PerturbationGraph,
    perturbation_set: GraphColoredVertices,
}
