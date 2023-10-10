use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_lib_param_bn::{ParameterId, VariableId};
use std::collections::HashMap;

/// Procedures for transforming Boolean networks so that they conform to our encoding.
///
/// In particular, here we have functions for normalizing the network before it is used in
/// a perturbed graph and then creating the "original" and "perturbed" network from the normalized
/// result.
mod _algo_network_transformations;
mod _impl_perturbation_graph;

/// Perturbation graph allows representing the *original* `SymbolicAsyncGraph` as well as
/// the async graph with perturbations encoded in parameters. Currently, we are "hacking"
/// the original symbolic implementation to do that. Please read all the documentation to
/// know what is actually allowed.
///
/// Rule of thumb: You can access the "inner" original and perturbed graph. Do not do that
/// unless you really have to! There is a lot of "safer" APIs on the `PerturbationGraph` itself.
/// Use those when possible. The only legitimate reason for accessing the inner graphs should
/// be to pass them into reachability and strong basin algorithms so that they can be used
/// for computing pre/post.
#[derive(Clone)]
pub struct PerturbationGraph {
    /// "Normal" unperturbed graph, but with the same encoding as the perturbed graph.
    original_graph: SymbolicAsyncGraph,
    /// Perturbed graph where each edge is also labelled with perturbations that enable it.
    perturbed_graph: SymbolicAsyncGraph,
    /// Variables which
    perturbable_vars: Vec<VariableId>,
    /// Obtain parameters that decide whether a specific variable is perturbed.
    perturbation_parameters: HashMap<VariableId, ParameterId>,
}
