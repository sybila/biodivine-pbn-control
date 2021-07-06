use crate::perturbation::PerturbationGraph;
use crate::perturbation::_algo_network_transformations::{
    make_original_network, make_perturbed_network, normalize_network,
};
use biodivine_lib_param_bn::biodivine_std::bitvector::{ArrayBitVector, BitVector};
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use biodivine_lib_param_bn::symbolic_async_graph::{
    GraphColoredVertices, GraphColors, SymbolicAsyncGraph, SymbolicContext,
};
use biodivine_lib_param_bn::{BooleanNetwork, ParameterId, VariableId, VariableIdIterator};
use std::collections::HashMap;

impl PerturbationGraph {

    pub fn new(network: &BooleanNetwork) -> PerturbationGraph {
        PerturbationGraph::with_restricted_variables(network, &network.variables().collect::<Vec<_>>())
    }

    /// Create a new perturbation graph for a given Boolean network.
    pub fn with_restricted_variables(network: &BooleanNetwork, perturb: &[VariableId]) -> PerturbationGraph {
        let normalized = normalize_network(network);

        let mut original_parameters = HashMap::new();
        let mut perturbed_parameters = HashMap::new();

        let original = make_original_network(&normalized, &mut original_parameters, perturb);
        let perturbed = make_perturbed_network(&normalized, &mut perturbed_parameters, perturb);

        assert_eq!(original_parameters, perturbed_parameters);

        PerturbationGraph {
            original_graph: SymbolicAsyncGraph::new(original).unwrap(),
            perturbed_graph: SymbolicAsyncGraph::new(perturbed).unwrap(),
            perturbation_parameters: original_parameters,
        }
    }

    pub fn as_original(&self) -> &SymbolicAsyncGraph {
        &self.original_graph
    }

    pub fn as_perturbed(&self) -> &SymbolicAsyncGraph {
        &self.perturbed_graph
    }

    pub fn as_symbolic_context(&self) -> &SymbolicContext {
        self.original_graph.symbolic_context()
    }

    pub fn variables(&self) -> VariableIdIterator {
        self.original_graph.as_network().variables()
    }

    pub fn get_perturbation_parameter(&self, variable: VariableId) -> Option<ParameterId> {
        self.perturbation_parameters.get(&variable).cloned()
    }

    /*
        WARNING: The unit color set in the perturbed graph is not correct! It enforces
        observability and for a regulation to be observable, the variable cannot be perturbed.
        So the unit set only contains one non-perturbed parametrisation.

        Consequently, we use the original graph where possible.
    */

    pub fn empty_colors(&self) -> &GraphColors {
        self.original_graph.empty_colors()
    }

    pub fn mk_empty_colors(&self) -> GraphColors {
        self.original_graph.mk_empty_colors()
    }

    pub fn empty_colored_vertices(&self) -> &GraphColoredVertices {
        self.original_graph.empty_vertices()
    }

    pub fn mk_empty_colored_vertices(&self) -> GraphColoredVertices {
        self.original_graph.mk_empty_vertices()
    }

    pub fn unit_colors(&self) -> &GraphColors {
        self.original_graph.unit_colors()
    }

    pub fn mk_unit_colors(&self) -> GraphColors {
        self.original_graph.mk_unit_colors()
    }

    pub fn unit_colored_vertices(&self) -> &GraphColoredVertices {
        self.original_graph.unit_colored_vertices()
    }

    pub fn mk_unit_colored_vertices(&self) -> GraphColoredVertices {
        self.original_graph.mk_unit_colored_vertices()
    }

    pub fn vertex(&self, state: &ArrayBitVector) -> GraphColoredVertices {
        self.original_graph.vertex(state)
    }

    pub fn fix_variable(&self, variable: VariableId, value: bool) -> GraphColoredVertices {
        self.original_graph.fix_network_variable(variable, value)
    }

    pub fn strong_basin(&self, target: &ArrayBitVector) -> GraphColoredVertices {
        let target_set = self.vertex(&target);
        let weak_basin = crate::aeon::reachability::backward(self.as_original(), &target_set);
        let strong_basin =
            crate::aeon::reachability::forward_closed(self.as_original(), &weak_basin);
        return strong_basin;
    }

    /// Return a subset of vertices and colors where the variable is perturbed to the given value.
    ///
    /// If no value is given, return vertices and colors where the variable is perturbed.
    ///
    /// If the value cannot be perturbed, return empty set.
    pub fn fix_perturbation(
        &self,
        variable: VariableId,
        value: Option<bool>,
    ) -> GraphColoredVertices {
        if let Some(is_perturbed) = self.perturbation_parameters.get(&variable) {
            let states = if let Some(value) = value {
                self.fix_variable(variable, value)
            } else {
                self.mk_unit_colored_vertices()
            };
            let bdd_is_perturbed = self
                .as_symbolic_context()
                .mk_uninterpreted_function_is_true(*is_perturbed, &[]);
            let colors_is_perturbed = self.unit_colors().copy(bdd_is_perturbed);
            states.intersect_colors(&colors_is_perturbed)
        } else {
            self.mk_empty_colored_vertices()
        }
    }

    /// Return a subset of colors for which the given `variable` is not perturbed.
    pub fn not_perturbed(&self, variable: VariableId) -> GraphColors {
        if let Some(is_perturbed) = self.perturbation_parameters.get(&variable) {
            let bdd_is_perturbed = self
                .as_symbolic_context()
                .mk_uninterpreted_function_is_true(*is_perturbed, &[]);
            self.unit_colors().copy(bdd_is_perturbed.not())
        } else {
            self.mk_unit_colors()
        }
    }

    /// Compute the subset of `target` to which a jump from `source` is possible using a perturbation.
    pub fn post_perturbation(
        &self,
        source: &ArrayBitVector,
        target: &GraphColoredVertices,
    ) -> GraphColoredVertices {
        /*
           The algorithm is based on iteratively removing the "negative results", i.e. the states
           which are not reachable by a perturbation. For this, we observe that a state is not
           reachable via a perturbation if there is a variable that is not perturbed and yet
           it is different from its value in source (since perturbation will not touch this value,
           it will not change and hence has to be the same).

           We thus specify a set of such unreachable states (per variable) and remove these states
           from the `target` set, together with all colours. The state-colour pairs that remain
           then represent the perturbations that can be applied to jump into target.

           (Alternatively, we could compute a "positive result", i.e.
               !is_perturbed[var] => source[var] = target[var]
           however, that would need some more custom BDD fiddling and is therefore easier to
           avoid at the moment)
        */
        let mut result = target.clone();
        for v in self.variables() {
            let value_in_source = source.get(usize::from(v));
            let mismatched_states = self.fix_variable(v, !value_in_source);
            let v_not_perturbed = self.not_perturbed(v);
            result = result.minus(&mismatched_states.intersect_colors(&v_not_perturbed));
        }

        result
    }
}
