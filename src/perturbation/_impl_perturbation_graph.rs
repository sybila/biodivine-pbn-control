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
        PerturbationGraph::with_restricted_variables(
            network,
            network.variables().collect::<Vec<_>>(),
        )
    }

    /// Create a new perturbation graph for a given Boolean network.
    pub fn with_restricted_variables(
        network: &BooleanNetwork,
        perturb: Vec<VariableId>,
    ) -> PerturbationGraph {
        // A network with all implicit parameters substituted for explicit parameters,
        // and with all observability constraints removed.
        let normalized = normalize_network(network);

        let mut original_parameters = HashMap::new();
        let mut perturbed_parameters = HashMap::new();

        // Two variants of the normalized network with unperturbed and perturbed dynamics.
        let original =
            make_original_network(&normalized, &mut original_parameters, perturb.clone());
        let perturbed =
            make_perturbed_network(&normalized, &mut perturbed_parameters, perturb.clone());

        assert_eq!(original_parameters, perturbed_parameters);

        // A graph based on the initial, non-normalized network.
        let basic_graph = SymbolicAsyncGraph::new(network).unwrap();

        // A symbolic context that should be valid for both the original and perturbed network
        // (these only differ in the usage of perturbation parameters in update functions).
        let perturbed_symbolic_context = SymbolicContext::new(&original).unwrap();

        // Transfer the BDD unit set from the non-normalized network to the "perturbed" context.
        // This should work, because the names of the implicit uninterpreted functions are the
        // same, and otherwise we have only added new parameters.
        //
        // This effectively transfers the static constraints from the original network into the
        // perturbed network, which we must do because we cannot apply them to implicit parameters
        // directly, and there are other problems with observability anyway.
        let perturbed_unit = perturbed_symbolic_context
            .transfer_from(basic_graph.unit_colored_vertices().as_bdd(), basic_graph.symbolic_context())
            .unwrap();

        PerturbationGraph {
            non_perturbable_graph: basic_graph,
            original_graph: SymbolicAsyncGraph::with_custom_context(
                &original,
                perturbed_symbolic_context.clone(),
                perturbed_unit.clone()
            ).unwrap(),
            perturbed_graph: SymbolicAsyncGraph::with_custom_context(
                &perturbed,
                perturbed_symbolic_context,
                perturbed_unit
            ).unwrap(),
            perturbable_vars: perturb.clone(),
            perturbation_parameters: original_parameters,
        }
    }

    pub fn as_non_perturbable(&self) -> &SymbolicAsyncGraph {
        &self.non_perturbable_graph
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
        self.original_graph.as_network().unwrap().variables()
    }

    pub fn perturbable_variables(&self) -> &Vec<VariableId> {
        &self.perturbable_vars
    }

    pub fn get_perturbation_parameter(&self, variable: VariableId) -> Option<ParameterId> {
        self.perturbation_parameters.get(&variable).cloned()
    }

    pub fn num_perturbation_parameters(&self) -> usize {
        self.perturbation_parameters.len()
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
        self.original_graph.empty_colored_vertices()
    }

    pub fn mk_empty_colored_vertices(&self) -> GraphColoredVertices {
        self.original_graph.mk_empty_colored_vertices()
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
        let target_set = self.vertex(target);
        let weak_basin = crate::aeon::reachability::backward(self.as_original(), &target_set, false);
        let strong_basin =
            crate::aeon::reachability::forward_closed(self.as_original(), &weak_basin, false);
        strong_basin
    }

    /// Return a subset of vertices and colors where the variable is perturbed to the given value.
    ///
    /// If no value is given, return vertices and colors where the variable is perturbed.
    ///
    /// If the value cannot be perturbed, return empty set.
    pub fn fix_perturbation(
        &self,
        variable: VariableId,
        value: Option<&bool>,
    ) -> GraphColoredVertices {
        if let Some(is_perturbed) = self.perturbation_parameters.get(&variable) {
            let states = if let Some(value) = value {
                self.fix_variable(variable, value.clone())
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

    pub fn build_colors_with_values(
        &self,
        bn: &BooleanNetwork,
        values: HashMap<String, bool>,
    ) -> GraphColors {
        let mut colors = self.as_original().mk_empty_colors();
        for v in bn.variables() {
            let function = bn.get_update_function(v);
            if !function.is_none() {
                continue;
            }
            println!("{:?}", bn.get_variable_name(v));

            let v_name = bn.get_variable_name(v);
            let value;
            if values.contains_key(v_name) {
                value = values.get(v_name).unwrap().clone();
            } else {
                value = false;
            }

            // let variable_parameter = bn.find_parameter((param_prefix + v_name).as_str()).unwrap();
            // assert_eq!(bn.get_parameter(variable_parameter).get_name().as_str(), format!("param_{}", v.as_str()));

            let colors_v_set;
            if value {
                let bdd = self
                    .as_symbolic_context()
                    .mk_implicit_function_is_true(v, &[]);
                colors_v_set = self.unit_colors().copy(bdd);
            } else {
                let bdd = self
                    .as_symbolic_context()
                    .mk_implicit_function_is_true(v, &[]);
                colors_v_set = self
                    .as_original()
                    .mk_unit_colors()
                    .minus(&self.unit_colors().copy(bdd));
            }

            colors = colors.minus(&colors_v_set);
        }
        colors
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

    pub fn create_perturbation_colors(
        &self,
        perturbation_size: usize,
        verbose: bool,
    ) -> GraphColors {
        // A map which gives us the symbolic variable of the perturbation parameter.
        let perturbation_bbd_vars_mapping =
            self.get_perturbation_bdd_mapping(self.perturbable_variables());
        let bdd_vars = self.as_symbolic_context().bdd_variable_set();
        // The list of symbolic variables of perturbation parameters.
        let perturbation_bdd_vars = Self::get_perturbation_bdd_vars(&perturbation_bbd_vars_mapping);

        let admissible_perturbations =
            crate::control::_symbolic_utils::mk_bdd_of_bound(bdd_vars, &perturbation_bdd_vars, perturbation_size);
        {
            let factor =
                2.0f64.powi(bdd_vars.num_vars() as i32 - perturbation_bdd_vars.len() as i32);
            if verbose {
                println!(
                    "[{}] >> Admissible fixed(Q) sets: {}",
                    perturbation_size,
                    admissible_perturbations.cardinality() / factor
                );
            }
        }
        let admissible_perturbations = self.empty_colors().copy(admissible_perturbations);
        admissible_perturbations
    }
}

#[cfg(test)]
mod tests {
    use biodivine_lib_param_bn::BooleanNetwork;
    use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
    use crate::perturbation::PerturbationGraph;

    #[test]
    pub fn test_unit_set_compatibility() {
        let network = BooleanNetwork::try_from(r#"
            # This network has:
            #   - implicit and explicit parameters
            #   - non-essential regulations
            #   - monotonicity constraints
            a -> b
            a -|? c
            b -> c
            c ->? a
            c -| b
            $c: f(a) | b
        "#).unwrap();

        // These two do not share a symbolic representation, but we should be able to transfer
        // colors between them, as long as the perturbation parameters are unconstrained.
        let stg = SymbolicAsyncGraph::new(&network).unwrap();
        let p_stg = PerturbationGraph::new(&network);

        let transferred = stg.transfer_from(&p_stg.mk_unit_colored_vertices(), p_stg.as_original());
        assert_eq!(stg.mk_unit_colored_vertices(), transferred.unwrap());

        let transferred = stg.transfer_from(&p_stg.as_original().mk_unit_colored_vertices(), p_stg.as_original());
        assert_eq!(stg.mk_unit_colored_vertices(), transferred.unwrap());

        let transferred = stg.transfer_from(&p_stg.as_perturbed().mk_unit_colored_vertices(), p_stg.as_perturbed());
        assert_eq!(stg.mk_unit_colored_vertices(), transferred.unwrap());

        let transferred = stg.transfer_from(&p_stg.as_non_perturbable().mk_unit_colored_vertices(), p_stg.as_non_perturbable());
        assert_eq!(stg.mk_unit_colored_vertices(), transferred.unwrap());
    }

}