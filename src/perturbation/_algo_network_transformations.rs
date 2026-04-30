use biodivine_lib_param_bn::{BooleanNetwork, FnUpdate, ParameterId, RegulatoryGraph, VariableId};
use std::collections::HashMap;
use std::convert::TryFrom;

/// Create a copy of the given `model`, but make every regulation non-observable, and add an
/// unspecified auto-regulation to each variable (that does not have an autoregulation already).
/// Additionally, convert every implicit update function to and explicit parameter.
///
/// This is a necessary pre-processing step before creating a `PerturbationGraph`.
pub fn normalize_network(network: &BooleanNetwork) -> BooleanNetwork {
    let mut network = network.clone();
    for var in network.variables() {
        if network.get_update_function(var).is_none() {
            // Create an explicit parameter to replace the implicit function.
            let regulators = network.regulators(var);
            let parameter = network
                .add_parameter(
                    format!("f_{}", network.get_variable_name(var)).as_str(),
                    u32::try_from(regulators.len()).unwrap(),
                )
                .unwrap();
            network
                .add_update_function(var, FnUpdate::mk_basic_param(parameter, &regulators))
                .unwrap();
        }
    }

    let network = network;

    // Copy graph, but with non-observable regulations.
    let mut result = RegulatoryGraph::new(
        network
            .variables()
            .map(|it| network.get_variable_name(it).clone())
            .collect(),
    );
    for regulation in network.as_graph().regulations() {
        // Make all self-regulations non-monotonic, since this can be influences
        // by the control transformation.
        let monotonicity = if regulation.get_regulator() == regulation.get_target() {
            None
        } else {
            regulation.get_monotonicity()
        };
        result
            .add_regulation(
                network.get_variable_name(regulation.get_regulator()),
                network.get_variable_name(regulation.get_target()),
                false,
                monotonicity,
            )
            .unwrap();
    }
    // Copy variables and regulations.
    for v in result.variables() {
        if result.find_regulation(v, v).is_none() {
            let name = result.get_variable_name(v).clone();
            result
                .add_regulation(name.as_str(), name.as_str(), false, None)
                .unwrap();
        }
    }

    let mut result = BooleanNetwork::new(result);

    // Copy parameters.
    for p in network.parameters() {
        let parameter = &network[p];
        result
            .add_parameter(parameter.get_name(), parameter.get_arity())
            .unwrap();
    }

    // Copy update functions.
    for v in result.variables() {
        // Technically, the models should have equivalent variable ids, so this should work.
        // Also, we already replaced all missing functions with parameters.
        let function = network.get_update_function(v).clone().unwrap();
        result.add_update_function(v, function).unwrap();
    }

    result
}

/// "Original" network contains the same parameters as perturbed network, but the parameters
/// actually do not matter. They are present only to ensure the two networks have the same
/// symbolic encoding.
///
/// The `perturb` parameter specifies which variables should be actually subject to perturbations.
pub fn make_original_network(
    network: &BooleanNetwork,
    perturbation_parameters: &mut HashMap<VariableId, ParameterId>,
    perturb: Vec<VariableId>,
) -> BooleanNetwork {
    let mut result = BooleanNetwork::new(network.as_graph().clone());

    // First, copy existing parameters from the original network:
    for old_id in network.parameters() {
        let parameter = network.get_parameter(old_id);
        let new_id = result
            .add_parameter(parameter.get_name(), parameter.get_arity())
            .unwrap();
        assert_eq!(
            new_id, old_id,
            "This should not happen. Encodings should be the same."
        );
    }

    // Then add control parameters and modify update functions:
    for v in network.variables().rev() {
        // We assume the function exists -- substituted implicit functions in normalization.
        let function = network.get_update_function(v).as_ref().unwrap();

        if perturb.contains(&v) {
            let v_perturbed = format!("{}_perturbed", network.get_variable_name(v));
            let parameter_id = result.add_parameter(v_perturbed.as_str(), 0).unwrap();
            perturbation_parameters.insert(v, parameter_id);

            // A little trick to avoid always cloning the value of fn_parameter...
            let fn_parameter = || FnUpdate::mk_param(parameter_id, &[]);

            // Set uncontrolled function to (v_perturbed || !v_perturbed) && f(...)
            // (The function has to *contain* the parameter to ensure the encoding is the same)
            let control_tautology = fn_parameter().or(FnUpdate::mk_not(fn_parameter()));
            let uncontrolled_function = control_tautology.and(function.clone());
            result
                .add_update_function(v, uncontrolled_function)
                .unwrap();
        } else {
            // If not perturbed, just copy what we have.
            result.add_update_function(v, function.clone()).unwrap();
        }
    }

    result
}

/// "Perturbed" network contains extra parameters which make it possible to disable perturbed
/// update functions (i.e. if the parameter is true, the network does not change that variable).
///
/// The `perturb` parameter specifies which variables should be actually subject to perturbations.
pub fn make_perturbed_network(
    network: &BooleanNetwork,
    perturbation_parameters: &mut HashMap<VariableId, ParameterId>,
    perturb: Vec<VariableId>,
) -> BooleanNetwork {
    let mut result = BooleanNetwork::new(network.as_graph().clone());

    // First, copy existing parameters from the original network:
    for old_id in network.parameters() {
        let parameter = network.get_parameter(old_id);
        let new_id = result
            .add_parameter(parameter.get_name(), parameter.get_arity())
            .unwrap();
        assert_eq!(
            new_id, old_id,
            "This should not happen. Encodings should be the same."
        );
    }

    // Then add control parameters and modify update functions:
    for v in network.variables().rev() {
        // We assume the function exists -- substituted implicit functions in normalization.
        let function = network.get_update_function(v).as_ref().unwrap();

        if perturb.contains(&v) {
            let v_perturbed = format!("{}_perturbed", network.get_variable_name(v));
            let parameter_id = result.add_parameter(v_perturbed.as_str(), 0).unwrap();
            perturbation_parameters.insert(v, parameter_id);
            // A little trick to avoid always cloning the value of fn_parameter...
            let fn_parameter = || FnUpdate::mk_param(parameter_id, &[]);

            // Set controlled function to (v_perturbed => v) && (!v_perturbed => f(...))
            let controlled_implies_v = fn_parameter().implies(FnUpdate::mk_var(v));
            let not_controlled_implies_f =
                FnUpdate::mk_not(fn_parameter()).implies(function.clone());
            let controlled_function = controlled_implies_v.and(not_controlled_implies_f);
            result.add_update_function(v, controlled_function).unwrap();
        } else {
            // If not perturbed, just copy what we have.
            result.add_update_function(v, function.clone()).unwrap();
        }
    }

    result
}
