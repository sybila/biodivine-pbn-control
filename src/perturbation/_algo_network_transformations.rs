use biodivine_lib_param_bn::{BooleanNetwork, FnUpdate, ParameterId, RegulatoryGraph, VariableId};
use std::collections::HashMap;
use std::convert::TryFrom;

/// Create a copy of the given `model`, but make every regulation non-observable, and add an
/// unspecified auto-regulation to each variable (that does not have an autoregulation already).
/// Additionally, convert every implicit update function to and explicit parameter.
///
/// This is a necessary pre-processing step before creating a `PerturbationGraph`.
pub fn normalize_network(network: &BooleanNetwork) -> BooleanNetwork {
    // Copy graph, but with non-observable regulations.
    let mut result = RegulatoryGraph::new(
        network
            .variables()
            .map(|it| network.get_variable_name(it).clone())
            .collect(),
    );
    for regulation in network.as_graph().regulations() {
        result
            .add_regulation(
                network.get_variable_name(regulation.get_regulator()),
                network.get_variable_name(regulation.get_target()),
                false,
                regulation.get_monotonicity(),
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
        // Technically, the models should have equivalent variable ids!
        if let Some(function) = network.get_update_function(v) {
            result.add_update_function(v, function.clone()).unwrap();
        } else {
            // Create an explicit parameter to replace the implicit function.
            let regulators = result.regulators(v);
            let parameter = result
                .add_parameter(
                    format!("update_{}", result.get_variable_name(v)).as_str(),
                    u32::try_from(regulators.len()).unwrap(),
                )
                .unwrap();
            result
                .add_update_function(v, FnUpdate::Param(parameter, regulators))
                .unwrap();
        }
    }

    result
}

/// "Original" network contains the same parameters as perturbed network, but the parameters
/// actually do not matter. They are present only to ensure the two networks have the same
/// symbolic encoding.
pub fn make_original_network(
    network: &BooleanNetwork,
    perturbation_parameters: &mut HashMap<VariableId, ParameterId>,
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
        let v_perturbed = format!("{}_perturbed", network.get_variable_name(v));
        let parameter_id = result.add_parameter(v_perturbed.as_str(), 0).unwrap();
        perturbation_parameters.insert(v, parameter_id);

        // We assume the function exists -- substituted implicit functions in normalization.
        let function = network.get_update_function(v).as_ref().unwrap();
        // A little trick to avoid always cloning the value of fn_parameter...
        let fn_parameter = || FnUpdate::mk_param(parameter_id, &[]);

        // Set uncontrolled function to (v_perturbed || !v_perturbed) && f(...)
        // (The function has to *contain* the parameter to ensure the encoding is the same)
        let control_tautology = fn_parameter().or(FnUpdate::mk_not(fn_parameter()));
        let uncontrolled_function = control_tautology.and(function.clone());
        result
            .add_update_function(v, uncontrolled_function)
            .unwrap();
    }

    result
}

/// "Perturbed" network contains extra parameters which make it possible to disable perturbed
/// update functions (i.e. if the parameter is true, the network does not change that variable).
///
/// The new parameters are saved into the provided hash map.
pub fn make_perturbed_network(
    network: &BooleanNetwork,
    perturbation_parameters: &mut HashMap<VariableId, ParameterId>,
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
        let v_perturbed = format!("{}_perturbed", network.get_variable_name(v));
        let parameter_id = result.add_parameter(v_perturbed.as_str(), 0).unwrap();
        perturbation_parameters.insert(v, parameter_id);

        // We assume the function exists -- substituted implicit functions in normalization.
        let function = network.get_update_function(v).as_ref().unwrap();
        // A little trick to avoid always cloning the value of fn_parameter...
        let fn_parameter = || FnUpdate::mk_param(parameter_id, &[]);

        // Set controlled function to (v_perturbed => v) && (!v_perturbed => f(...))
        let controlled_implies_v = fn_parameter().implies(FnUpdate::mk_var(v));
        let not_controlled_implies_f = FnUpdate::mk_not(fn_parameter()).implies(function.clone());
        let controlled_function = controlled_implies_v.and(not_controlled_implies_f);
        result.add_update_function(v, controlled_function).unwrap();
    }

    result
}