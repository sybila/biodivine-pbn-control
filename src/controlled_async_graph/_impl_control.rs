use biodivine_lib_param_bn::{VariableId, BooleanNetwork};
use biodivine_lib_param_bn::async_graph::{AsyncGraph, DefaultEdgeParams, AsyncGraphEdgeParams};
use std::collections::HashMap;
use crate::controlled_async_graph::{ControlledAsyncGraph, Bwd, ControlVariable};
use biodivine_lib_param_bn::bdd_params::{BddParams, BddParameterEncoder};
use biodivine_lib_std::IdState;
use biodivine_aeon_server::scc::algo_par_reach::guarded_reach;
use std::sync::atomic::AtomicBool;
use biodivine_aeon_server::scc::{ProgressTracker, StateSet};
use crate::strong_basin::_algo_sb_parallel_fixed_point::find_strong_basin;
use biodivine_lib_std::param_graph::{Params};
use biodivine_lib_bdd::{BddVariableSetBuilder};
use std::borrow::Borrow;
use std::fs::read_to_string;
use crate::strong_basin::_algo_utils::get_all_params_with_attractor;

impl ControlledAsyncGraph {
    /// Create a new `AsyncGraph` from the given `BooleanNetwork`.
    pub fn new(network: BooleanNetwork) -> ControlledAsyncGraph {
        let bn = network.clone();
        let mut vars = BddVariableSetBuilder::new();
        let mut control_vars = HashMap::new();

        for vid in bn.graph().variable_ids() {
            let v = bn.graph().get_variable(vid);
            let bdd_name = format!("{}_is_controlled", v.get_name());
            let is_controlled = vars.make_variable(&bdd_name);
            let bdd_name2 = format!("{}_control_value", v.get_name());
            let control_val = vars.make_variable(&bdd_name2);
            control_vars.insert(vid,ControlVariable {
                original_name: v.get_name().clone(),
                is_controlled_variable: is_controlled,
                control_value_variable: control_val,
            });
        }

        let encoder = BddParameterEncoder::new_with_custom_variables(&bn, vars);
        let edges = DefaultEdgeParams::new_with_custom_encoder(bn.clone(), encoder.clone()).unwrap();
        let g = AsyncGraph::new_with_edges(edges).unwrap();

        return ControlledAsyncGraph {
            network: bn.clone(),
            encoder: encoder.clone(),
            graph: g,
            controlVariables: control_vars,
        }
    }

    /// Return the total number of states in this graph.
    pub fn num_states(&self) -> usize {
        return self.graph.num_states()
    }

    /// Return an empty parameter set.
    pub fn empty_params(&self) -> &BddParams { return &self.graph.empty_params(); }

    /// Return a unit parameter set, i.e. all parameters that satisfy all static conditions.
    pub fn unit_params(&self) -> &BddParams { return &self.graph.unit_params(); }

    /// Compute the parameter set which enables the value of `variable` to be flipped
    /// in the given `state`.
    ///
    /// Note that this set is not a subset of the `unit_set`! It is really just the flip condition.
    pub fn edge_params(&self, state: IdState, variable: VariableId) -> BddParams {
        let default_edge_params = self.graph.edges().edge_params(state, variable);
        let mut allowed_control_params =  self.unit_params().clone();
        //Ensure, that in edge params:;
        // * variableId_is_controlled = FALSE && variableId_control_value = FALSE
        // * for all other variables V:
        // ** V_is_controlled=FALSE && V_control_value=FALSE
        // ** V_is_controlled=TRUE && V_control_value=state.V
        for v in self.network.graph().variable_ids() {
            let controlled_var = self.controlVariables.get(&v).unwrap();
            if v == variable {
                let i_c = self.encoder.bdd_variables.mk_not_var(controlled_var.is_controlled_variable);
                let c_v = self.encoder.bdd_variables.mk_not_var(controlled_var.control_value_variable);
                let not_controlled = i_c.and(&c_v);
                allowed_control_params = allowed_control_params.intersect(&BddParams::from(not_controlled));
            } else {
                let i_c1 = self.encoder.bdd_variables.mk_not_var(controlled_var.is_controlled_variable);
                let c_v1 = self.encoder.bdd_variables.mk_not_var(controlled_var.control_value_variable);
                let i_c2 = self.encoder.bdd_variables.mk_var(controlled_var.is_controlled_variable);
                let c_v2;
                if state.get_bit(v.into()) {
                    c_v2 = self.encoder.bdd_variables.mk_var(controlled_var.control_value_variable);
                } else {
                    c_v2 = self.encoder.bdd_variables.mk_not_var(controlled_var.control_value_variable);
                }
                let not_controlled = i_c1.and(&c_v1);
                let controlled_correctly = i_c2.and(&c_v2);
                let viable_params = not_controlled.or(&controlled_correctly);

                allowed_control_params = allowed_control_params.intersect(&BddParams::from(viable_params));
            }
        }

        return default_edge_params.intersect(&allowed_control_params);
    }

    pub fn make_witness(&self, params: &BddParams) -> BooleanNetwork {
        return self.graph.make_witness(params);
    }

    pub fn find_permanent_control(&self, source: IdState, target: &StateSet) -> HashMap<IdState, BddParams> {
        let mut controls: HashMap<IdState, BddParams> = HashMap::new();
        let no_guard = StateSet::new_with_initial(self.num_states(), self.unit_params());

        let b = Bwd { graph: self };

        let backward_reach = guarded_reach(&b, target, &no_guard, &AtomicBool::new(false), &ProgressTracker::new(&self.graph));
        let mut weakBasin = HashMap::new();
        for (n, p) in backward_reach.iter() {
            weakBasin.insert(n, p.clone());
        }

        let strongBasins = find_strong_basin(self, target, self.unit_params());
        for (state, params) in weakBasin {
            let params_t = strongBasins.get(&state);
            if params_t.is_some() {
                let mut params = params_t.unwrap().clone();
                if !params.is_empty() {
                    // Get strong basin of the control
                    let mut control_params = self.encoder.bdd_variables.mk_true();
                    for v in self.network.graph().variable_ids() {
                        let controlled_var = self.controlVariables.get(&v).unwrap();
                        if state.get_bit(v.into()) != source.get_bit(v.into()) {
                            let i_c = self.encoder.bdd_variables.mk_var(controlled_var.is_controlled_variable);
                            let c_v;
                            if state.get_bit(v.into()) {
                                c_v = self.encoder.bdd_variables.mk_var(controlled_var.control_value_variable);
                            } else {
                                c_v = self.encoder.bdd_variables.mk_not_var(controlled_var.control_value_variable);
                            }
                            let is_controlled = i_c.and(&c_v);
                            control_params = control_params.and(&is_controlled);
                        } else {
                            let i_c = self.encoder.bdd_variables.mk_not_var(controlled_var.is_controlled_variable);
                            let c_v = self.encoder.bdd_variables.mk_not_var(controlled_var.control_value_variable);
                            let not_controlled = i_c.and(&c_v);
                            control_params = control_params.and(&not_controlled);
                        }
                    }
                    params = params.intersect(&BddParams::from(control_params));
                    if !params.is_empty() {
                        controls.insert(state, params.clone());
                    }
                }
            }
        }

        return controls
    }

    pub fn find_temporary_control(&self, source: IdState, target: &StateSet) -> HashMap<IdState, BddParams> {
        let mut controls: HashMap<IdState, BddParams> = HashMap::new();
        let no_guard = StateSet::new_with_initial(self.num_states(), self.unit_params());

        let b = Bwd { graph: self };

        let backward_reach = guarded_reach(&b, target, &no_guard, &AtomicBool::new(false), &ProgressTracker::new(&self.graph));
        let mut weak_basin = HashMap::new();
        for (n, p) in backward_reach.iter() {
            weak_basin.insert(n, p.clone());
        }

        let mut without_control = self.encoder.bdd_variables.mk_true();
        for v in self.network.graph().variable_ids() {
            let controlled_var = self.controlVariables.get(&v).unwrap();
            let i_c = self.encoder.bdd_variables.mk_not_var(controlled_var.is_controlled_variable);
            let c_v = self.encoder.bdd_variables.mk_not_var(controlled_var.control_value_variable);
            let not_controlled = i_c.and(&c_v);
            without_control = without_control.and(&not_controlled);
        }

        let strong_basin = find_strong_basin(self, target, &self.unit_params().intersect(&BddParams::from(without_control)));

        let strong_basin_seed = &StateSet::new_with_fun(self.graph.num_states(), |s| if strong_basin.contains_key(&s) { Some(BddParams::from(self.encoder.bdd_variables.mk_true())) } else { None });
        let extendedStrongBasin = find_strong_basin(self, strong_basin_seed, self.unit_params());

        for (state, params) in weak_basin {
            let params_t = extendedStrongBasin.get(&state);
            if params_t.is_some() {
                let mut params = params_t.unwrap().clone();
                if !params.is_empty() {
                    let mut control_params = self.encoder.bdd_variables.mk_true();
                    for v in self.network.graph().variable_ids() {
                        let controlled_var = self.controlVariables.get(&v).unwrap();
                        if state.get_bit(v.into()) != source.get_bit(v.into()) {
                            let i_c = self.encoder.bdd_variables.mk_var(controlled_var.is_controlled_variable);
                            let c_v;
                            if state.get_bit(v.into()) {
                                c_v = self.encoder.bdd_variables.mk_var(controlled_var.control_value_variable);
                            } else {
                                c_v = self.encoder.bdd_variables.mk_not_var(controlled_var.control_value_variable);
                            }
                            let is_controlled = i_c.and(&c_v);
                            control_params = control_params.and(&is_controlled);
                        } else {
                            let i_c = self.encoder.bdd_variables.mk_not_var(controlled_var.is_controlled_variable);
                            let c_v = self.encoder.bdd_variables.mk_not_var(controlled_var.control_value_variable);
                            let not_controlled = i_c.and(&c_v);
                            control_params = control_params.and(&not_controlled);
                        }
                    }
                    params = params.intersect(&BddParams::from(control_params));
                    if !params.is_empty() {
                        controls.insert(state, params.clone());
                    }
                }
            }
        }

        return controls
    }
}

