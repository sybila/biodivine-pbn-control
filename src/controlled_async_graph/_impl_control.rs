use biodivine_lib_param_bn::{VariableId, BooleanNetwork, VariableIdIterator};
use biodivine_lib_param_bn::async_graph::{AsyncGraph, DefaultEdgeParams};
use std::collections::HashMap;
use crate::controlled_async_graph::{ControlledAsyncGraph, Fwd, FwdIterator, Bwd, BwdIterator};
use biodivine_lib_param_bn::bdd_params::{BddParams, BddParameterEncoder};
use biodivine_lib_std::IdState;
use std::env::var;
use biodivine_aeon_server::scc::algo_reach::reach;
use biodivine_aeon_server::scc::algo_par_reach::guarded_reach;
use std::sync::atomic::AtomicBool;
use biodivine_aeon_server::scc::{ProgressTracker, StateSet};
use crate::strong_basin::_algo_sb_parallel_fixed_point::find_strong_basin;
use biodivine_lib_std::param_graph::{Params, EvolutionOperator, InvertibleEvolutionOperator};
use biodivine_lib_bdd::BddVariableSetBuilder;
use crate::async_graph_with_control::AsyncGraphWithControl;

impl ControlledAsyncGraph {
    /// Create a new `AsyncGraph` from the given `BooleanNetwork`.
    pub fn new(network: BooleanNetwork) -> ControlledAsyncGraph {
        let bn = network.clone();
        let mut vars = BddVariableSetBuilder::new();

        for vid in bn.graph().variable_ids() {
            let v = bn.graph().get_variable(vid);
            let bdd_name = format!("{}_is_controlled", v.get_name());
            bdd.make_variable(&bdd_name);
            let bdd_name2 = format!("{}_control_value", v.get_name());
            bdd.make_variable(&bdd_name2);
        }

        let encoder = BddParameterEncoder::new_with_custom_variables(&n, vars);
        let edges = DefaultEdgeParams::new_with_custom_encoder(n, encoder).unwrap();
        let g = AsyncGraph::new_with_edges(edges).unwrap();
        return ControlledAsyncGraph {
            network: n,
            graph: g,
        }
    }

    /// Return the total number of states in this graph.
    pub fn num_states(&self) -> usize {
        return self.graph.num_states()
    }

    /// Return an empty parameter set.
    pub fn empty_params(&self) -> BddParams { return self.graph.empty_params() }

    /// Return a unit parameter set, i.e. all parameters that satisfy all static conditions.
    pub fn unit_params(&self) -> BddParams {
        return self.graph.unit_params().clone();
    }

    /// Compute the parameter set which enables the value of `variable` to be flipped
    /// in the given `state`.
    ///
    /// Note that this set is not a subset of the `unit_set`! It is really just the flip condition.
    pub fn edge_params(&self, state: IdState, variable: VariableId) -> BddParams {
        let default_edge_params = self.graph.edge_params(state, variable);
        //TODO ensure, that in edge params:
        // * variableId_is_controlled = FALSE && variableId_control_value = FALSE
        // * for all other variables V:
        // ** V_is_controlled=FALSE && V_control_value=FALSE
        // ** V_is_controlled=TRUE && V_control_value=state.V
        return default_edge_params;
    }

    pub fn make_witness(&self, params: &BddParams) -> BooleanNetwork {
        return self.graph.make_witness(params);
    }

    pub fn find_permanent_control(&mut self, source: IdState, target: &StateSet) -> HashMap<IdState, BddParams> {
        let mut currentParams = self.graph.unit_params();
        let bdd = currentParams.into_bdd();

        let mut controls: HashMap<IdState, BddParams> = HashMap::new();
        let no_guard = StateSet::new_with_initial(self.num_states(), self.unit_params());

        let b = Bwd { graph: self };

        let backward_reach = guarded_reach(&b, target, &no_guard, &AtomicBool::new(false), &ProgressTracker::new(&self.graph));
        let mut weakBasin = HashMap::new();
        for (n, p) in backward_reach.iter() {
            weakBasin.insert(n, p.clone());
        }

        let strongBasings = find_strong_basin(self, target);
        for (state, params) in weakBasin {
            let params_t = permanentStrongBasin.get(&state);
            if params_t.is_some() {
                let params = params_t.unwrap();
                if !params.is_empty() {
                    //TODO filter params so that there are only params where is control set to Hamming difference of state and target, insert these params
                    controls.insert(state, params.clone());
                }
            }
        }

        return controls
    }

    pub fn find_temporary_control(&mut self, source: IdState, target: &StateSet) -> HashMap<IdState, BddParams> {
        let mut controls: HashMap<IdState, BddParams> = HashMap::new();
        let no_guard = StateSet::new_with_initial(self.num_states(), self.unit_params());

        let b = Bwd { graph: self };

        let backward_reach = guarded_reach(&b, target, &no_guard, &AtomicBool::new(false), &ProgressTracker::new(&self.graph));
        let mut weakBasin = HashMap::new();
        for (n, p) in backward_reach.iter() {
            weakBasin.insert(n, p.clone());
        }

        let strongBasin = find_strong_basin(self, target);
        let strongBasinSeed = &StateSet::new_with_fun(graph.num_states(), |s| if strongBasin.contains_key(&s) { Some(strongBasin.get(&s).unwrap().clone()) } else { None });
        let extendedStrongBasin = find_strong_basin(self, strongBasinSeed);

        for (state, params) in weakBasin {
            self.set_controls(source, state);
            let params_t = extendedStrongBasin.get(&state);
            if params_t.is_some() {
                let params = params_t.unwrap();
                if !params.is_empty() {
                    //TODO filter params so that there are only params where is control set to Hamming difference of state and target, insert these params
                    controls.insert(state, params.clone());
                }
            }
        }

        return controls
    }
}

