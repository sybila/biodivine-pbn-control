use biodivine_lib_param_bn::{VariableId, BooleanNetwork, VariableIdIterator};
use biodivine_lib_param_bn::async_graph::AsyncGraph;
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

impl ControlledAsyncGraph {
    /// Create a new `AsyncGraph` from the given `BooleanNetwork`.
    pub fn new(network: BooleanNetwork) -> ControlledAsyncGraph {
        let n = network.clone();
        let g = AsyncGraph::new(network);
           ControlledAsyncGraph {
               network: n,
               graph: g.unwrap(),
               controls: HashMap::new(),
           }
    }

    pub fn set_controls(&mut self, from: IdState, to: IdState) {
        self.controls.clear();

        for v in self.network.graph().variable_ids() {
            if from.get_bit(v.into()) != to.get_bit(v.into()) {
                self.controls.insert(v, to.get_bit(v.into()));
            }
        }
    }

    /// Return the total number of states in this graph.
    pub fn num_states(&self) -> usize {
        return self.graph.num_states()
    }

    /// Return an empty parameter set.
    pub fn empty_params(&self) -> BddParams {
        return self.graph.empty_params()
    }

    /// Return a unit parameter set, i.e. all parameters that satisfy all static conditions.
    pub fn unit_params(&self) -> &BddParams {
        return &self.graph.unit_params();
    }

    /// Compute the parameter set which enables the value of `variable` to be flipped
    /// in the given `state`.
    ///
    /// Note that this set is not a subset of the `unit_set`! It is really just the flip condition.
    pub fn edge_params(&self, state: IdState, variable: VariableId) -> BddParams {
        let control = self.controls.get(&variable);
        if control != None {
            // Controlled variable can not be flipped
            return self.empty_params()
        }

        for v in self.network.graph().variable_ids() {
            let c = self.controls.get(&v);
            if c.is_some() && state.get_bit(variable.into()) != *c.unwrap() {
                // State is not valid as it does not fulfill the control condition
                return self.empty_params()
            }
        }

        return self.graph.edge_params(state, variable)
    }

    pub fn make_witness(&self, params: &BddParams) -> BooleanNetwork {
        return self.graph.make_witness(params);
    }

    pub fn find_permanent_control(&mut self, source: IdState, target: &StateSet) -> HashMap<IdState, BddParams> {
        let mut controls: HashMap<IdState, BddParams> = HashMap::new();
        let no_guard = StateSet::new_with_initial(self.num_states(), self.unit_params());

        let b = Bwd { graph: self };

        let backward_reach = guarded_reach(&b, target, &no_guard, &AtomicBool::new(false), &ProgressTracker::new(&self.graph));
        let mut weakBasin = HashMap::new();
        for (n, p) in backward_reach.iter() {
            weakBasin.insert(n, p.clone());
        }

        for (state, params) in weakBasin {
            self.set_controls(source, state);
            let permanentStrongBasin = find_strong_basin(self, target);
            let params_t = permanentStrongBasin.get(&state);
            if params_t.is_some() {
                let params = params_t.unwrap();
                if !params.is_empty() {
                    controls.insert(state, params.clone());
                }
            }
        }

        return controls
    }
}
