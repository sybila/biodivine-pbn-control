use biodivine_lib_param_bn::{VariableId, BooleanNetwork};
use biodivine_lib_param_bn::async_graph::{AsyncGraph, AsyncGraphEdgeParams};
use std::collections::HashMap;
use crate::async_graph_with_control::{AsyncGraphWithControl, Bwd};
use biodivine_lib_param_bn::bdd_params::{BddParams};
use biodivine_lib_std::IdState;
use biodivine_aeon_server::scc::algo_par_reach::guarded_reach;
use std::sync::atomic::AtomicBool;
use biodivine_aeon_server::scc::{ProgressTracker, StateSet};
use crate::strong_basin_control_experimental::_algo_sb_parallel_fixed_point::find_strong_basin;
use biodivine_lib_std::param_graph::{Params};

impl AsyncGraphWithControl {
    /// Create a new `AsyncGraph` from the given `BooleanNetwork`.
    pub fn new(network: BooleanNetwork) -> AsyncGraphWithControl {

        let g = AsyncGraph::new(network.clone());
           AsyncGraphWithControl {
               network: network.clone(),
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
    pub fn empty_params(&self) -> &BddParams { return self.graph.empty_params(); }

    /// Return a unit parameter set, i.e. all parameters that satisfy all static conditions.
    pub fn unit_params(&self) -> &BddParams {
        return self.graph.unit_params();
    }

    /// Compute the parameter set which enables the value of `variable` to be flipped
    /// in the given `state`.
    ///
    /// Note that this set is not a subset of the `unit_set`! It is really just the flip condition.
    pub fn edge_params(&self, state: IdState, variable: VariableId) -> BddParams {
        let control = self.controls.get(&variable);
        if control != None {
            // Controlled variable can not be flipped
            return self.empty_params().clone();
        }

        for v in self.network.graph().variable_ids() {
            match self.controls.get(&v) {
                Some(c) if *c != state.get_bit(v.into()) => {
                    // State is not valid as it does not fulfill the control condition
                    return self.empty_params().clone()
                },
                _ => ()
            }
        }

        return self.graph.edges().edge_params(state, variable)
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

        self.set_controls(source, source);
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
        let strongBasinSeed = &StateSet::new_with_fun(self.graph.num_states(), |s| if strongBasin.contains_key(&s) { Some(strongBasin.get(&s).unwrap().clone()) } else { None });

        for (state, params) in weakBasin {
            self.set_controls(source, state);
            let extendedStrongBasin = find_strong_basin(self, strongBasinSeed);
            let params_t = extendedStrongBasin.get(&state);
            if params_t.is_some() {
                let params = params_t.unwrap();
                if !params.is_empty() {
                    controls.insert(state, params.clone());
                }
            }
        }

        self.set_controls(source, source);
        return controls
    }
}

