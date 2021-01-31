mod _impl_control;
mod _impl_evolution_operators;

use biodivine_lib_param_bn::{VariableId, BooleanNetwork, VariableIdIterator};
use std::collections::HashMap;
use biodivine_lib_param_bn::async_graph::{AsyncGraph, DefaultEdgeParams};
use biodivine_lib_std::param_graph::{InvertibleGraph, Graph};
use biodivine_lib_std::{IdStateRange, IdState};
use biodivine_lib_param_bn::bdd_params::BddParams;

pub struct AsyncGraphWithControl {
    pub graph: AsyncGraph<DefaultEdgeParams>,
    network: BooleanNetwork,
    controls: HashMap<VariableId, bool>
}

/// A forward `EvolutionOperator` of the `ControlledAsyncGraph`.
pub struct Fwd<'a> {
    graph: &'a AsyncGraphWithControl,
}

/// A backward `EvolutionOperator` of the `ControlledAsyncGraph`.
pub struct Bwd<'a> {
    graph: &'a AsyncGraphWithControl,
}

/// An iterator over successors of a state in the `ControlledAsyncGraph`.
pub struct FwdIterator<'a> {
    graph: &'a AsyncGraphWithControl,
    state: IdState,
    variables: VariableIdIterator,
}

/// An iterator over predecessors of a state in the `ControlledAsyncGraph`.
pub struct BwdIterator<'a> {
    graph: &'a AsyncGraphWithControl,
    state: IdState,
    variables: VariableIdIterator,
}

impl<'a> Graph for &'a AsyncGraphWithControl {
    type State = IdState;
    type Params = BddParams;
    type States = IdStateRange;
    type FwdEdges = Fwd<'a>;
    type BwdEdges = Bwd<'a>;

    fn states(&self) -> Self::States {
        return IdStateRange::new(self.num_states());
    }

    fn fwd(&self) -> Self::FwdEdges {
        return Fwd { graph: self };
    }

    fn bwd(&self) -> Self::BwdEdges {
        return Bwd { graph: self };
    }
}

impl<'a> InvertibleGraph for &'a AsyncGraphWithControl {
    type FwdEdges = Fwd<'a>;
    type BwdEdges = Bwd<'a>;
}
