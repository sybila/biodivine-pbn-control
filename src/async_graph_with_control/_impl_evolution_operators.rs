use biodivine_lib_std::param_graph::{EvolutionOperator, InvertibleEvolutionOperator};
use crate::async_graph_with_control::{Fwd, FwdIterator, Bwd, BwdIterator};
use biodivine_lib_std::IdState;
use biodivine_lib_param_bn::bdd_params::BddParams;

impl<'a> EvolutionOperator for Fwd<'a> {
    type State = IdState;
    type Params = BddParams;
    type Iterator = FwdIterator<'a>;

    fn step(&self, current: IdState) -> Self::Iterator {
        return FwdIterator {
            graph: self.graph,
            variables: self.graph.network.graph().variable_ids(),
            state: current,
        };
    }
}

impl<'a> InvertibleEvolutionOperator for Fwd<'a> {
    type InvertedOperator = Bwd<'a>;

    fn invert(&self) -> Self::InvertedOperator {
        return Bwd { graph: self.graph };
    }
}

impl Iterator for FwdIterator<'_> {
    type Item = (IdState, BddParams);

    fn next(&mut self) -> Option<Self::Item> {
        return if let Some(var) = self.variables.next() {
            let target = self.state.flip_bit(var.into());
            let edge_params = self.graph.edge_params(self.state, var);
            Some((target, edge_params))
        } else {
            None
        };
    }
}

impl<'a> EvolutionOperator for Bwd<'a> {
    type State = IdState;
    type Params = BddParams;
    type Iterator = BwdIterator<'a>;

    fn step(&self, current: IdState) -> Self::Iterator {
        return BwdIterator {
            graph: self.graph,
            variables: self.graph.network.graph().variable_ids(),
            state: current,
        };
    }
}

impl<'a> InvertibleEvolutionOperator for Bwd<'a> {
    type InvertedOperator = Fwd<'a>;

    fn invert(&self) -> Self::InvertedOperator {
        return Fwd { graph: self.graph };
    }
}

impl Iterator for BwdIterator<'_> {
    type Item = (IdState, BddParams);

    fn next(&mut self) -> Option<Self::Item> {
        return if let Some(var) = self.variables.next() {
            let source = self.state.flip_bit(var.into());
            let edge_params = self.graph.edge_params(source, var);
            Some((source, edge_params))
        } else {
            None
        };
    }
}