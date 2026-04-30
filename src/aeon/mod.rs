/// Xie-Beerel TSCC algorithm enhanced with TGR as preprocessing.
pub mod attractors;
//utils for working with phenotypes
pub mod phentoype;
/// Reachability algorithms that use saturation for improved efficiency.
pub mod reachability;
/// Transition guided reduction quickly eliminates most non-attractor states in a graph.
mod tgr;

#[cfg(test)]
pub mod config;
