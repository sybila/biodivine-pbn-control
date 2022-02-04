pub mod control;

/// "Standard" algorithms mostly adapted from Aeon.
pub mod aeon;

/// A perturbed version of the symbolic asynchronous graph.  We keep it in a separate
/// module because the implementation is very fragile and we want to hide as much of the
/// API as possible even from ourselves!
pub mod perturbation;
