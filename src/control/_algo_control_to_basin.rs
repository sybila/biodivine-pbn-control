use std::collections::{HashMap, HashSet};
use biodivine_lib_std::IdState;
use biodivine_lib_param_bn::bdd_params::BddParams;
use biodivine_lib_bdd::Bdd;

pub fn find_smallest_control_to_basin(source:&IdState, basin: &HashMap<IdState, BddParams>) -> Vec<(IdState, BddParams)> {
    let mut smallest_vec: Vec<(IdState, BddParams)> = basin.iter().map(|(k, val)| (k.clone(), val.clone())).collect();
    smallest_vec.sort_by(|a, b| control_dist(source, &a.0).cmp(&control_dist(source, &b.0)));
    return smallest_vec;
}

pub fn find_robust_control_to_basin(basin: &HashMap<IdState, BddParams>) -> Vec<(IdState, BddParams)> {
    let mut robust_vec: Vec<(IdState, BddParams)> = basin.iter().map(|(k, val)| (k.clone(), val.clone())).collect();
    robust_vec.sort_by(|a, b| a.1.cardinality().partial_cmp(&b.1.cardinality()).unwrap() ) ;
    return robust_vec;
}

pub fn find_best_control_to_basin(source:&IdState, basin: &HashSet<IdState>) -> IdState {
    if basin.is_empty() {
        return source.clone();
    }

    let mut b_iter = basin.iter();
    let mut best_state = b_iter.next().unwrap();
    let mut best_dist = control_dist(source, best_state);

    for i in b_iter {
        let dist = control_dist(source, i);
        if dist < best_dist {
            best_state = i;
            best_dist = dist;
        }
    }

    return best_state.clone();
}

pub fn control_dist(source:&IdState, target:&IdState) -> u32 {
    let s: usize = source.clone().into();
    let t: usize = target.clone().into();
    let res = s ^ t;
    return res.count_ones();
}