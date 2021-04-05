use std::fs;
use biodivine_lib_param_bn::symbolic_async_graph::SymbolicAsyncGraph;
use biodivine_lib_param_bn::BooleanNetwork;
use std::convert::TryFrom;
use biodivine_lib_param_bn::biodivine_std::traits::Set;
use std::collections::HashMap;
use std::time::Instant;
use biodivine_pbn_control::control::_algo_utils::find_attractors;
use biodivine_pbn_control::control::_impl_temporary_control::TemporaryControl;

fn main() {
    for m in ["myeloid", "cardiac", "erbb", "tumour", "mapk", "hgf"].iter() {
    //for m in ["myeloid", "cardiac"].iter() {
    let witness_str: &str = &fs::read_to_string(format!("models/{}_witness.aeon", m)).unwrap();
        let model = BooleanNetwork::try_from(witness_str).unwrap();
        let graph = &SymbolicAsyncGraph::new(model.clone()).unwrap();


        println!("Model: {}", m);
        println!("# variables: {}", model.clone().variables().len());
        //println!("{}", graph.unit_colors().approx_cardinality());

        let attractors = &find_attractors(graph);

        let all = attractors.iter().fold(graph.mk_empty_vertices(), |a, b| a.union(b));
        println!("# attractors: {}", attractors.clone().len());
        println!("# attractor vertices: {}", all.vertices().approx_cardinality());

        let mut vertices = Vec::new();
        for a in attractors {
            vertices.push(a.pick_vertex());
        }

        for (t_i, t) in vertices.clone().iter().enumerate() {
            for (s_i, s) in vertices.clone().iter().enumerate() {
                if t_i == s_i {
                    continue
                }
                for source in s.vertices().materialize().iter() {
                    for target in t.vertices().materialize().iter() {
                        let begin = Instant::now();
                        let _ = TemporaryControl::new(BooleanNetwork::try_from(witness_str).unwrap(), &source, &target);
                        let elapsed = begin.elapsed().as_millis();
                        println!("From attractor {} to attractor {} elapsed {} ms", s_i, t_i, elapsed);
                    }
                }

            }
        }
    }
}
