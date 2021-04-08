fn main() {}
/*
fn main() {
    for m in ["myeloid", "cardiac", "erbb", "tumour", "mapk", "hgf"].iter() {
        let witness_str: &str = &fs::read_to_string(format!("models/{}_witness.aeon", m)).unwrap();
        let witness_model = BooleanNetwork::try_from(witness_str).unwrap();
        let witness_graph = &SymbolicAsyncGraph::new(witness_model.clone()).unwrap();

        println!("Model: {}", m);
        println!("# variables: {}", witness_model.clone().variables().len());

        let attractors = &find_attractors(witness_graph);
        println!("# attractors: {}", attractors.clone().len());
        let mut vertices = Vec::new();
        for a in attractors {
            vertices.push(a.pick_vertex());
        }

        let param_str: &str = &fs::read_to_string(format!("models/{}_4unknown.aeon", m)).unwrap();
        let param_model = BooleanNetwork::try_from(param_str).unwrap();
        let param_graph = &SymbolicAsyncGraph::new(param_model.clone()).unwrap();
        println!("# params: {}", param_graph.unit_colors().approx_cardinality());

        println!("WITNESS TEMPORARY");

        for (t_i, t) in vertices.clone().iter().enumerate() {
            if t_i > 3 {
                break;
            }
            let begin = Instant::now();
            for (s_i, s) in vertices.clone().iter().enumerate() {
                if s_i > 3 {
                    break;
                }
                if t_i == s_i {
                    continue;
                }
                for source in s.vertices().materialize().iter() {
                    for target in t.vertices().materialize().iter() {
                        let _ = TemporaryControl::new(BooleanNetwork::try_from(witness_str).unwrap(), &source, &target);
                    }
                }
            }
            let elapsed = begin.elapsed().as_millis();
            println!("To attractor {} elapsed {} ms", t_i, elapsed);
        }

        println!("WITNESS PERMANENT");

        for (t_i, t) in vertices.clone().iter().enumerate() {
            if t_i > 3 {
                break;
            }
            let begin = Instant::now();
            for (s_i, s) in vertices.clone().iter().enumerate() {
                if s_i > 3 {
                    break;
                }
                if t_i == s_i {
                    continue;
                }
                for source in s.vertices().materialize().iter() {
                    for target in t.vertices().materialize().iter() {
                        let _ = PermanentControl::new(BooleanNetwork::try_from(witness_str).unwrap(), &source, &target);
                    }
                }
            }
            let elapsed = begin.elapsed().as_millis();
            println!("To attractor {} elapsed {} ms", t_i, elapsed);
        }

        println!("4 PARAMS TEMPORARY");
        for (t_i, t) in vertices.clone().iter().enumerate() {
            if t_i > 3 {
                break;
            }
            let begin = Instant::now();
            for (s_i, s) in vertices.clone().iter().enumerate() {
                if s_i > 3 {
                    break;
                }
                if t_i == s_i {
                    continue;
                }
                for source in s.vertices().materialize().iter() {
                    for target in t.vertices().materialize().iter() {
                        let _ = TemporaryControl::new(BooleanNetwork::try_from(param_str).unwrap(), &source, &target);
                    }
                }
            }
            let elapsed = begin.elapsed().as_millis();
            println!("To attractor {} elapsed {} ms", t_i, elapsed);
        }

        println!("4 PARAMS PERMANENT");
        for (t_i, t) in vertices.clone().iter().enumerate() {
            if t_i > 3 {
                break;
            }
            let begin = Instant::now();
            for (s_i, s) in vertices.clone().iter().enumerate() {
                if s_i > 3 {
                    break;
                }
                if t_i == s_i {
                    continue;
                }
                for source in s.vertices().materialize().iter() {
                    for target in t.vertices().materialize().iter() {
                        let _ = PermanentControl::new(BooleanNetwork::try_from(param_str).unwrap(), &source, &target);
                    }
                }
            }
            let elapsed = begin.elapsed().as_millis();
            println!("To attractor {} elapsed {} ms", t_i, elapsed);
        }
    }
}
*/
