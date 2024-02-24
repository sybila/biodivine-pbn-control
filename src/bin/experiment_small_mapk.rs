
use biodivine_lib_param_bn::BooleanNetwork;
use biodivine_pbn_control::aeon::phentoype::build_phenotype;
use biodivine_pbn_control::perturbation::PerturbationGraph;
use biodivine_pbn_control::control::{ControlMap, PhenotypeOscillationType};
use std::borrow::Borrow;
use std::collections::HashMap;
use std::hash::Hash;

fn main() {
    let model_string = std::fs::read_to_string("./model.aeon").unwrap();
    let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    let perturbation_graph = PerturbationGraph::new(&bn);
    let mut phenotype_vals = HashMap::new();
    phenotype_vals.insert("v_Apoptosis", true);
    phenotype_vals.insert("v_Growth_Arrest", true);
    phenotype_vals.insert("v_Proliferation", false);

    let phenotype = build_phenotype(perturbation_graph.as_perturbed(), phenotype_vals);

    let map = perturbation_graph.ceiled_phenotype_permanent_control(
        phenotype.clone(),
        1,
        PhenotypeOscillationType::Forbidden,
        false,
        false,
    );
    let perturbations = map.working_perturbations(0.1, false);

    // let model_string = std::fs::read_to_string( "./models/small_mapk.aeon").unwrap();
    // let bn = BooleanNetwork::try_from(model_string.as_str()).unwrap();
    //
    // let mut controllable_vars = Vec::new();
    // let mut controllable_vars_size = 0;
    // for v in bn.clone().variables() {
    //     if bn.get_variable_name(v).as_str() != "v_ERK" {
    //         controllable_vars.push(v);
    //         controllable_vars_size+=1;
    //     }
    // }
    //
    // let perturbation_graph = PerturbationGraph::with_restricted_variables(&bn, controllable_vars);
    // let normal_graph = SymbolicAsyncGraph::new(bn.clone()).unwrap();
    //
    // let mut phenotype_vals = HashMap::new();
    // phenotype_vals.insert("v_ERK", true);

    // let phenotype = build_phenotype(perturbation_graph.as_perturbed(),
    //                                 phenotype_vals);

    // Print model statistics:
    // let model_variables = bn.num_vars();
    // assert!(i32::try_from(model_variables).is_ok());

    // All colors considered by the perturbation graph
    // let all_colors = perturbation_graph.unit_colors().approx_cardinality();
    // The (combinatorial) portion of colours that appear due to perturbation parameters.
    // let perturbation_colors = 2.0f64.powi(controllable_vars_size);
    // The (combinatorial) portion of colours that are carried over from the original model.
    // let model_colors = all_colors / perturbation_colors;

    // println!("Original model colors {}", normal_graph.unit_colors().approx_cardinality());
    // println!("Original model colors BDD vars {}", normal_graph.unit_colors().as_bdd().num_vars());
    // println!("Perturbed model colors {}", perturbation_graph.as_perturbed().unit_colors().approx_cardinality());
    // println!("As original model colors {}", perturbation_graph.as_original().unit_colors().approx_cardinality());
    // println!("Perturbed model colors BDD vars {}", perturbation_graph.unit_colors().as_bdd().num_vars());
    // println!("Variables: {}", model_variables);
    // println!("Controllable variables: {}", controllable_vars_size);
    // println!("Uncertainty colors: {}", model_colors);
    // println!("Perturbation colors: {}", perturbation_colors);
    //
    // println!(">>>>>>>>>>>>> All colors, do not restrict control call in any way");
    // println!("All colors: {}", all_colors);

    // let map = perturbation_graph.ceiled_phenotype_permanent_control(phenotype.clone(), 1, PhenotypeOscillationType::Required, false, false);
    // let perturbations = map.working_perturbations(0.1,  false);
    //
    // // Kinda OK
    // println!("Set working perturbations of EMPTY: {:?}", perturbations[0].1.approx_cardinality());
    //
    // // TODO: Why the below doesn't work? shouldn't be BDD vars removed thanks to projections in the previous step? (We copied colors INTO normal_grah unit colors
    // // println!("Set working perturbations of EMPTY: {:?}", perturbations[0].1.intersect(normal_graph.unit_colors()));
    //
    // // OK
    // println!("BDD num vars: {:?}", perturbations[0].1.as_bdd().num_vars());
    //
    // // Kinda OK
    // println!("Set working perturbations of {:?}: {:?}",perturbations[2].0, perturbations[2].1.approx_cardinality());
    //
    // // OK
    // println!("BDD num vars of {:?}: {:?}", perturbations[2].0, perturbations[2].1.as_bdd().num_vars());
    //
    // // Kinda OK
    // println!("Single working perturbations of EMPTY: {:?}", map.perturbation_working_colors( &HashMap::new()).approx_cardinality());
    // println!("BDD num vars: {:?}", map.perturbation_working_colors( &HashMap::new()).as_bdd().num_vars());
    //
    // // println!("Single working perturbations of EMPTY: {:?}", map.perturbation_working_colors(relevant_colors.clone(), &HashMap::new()).to_dot_string(perturbation_graph.as_symbolic_context()));
    // // TODO: Following doesn't work, because params of unperturbed variables were not removed using universal projection (see BDD above)
    // // println!("Single working perturbations of EMPTY (translated): {:?}", normal_graph.transfer_colors_from(&map.perturbation_working_colors(relevant_colors.clone(), &HashMap::new()), &perturbation_graph.as_original()).unwrap().intersect(&normal_graph.unit_colors()).approx_cardinality());
    //
    // // Kinda OK
    // println!("Single working perturbations: {:?}", map.perturbation_working_colors(&HashMap::from([("v_SPRY".to_string(), true)])).approx_cardinality());

    // // Same as above
    // // println!("Single working perturbations: {:?}", normal_graph.transfer_colors_from(&map.perturbation_working_colors(&HashMap::from([("v_SPRY".to_string(), true)])), &perturbation_graph.as_original()).unwrap().intersect(&normal_graph.unit_colors()).approx_cardinality());
    //
    // println!(">>>>>>>>>>>>> Original model colors translated into perturbed ones");
    // // 888 * 2024 -> OK
    // println!("Relevant colors: {}", relevant_colors.approx_cardinality());
    //
    // let mut map2 = perturbation_graph.ceiled_phenotype_permanent_control(phenotype.clone(), 1, PhenotypeOscillationType::Allowed, true, false);
    // let mut perturbations2 = map2.working_perturbations(0.1, false);
    //
    // // All OK
    // println!("Set working perturbations of EMPTY: {:?}", perturbations2[0].1.approx_cardinality());
    // println!("Set working perturbations of {:?}: {:?}",perturbations2[2].0, perturbations2[2].1.approx_cardinality());
    // println!("Single working perturbations: {:?}", map2.perturbation_working_colors(&HashMap::new()).approx_cardinality());
    // println!("Single working perturbations: {:?}", map2.perturbation_working_colors( &HashMap::from([("v_SPRY".to_string(), true)])).approx_cardinality());
    //
    // println!(">>>>>>>>>>>>> Original model colors with crossjoin");
    //
    // let irrelevant_colors = perturbation_graph.as_perturbed().unit_colors().minus(&relevant_colors).intersect(perturbation_graph.as_original().unit_colors());
    // let relevant_colors = perturbation_graph.as_perturbed().unit_colors().minus(&irrelevant_colors);
    //
    // println!("Relevant colors: {}", relevant_colors.approx_cardinality());
    //
    // map2 = perturbation_graph.ceiled_phenotype_permanent_control(phenotype, 1, PhenotypeOscillationType::Allowed, true, false);
    // perturbations2 = map2.working_perturbations(0.1, false);
    //
    // // TODO: All OK, I was just expecting, that maybe less colors will work, wondering, what I have computed here :D
    // println!("Set working perturbations of EMPTY: {:?}", perturbations2[0].1.approx_cardinality());
    // println!("Set working perturbations of {:?}: {:?}",perturbations2[2].0, perturbations2[2].1.approx_cardinality());
    // println!("Single working perturbations: {:?}", map2.perturbation_working_colors( &HashMap::new()).approx_cardinality());
    // println!("Single working perturbations: {:?}", map2.perturbation_working_colors( &HashMap::from([("v_SPRY".to_string(), true)])).approx_cardinality());
}
