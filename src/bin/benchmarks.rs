use std::fs;
use biodivine_lib_param_bn::async_graph::AsyncGraph;
use biodivine_pbn_control::strong_basin::_algo_utils::{find_attractors, get_all_params_with_attractor};
use clap::{Arg, App, ArgMatches};
use biodivine_pbn_control::strong_basin::_algo_sb_parallel_fixed_point::find_strong_basin;
use biodivine_pbn_control::controlled_async_graph::ControlledAsyncGraph;
use biodivine_aeon_server::scc::StateSet;
use std::time::{Instant, Duration, SystemTime};
use serde_json::json;
use std::path::Path;
use biodivine_lib_std::{IdState};
use biodivine_lib_param_bn::BooleanNetwork;
use std::convert::TryFrom;
use std::fs::{File, OpenOptions};
use std::io::prelude::*;


const STRONG_BASIN_ARG: &str = "strong-basin";
const TEMP_CONTROL_ARG: &str = "temporary-control";
const PERM_CONTROL_ARG: &str = "permanent-control";
const CONFIG_ARG: &str = "config-path";
const MODELS_DIR_ARG: &str = "models-path";
const OUT_FILE_ARG: &str = "out-file";
const TAG_ARG: &str = "tag";


fn main() {
    let matches = App::new("Parametrised Boolean networks control benchmarks")
        .version("0.1")
        .author("Eva Smijakova <xsmijak1@fi.muni.cz>")
        .about("Measures performance of current parametrised Boolean network control implementation.")
        .arg(Arg::with_name(STRONG_BASIN_ARG)
            .short("sb")
            .long(STRONG_BASIN_ARG)
            .help("Script will measure performance of strong basin"))
        .arg(Arg::with_name(TEMP_CONTROL_ARG)
             .short("tc")
             .long(TEMP_CONTROL_ARG)
             .help("Script will measure performance of temporary control"))
        .arg(Arg::with_name(PERM_CONTROL_ARG)
             .short("pc")
             .long(PERM_CONTROL_ARG)
             .help("Script will measure performance of permanent control"))
        .arg(Arg::with_name(MODELS_DIR_ARG)
            .short("md")
            .long(MODELS_DIR_ARG)
            .help("Path to the directory containing benchmark models")
            .takes_value(true)
            .value_name("DIRECTORY"))
        .arg(Arg::with_name(CONFIG_ARG)
            .short("c")
            .long(CONFIG_ARG)
            .help("Path to the configuration file containing information about models")
            .takes_value(true)
            .value_name("FILE"))
        .arg(Arg::with_name(OUT_FILE_ARG)
            .short("of")
            .long(OUT_FILE_ARG)
            .help("Name of the output file")
            .takes_value(true)
            .value_name("FILE"))
        .arg(Arg::with_name(TAG_ARG)
            .short("t")
            .long(TAG_ARG)
            .help("Tag of the current implementation version")
            .value_name("TAG")
            .takes_value(true)
            .required(true))
        .get_matches();

    let measure_basin = matches.is_present(STRONG_BASIN_ARG);
    let measure_temp =  matches.is_present(TEMP_CONTROL_ARG);
    let measure_perm =  matches.is_present(PERM_CONTROL_ARG);

    if !measure_basin && !measure_temp && !measure_perm {
        println!("No measurements to be done, please provide at least one of: --{}, --{}, --{}", STRONG_BASIN_ARG, TEMP_CONTROL_ARG, PERM_CONTROL_ARG)
    }

    let models_path = matches.value_of(MODELS_DIR_ARG).unwrap_or("../../models");
    let config_file = matches.value_of(CONFIG_ARG).unwrap_or("benchmark_config.json");
    let out_file = matches.value_of(OUT_FILE_ARG).unwrap_or("benchmark_results.csv");
    let tag = matches.value_of(TAG_ARG).unwrap();

    create_csv(out_file);
    let config_str: &str = &fs::read_to_string(config_file).unwrap();

    println!("{}", config_str);
    println!("{:?}", json!(config_str).as_object());

    for (model_name, attractors) in json!(config_str).as_object().unwrap() {
        let path = Path::new(models_path).join(model_name);
        let content: &str = &fs::read_to_string(path).unwrap();
        let model = BooleanNetwork::try_from(content).unwrap();
        let graph = &ControlledAsyncGraph::new(model);

        for a in attractors.as_array().unwrap() {
            let target = IdState::from(a.as_u64().unwrap() as usize);
            let target_seed = &StateSet::new_with_fun(graph.num_states(), |s| if s.eq(&target) { Some(graph.unit_params().clone()) } else { None });

            let relevant_params = get_all_params_with_attractor(graph, target);
            let relevant_params_cardinality =  relevant_params.cardinality();

            let begin = Instant::now();
            let basin = find_strong_basin(graph, target_seed, graph.unit_params());
            println!("Strong basin computation time (ms): {:?}", begin.elapsed().as_millis());
            if measure_basin {
                write_to_csv(out_file, 4, "SB", begin.elapsed().as_millis(),
                             model_name,target, target, graph.num_states(),
                                relevant_params_cardinality, basin.len(), tag)
            }

            for b in attractors.as_array().unwrap() {
                let source = IdState::from(b.as_u64().unwrap() as usize);

                if measure_temp {
                    let begin = Instant::now();
                    let temporary = graph.find_temporary_control(source, target_seed);
                    println!("Temporary control can be done in {} ways.", temporary.len());
                    println!("Temporary control computation time (ms): {:?}", begin.elapsed().as_millis());

                    write_to_csv(out_file, 4, "TEMP", begin.elapsed().as_millis(),
                                 model_name,source, target, graph.num_states(),
                                 relevant_params_cardinality, basin.len(), tag)
                }

                if measure_perm {
                    let begin = Instant::now();
                    let permanent = graph.find_permanent_control(source, target_seed);
                    println!("Permanent control can be done in {} ways.", permanent.len());
                    println!("Permanent control computation time (ms): {:?}", begin.elapsed().as_millis());

                    write_to_csv(out_file, 4, "PERM", begin.elapsed().as_millis(),
                                 model_name,source, target, graph.num_states(),
                                 relevant_params_cardinality, basin.len(), tag)
                }

                println!();
            }
        }
    }
}

fn create_csv(file_path: &str) {
    if !Path::new(file_path).exists() {
        let mut file = OpenOptions::new()
            .write(true)
            .create(true)
            .open(file_path)
            .unwrap();
        if let Err(e) = writeln!(file, "time,parallelism,method,duration,model,source,target,model_size,parameters_size,sb_size,tag") {
            eprintln!("Couldn't write to file: {}", e);
        }
    }
}

fn write_to_csv(file_path: &str, parallelism: i64, method: &str, duration: u128, model_name: &str,
                source: IdState, target: IdState, model_size: usize, param_size: f64, sb_size: usize, tag: &str) {
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(file_path)
        .unwrap();

    if let Err(e) = writeln!(file, "{:?},{},{},{:?},{},{},{},{},{},{},{}",
                             SystemTime::now(), parallelism, method, duration, model_name, source,
                             target, model_size, param_size, sb_size, tag) {
        eprintln!("Couldn't write to file: {}", e);
    }
}
