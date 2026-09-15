mod ising;
mod ising_reader;
mod anneal;
use crate::ising_reader::*;
use crate::anneal::*;
use crate::ising::*;
use clap::{Parser,Subcommand};
use std::fs::File;
use std::path::PathBuf;
use std::io::{Write,BufWriter};
use bitvec::prelude::*;
use std::process;

#[derive(Parser)]
#[command(about = "Implementation of High Temperature Simulated \
Annealing Ising Problem with Approximations for Convergence Rate")]
struct Args{
  input_file: PathBuf,
  #[command(subcommand)]
  mode: Mode,
}

#[derive(Subcommand)]
enum Mode{
  #[command(about = "Infinite Temperature Experiments")]
  Inf{
    #[arg(short = 's', long = "stationary_file",
      default_value = "stationary_hits_infinite.txt")]
    stationary_file_name:PathBuf,
    #[arg(short = 'n', long = "max_steps", default_value="100")]
    max_steps:usize,
    #[arg(short = 'r', long = "runs", default_value = "100")]
    runs:usize,
  },
  #[command(about = "Theoretical Convergence Rate, High Temperature")]
  HighConv{
    #[arg(short = 'T', long = "temperature_override", default_value = None)]
    temperature_override:Option<f64>,
  },
  #[command(about = "Finite Temperature Experiments")]
  FiniteT{
    #[arg(short = 's', long = "stationary_file",
      default_value = "stationary_hits_finite.txt")]
    stationary_file_name:PathBuf,
    #[arg(short = 'n', long = "max_steps", default_value="100")]
    max_steps:usize,
    #[arg(short = 'r', long = "runs", default_value = "100")]
    runs:usize,
    #[arg(short = 'T', long = "temperature_override", default_value="None")]
    temperature_override: Option<f64>,
  }
}

fn finite_temperature_option<I:Ising>(
  ising_model: &mut I,
  max_steps: usize,
  runs: usize,
  stationary_file_name: PathBuf,
  temperature:f64
  ){
  let stationary_file = File::create(stationary_file_name).expect("main.rs\
Error:could not create new stationary file");
  let mut stationary_out = BufWriter::new(stationary_file);
  let all_zeros = bitvec![usize, LocalBits; 0; ising_model.n_points()];
  for _ in 0..runs{
  let hitting_time = stationary_finite_temperature(
    ising_model,
    &all_zeros,
    max_steps,
    temperature
  );
    match hitting_time{
      Some(steps) => writeln!(stationary_out, "{}", steps).unwrap(),
      None => writeln!(stationary_out, "NA").unwrap()
    }
    ising_model.init_spins();
    ising_model.init_cost();
  }
}

fn infinite_temperature_option<I:Ising>(
  ising_model: &mut I,
  max_steps:usize,
  runs:usize,
  stationary_file_name:PathBuf
  ){
  let stationary_file = File::create(stationary_file_name).expect("main.rs\
Error:could not create new stationary file");
  let mut stationary_out = BufWriter::new(stationary_file);
  let all_zeros = bitvec![usize, LocalBits; 0; ising_model.n_points()];
  for _ in 0..runs{
    let hitting_time = stationary_infinite_temperature(
      ising_model,
      &all_zeros,
      max_steps
      );
    match hitting_time{
      Some(steps) => writeln!(stationary_out, "{}", steps).unwrap(),
      None => writeln!(stationary_out, "NA").unwrap()
    }
    ising_model.init_spins();
    ising_model.init_cost();
  }
}

fn main(){
  let args = Args::parse();
  let input_file = args.input_file;
  let (mut ising_model, temperature_from_file)=from_ising_file_disjoint_simple(
    input_file
  );
  match args.mode{
    Mode::Inf{stationary_file_name, max_steps, runs} =>{
      infinite_temperature_option(
        &mut ising_model,
        max_steps,
        runs,
        stationary_file_name)
      }
    Mode::FiniteT{
      stationary_file_name,
      max_steps,
      runs,
      temperature_override} => {
      let temperature:f64;
      match temperature_override{
        Some(temperature_opt) => { 
          temperature = temperature_opt;
        }
        None => {
          temperature = temperature_from_file;
        }
      }
      finite_temperature_option(
        &mut ising_model,
        max_steps,
        runs,
        stationary_file_name,
        temperature);
    }
    Mode::HighConv{temperature_override} => {
      let temperature:f64;
      match temperature_override{
        Some(temperature_opt) => { 
          temperature = temperature_opt;
        }
        None => {
          temperature = temperature_from_file;
        }
      }
      let (sign, log_sum) = theoretical_perturbation_naive(
        &mut ising_model,
        &None,
        temperature);
      println!("log absolute perturbation:{}", log_sum);
      println!("sign:{}", if sign {"-"} else {"+"});
    }
  }
}
