mod ising;
mod ising_reader;
mod anneal;
use crate::ising_reader::*;
use crate::anneal::*;
use clap::Parser;
use std::fs::File;
use std::path::PathBuf;
use std::io::{Write,BufWriter};

#[derive(Parser)]
#[command(about = "Implementation of High Temperature Simulated
  Annealing Ising Problem with Approximations for Convergence Rate")]
struct Args{
  input_file: PathBuf,
  #[arg(short = 'n', long = "max_steps", default_value="100")]
  max_steps:usize,
  #[arg(short = 'r', long = "runs", default_value = "100")]
  runs:usize,
  #[arg(short = 's', long = "stationary_file",
    default_value = "stationary_hits.txt")]
  stationary_file:PathBuf,
}

fn main(){
  let args = Args::parse();
  let input_file = args.input_file;
  let runs = args.runs;
  let max_steps = args.max_steps;
  let stationary_file = File::create(args.stationary_file).unwrap();
  let mut stationary_out  = BufWriter::new(stationary_file);

  let (mut ising_model, temperature)=from_ising_file_disjoint_simple(
    input_file
  );
  let all_zeros=vec![false; ising_model.n_points];
  for _ in 0..runs{
    let hitting_time = stationary_infinite_temperature(
      &mut ising_model,
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
