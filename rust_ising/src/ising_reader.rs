use std::io::{BufReader, BufRead};
use std::path::Path;
use crate::ising::*;
use std::fs::File;
use bitvec::prelude::*;

 macro_rules! local_index { 
  ($rel_index:expr) => 
  {(
    if $rel_index % 2 == 0 { $rel_index + 1 } else { $rel_index - 1 }
    )}
 }


pub fn from_ising_file_disjoint_simple(path: impl AsRef<Path>) -> 
(IsingDisjoint, f64)
{
  let file = File::open(path).unwrap();
  let reader = BufReader::new(file);
  let mut lines = reader.lines();
  let mut in_section = false;
  for line in lines.by_ref(){
    let line = line.unwrap();
    if line.starts_with("$temp"){
      in_section = true;
      break
    }
  }
  assert!(in_section, "Could not find temperature");
  let line = lines.by_ref()
    .next()
    .expect("Could not find dimension and sizes")
    .expect("error reading line");
  let mut parser = line.split_whitespace();
  let temp:f64 = parser
    .next().unwrap()
    .parse().unwrap();

  for line in lines.by_ref(){
    let line = line.unwrap();
    if line.starts_with("$dim_sizes"){
      in_section = true;
      break
    }
  }
  assert!(in_section, "Could not find dimension and sizes");
  let line = lines.by_ref()
    .next()
    .expect("Could not find dimension and sizes")
    .expect("error reading line");
  let mut parser = line.split_whitespace();
  let dim:usize = parser
    .next().unwrap()
    .parse().unwrap();
  assert!(dim > 0, "Dimension cannot be 0, that is nonsensical");
  let mut sizes = Vec::<usize>::new();
  for _i in 0..dim{
    let length:usize = parser
      .next().expect("Error: Dimension and Sizes mismatch")
      .parse().expect("Error parsing sizes");
    sizes.push(length);
  }
  let mut ising_instance = IsingDisjoint::new(dim, sizes);
  in_section = false;
  for line in lines.by_ref(){
    let line = line.unwrap();
    if line.starts_with("$edge_weights_start") {
      in_section = true;
      break;
    }
  }
  assert!(in_section, "Could not find edge_weights");
  let mut weights = vec![0.0; ising_instance.deg*ising_instance.n_points];
  /* Difficult command time! */
  for node in 0..ising_instance.n_points{
    let mut i = 0;
    for neighbour in ising_instance.neighbours(node){
      if neighbour > node {
        let line = lines.by_ref()
          .next()
          .expect("Insufficient number of weights")
          .expect("Error reading line");
        let weight:f64 = line.split_whitespace()
          .next()
          .expect("Blank weight")
          .parse()
          .expect("Error reading line");
        weights[ising_instance.deg*node + i] = weight;
        weights[ising_instance.deg*neighbour + local_index!(i)] = weight;
      }
      i += 1;
    }
  }
  ising_instance.set_edges(weights);
  in_section=false;
  for line in lines.by_ref(){
    let line = line.unwrap();
    if line.starts_with("$mu") {
      in_section = true;
      break;
    }
  }
  assert!(in_section, "Could not find magnetic moment mu");
  ising_instance.set_mu(
    lines.by_ref()
      .next()
      .expect("Could not find magnetic moment mu")
      .expect("Error reading line below mu")
      .split_whitespace()
      .next()
      .expect("Blank magnetic moment")
      .parse::<f64>()
      .expect("Error reading magnetic moment mu")
  );
  in_section=false;
  for line in lines.by_ref(){
    let line = line.unwrap();
    if line.starts_with("$external_magnetic_field") {
      in_section = true;
      break;
    }
  }
  assert!(in_section, "Could not find external magnetic field section");
  ising_instance.set_magnetic_field(
    (0..ising_instance.n_points).map(|_| 
      lines.by_ref()
        .next()
        .expect("Error: Not enough entries in the external magnetic field")
        .expect("Error reading external magnetic field entries")
        .split_whitespace()
        .next()
        .expect("Blank magnetic field entry")
        .parse::<f64>()
        .expect("Error reading magnetic field entry")
    ).collect()
  );
  in_section=false;
  for line in lines.by_ref(){
    let line = line.unwrap();
    if line.starts_with("$starting_configs") {
      in_section = true;
      break;
    }
  }
  let starting_configs: Option<Vec<StartingConfig>>;
  if in_section {
    let mut starting_configs_base = Vec::<StartingConfig>::new();
    let hex_length = ising_instance.n_points.div_ceil(8);
    while let Some(line) = lines.next() {
      let line = line.expect("Error reading line");
      if line.trim().is_empty() {
        break;
      }
      let mut line = line.split_whitespace();
      let mut hex = Vec::<u8>::new();
      for _ in 0..hex_length{
        hex.push(
          u8::from_str_radix(
            line.next()
              .expect(
 "Error: Not enough entries for hex code of starting config"),
            16
          )
          .expect("Error: No integer detected")
        );
        }
      let chunk_size = std::mem::size_of::<usize>();
      let config_usize:Vec<usize> = hex.chunks(chunk_size)
        .map(|chunk| {
          chunk.iter()
            .enumerate()
            .fold(0usize, |acc, (i, &byte)|{
              acc | ((byte as usize) << i*8)
            })
          })
        .collect();
      let starting_config = usize_to_config(config_usize, 
        ising_instance.n_points);
      let prob = line.next()
        .expect("Could not find probability associated with starting 
          configuration")
        .parse::<f64>()
        .expect("Probability for starting configuration not a float");
      //DEBUG 
      #[cfg(debug_assertions)]
      {
        if prob > 1.0{
          println!("Warning: Weight greater than 1.");
          }
      }
      starting_configs_base.push(StartingConfig{
        config:starting_config, weight:prob
      })
      }
    starting_configs=Some(starting_configs_base)
  }
  else {
    println!("Warning: Starting Configurations not found. Defaulting \
to random Starting Configurations");
    starting_configs = None;
  }
  ising_instance.set_starting_configs(starting_configs);
  ising_instance.init_spins();
  ising_instance.init_cost();
  (ising_instance, temp)
}

fn usize_to_config(uvec:Vec<usize>, n_points:usize)->BitVec{
    let mut bit_vec = BitVec::from_vec(uvec);
    bit_vec.truncate(n_points);
    bit_vec
}


