use crate::ising::*;
use bitvec::prelude::*;

pub fn stationary_infinite_temperature<I:Ising>(
  ising: &mut I, 
  ground_state: &BitVec,
  max_steps:usize) -> Option<usize>{
  let mut hitting_time = 0;
  while hitting_time < max_steps{
    if ground_state == ising.config(){
      return Some(hitting_time)
    }
    ising.inf_anneal();
    hitting_time += 1;
  }
  if ground_state == ising.config(){
    return Some(hitting_time)
  }
  None
}

pub fn theoretical_perturbation_naive<I:Ising>(
  ising: &mut I,
  ground_state_opt: &Option<BitVec>,
  temperature: f64
  )->(bool, f64){
  /* since $N/2^{2N+1}$ is quite a tiny number, *
   * we will to change the method.  */
  let mut sum = 0.0;
  let mut sign = false;
  for node in 0..ising.n_points(){
    let ground_state_cluster:BitVec;
    let ground_state_cluster_opt: Option<&Bits>;
    let ground_state_center:Option<bool>;
    match ground_state_opt{
      Some(ground_state) => {
        let neighbours = ising.neighbours(node);
        ground_state_cluster = (0..neighbours.len())
          .map(|i| ground_state[neighbours[i]])
          .collect();
        ground_state_cluster_opt = Some(ground_state_cluster.as_bitslice());
        ground_state_center = Some(ground_state[node]);
      }
      None =>{
      ground_state_cluster_opt = None;
      ground_state_center = None;
      }
    }
    let mut node_sum = 0.0;
    for cluster_index in 0 .. 1<<ising.deg(){
      let cluster = &cluster_index.view_bits::<LocalBits>()[..ising.deg()]; 
      node_sum += ising.cluster_cost(
        temperature,
        node,
        cluster,
        ground_state_cluster_opt,
        ground_state_center);
    }
    #[cfg(debug_assertions)]
    println!("node: {}: sum = {}", node, node_sum);
    sum += node_sum;
  }
  if sum < 0.0{
    sign = true;
    sum = -1.0 * sum;
  }
  (
    sign,
    sum.log2()-((ising.n_points() + 1+ising.deg()) as f64) +
      (ising.n_points() as f64).log2()
  )
}
pub fn stationary_finite_temperature<I:Ising>(
  ising: &mut I,
  ground_state: &BitVec,
  max_steps:usize,
  temperature:f64) -> Option<usize>
{
  let mut hitting_time = 0;
  while hitting_time < max_steps{
    if ground_state == ising.config(){
      return Some(hitting_time)
    }
    ising.anneal(temperature);
    hitting_time += 1;
  }
  if ground_state == ising.config(){
    return Some(hitting_time)
  }
  None
}
