use crate::ising::*;

pub fn stationary_infinite_temperature<I:Ising>(
  ising: &mut I, 
  ground_state: &Vec<bool>,
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
