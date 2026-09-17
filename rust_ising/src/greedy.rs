use crate::ising::*;
use bitvec::prelude::*;
pub fn greedy<I:Ising>(
  ising:&mut I, 
  inner_sweeps:usize, 
  outer_sweeps:usize
  ) -> BitVec
{
  let mut results = Vec::new();
  for _ in 0..outer_sweeps{
    ising.init_spins();
    ising.init_cost();
    let mut new_cost = ising.cost();
    for _ in 0..inner_sweeps{
      let mut changed=false;
      
      for node in 0..ising.n_points() {
        if ising.cost_diff(node) < 0.0 {
          ising.flip_spin(node);
          changed=true;
        }
      }
      new_cost = ising.cost();
      if !changed{
        break
      }
    }
    results.push((ising.config().clone(), new_cost));
  }
  results.into_iter()
    .min_by(|(_,cost_a), (_, cost_b)| cost_a.total_cmp(cost_b))
    .unwrap()
    .0
}
