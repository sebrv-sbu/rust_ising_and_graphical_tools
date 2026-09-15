use std::env;
use std::path::Path;
use std::fs::write;
fn binomial(n:usize, k:usize) -> usize{
  if k > n { return 0; }
  let k = if k>n-k {n-k} else {k};
  let mut result = 1usize;
  let mut i =0;
  while i < k{
    result = result * (n-i)/(i+1);
    i+=1;
  }
  result
}

fn powi_const(x:f64, n:i32) -> f64{
  let mut result = 1.0;
  let mut applications = 0;
  while n > applications{
    result *= x;
    applications += 1;
  }
  result
}

fn integral_weight_precomputed(m:usize, d:usize) -> f64{
  let n = 2 * d -m;
  let mut sum = 0.0;
  let mut k = 0;
  while k <= n{
    let summand = binomial(n,k) as f64 * powi_const(2.0, (n-k) as i32)
      / (m+k+1) as f64;
    if k % 2 == 0{
      sum += summand;
    } else {
      sum -= summand;
    }
    k+=1;
  }
  sum
}

fn make_integral_weights(n:usize, d:usize)->Vec<f64> {
  (0..n).map(|m| integral_weight_precomputed(m,d)).collect()
}

fn main(){
  let dim = env::var("ISING_DIM")
    .unwrap_or_else(|_| "2".to_string())
    .parse::<usize>()
    .expect("ISING_DIM must be a positive integer");
  let out_dir = env::var("OUT_DIR").unwrap();
  let dest = std::path::Path::new(&out_dir).join("dim.rs");
  
  std::fs::write(
    &dest,
    format!("pub const D: usize={dim};\n"),
  ).unwrap();
  
  let n = 2 * dim + 1;
  let weights = make_integral_weights(n, dim);
  let weights_src = weights.iter()
    .map(|w| format!("{:?}", w))
    .collect::<Vec<_>>()
    .join(", ");
  write(Path::new(&out_dir).join("weights.rs"),
    format!("pub const INTEGRAL_WEIGHTS: [f64; {}] = [{weights_src}];\n", n))
    .unwrap();
  println!("cargo:rerun-if-env-changed=ISING_DIM");
}
