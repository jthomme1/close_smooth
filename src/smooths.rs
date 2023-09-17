use std::vec::Vec;
use std::thread;
use crate::composite::Composite;
//use once_cell::sync::Lazy;
use super::PRIMES;
use rayon::prelude::*;

//static NUM_THREADS: Lazy<usize> = Lazy::new(|| thread::available_parallelism().unwrap().get());

pub struct Smooths {
    pub b_ind: usize,
    pub lower_bound: u128,
    pub upper_bound: u128,
    // each exponent per prime
    pub state: Vec<Vec<Composite>>,
    pub smooths: Vec<u128>,
}

impl Smooths {
    pub fn new(lower_bound: u128, upper_bound: u128, b_ind: usize, state: &Vec<Vec<Composite>>) -> Self {
        Smooths {
            b_ind: b_ind,
            lower_bound: lower_bound,
            upper_bound: upper_bound,
            state: state.clone(),
            smooths: vec![],
        }
    }

    pub fn add_smooths(&mut self) {
        let lower_bound = self.lower_bound;
        let upper_bound = self.upper_bound;
        let b_ind = self.b_ind;
        let state_c = &self.state;
        let generate_some = |i: usize, e: u32| {
            let mut new_smooths = vec![];
            let ec = usize::try_from(e).unwrap();
            let mut state = {
                if state_c.len() == 0 {
                    Composite::new(i, e)
                } else {
                    state_c[i][ec-1].clone()
                }
            };
            let mut stop = Composite::new(i, e+1);
            if stop.value > upper_bound {
                stop = Composite::new(i, 0);
            }
            let mut c: usize = 0;
            //let mut flag = false;
            while state != stop && c != 20000 {
                //flag = true;
                if state.value >= lower_bound {
                    new_smooths.push(state.value);
                }
                state.inc_vec_with_bound(upper_bound);
                c += 1;
            }
            /*if flag && state == stop {
                println!("Generated all {}-smooth with exponent {e} numbers up to {}", PRIMES[i], upper_bound);
            } else {
                println!("ind: {i}, {:?}", state);
            }*/
            (new_smooths, state)
        };
        let (mut new_smooths, new_state): (Vec<u128>, Vec<Vec<Composite>>) = thread::scope(|s| {
            let mut handles = vec![];
            for i in 0..=b_ind {
                let mut inner_handles = vec![];
                let mut e = 1;
                let p = u128::try_from(PRIMES[i]).unwrap();
                let mut cur = p;
                loop {
                    let h = s.spawn(move || generate_some(i, e));
                    inner_handles.push(h);
                    if cur > upper_bound/p {
                        break;
                    }
                    cur *= p;
                    e += 1;
                }
                handles.push(inner_handles);
            }
            let mut all_smooth_vec: Vec<Vec<u128>> = vec![];
            let mut new_state: Vec<Vec<Composite>> = vec![];
            for ih in handles {
                let (smooths, states): (Vec<Vec<u128>>, Vec<Composite>) = ih.into_iter()
                  .map(|h| h.join().unwrap())
                  .unzip();
                new_state.push(states);
                all_smooth_vec.push(smooths.concat());
            }
            (all_smooth_vec.concat(), new_state)
        });
        new_smooths.par_sort_unstable();
        self.smooths.append(&mut new_smooths);
        self.smooths.par_sort_unstable();
        self.state = new_state;
    }

    pub fn len(&self) -> usize {
        self.smooths.len()
    }
}
