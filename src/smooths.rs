use std::vec::Vec;
use std::thread;
use crate::composite::Composite;
use once_cell::sync::Lazy;
use super::PRIMES;
use rayon::prelude::*;

static NUM_THREADS: Lazy<usize> = Lazy::new(|| thread::available_parallelism().unwrap().get());

pub struct Smooths {
    pub b_ind: usize,
    pub lower_bound: u128,
    pub upper_bound: u128,
    pub smooths: Vec<u128>,
}

impl Smooths {
    pub fn new(lower_bound: u128, upper_bound: u128, b_ind: usize) -> Self {
        let mut ret = Smooths {
            b_ind: b_ind,
            lower_bound: lower_bound,
            upper_bound: upper_bound,
            smooths: vec![],
        };
        for i in 0..=b_ind {
            ret.add_prime(i);
            println!("Done with {}", PRIMES[i]);
        }
        ret
    }

    fn add_prime(&mut self, ind: usize) {
        let lower_bound = self.lower_bound;
        let upper_bound = self.upper_bound;
        let generate_with_fixed = |start_off: usize| {
            let mut c = Composite::new(ind, 1);

            if !c.inc_vec_by_n_with_bound(start_off, upper_bound) {
                return vec![];
            }

            let mut new_smooths: Vec<u128> = vec![];

            loop {
                if c.value >= lower_bound {
                    new_smooths.push(c.value);
                }
                if !c.inc_vec_by_n_with_bound(*NUM_THREADS, upper_bound) {
                    break;
                }
            }
            new_smooths.par_sort_unstable();
            new_smooths
        };
        let mut new_smooths = thread::scope(|s| {
            let mut handles = vec![];
            for i in 0..*NUM_THREADS {
                let h = s.spawn(move || generate_with_fixed(i));
                handles.push(h);
            }
            handles.into_iter()
                .map(|h| h.join().unwrap())
                .collect::<Vec<Vec<u128>>>()
                .concat()
        });
        new_smooths.par_sort_unstable();
        self.smooths.append(&mut new_smooths);
        self.smooths.par_sort_unstable();
    }
}

