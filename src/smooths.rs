use std::vec::Vec;
use std::thread;
use crate::composite::Composite;
use once_cell::sync::Lazy;
use super::PRIMES;
use rayon::prelude::*;

static NUM_THREADS: Lazy<usize> = Lazy::new(|| thread::available_parallelism().unwrap().get());

pub struct Smooths {
    pub b_ind: usize,
    pub upper_bound: u128,
    pub smooths: Vec<u128>,
}

impl Smooths {
    pub fn new(upper_bound: u128, b_ind: usize) -> Self {
        let mut ret = Smooths {
            b_ind: b_ind,
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
        let upper_bound = self.upper_bound;
        let generate_with_fixed = |start_off: usize| {
            let mut c = Composite::new(ind, 1);

            if !c.inc_vec_by_n_with_bound(start_off, upper_bound) {
                return vec![];
            }

            let mut new_smooths: Vec<u128> = vec![];

            loop {
                new_smooths.push(c.value);
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

    //TODO: CHECK AGAIN
    pub fn find_ind_gt(&self, b: u128) -> Option<usize> {
        if self.smooths.len() == 0 || self.smooths[self.smooths.len()-1] <= b {
            return None;
        }
        let ind = match self.smooths.binary_search(&b) {
            Ok(x) => x+1,
            Err(x) => x,
        };
        assert!(self.smooths[ind] > b);
        Some(ind)
    }

    pub fn find_ind_le(&self, b: u128) -> Option<usize> {
        if self.smooths.len() == 0 || self.smooths[0] > b {
            return None;
        }
        let ind = match self.smooths.binary_search(&b) {
            Ok(x) => x,
            Err(x) => x-1,
        };
        assert!(self.smooths[ind] <= b);
        Some(ind)
    }

}

