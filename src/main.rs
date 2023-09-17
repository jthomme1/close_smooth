use crate::{smooths::Smooths, composite::Composite};
use std::vec::Vec;
use once_cell::sync::Lazy;
use primal;
use std::env;
use std::thread;
use rayon::prelude::*;
//use std::collections::BTreeSet;
use rand::{self, rngs::ThreadRng, seq::SliceRandom};
use std::cmp::{min, max};
//use average::Skewness;
use num::One;
use num::bigint::BigUint;
use num_prime::RandPrime;

pub mod composite;
pub mod smooths;

// this should suffice for now
static PRIME_BOUND: usize = 1<<20;
static PRIMES: Lazy<Vec<u128>> = Lazy::new(||
                                           primal::Sieve::new(PRIME_BOUND)
                                           .primes_from(0)
                                           .map(|x| u128::try_from(x).unwrap())
                                           .collect());
static NUM_THREADS: Lazy<usize> = Lazy::new(|| thread::available_parallelism().unwrap().get());

#[allow(non_snake_case)]
fn one_smooth_in(a: BigUint, b: BigUint, B: u64) {
    assert!(a < b);
    assert!(a > (B*B).into());
    assert!(B*25 > a.bits()*a.bits());
    // rounding down is fine, because rounding up would result in a product outside the bounds
    // also, rounding down we don't lose anything because we're considering integers
    let ub_x1 = b.sqrt();
    //println!("upper_bound_smaller_factor: {ub_x1}");
    // we need to round up because otherwise we could get a product outside of the interval
    // also, we don't lose anything rounding up because we're considering integers
    let mut lb_x2 = a.sqrt();
    if lb_x2.clone() * lb_x2.clone() < b {
        lb_x2 += 1u32;
    }
    //println!("lower_bound_bigger_factor: {lb_x2}");

    let ub_x2 = (B*b.clone()).sqrt();
    let mut lb_x1 = (a.clone()/B).sqrt();
    if lb_x1.clone() * lb_x1.clone() < a.clone()/B {
        lb_x1 += 1u32;
    }
    let mut smooths_per_level: Vec<Vec<BigUint>> = vec![];
    let mut base_case_state: Vec<Vec<Composite>> = vec![];
    loop {
        some_smooths_in(&lb_x1, &ub_x2, B, &mut smooths_per_level, &mut base_case_state, 0);
        let smaller_smooths = &smooths_per_level[0];
        let ub_x2_ind = smaller_smooths.len();
        let lb_x1_ind = 0;

        let ub_x1_ind = match smaller_smooths.binary_search(&ub_x1) {
            Ok(i) => i+1,
            Err(i) => i,
        };
        let lb_x2_ind = match smaller_smooths.binary_search(&lb_x2) {
            Ok(i) => i,
            Err(i) => i,
        };

        let x1_range = &smaller_smooths[lb_x1_ind..ub_x1_ind];
        let x2_range = &smaller_smooths[lb_x2_ind..ub_x2_ind];

        //let mut rng = &mut rand::thread_rng();

        let found = x1_range.into_par_iter().find_any(|&x1| {
            let up = b.clone()/x1;
            let mut low = a.clone()/x1;
            if low.clone() * x1 < a {
                low += 1u32;
            }
            low = max(low, x1.clone());
            //println!("For {small}, the interval for the bigger number is: [{low}, {up}]");
            // up_ind is not inclusive
            let up_ind = match x2_range.binary_search(&up) {
                Ok(i) => i+1,
                Err(i) => i,
            };
            // low_ind is inclusive
            let low_ind = match x2_range.binary_search(&low) {
                Ok(i) => i,
                Err(i) => i,
            };
            let elem_x2_range = &x2_range[low_ind..up_ind];
            for x2 in elem_x2_range {
                let smooth = x1 * x2;
                assert!(smooth >= a && smooth <= b);
                println!("Found {smooth}");
                return true;
            }
            return false;
        });
        match found {
            Some(_) => break,
            None => continue,
        };
    }
}

#[allow(non_snake_case)]
fn some_smooths_in(a: &BigUint, b: &BigUint, B: u64, smooths_per_level: &mut Vec<Vec<BigUint>>, base_case_state: &mut Vec<Vec<Composite>>, level: usize) {
    assert!(a < b);
    // base case
    if *a < (B*B).into() {
        if smooths_per_level.len() <= level {
            smooths_per_level.push(vec![]);
        }
        println!("Entering base case for interval: [{a}, {b}]");
        let ind: usize = PRIMES.binary_search(&(B.into())).unwrap();
        let mut smooths = Smooths::new(a.to_u64_digits()[0] as u128, b.to_u64_digits()[0] as u128, ind, base_case_state);
        smooths.add_smooths();
        println!("Exiting from base case with {} nrs for interval: [{a}, {b}]", smooths.len());
        *base_case_state = smooths.state;
        smooths_per_level[level].append(&mut smooths.smooths.iter().map(|&x| BigUint::from(x)).collect());
        smooths_per_level[level].par_sort_unstable();
        return;
    }
    if smooths_per_level.len() <= level {
        smooths_per_level.push(vec![]);
    }
    let start_num = smooths_per_level[level].len();
    println!("Considering interval [{a}, {b}] with currently {start_num} numbers");
    // rounding down is fine, because rounding up would result in a product outside the bounds
    // also, rounding down we don't lose anything because we're considering integers
    let ub_x1 = b.sqrt();
    //println!("upper_bound_smaller_factor: {ub_x1}");
    // we need to round up because otherwise we could get a product outside of the interval
    // also, we don't lose anything rounding up because we're considering integers
    let mut lb_x2 = a.sqrt();
    if lb_x2.clone() * lb_x2.clone() < *b {
        lb_x2 += 1u32;
    }
    //println!("lower_bound_bigger_factor: {lb_x2}");

    let ub_x2 = (B*b.clone()).sqrt();
    let mut lb_x1 = (a.clone()/B).sqrt();
    if lb_x1.clone() * lb_x1.clone() < a.clone()/B {
        lb_x1 += 1u32;
    }
    some_smooths_in(&lb_x1, &ub_x2, B, smooths_per_level, base_case_state, level+1);
    let smaller_smooths = &smooths_per_level[level+1];

    let ub_x2_ind = smaller_smooths.len();
    let lb_x1_ind = 0;

    let ub_x1_ind = match smaller_smooths.binary_search(&ub_x1) {
        Ok(i) => i+1,
        Err(i) => i,
    };
    let lb_x2_ind = match smaller_smooths.binary_search(&lb_x2) {
        Ok(i) => i,
        Err(i) => i,
    };
    /*
    println!("x1_range: [{lb_x1}, {ub_x1}]");
    println!("The first smooth number in there: {}", smaller_smooths[lb_x1_ind]);
    println!("The first smooth number after: {}", smaller_smooths[ub_x1_ind]);
    println!("x2_range [{lb_x2}, {ub_x2}]");
    println!("The first smooth number in there: {}", smaller_smooths[lb_x2_ind]);
    println!("The first smooth number after: {}", smaller_smooths[ub_x2_ind]);
    */
    let x1_range = &smaller_smooths[lb_x1_ind..ub_x1_ind];
    let x2_range = &smaller_smooths[lb_x2_ind..ub_x2_ind];


    let smooth_combs = || {
        let mut new_smooths = vec![];
        let mut rng = &mut rand::thread_rng();

        //let mut x1_hits: Vec<u32> = vec![];
        for x1 in x1_range.choose_multiple(&mut rng, min(10000, ub_x1_ind - lb_x1_ind)) {
            let up = b.clone()/x1;
            let mut low = a.clone()/x1;
            if low.clone() * x1 < *a {
                low += 1u32;
            }
            low = max(low, x1.clone());
            //println!("For {small}, the interval for the bigger number is: [{low}, {up}]");
            // up_ind is not inclusive
            let up_ind = match x2_range.binary_search(&up) {
                Ok(i) => i+1,
                Err(i) => i,
            };
            // low_ind is inclusive
            let low_ind = match x2_range.binary_search(&low) {
                Ok(i) => i,
                Err(i) => i,
            };
            //x1_hits.push(u32::try_from(if low_ind >= up_ind {0} else {up_ind - low_ind}).unwrap());
            let elem_x2_range = &x2_range[low_ind..up_ind];
            //println!("x1: {x1}, x2 range: {elem_x2_range:?}");
            for x2 in elem_x2_range.choose_multiple(&mut rng, min(1000, up_ind - low_ind)) {
                let smooth = x1 * x2;
                //println!("smooth: {smooth}, x1: {x1}, x2: {x2}");
                assert!(smooth >= *a && smooth <= *b);
                /*match smooth_set.get_mut(&smooth) {
                    Some(v) => *v += 1,//.push((*x1, *x2)),
                    None => {smooth_set.insert(smooth, 1); ()},//vec![(*x1, *x2)]); ()},
                };*/
                new_smooths.push(smooth);
            }
        }
        new_smooths.sort_unstable();
        new_smooths
    };

    let new_smooths: Vec<BigUint> = thread::scope(|s| {
        let mut handles = vec![];
        for _ in 0..*NUM_THREADS{
            let h = s.spawn(move || smooth_combs());
            handles.push(h);
        }
        handles.into_iter()
            .map(|h| h.join().unwrap())
            .collect::<Vec<Vec<BigUint>>>()
            .concat()
    });
    //let x1_skew: Skewness = x1_hits.into_iter().map(f64::from).collect();
    //let mult_skew: Skewness = smooth_set.values().cloned().map(|x| f64::try_from(x).unwrap()).collect();
    //println!("x1s: {}, mean: {}, variance: {}, skewness: {}", x1_skew.len(), x1_skew.mean(), x1_skew.sample_variance(), x1_skew.skewness());
    //println!("multiplicities: {}, mean: {}, variance: {}, skewness: {}", smooth_set.values().len(), mult_skew.mean(), mult_skew.sample_variance(), mult_skew.skewness());
    //println!("{:?}", smooth_set);
    //let vals = smooth_set.values().cloned().collect::<Vec<usize>>();
    //let mut val_map = HashMap::new();
    //for val in vals.iter() {
    //    match val_map.get(val) {
    //        Some(x) => val_map.insert(val, x+1),
    //        None => val_map.insert(val, 1),
    //    };
    //}
    //println!("{:?}", val_map);
    smooths_per_level[level] = new_smooths;
    smooths_per_level[level].par_sort_unstable();
    smooths_per_level[level].dedup();
    let end_num = smooths_per_level[level].len();
    println!("Exiting from level {level} with {} nrs for interval: [{a}, {b}]", end_num);
}


#[allow(non_snake_case)]
fn main() {
    let args: Vec<String> = env::args().collect();
    let a: BigUint;
    let b: BigUint;
    let mut B: u64;
    if args.len() == 2 {
        let bit_size = usize::from_str_radix(&args[1], 10).unwrap();
        let mut rng = rand::thread_rng();
        let n: BigUint = <ThreadRng as RandPrime<BigUint>>::gen_prime_exact(&mut rng, bit_size, None).into();
        let four_n: BigUint = 4u32*n.clone();
        B = n.bits()*n.bits()/10;
        a = n.clone()-(four_n.clone()).sqrt()+BigUint::one();
        b = n.clone()+(four_n.clone()).sqrt()+BigUint::one();
    } else if args.len() == 3 {
        let n: BigUint = args[1].parse::<BigUint>().unwrap();
        let four_n: BigUint = 4u32*n.clone();
        B = u64::from_str_radix(&args[2], 10).unwrap();
        a = n.clone()-(four_n.clone()).sqrt()+BigUint::one();
        b = n.clone()+(four_n.clone()).sqrt()+BigUint::one();
    } else if args.len() == 4 {
        a = args[1].parse::<BigUint>().unwrap();
        b = args[2].parse::<BigUint>().unwrap();
        B = u64::from_str_radix(&args[3], 10).unwrap();
    } else {
        println!("Please supply 1, 2 or 3 arguments (bits of n; n, B; lower_bound, upper_bound, B)");
        return;
    }
    B = {match PRIMES.binary_search(&(B.into())) {
        Ok(i) => PRIMES[i],
        Err(i) => PRIMES[i],
    }}.try_into().unwrap();
    println!("Interval: [{a} to {b}], B: {B}");
    one_smooth_in(a, b, B);


}
