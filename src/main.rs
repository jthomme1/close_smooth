use crate::smooths::Smooths;
use std::vec::Vec;
use once_cell::sync::Lazy;
use primal;
use std::env;
use integer_sqrt::IntegerSquareRoot;
use rayon::prelude::*;
use std::collections::HashMap;
//use rand::seq::SliceRandom;
use std::cmp::{/*min, */max};
//use average::Skewness;

pub mod composite;
pub mod smooths;

// this should suffice for now
static PRIME_BOUND: usize = 1<<16;
static PRIMES: Lazy<Vec<u128>> = Lazy::new(||
                                           primal::Sieve::new(PRIME_BOUND)
                                           .primes_from(0)
                                           .collect::<Vec<usize>>()
                                           .into_iter()
                                           .map(|x| u128::try_from(x).unwrap())
                                           .collect());

#[allow(non_snake_case)]
fn smooths_in(a: u128, b: u128, B: u128) -> Vec<u128> {
    println!("Recursing into interval: [{a}, {b}]");
    assert!(a < b);
    // base case
    if a < B*B {
        println!("Entering base case for interval: [{a}, {b}]");
        let ind: usize = PRIMES.binary_search(&B).unwrap();
        let smooths = Smooths::new(a, b, ind);
        println!("Returning {} nrs from base case for interval: [{a}, {b}]", smooths.len());
        return smooths.to_u128();
    }
    // rounding down is fine, because rounding up would result in a product outside the bounds
    // also, rounding down we don't lose anything because we're considering integers
    let ub_x1 = b.integer_sqrt();
    //println!("upper_bound_smaller_factor: {ub_x1}");
    // we need to round up because otherwise we could get a product outside of the interval
    // also, we don't lose anything rounding up because we're considering integers
    let mut lb_x2 = a.integer_sqrt();
    if lb_x2 * lb_x2 < b {
        lb_x2 += 1;
    }
    //println!("lower_bound_bigger_factor: {lb_x2}");

    let ub_x2 = (B*b).integer_sqrt();
    let mut lb_x1 = (a/B).integer_sqrt();
    if lb_x1 * lb_x1 < a/B {
        lb_x1 += 1;
    }
    let smaller_smooths = smooths_in(lb_x1, ub_x2, B);

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

    let mut smooth_set = HashMap::new();

    //let mut rng = &mut rand::thread_rng();

    //let mut x1_hits: Vec<u32> = vec![];

    for x1 in x1_range {//.choose_multiple(&mut rng, (ub_x1_ind - lb_x1_ind)/usize::try_from(B.integer_sqrt()).unwrap()) {
        let up = b/x1;
        let mut low = a/x1;
        if low * x1 < a {
            low += 1;
        }
        low = max(low, *x1);
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
        for x2 in elem_x2_range {//.choose_multiple(&mut rng, min(100000, up_ind - low_ind)) {
            let smooth = x1 * x2;
            //println!("smooth: {smooth}, x1: {x1}, x2: {x2}");
            assert!(smooth >= a && smooth <= b);
            match smooth_set.get_mut(&smooth) {
                Some(v) => *v += 1,//.push((*x1, *x2)),
                None => {smooth_set.insert(smooth, 1); ()},//vec![(*x1, *x2)]); ()},
            };
        }
    }
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
    let mut smooths = smooth_set.keys().into_iter().map(|x| *x).collect::<Vec<u128>>();
    smooths.par_sort_unstable();
    println!("Returning {} nrs for interval: [{a}, {b}]", smooths.len());
    smooths
}


#[allow(non_snake_case)]
fn main() {
    let args: Vec<String> = env::args().collect();
    let a: u128;
    let b: u128;
    let B: u128;
    if args.len() == 3 {
        let n = u128::from_str_radix(&args[1], 10).unwrap();
        B = u128::from_str_radix(&args[2], 10).unwrap();
        a = n-(4*n).integer_sqrt()+1;
        b = n+(4*n).integer_sqrt()+1;
    } else if args.len() == 4 {
        a = u128::from_str_radix(&args[1], 10).unwrap();
        b = u128::from_str_radix(&args[2], 10).unwrap();
        B = u128::from_str_radix(&args[3], 10).unwrap();
    } else {
        println!("Please supply 2 or 3 arguments (either n b or lower_bound upper_bound b)");
        return;
    }
    println!("{} primes generated.", PRIMES.len());
    println!("Interval: [{a} to {b}]");
    println!("{:?}", smooths_in(a, b, B));


}
