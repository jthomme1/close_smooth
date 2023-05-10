use crate::smooths::Smooths;
use std::vec::Vec;
use once_cell::sync::Lazy;
use primal;
use std::env;
use integer_sqrt::IntegerSquareRoot;
use rayon::prelude::*;
use std::collections::HashSet;

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

fn smooths_in(a: u128, b: u128, B: u128) -> Vec<u128> {
    println!("Recursing into interval: [{a}, {b}]");
    assert!(a < b);
    // base case
    if a < B*B {
        println!("Entering base case for interval: [{a}, {b}]");
        let ind: usize = PRIMES.binary_search(&B).unwrap();
        let smooths = Smooths::new(a, b, ind);
        println!("Returning {} nrs from base case for interval: [{a}, {b}]", smooths.smooths.len());
        return smooths.smooths;
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
    let mut lb_x1= (a/B).integer_sqrt();
    if lb_x1 * lb_x1 < a/B {
        lb_x1 += 1;
    }
    let smaller_smooths = smooths_in(lb_x1, ub_x2, B);

    let ub_x2_ind = smaller_smooths.len()-1;
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

    let mut smooth_set = HashSet::new();

    for x1 in x1_range {
        let up = b/x1;
        let mut low = a/x1;
        if low * x1 < a {
            low += 1;
        }
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
        //println!("x2_range starting from {low_ind} to {up_ind}");
        for x2_ind in low_ind..up_ind {
            let x2 = x2_range[x2_ind];
            let smooth = x1 * x2;
            //println!("smooth: {smooth}, x1: {x1}, x2: {x2}");
            assert!(smooth >= a && smooth <= b);
            smooth_set.insert(smooth);
        }
    }
    let mut smooths = smooth_set.into_iter().collect::<Vec<u128>>();
    smooths.par_sort_unstable();
    println!("Returning {} nrs for interval: [{a}, {b}]", smooths.len());
    smooths
}

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
