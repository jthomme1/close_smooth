use crate::smooths::Smooths;
use std::vec::Vec;
use once_cell::sync::Lazy;
use primal;
use std::env;
use integer_sqrt::IntegerSquareRoot;

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

fn main() {
    let args: Vec<String> = env::args().collect();
    let interval_upper_bound: u128;
    let interval_lower_bound: u128;
    let b: u128;
    if args.len() == 3 {
        let n = u128::from_str_radix(&args[1], 10).unwrap();
        b = u128::from_str_radix(&args[2], 10).unwrap();
        interval_upper_bound = n+(4*n).integer_sqrt()+1;
        interval_lower_bound = n-(4*n).integer_sqrt()+1;
    } else if args.len() == 4 {
        interval_lower_bound = u128::from_str_radix(&args[1], 10).unwrap();
        interval_upper_bound = u128::from_str_radix(&args[2], 10).unwrap();
        b = u128::from_str_radix(&args[3], 10).unwrap();
    } else {
        println!("Please supply 2 or 3 arguments (either n b or lower_bound upper_bound b)");
        return;
    }
    println!("{} primes generated.", PRIMES.len());
    println!("Interval: [{interval_lower_bound} to {interval_upper_bound}]");
    let ind: usize = PRIMES.binary_search(&b).unwrap();

    let smooths_upper_bound = (b*b*interval_upper_bound).integer_sqrt();
    let smooths_lower_bound = (interval_lower_bound/(b*b)).integer_sqrt();
    println!("Considering smooths in [{smooths_lower_bound}, {smooths_upper_bound}]");

    let smooths = Smooths::new(smooths_upper_bound, ind);

    let mut upper_bound_smaller_factor = interval_upper_bound.integer_sqrt();
    if upper_bound_smaller_factor * upper_bound_smaller_factor < interval_upper_bound {
        upper_bound_smaller_factor += 1;
    }
    //println!("upper_bound_smaller_factor: {upper_bound_smaller_factor }");
    let lower_bound_bigger_factor = interval_lower_bound.integer_sqrt();
    //println!("lower_bound_bigger_factor: {lower_bound_bigger_factor}");

    let sub_ind = smooths.smooths.len()-1;
    let slb_ind = match smooths.smooths.binary_search(&smooths_lower_bound) {
        Ok(i) => i,
        Err(i) => i,
    };

    let ub_small = smooths.find_ind_le(upper_bound_smaller_factor).unwrap();
    let lb_big = match smooths.smooths.binary_search(&lower_bound_bigger_factor) {
        Ok(i) => i,
        Err(i) => i,
    };


    //println!("Lower part [{smooths_lower_bound}, {upper_bound_smaller_factor}]");
    //println!("The first smooth number in there: {}", smooths.smooths[slb_ind]);
    //println!("The last smooth number in there: {}", smooths.smooths[ub_small]);
    //println!("Upper part [{lower_bound_bigger_factor}, {smooths_upper_bound}]");
    //println!("The first smooth number in there: {}", smooths.smooths[lb_big]);
    //println!("The last smooth number in there: {}", smooths.smooths[sub_ind]);
    let upper_part = &smooths.smooths[lb_big..=sub_ind];
    let lower_part = &smooths.smooths[slb_ind..=ub_small];

    for small in lower_part {
        let up = interval_upper_bound/small;
        let low = interval_lower_bound/small;
        //println!("For {small}, the interval for the bigger number is: [{low}, {up}]");
        // up_ind is not inclusive
        let up_ind = match upper_part.binary_search(&up) {
            Ok(i) => {
                println!("Found {b}-smooth number {} = {}*{}", upper_part[i]*small, upper_part[i], small);
                return;
            },
            Err(i) => i,
        };
        // low_ind is inclusive
        let low_ind = match upper_part.binary_search(&low) {
            Ok(i) => {
                println!("Found {b}-smooth number {} = {}*{}", upper_part[i]*small, upper_part[i], small);
                return;
            },
            Err(i) => i,
        };
        //println!("Checking interval from {low_ind} to {up_ind}");
        for big in low_ind..up_ind {
            println!("Found {b}-smooth number {} = {}*{}", upper_part[big]*small, upper_part[big], small);
            return;
        }
    }
    println!("No {b}-smooth number in [{interval_lower_bound}, {interval_upper_bound}]");
}
