use std::{io::{BufReader, Write, BufRead}, env, fs::File};

fn main() {
    let args = env::args().collect::<Vec<_>>();
    let buf = BufReader::new(File::open(&args[1]).unwrap());
    let mut out = File::create(&args[2]).unwrap();
    
    let mut cur_chr = String::new();
    let mut pos: u64 = 1;
    for line in buf.lines() {
        let parts = line.unwrap().split(" ").map(|s| s.to_string()).collect::<Vec<_>>();
        
        // bed is 0-indexed and wiggle is 1-indexed, so we need to add 1 here
        let start = parts[1].parse::<u64>().unwrap() + 1;
        // bed is end-inclusive and wiggle is end-inclusive, so no need to add 1 here
        let mut end = parts[2].parse::<u64>().unwrap();
        // for insertion we will have start > end after above transformation - account for that
        if start > end {
            end += 1;
        }
        
        
        if cur_chr.is_empty() || cur_chr != parts[0] {
            write!(out, "fixedStep  chrom={} start=1 step=1\n", parts[0]).unwrap();
            
            pos = 1;
            cur_chr = parts[0].clone();
        }
        
        while pos <= end {
            if pos < start {
                write!(out, "0\n").unwrap();
            }
            else {
                write!(out, "{}\n", parts[7]).unwrap();
            }
            
            pos += 1;
        }
    }
}