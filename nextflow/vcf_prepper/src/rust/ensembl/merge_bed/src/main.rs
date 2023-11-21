use std::{io::{BufReader, BufRead, Write}, fs::File, env, collections::{HashSet}};

fn main() {
    // read cli arguments
    let args = env::args().collect::<Vec<_>>();
    let mut out = File::create(&args[1]).unwrap();

    let mut current_ids = HashSet::new();
    let mut file_counter = 2;
    while file_counter < args.len() {
        let reader = BufReader::new(File::open(&args[file_counter]).unwrap());
        for line in reader.lines() {
            
            let parts = line.unwrap().split(" ").map(|s| s.to_string()).collect::<Vec<_>>();
            
            if !current_ids.contains(&parts[3]) {
                write!(out, "{} {} {} {} {} {} {} {} {}\n",
                    parts[0], 
                    parts[1],
                    parts[2],
                    parts[3], 
                    parts[4], 
                    parts[5],
                    parts[6],
                    parts[7],
                    parts[8]
                ).unwrap();
                
                current_ids.insert(parts[3].clone());
            }
        }
        
        file_counter += 1;
    }
}