use std::{io::{BufReader, BufRead, Write}, fs::File, env, collections::{HashSet, HashMap}};


struct Line {
    chromosome: String,
    start: u64,
    end: u64,
    id: String,
    variety: String,
    reference: String,
    alts: HashSet<String>,
    group: u8,
    severity: String,
    severity_rank: u8
}

impl Line {
    fn compatible(&self, other: &Line) -> bool {
        self.chromosome == other.chromosome &&
        self.start == other.start &&
        self.variety == other.variety &&
        self.reference == other.reference
    }
    
    fn redundant(&self, other: &Line) -> bool {
        self.id == other.id && 
        self.variety != other.variety 
    }
    
    fn merge(&mut self, mut more: Option<Line>, out: &mut File) {
        // merge new line if not empty (and a Line instance)
        if let Some(ref mut more) = more {
            if self.compatible(more) {
                self.alts.extend(more.alts.clone());
                if more.severity_rank < self.severity_rank {
                    self.end = more.end;
                    self.variety = more.variety.clone();
                    self.group = more.group;
                    self.severity = more.severity.to_string();
                    self.severity_rank = more.severity_rank;
                }
                return;
            }
            
            // if somehow with same rs id we have different variety of variant we skip the later ones
            if self.redundant(more) {
                return
            }
        }
        
        // if new Line is not compatible with the current one it is a new variant
        // print out the current line
        if self.alts.len() > 0 {
            let alts = Vec::from_iter(self.alts.clone());
            write!(out, "{} {} {} {} {} {} {} {} {}\n",
                self.chromosome, self.start, self.end,
                self.id, self.variety, self.reference,
                alts.join("/"), self.group, self.severity
            ).unwrap();
        }
        
        // make the new Line as the current one
        if let Some(more) = more {
            *self = more;
        }
    }
}

fn main() {
    // read cli arguments
    let args = env::args().collect::<Vec<_>>();
    let reader = BufReader::new(File::open(&args[1]).unwrap());
    let mut out = File::create(&args[2]).unwrap();
    let json = std::fs::read_to_string(&args[3]).unwrap();
    
    let severity = {
        serde_json::from_str::<HashMap<String, String>>(&json).unwrap()
    };
    
    // dummy initial value for the object to read line from vcf
    // this line is guranteed to not get printed as alt.len == 0
    let mut lines = Line {
        chromosome: "".to_string(),
        start: 1,
        end: 0,
        id: "".to_string(),
        variety: "".to_string(),
        reference: "".to_string(),
        alts: HashSet::new(),
        group: 0,
        severity: "".to_string(),
        severity_rank: 255
    };
    for line in reader.lines() {
        
        let parts = line.unwrap().split(" ").map(|s| s.to_string()).collect::<Vec<_>>();
        let severity_rank = (*severity.get(&parts[8]).unwrap_or(&String::from("0"))).parse::<u8>().unwrap();
        
        let more = Line {
            chromosome: parts[0].clone(),
            start: parts[1].parse::<u64>().unwrap(),
            end: parts[2].parse::<u64>().unwrap(),
            id: parts[3].clone(),
            variety: parts[4].clone(),
            reference: parts[5].clone(),
            alts: HashSet::from_iter(parts[6].split("/").map(|s| s.to_string())),
            group: parts[7].parse::<u8>().unwrap(),
            severity: parts[8].clone(),
            severity_rank: severity_rank
        };
                
        lines.merge(Some(more), &mut out);
        
    }
    
    lines.merge(None, &mut out);
}