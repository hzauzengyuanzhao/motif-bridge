//! meme2homer - Convert MEME motif format to HOMER motif format
//!
//! Build (Cargo, recommended - includes .gz support):
//!   cd rust_scripts && cargo build --release
//!   # Binaries: rust_scripts/target/release/meme2homer  homer2meme
//!
//! Install to PATH:
//!   cd rust_scripts && cargo install --path .
//!
//! Usage:
//!   meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer
//!   meme2homer -i motifs.meme.gz -j JASPAR2026 > motifs.homer
//!   zcat motifs.meme.gz | meme2homer -i -

use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::process;
use flate2::read::MultiGzDecoder;

// ---------------------------------------------------------------------------
// Argument parsing
// ---------------------------------------------------------------------------

struct Args {
    input: String,
    db: String,
    motif_name: String,
    extract: String,
    bg: f64,
    t_offset: f64,
}

macro_rules! next_arg {
    ($argv:expr, $i:expr, $flag:expr) => {{
        $i += 1;
        if $i >= $argv.len() {
            eprintln!("Error: {} requires an argument.", $flag);
            process::exit(1);
        }
        $argv[$i].clone()
    }};
}

fn parse_args() -> Args {
    let argv: Vec<String> = env::args().collect();
    let mut input = String::new();
    let mut db = String::from("NA");
    let mut motif_name = String::new();
    let mut extract = String::new();
    let mut bg = 0.25_f64;
    let mut t_offset = 4.0_f64;
    let mut i = 1usize;
    while i < argv.len() {
        match argv[i].as_str() {
            "-i" => { input = next_arg!(argv, i, "-i"); }
            "-j" => { db = next_arg!(argv, i, "-j"); }
            "-k" => { motif_name = next_arg!(argv, i, "-k"); }
            "-e" => { extract = next_arg!(argv, i, "-e"); }
            "-b" => {
                let val = next_arg!(argv, i, "-b");
                bg = val.parse().unwrap_or_else(|_| {
                    eprintln!("Error: invalid -b value: {}", val);
                    process::exit(1);
                });
                if bg <= 0.0 || bg > 1.0 {
                    eprintln!("Error: -b must be in (0, 1], got {}", bg);
                    process::exit(1);
                }
            }
            "-t" => {
                let val = next_arg!(argv, i, "-t");
                t_offset = val.parse().unwrap_or_else(|_| {
                    eprintln!("Error: invalid -t value: {}", val);
                    process::exit(1);
                });
            }
            "-h" | "--help" => { usage(); process::exit(0); }
            other => { eprintln!("Unknown option: {}", other); usage(); process::exit(1); }
        }
        i += 1;
    }
    if input.is_empty() {
        eprintln!("Error: -i <input_file> is required.");
        usage();
        process::exit(1);
    }
    Args { input, db, motif_name, extract, bg, t_offset }
}

fn usage() {
    eprintln!(r#"Usage: meme2homer -i <input_file> [OPTIONS]

Convert MEME format to HOMER motif format.

Options:
    -i <file>    Input MEME file (.meme or .meme.gz), or '-' for stdin
    -j <string>  Database name (default: NA)
    -k <string>  Override motif name
    -e <string>  Extract only specified motif by id or name
    -b <float>   Background probability in (0, 1] (default: 0.25)
    -t <float>   Threshold offset in log2 bits (default: 4.0)
    -h           Show this help

Examples:
    meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer
    meme2homer -i motifs.meme.gz -j JASPAR2026 > motifs.homer
    meme2homer -i motifs.meme -b 0.25 -t 6
"#);
}

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

struct Motif {
    id: String,
    description: String,
    matrix: Vec<Vec<f64>>,
}

impl Motif {
    fn calculate_score(&self, bg: f64, t_offset: f64) -> f64 {
        let raw: f64 = self.matrix.iter().map(|row| {
            let max_p = row.iter().cloned().fold(0.0_f64, f64::max);
            if max_p > 0.0 { (max_p / bg).log2() } else { 0.0 }
        }).sum();
        (raw - t_offset).max(0.0)
    }

    fn print_homer(&self, bg: f64, t_offset: f64) {
        let score = self.calculate_score(bg, t_offset);
        println!(">{}\t{}\t{:.6}\t0\t0\t0", self.id, self.description, score);
        for row in &self.matrix {
            let line: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
            println!("{}", line.join("\t"));
        }
    }
}

// ---------------------------------------------------------------------------
// Core parsing
// ---------------------------------------------------------------------------

fn parse_and_convert<R: BufRead>(reader: R, args: &Args) {
    let mut in_motif  = false;
    let mut in_matrix = false;
    let mut motif_id  = String::new();
    let mut description = String::new();
    let mut matrix: Vec<Vec<f64>> = Vec::new();

    for line_result in reader.lines() {
        let line = line_result.expect("Failed to read line");
        let trimmed = line.trim();

        if let Some(rest) = trimmed.strip_prefix("MOTIF") {
            if in_motif && !matrix.is_empty() {
                Motif {
                    id: motif_id.clone(),
                    description: description.clone(),
                    matrix: std::mem::take(&mut matrix),
                }.print_homer(args.bg, args.t_offset);
            }
            in_matrix = false;

            let parts: Vec<&str> = rest.split_whitespace().collect();
            if parts.is_empty() { in_motif = false; matrix.clear(); continue; }
            let id = parts[0].to_string();
            let original_name = if parts.len() > 1 { parts[1..].join(" ") } else { id.clone() };

            if !args.extract.is_empty() && id != args.extract && original_name != args.extract {
                in_motif = false;
                continue;
            }
            motif_id = id;
            description = if !args.motif_name.is_empty() {
                format!("{}/{}", args.motif_name, args.db)
            } else {
                format!("{}/{}", original_name, args.db)
            };
            in_motif = true;
            continue;
        }

        if !in_motif { continue; }

        if trimmed.starts_with("URL ") || trimmed == "URL" { continue; }

        if trimmed.starts_with("//") {
            if !matrix.is_empty() {
                Motif {
                    id: motif_id.clone(),
                    description: description.clone(),
                    matrix: std::mem::take(&mut matrix),
                }.print_homer(args.bg, args.t_offset);
            }
            in_motif  = false;
            in_matrix = false;
            continue;
        }

        if trimmed.starts_with("letter-probability matrix:") {
            in_matrix = true;
            continue;
        }

        if in_matrix {
            let first = trimmed.chars().next();
            if !matches!(first, Some('0'..='9') | Some('.')) { continue; }
            let row: Result<Vec<f64>, _> = trimmed
                .split_whitespace()
                .map(|s| s.parse::<f64>())
                .collect();
            if let Ok(values) = row {
                if values.len() == 4 {
                    matrix.push(values);
                } else if !values.is_empty() {
                    eprintln!("Warning: skipping malformed matrix row with {} cols (expected 4)", values.len());
                }
            }
        }
    }

    if in_motif && !matrix.is_empty() {
        Motif { id: motif_id, description, matrix }.print_homer(args.bg, args.t_offset);
    }
}

fn main() {
    let args = parse_args();
    if args.input == "-" {
        parse_and_convert(BufReader::new(io::stdin().lock()), &args);
    } else if args.input.ends_with(".gz") {
        let file = File::open(&args.input)
            .unwrap_or_else(|e| { eprintln!("Cannot open {}: {}", args.input, e); process::exit(1); });
        parse_and_convert(BufReader::new(MultiGzDecoder::new(file)), &args);
    } else {
        let file = File::open(&args.input)
            .unwrap_or_else(|e| { eprintln!("Cannot open {}: {}", args.input, e); process::exit(1); });
        parse_and_convert(BufReader::new(file), &args);
    }
}
