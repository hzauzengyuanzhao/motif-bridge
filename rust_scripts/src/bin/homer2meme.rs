//! homer2meme - Convert HOMER motif format to MEME motif format
//!
//! Build (Cargo, recommended - includes .gz support):
//!   cd rust_scripts && cargo build --release
//!   # Binaries: rust_scripts/target/release/meme2homer  homer2meme
//!
//! Install to PATH:
//!   cd rust_scripts && cargo install --path .
//!
//! Usage:
//!   homer2meme -i motifs.homer > motifs.meme
//!   homer2meme -i motifs.homer.gz > motifs.meme
//!   cat motifs.homer | homer2meme -i -

use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::process;
use flate2::read::MultiGzDecoder;

struct Args {
    input: String,
    extract: String,
    pseudocount: f64,
}

fn parse_args() -> Args {
    let argv: Vec<String> = env::args().collect();
    let mut input = String::new();
    let mut extract = String::new();
    let mut pseudocount = 0.01_f64;
    let mut i = 1usize;
    while i < argv.len() {
        match argv[i].as_str() {
            "-i" => {
                i += 1;
                if i >= argv.len() { eprintln!("Error: -i requires a value"); process::exit(1); }
                input = argv[i].clone();
            }
            "-e" => {
                i += 1;
                if i >= argv.len() { eprintln!("Error: -e requires a value"); process::exit(1); }
                extract = argv[i].clone();
            }
            "-a" => {
                i += 1;
                if i >= argv.len() { eprintln!("Error: -a requires a value"); process::exit(1); }
                let raw = argv[i].clone();
                pseudocount = raw.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid -a value: {}", raw);
                    process::exit(1);
                });
                if pseudocount <= 0.0 {
                    eprintln!("Error: -a must be > 0, got {}", pseudocount);
                    process::exit(1);
                }
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
    Args { input, extract, pseudocount }
}

fn usage() {
    eprintln!(r#"Usage: homer2meme -i <input_file> [OPTIONS]

Convert HOMER motif format to MEME format.

Options:
    -i <file>    Input HOMER motif file (.homer, .motif, or .gz), or '-' for stdin
    -e <string>  Extract only specified motif by id or description
    -a <float>   Pseudocount for log-odds -> probability conversion (default: 0.01)
    -h           Show this help

Examples:
    homer2meme -i motifs.homer > motifs.meme
    homer2meme -i motifs.homer.gz > motifs.meme
    homer2meme -i motifs.homer -e "CTCF/Jaspar"
    cat motifs.homer | homer2meme -i -
"#);
}

fn is_logodds(row: &[f64]) -> bool {
    let s: f64 = row.iter().sum();
    !(0.98..=1.02).contains(&s)
}

fn logodds_to_prob(row: &[f64], pseudocount: f64) -> Vec<f64> {
    let background = 0.25_f64;
    let raw: Vec<f64> = row.iter().map(|&v| 2.0_f64.powf(v) * background).collect();
    let total: f64 = raw.iter().sum::<f64>() + pseudocount * raw.len() as f64;
    raw.iter().map(|&v| (v + pseudocount) / total).collect()
}

fn print_meme_header() {
    println!("MEME version 4");
    println!();
    println!("ALPHABET= ACGT");
    println!();
    println!("strands: + -");
    println!();
    println!("Background letter frequencies");
    println!("A 0.25 C 0.25 G 0.25 T 0.25");
    println!();
}

fn print_meme_motif(id: &str, desc: &str, matrix: &[Vec<f64>]) {
    let width = matrix.len();
    println!("MOTIF {} {}", id, desc);
    println!();
    println!("letter-probability matrix: alength= 4 w= {} nsites= 20 E= 0", width);
    for row in matrix {
        let line: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
        println!("  {}", line.join("  "));
    }
    println!();
}

fn parse_and_convert<R: BufRead>(reader: R, args: &Args) {
    let mut header_printed = false;
    let mut in_motif = false;
    let mut motif_id = String::new();
    let mut description = String::new();
    let mut matrix: Vec<Vec<f64>> = Vec::new();

    for line_result in reader.lines() {
        let line = line_result.expect("Failed to read line");
        let trimmed = line.trim();

        if trimmed.is_empty() { continue; }

        if let Some(rest) = trimmed.strip_prefix('>') {
            if in_motif && !matrix.is_empty() {
                if !header_printed { print_meme_header(); header_printed = true; }
                print_meme_motif(&motif_id, &description, &matrix);
            }
            matrix.clear();

            let parts: Vec<&str> = rest.splitn(6, '\t').collect();
            let mid = parts.first().copied().unwrap_or("motif").to_string();
            let raw_desc = parts.get(1).copied().unwrap_or("").to_string();
            let final_desc = if raw_desc.is_empty() { mid.clone() } else { raw_desc };

            if !args.extract.is_empty() && mid != args.extract && final_desc != args.extract {
                in_motif = false;
                continue;
            }
            motif_id = mid;
            description = final_desc;
            in_motif = true;
            continue;
        }

        if !in_motif { continue; }

        let row_result: Result<Vec<f64>, _> = trimmed.split_whitespace()
            .map(|s| s.parse::<f64>()).collect();
        if let Ok(mut row) = row_result {
            if !row.is_empty() {
                if row.len() != 4 {
                    eprintln!("Warning: skipping malformed matrix row with {} cols (expected 4): {}", row.len(), trimmed);
                    continue;
                }
                if is_logodds(&row) { row = logodds_to_prob(&row, args.pseudocount); }
                matrix.push(row);
            }
        }
    }

    if in_motif && !matrix.is_empty() {
        if !header_printed { print_meme_header(); }
        print_meme_motif(&motif_id, &description, &matrix);
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
