use std::io::stdin;
use std::env;


// define constants established from biological principles
const NUCLEOSOME_SIZE: u16 = 147;
const MAX_GAP_LENGTH:  u16 = 145; // 2bp less than the nucleosome size

// define a struct to hold and pass operating parameters
struct Config {
    min_insert_size: u16,
    max_mapped_insert_size: u16,
    n_insert_sizes: usize, // number of insert sizes to consider
    gap_lengths: Vec<u16>,
}
impl Config {
    fn new() -> Self {

        // set the expected range of insert sizes
        // MAX_MAPPED_INSERT_SIZE may be smaller than MAX_INSERT_SIZE (and usually is)
        let min_insert_size: u16 = env::var("MIN_INSERT_SIZE").unwrap().parse().expect("MIN_INSERT_SIZE must be u16");
        let max_mapped_insert_size: u16 = env::var("MAX_MAPPED_INSERT_SIZE").unwrap().parse().expect("MAX_MAPPED_INSERT_SIZE must be u16");

        // set the predetermined gap lengths, i.e., the allowed size of the nucleosome-free regions
        // smallest allowed sizes it 10 bp larger than the minimum insert size to allow some potential for sub-nucleosomal fragments
        let gap_lengths: Vec<u16> = ((min_insert_size + 10)..=MAX_GAP_LENGTH).collect();
        Config {
            min_insert_size,
            max_mapped_insert_size,
            n_insert_sizes: (max_mapped_insert_size - min_insert_size + 1) as usize,
            gap_lengths,
        }
    }
}

fn main() {
    let config = Config::new();

    // TODO: presumably should just use ndarray now that we have a large number of insert size rows

    // instantiate and allocate a set of four large vectors
    // each vector represents one insert size from MIN_INSERT_SIZE..=MAX_MAPPED_INSERT_SIZE
    // each slot in a vector represent one position on the working chromosome
    // due to memory constraints, it is not possible hold an entire genome or chromosome in memory
    // therefore, we use a vector of a million slots each and rotate it as we read the input
    // include an extra max_mapped_insert_size slots to hold the overrunning inserts
    let mut insert_map: Vec<Vec<u32>> = Vec::with_capacity(config.n_insert_sizes);
    for _ in 0..config.n_insert_sizes {
        insert_map.push(vec![0; 1e6 as usize + config.max_mapped_insert_size as usize]);
    }

    // open a file stream from STDIN
    // input column format is  chrom, start0, end1, i.e., BED3 sorted by chrom, start0
    let mut input = stdin().lines().peekable();

    // peek at the first line to determine the working chromosome
    let first_line = input.peek().unwrap().unwrap();
    let mut working_chrom: String = first_line.split('\t').next().unwrap().to_string();
    

}


fn load_inserts_chunk(&mut input) {


    // read the file line by line
    for line in input {
        let parts: Vec<&str> = line.unwrap().split('\t').collect();
        let chrom = parts[0].to_string();
        let start: u32    = parts[1].parse().unwrap();
        let end: u32      = parts[2].parse().unwrap();

        // if the chromosome has changed, process the previous chromosome
        if chrom != working_chrom {
            process_chrom();
            // reset the working chromosome
            working_chrom = chrom;
        }

        // increment
    }
}

fn process_chrom(){

}
