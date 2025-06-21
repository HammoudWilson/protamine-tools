//! This program reads insert endpoints from standard input, stores the 
//! counts in a coordinate map, and processes them by chromosome.
//! 
//! The input is expected to be in pre-sorted BED3 format. The inserts  
//! in the input stream may be from one sample or a set of samples as 
//! dictated by the user's needs.
//!  
//! Insert endpoints are counted in two maps, one for starts, one for ends.
//! After loading a chromosome, it is interrogated for insert patterns
//! consistent with two positioned nucleosomes flanking an intervening
//! nucleosome-free region (NFR), i.e., gap. A series of allowed sizes for  
//! the central NFR are interrogated and the best score and size are kept for 
//! every other position on the chromosome (to help with computational speed).
//! 
//! Scoring is achieved by moving a set of masks, one for each gap length, 
//! over the insert endpoint count maps. The masks multiply endpoint counts 
//! by +1 in regions expected to have insert starts or ends, by -1 in nucleocome 
//! regions expected to be flanked by (but not to contain) insert endpoints, and 
//! by 0 in regions that do not contribute to the score. Specifically, insert
//! ends are not scored on the leftmost part of the mask, and insert starts
//! are not scored on the rightmost part of the mask, to focus scores on inserts
//! between and flanking the positioned dinucleosomes.
//! 
//! Regions without positioned nucleosomes are expected to yield low insert counts,
//! reflecting a lower overall accessibility of the chromatin in those regions,
//! with scores that range around zero given a random likelihood of non-positioned
//! insert endpoints being in positive vs. negatively scored positions.
//! 
//! In constrast, regions with positioned nucleosomes are expected to yield
//! high insert counts, reflecting a higher overall accessibility of the chromatin 
//! in those regions. Moreover, scores are expected to be strongly positive
//! when a scoring mask of the appropriate size is centered on the NFR and thus the 
//! counts and mask are in phase and reinforcing, and negative when the scoring mask 
//! is centered on a positioned nucleosome or of the wrong size, leading to 
//! peaks in the score matrix within central NFRs.

use std::io::stdin;
use ab_initio::InsertMap;

// load inserts into endoint buffers
// process intermittently by chromosome
fn main() {

    // instantiate maps to hold insert endpoint counts
    let mut insert_map = InsertMap::new();

    // open a file stream from STDIN
    // input column format is pre-sorted chrom, start0, end1, i.e., BED3
    let input = stdin().lines();

    // read the file line by line, adding each insert endpoint to the map
    // the insert map handles data processing by chromosome
    for line in input {
        let line = line.unwrap();
        let parts: Vec<&str> = line.split('\t').collect();
        let chrom = parts[0].to_string();
        let start0: usize = parts[1].parse().unwrap();
        let end1:   usize = parts[2].parse().unwrap();
        insert_map.add_insert(chrom, start0, end1);
    }

    // process the last chromosome
    insert_map.process_chrom();
}
