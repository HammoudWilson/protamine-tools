//! This program reads insert endpoints from stage-level insert_bgz files, 
//! stores the counts in a coordinate map, and searches for positioned 
//! dinucleosomes ab initio, i.e., with no other guidance than the inserts 
//! themselves. Work is done for one chromosome provided with stage 
//! information as command line arguments.
//!  
//! Insert endpoints are counted in two maps per stage, one for insert endpoints, 
//! one for insert centers of fragments of sizes consistent with them spanning a 
//! single nucleosome. After loading the chromosome, it is interrogated for insert 
//! patterns consistent with two positioned nucleosomes flanking an intervening
//! nucleosome-free region (NFR), i.e., gap. A series of allowed sizes for 
//! the central NFR are interrogated and the best score and size are kept for 
//! every third position on the chromosome (for speed).
//! 
//! Scoring is achieved by moving a set of masks, one for each gap length, 
//! over the insert count maps. The masks retain and sum endpoint counts in NFR 
//! regions expected to have insert starts or ends, and insert center counts 
//! within the putative spans of the two positioned nucleosomes. Insert centers
//! are counted twice, once for each endpoint, for consistent weighting. Increasing
//! gap lengths are linearly penalized to counteract expansion to larger and larger
//! gaps simply because they capture more inserts, thus helping ensure that 
//! nucleosomes are taken into account in the scoring.
//! 
//! Regions without positioned nucleosomes are expected to yield low insert counts,
//! reflecting a lower overall accessibility of the chromatin in those regions.
//! 
//! In constrast, regions with positioned nucleosomes are expected to yield
//! high insert counts, reflecting a higher overall accessibility. Moreover, scores 
//! are expected to be maximal when a scoring mask of the appropriate size is 
//! centered on the NFR and thus the counts and mask are in phase and reinforcing, 
//! i.e., such that nucleosomal inserts get scored at both their endpoints and centers.
//! 
//! Searhing is done independently by stage. The process is the same for all stages,
//! except that the minimum scoring threshold is scaled to the total number of inserts 
//! for each stage such that stages with more data must have a higher score to call 
//! a region as a positioned dinucleosome.
//! 
//! Finally, for each candidate dinucleosome discovered for each stage, a score
//! is calculated for the same span in all other stages for cross-comparison.
//! These and other final scores are reported without gap extension penalties,
//! which were only used to optimize the initial dinucleosome search.

use std::env;
use std::collections::HashMap;
use ab_initio::InsertMap;

// load and process data
fn main() -> Result<(), Box<dyn std::error::Error>> {

    // read command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: ab_initio <chrom> <chrom_length> <stages> <stage_insert_files>");
        std::process::exit(1);
    }
    let chrom = &args[1];
    let chrom_length: usize = args[2].parse()?;
    let stages = &args[3];
    let stage_insert_files = &args[4];

    // parse stages and stage insert file paths
    let stages: Vec<&str> = stages.split(',').collect();
    let stage_insert_file_paths: Vec<&str> = stage_insert_files.split(',').collect();
    let mut stage_insert_files: HashMap<String, String> = HashMap::new();
    for (stage, path) in stages.iter().zip(stage_insert_file_paths.iter()) {
        stage_insert_files.insert(stage.to_string(), path.to_string());
    }

    // instantiate maps to hold insert endpoint counts
    let mut insert_map = InsertMap::new(chrom_length, stages);

    // load inserts from files for all stages
    insert_map.load_inserts_by_stage(&chrom, &stage_insert_files)?;

    // use insert counts to set stage-specific score thresholds
    insert_map.set_min_scores();

    // find positioned dinucleosomes independently by stage, print to STDOUT
    insert_map.find_dinucs_by_stage(&chrom);

    Ok(())
}
