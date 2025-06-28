// dependencies
use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead, Write, BufWriter, stdout};
use std::collections::HashMap;
use rayon::prelude::*;
use rayon::join;
use wide::i16x8; // must use i16x8 as u16x8 doesn't support all desired functions with stable Rust

// define the dinucleosome scoring algorithm, i.e., scoring masks
const FLANK_LEN:    i16   = 75;  // the length of the NFRs flanking the positioned nucleosomes
const NUC_LEN:      i16   = 147; // nucleosome span in bp, i16 to match the counts and masks
const MAX_NUC_LEN:  usize = 320; // upper size limit of fragments counted as nuc_centers (WAS UNLIMITED)
const MIN_GAP_LEN:  i16   = 51;  // minimum allowed central NFR length
const MAX_GAP_LEN:  i16   = 351; // maximum allowed central NFR length
const MAX_J0:       usize = MAX_GAP_LEN as usize; // same as  MAX_GAP_LEN, but as usize for indexing
const GAP_LEN_STEP: usize = 4;   // step size for the central NFR gap lengths, every other odd number is used
const MIN_HALF_GAP: i16   = (MIN_GAP_LEN - 1) / 2; // used for dinuc collision avoidance
const BLOCK_HALF_LEN: i16 = NUC_LEN + MIN_HALF_GAP; 
const COLLISION_TOLERANCE: i16 = 10; // some slop to allow for imprecise nucleosome positioning

// define the observed inserts and chromosome coordinate search
const CHROM_POS_STEP: usize = 3; // step size for the chromosome positions, only every third position is processed to reduce loop iterations

// define parallelization implementation
const MIN_PAR_ITER_LEN:     usize = 25_000;   // minimum length of a parallel iterator to use parallel processing of chromosome positions
const N_LANES:              usize = 8;         // number of lanes in the SIMD vector, i.e., 8 for i16x8
const MIN_DINUC_SEARCH_LEN: usize = 1_000_000; // minimum length of a slice to use parallel search for positioned dinucleosomes

// define scoring parameters
const NFR_WEIGHT: i16 = 1;   // weight of the NFR endpoints in scoring masks
const NUC_WEIGHT: i16 = 2;   // weight of the nucleosome centers in scoring masks
const MIN_SCORE:  i16 = 150; // minimum score to report a positioned dinucleosome for the stage with the lowest number of inserts (WAS 200)
const EARLY_REJECTION_SCORE: i16 = -999; // early rejection score for positions with too few endpoints that weren't subject to full mask processing
const GAP_LEN_PENALTY_BASE:  f64 = 1.0; // penalty for extending the gap length by one base, used to disfavor longer gaps with few additional endpoints (WAS 0.5)

// struct to hold insert endpoint counts for one stage
struct StageInsertMap {
    // vectors of insert counts
    // one vector for counting base positions of each insert endpoint, for scoring nucleosome-free (accessible) regions
    // one vector for counting base positions of each insert centers >NUC_LEN, for scoring positioned nucleosome regions
    // each vector element represents one coordinate position on the working chromosome
    // each count is held as i16 to support max count of 32,767, more than enough for ATAC-seq (see note above why i16 is used, not u16)
    pub endpoints:   Vec<i16>,
    pub nuc_centers: Vec<i16>,

    // the range of endpoint positions encountered on the current working chromosome for this stage
    pub min_start0:  usize,
    pub max_end1:    usize, 

    // the total number of inserts for this stage+chromosome, used to scale the stage-specific minimum score threshold
    pub n_inserts: usize,
}
impl StageInsertMap {
    pub fn new(chrom_length: usize) -> Self {
        StageInsertMap {
            endpoints:      vec![0; chrom_length],
            nuc_centers:    vec![0; chrom_length],
            min_start0:     chrom_length,
            max_end1:       0,
            n_inserts:  0,
        }
    }
}

// struct to assemble chains of overlapping dinucleosomes
struct NucChain {
    start0:       u32,      // start position of the overall chain span, including FLANK_LEN
    end1:         u32,      // end position of the overall chain span
    nuc_starts0:  Vec<u32>, // two or more 0-referenced start positions of the nucleosomes in the chain
    merge_types:  Vec<u8>,  // types of the nucleosome merges, 0 for flank overlap, 1 for nucleosome overlap
    dinuc_scores: Vec<i16>, // the original scores of all the merged dinucleosomes
}
impl NucChain {

    // start a new nucleosome chain from a single input dinucleosome
    pub fn new(i0: u32, score: i16, gap_len: i16) -> Self {
        let mut chain = NucChain{
            start0: 0,
            end1: 0,
            nuc_starts0:  vec![0; 10000], // more nucleosomes than will ever be required in a chain
            merge_types:  vec![0; 10000],
            dinuc_scores: vec![0; 10000],
        };
        chain.reset(i0, score, gap_len);
        chain
    }

    // reset the chain by jumping to the next non-overlapping dinucleosome without reallocation
    pub fn reset(&mut self, i0: u32, score: i16, gap_len: i16) {
        let half_gap_len = ((gap_len - 1) / 2) as u32;
        let half_len = (FLANK_LEN + NUC_LEN) as u32 + half_gap_len;

        if half_len > i0 {
            panic!("RESET i0 = {}, half_len = {}", i0, half_len);
        }
        if half_gap_len + NUC_LEN as u32 > i0 {
            panic!("RESET i0 = {}, half_gap_len = {}, NUC_LEN = {}", i0, half_gap_len, NUC_LEN);
        }

        self.start0 = i0 - half_len;
        self.end1 = i0 + half_len + 1;
        self.nuc_starts0[0] = i0 - half_gap_len - NUC_LEN as u32;
        self.nuc_starts0[1] = i0 + half_gap_len + 1;
        self.nuc_starts0.truncate(2);
        self.merge_types.truncate(0);
        self.dinuc_scores[0] = score;
        self.dinuc_scores.truncate(1);
    }

    // attempt to merge an overlapping dinucleosome into the current working chain
    // return true if the dinucleosome was merged, false if no overlap
    pub fn merge_overlap(&mut self, i0: u32, score: i16, gap_len: i16) -> bool {
        let half_gap_len = ((gap_len - 1) / 2) as u32;
        let half_len = (FLANK_LEN + NUC_LEN) as u32 + half_gap_len;

        if half_len > i0 {
            panic!("MERGE i0 = {}, half_len = {}", i0, half_len);
        }
        if half_gap_len + NUC_LEN as u32 > i0 {
            panic!("MERGE i0 = {}, half_gap_len = {}, NUC_LEN = {}", i0, half_gap_len, NUC_LEN);
        }

        let nuc_start0_left  = *self.nuc_starts0.last().unwrap();  // rightmost nucleosome start in the working chain
        let nuc_start0_right = i0 - half_gap_len - NUC_LEN as u32; // leftmost  nucleosome start in the query dinucleosome

        // handle overlapping nucleosomes by taking a weighted average of their start positions
        // modifies the last nucleosome in the chain, then adds one new nucleosome
        if nuc_start0_right < nuc_start0_left + NUC_LEN as u32 {
            self.end1 = i0 + half_len + 1;
            let i0_l = nuc_start0_left  as f64;
            let i0_r = nuc_start0_right as f64;
            let s_l  = *self.dinuc_scores.last().unwrap() as f64;
            let s_r  = score as f64;
            self.nuc_starts0.pop();
            self.nuc_starts0.extend(vec![
                ((i0_l * s_l + i0_r * s_r) / (s_l + s_r)) as u32, // not necessarily a multiple of CHROM_POS_STEP
                i0 + half_gap_len + 1,
            ]);
            self.merge_types.push(1); // 1 for nucleosome overlap
            self.dinuc_scores.push(score);
            true

        // handle overlapping flanks by fusing them into a new gap
        // adds two new nucleosomes to the chain, one on each side of the new merged gap
        } else if nuc_start0_right <= nuc_start0_left + (NUC_LEN + MAX_GAP_LEN) as u32 {
            self.end1 = i0 + half_len + 1;
            self.nuc_starts0.extend(vec![
                nuc_start0_right, 
                i0 + half_gap_len + 1,
            ]);
            self.merge_types.push(0); // 0 for flank overlap, i.e. without a nucleosome merge
            self.dinuc_scores.push(score);
            true

        // otherwise return false if no overlap
        } else {
            false
        }
    }

    // evaluate scores for all stages for a completed chain, i.e., one with no pending overlaps
    // print chain to STDOUT for final table assembly in R
    pub fn emit(&self, chrom: &str, stage: &str, insert_map: &InsertMap, writer: &mut BufWriter<impl Write>) {

        // create custom masks for the exact chain
        let chain_length = (self.end1 - self.start0) as usize;
        let mut nfr_mask = vec![NFR_WEIGHT as u32; chain_length];
        let mut nuc_mask = vec![0_u32;             chain_length];
        for nuc_start0 in self.nuc_starts0.iter() {
            let nuc_start_idx = (*nuc_start0 - self.start0) as usize;
            let nuc_end_idx   = nuc_start_idx + NUC_LEN as usize;
            nfr_mask[nuc_start_idx..nuc_end_idx].fill(0);
            nuc_mask[nuc_start_idx..nuc_end_idx].fill(NUC_WEIGHT as u32);
        }

        // collect the scores for the chain across all stages, as well as the gap and nuc scores for the index stage
        let (stage_scores, stage_counts, index_nfr_score, index_nuc_score) = insert_map.collect_stage_scores(
            stage, self.start0 as usize, chain_length, 
            nfr_mask, nuc_mask
        );

        // print the chain to STDOUT
        let merge_types = if self.merge_types.is_empty() {
            "NA".to_string()
        } else {
            self.merge_types.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
        };
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
            chrom, // BED3 span of the overall chain
            self.start0,
            self.end1,
            stage,
            self.nuc_starts0.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
            merge_types,
            self.dinuc_scores.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
            index_nfr_score,
            index_nuc_score,
            stage_scores,
            stage_counts,
        ).unwrap(); // do not expect any errors here
    }
}

// struct to hold and process insert counts for positioned nucleosome discovery
pub struct InsertMap {

    // the insert endpoint counts for each stage
    data:   HashMap<String, StageInsertMap>,
    stages: Vec<String>, // the stages for which insert maps are created, for maintaining strict stage ordering

    // central NFR, i.e., gap lengths used in analysis
    gap_lens: Vec<i16>,
    odd_j0s:  Vec<usize>, // same as gap_lens, as usize for indexing

    // pre-calculated masks for scoring inserts
    // one mask for each NFR gap length from 0 to MAX_GAP_LEN in the outer vector (only some of these are used)
    // each inner mask is a variable-length vector of roughly 600-700 bases
    // masks of 0|1 i16 values are multiplied by windows moved along the endpoints and nuc_centers vectors
    endpoint_masks:   Vec<Vec<i16>>, // pre-calculated masks for NFR endpoint counting in scores
    nuc_center_masks: Vec<Vec<i16>>, // pre-calculated masks for nucleosome center counting

    // the half-widths of the masks, i.e., the total number of positions on each side of the central base, by gap_len
    half_mask_widths:    Vec<usize>,
    max_half_mask_width: usize, // used to determine the working size of the chromosome score matrix

    // combined mask lengths over all scored spans, by gap_len
    mask_lens: Vec<usize>,

    // penalties for longer gaps, i.e., disfavor gap extensions that add few additional endpoints, by gap_len
    gap_len_penalties: Vec<i16>, 

    // the first 0-referenced index past the last full SIMD lane, by gap_len
    mask_tail_i0s: Vec<usize>,

    // the stage-specific score thresholds
    min_scores: HashMap<String, i16>,
}
impl InsertMap {

    // instantiate a new InsertMap
    pub fn new(chrom_length: usize, stages: Vec<&str>) -> Self {
        eprintln!("    instantiating insert map");
        let half_mask_widths: Vec<usize> = (0..=MAX_GAP_LEN).into_iter()
            .map(|gap_len| {
                (FLANK_LEN + NUC_LEN + gap_len.saturating_sub(1) / 2) as usize
            })
            .collect();
        InsertMap {
            data: stages.iter().map(|stage| {
                (stage.to_string(), StageInsertMap::new(chrom_length))
            }).collect(),
            stages: stages.iter().map(|s| s.to_string()).collect(),
            gap_lens: (MIN_GAP_LEN..=MAX_GAP_LEN).step_by(GAP_LEN_STEP).collect(),
            odd_j0s:  (MIN_GAP_LEN..=MAX_GAP_LEN).step_by(GAP_LEN_STEP).map(|j| j as usize).collect(),
            endpoint_masks:   Self::get_endpoint_masks(),
            nuc_center_masks: Self::get_nuc_center_masks(),
            half_mask_widths:    half_mask_widths.clone(),
            max_half_mask_width: half_mask_widths[MAX_GAP_LEN as usize],
            mask_lens: (0..=MAX_GAP_LEN).into_iter().map(|gap_len| {
                ((FLANK_LEN + NUC_LEN) * 2 + gap_len) as usize
            }).collect(),
            gap_len_penalties: (0..=MAX_GAP_LEN).into_iter().map(|gap_len| {
                ((gap_len - MIN_GAP_LEN) as f64 * GAP_LEN_PENALTY_BASE) as i16
            }).collect(),
            mask_tail_i0s: (0..=MAX_GAP_LEN).into_iter().map(|gap_len| {
                let mask_len = ((FLANK_LEN + NUC_LEN) * 2 + gap_len) as usize;
                mask_len / N_LANES * N_LANES
            }).collect(),
            min_scores: stages.iter().map(|stage| {
                (stage.to_string(), 0) // initialized to 0, will be set later
            }).collect(),
        }
    }

    // insert endpoint and center masks are multiplied by a query span of the chromosome insert counts
    // they count:
    //      insert endpoints with a weight of 1 per endpoint in putative NFR accessible regions
    //      insert centers   with a weight of 2 per center (1 per endpoint) for putative positioned-nucleosome spanning inserts
    //          thus, inserts that span a nucleocome may be counted with a net weight of 4 (1 for each endpoint, 2 for the center)
    fn get_endpoint_masks() -> Vec<Vec<i16>> {
        let mut masks: Vec<Vec<i16>> = Vec::new();
        (0..=MAX_GAP_LEN).into_iter().for_each(|gap_len| {
            let mask: Vec<i16> = [
                vec![NFR_WEIGHT; FLANK_LEN as usize],
                vec![0;          NUC_LEN   as usize],
                vec![NFR_WEIGHT; gap_len   as usize],
                vec![0;          NUC_LEN   as usize],
                vec![NFR_WEIGHT; FLANK_LEN as usize],
            ].concat();
            masks.push(mask);
        });
        masks
    }
    fn get_nuc_center_masks() -> Vec<Vec<i16>> {
        let mut masks: Vec<Vec<i16>> = Vec::new();
        (0..=MAX_GAP_LEN).into_iter().for_each(|gap_len| {
            let mask: Vec<i16> = [
                vec![0;          FLANK_LEN as usize],
                vec![NUC_WEIGHT; NUC_LEN   as usize],
                vec![0;          gap_len   as usize],
                vec![NUC_WEIGHT; NUC_LEN   as usize],
                vec![0;          FLANK_LEN as usize],
            ].concat();
            masks.push(mask);
        });
        masks
    }

    // load the inserts for the chromosome from each stage file
    pub fn load_inserts_by_stage(
        &mut self, chrom: &str, stage_insert_files: &HashMap<String, String>
    ) -> Result<(), Box<dyn std::error::Error>> {
        eprintln!("    loading inserts");
        for stage in self.stages.iter() {
            eprintln!("      stage {}", stage);
            let file_path = stage_insert_files.get(stage).unwrap();
            let stage_insert_map = self.data.get_mut(stage).unwrap();
            let mut cmd = Command::new("tabix")
                .arg(file_path)

                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                .arg(chrom)
                // .arg(format!("{}:25000000-50000000", chrom))

                .stdout(Stdio::piped())
                .spawn()?;
            let stdout = cmd.stdout.take().unwrap();
            let reader = BufReader::new(stdout);
            for line in reader.lines() {
                let line = line?;
                let parts: Vec<&str> = line.split('\t').collect();
                let start0: usize = parts[1].parse()?;
                let end1:   usize = parts[2].parse()?;
                stage_insert_map.endpoints[start0]   += 1;
                stage_insert_map.endpoints[end1 - 1] += 1;
                let insert_size = end1 - start0;
                if insert_size >  NUC_LEN as usize &&
                   insert_size <= MAX_NUC_LEN {
                    let center0 = start0 + insert_size / 2;
                    stage_insert_map.nuc_centers[center0] += 1; // nucleosome insert centers counted as 1, weighted as 2 in the masks
                }
                if start0 < stage_insert_map.min_start0 { stage_insert_map.min_start0 = start0; }
                if end1   > stage_insert_map.max_end1   { stage_insert_map.max_end1 = end1; }
                stage_insert_map.n_inserts += 1;
            }
        }
        Ok(())
    }

    // set the score thresholds for each stage based on the total number of inserts
    // samples with more inserts must have a higher threshold to be called a positioned dinucleosome
    // at present, only uses insert counts on the working chromosome; if insufficient, need to provide genome insert counts as argument
    pub fn set_min_scores(&mut self) {
        eprintln!("    setting score thresholds");
        let min_n_inserts = self.data.values().map(|m| m.n_inserts).min().unwrap() as f64;
        for (stage, stage_insert_map) in self.data.iter_mut() {
            let stage_min_score = (MIN_SCORE as f64 * stage_insert_map.n_inserts as f64 / min_n_inserts) as i16;
            self.min_scores.insert(stage.clone(), stage_min_score);
        }
    }

    // move a scoring mask over the mapped counts for all central NFR gap lengths in use
    // establish a score for the evidence that each position is a nucleosome-free region at each gap length
    // choose the best score and its corresponding gap length at each chromosome position
    // extract non-overlapping local maxima above a threshold score as positioned dinucleosomes
    // cross-compare scores for all candidate dinucleosomes to all other stages
    // i is the index for the chromosome position, j is the index for the central NFR gap length
    pub fn find_dinucs_by_stage(&self, chrom: &str) {
        eprintln!("    finding positioned nucleosomes");
        let stdout = stdout();
        let mut writer = BufWriter::new(stdout.lock());
        for stage in self.stages.iter() {
            eprintln!("      stage {}", stage);
            let stage_insert_map = self.data.get(stage).unwrap();
            let stage_min_score = *self.min_scores.get(stage).unwrap();

            // set an iterator over the range of chromosome positions with data
            // ensure that all stages interrogate the same chromosome positions at multiples of CHROM_POS_STEP
            let min_chrom_i0 = (stage_insert_map.min_start0 + self.max_half_mask_width) / CHROM_POS_STEP *  CHROM_POS_STEP;
            let max_chrom_i1 = (stage_insert_map.max_end1   - self.max_half_mask_width) / CHROM_POS_STEP *  CHROM_POS_STEP;
            let i_par_iter = (min_chrom_i0..max_chrom_i1).into_par_iter()
                .step_by(CHROM_POS_STEP) // every nth coordinate position to reduce the number of iterations
                .with_min_len(MIN_PAR_ITER_LEN);

            // calculate the best gap-length-score at each chromosome position when sufficient inserts are present
            let best_scores: Vec<(i16, i16)> = i_par_iter.map(|i0| {
                if self.max_possible_score(stage_insert_map, i0) < stage_min_score { 
                    return (EARLY_REJECTION_SCORE, MAX_GAP_LEN);  // short circuit most positions with too-few endpoints
                }
                let gap_scores = self.odd_j0s.iter().map(|j0| {
                    self.dot_product(stage_insert_map, i0, *j0)
                }).collect::<Vec<i16>>();
                gap_scores.into_iter()
                    .zip(self.gap_lens.iter().cloned())
                    .max_by_key(|(score, _)| *score)
                    .unwrap() // cannot be empty
            }).collect();

            // recursively find positioned dinucleosomes in the best scores
            // the net algorithm optimizes local combinations of NFR center position and gap length
            let mut dinuc_centers: Vec<(u32, i16, i16)> = Vec::new();
            let length = best_scores.len();
            find_positioned_dinuc(
                &best_scores, 
                0, 
                length, 
                &mut dinuc_centers,
                min_chrom_i0,
                stage_min_score
            );
            dinuc_centers.par_sort_unstable_by_key(|&(i0, _, _)| i0); // sort by chromosome position of the gap center

            // collapse overlapping dinucleosomes into nucleosome chains and print the results
            let mut dinuc_iter = dinuc_centers.into_iter();
            let mut chain = dinuc_iter.next().map(|(i0, score, gap_len)| {
                NucChain::new(i0, score, gap_len) // initialize the working chain, allocated once
            }).unwrap();
            for (i0, score, gap_len) in dinuc_iter {
                // check for two types of overlap between the working chain and the next dinucleosome
                // if overlap found, merge the dinucleosome into the working chain
                // otherwise emit the chain and start a new one
                if !chain.merge_overlap(i0, score, gap_len) {
                    chain.emit(chrom, stage, self, &mut writer);
                    chain.reset(i0, score, gap_len);
                }
            }
            chain.emit(chrom, stage, self, &mut writer); // commit the last working chain

            // // cross-correlate scores of all dinucs found in this stage against all other stages
            // // print results for final table assembly in R
            // for (i0, score, gap_len) in dinuc_centers {
            //     let j0 = gap_len as usize;
            //     let gap_start_i0 = i0 as usize - (j0 - 1) / 2;
            //     let gap_score = parse_score_part( // the portion of the score contributed by the endpoints in the central NFR gap
            //         &stage_insert_map.endpoints, 
            //         gap_start_i0, 
            //         j0, 
            //         NFR_WEIGHT
            //     );
            //     let nuc_score = parse_score_part( // the portion of the score contributed by the nucleosome centers flanking the central NFR gap
            //         &stage_insert_map.nuc_centers, 
            //         gap_start_i0 - NUC_LEN as usize, // putative nucleosome center on the left of the gap
            //         NUC_LEN as usize, 
            //         NUC_WEIGHT
            //     ) + parse_score_part(
            //         &stage_insert_map.nuc_centers, 
            //         gap_start_i0 + j0, // putative nucleosome center on the right of the gap
            //         NUC_LEN as usize, 
            //         NUC_WEIGHT
            //     );
            //     write!(
            //         &mut writer,
            //         "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
            //         chrom, 
            //         i0 - self.half_mask_widths[j0] as u32,
            //         i0 + self.half_mask_widths[j0] as u32 + 1,

            //         i0, // gap center
            //         gap_len,

            //         stage,

            //         score + self.gap_len_penalties[j0], // report the actual score without the gap extension penalty used for optimizing fits
            //         gap_score,
            //         nuc_score,
            //         self.collect_stage_scores(stage, score, i0 as usize, j0)

            //     ).unwrap(); // do not expect any errors here
            // }
        }
    }

    // use a simple endpoint+nuc_center count to estimate the maximum possible score at a given position
    // assume the largest possible mask, i.e., the largest gap length
    fn max_possible_score(&self, stage_insert_map: &StageInsertMap, i0: usize) -> i16 {
        let min_chrom_i0 = i0 - self.half_mask_widths[MAX_J0];
        let mut sum = i16x8::splat(0);
        for mask_i0 in (0..self.mask_tail_i0s[MAX_J0]).step_by(N_LANES) {
            let chrom_mask_i0 = min_chrom_i0 + mask_i0;
            sum += i16x8::from_slice_unaligned(&stage_insert_map.endpoints[  chrom_mask_i0..chrom_mask_i0 + N_LANES]) +
                   i16x8::from_slice_unaligned(&stage_insert_map.nuc_centers[chrom_mask_i0..chrom_mask_i0 + N_LANES]);
        }
        sum.as_array_ref().iter().sum() // don't assess the tail here, it is too small a fraction of the total to matter
    }

    // use SIMD  for efficient dot product calculation of counts * mask
    fn dot_product(&self, stage_insert_map: &StageInsertMap, i0: usize, j0: usize) -> i16 {
        let min_chrom_i0 = i0 - self.half_mask_widths[j0];

        // calculate the SIMD dot product on the full chunks of N_LANES
        let mut sum = i16x8::splat(0);
        for mask_i0 in (0..self.mask_tail_i0s[j0]).step_by(N_LANES) {
            let chrom_mask_i0 = min_chrom_i0 + mask_i0;
            sum += i16x8::from_slice_unaligned(&stage_insert_map.endpoints[  chrom_mask_i0..chrom_mask_i0 + N_LANES]) * 
                   i16x8::from_slice_unaligned(&self.endpoint_masks[j0][           mask_i0..mask_i0       + N_LANES]) +
                   i16x8::from_slice_unaligned(&stage_insert_map.nuc_centers[chrom_mask_i0..chrom_mask_i0 + N_LANES]) * 
                   i16x8::from_slice_unaligned(&self.nuc_center_masks[j0][         mask_i0..mask_i0       + N_LANES]);
        }

        // calculate the horizontal sum of the SIMD lanes
        // the nature of the data ensures that the sum fits in i16
        let mut dot_product: i16 = sum.as_array_ref().iter().sum();

        // handle any remaining elements that did not fit into a full chunk
        if self.mask_tail_i0s[j0] < self.mask_lens[j0] {
            for mask_i0 in self.mask_tail_i0s[j0]..self.mask_lens[j0] {
                let chrom_mask_i0 = min_chrom_i0 + mask_i0;
                dot_product += stage_insert_map.endpoints[  chrom_mask_i0] * self.endpoint_masks[  j0][mask_i0] +
                               stage_insert_map.nuc_centers[chrom_mask_i0] * self.nuc_center_masks[j0][mask_i0];
            }
        }

        // return the final dot product result
        // add a penalty against overly long gap lengths, i.e., disfavor gap extension to add few additional endpoints
        dot_product - self.gap_len_penalties[j0] // the penalty is 0 for the shortest gap length, otherwise GAP_LEN / 2
    }

    // // collect the scores for a given candidate dinucleosome across all stages
    // fn collect_stage_scores(&self, index_stage: &str, index_stage_score: i16, i0: usize, j0: usize) -> String {
    //     let mut scores: Vec<String> = Vec::new();
    //     for stage in self.stages.iter() {
    //         let stage_score = if stage == index_stage {
    //             index_stage_score // use the previously score calculated for the index stage
    //         } else {
    //             let stage_insert_map = self.data.get(stage).unwrap();
    //             self.dot_product(stage_insert_map, i0, j0) // otherwise calculate the score for the other stage
    //         };
    //         let final_score = stage_score + self.gap_len_penalties[j0]; // remove the gap extension penalty for best score linearity
    //         scores.push(final_score.to_string());
    //     }
    //     scores.join("\t")
    // }

    // collect the scores for a given candidate dinucleosome across all stages
    fn collect_stage_scores(
        &self, index_stage: &str, start0: usize, chain_length: usize, 
        nfr_mask: Vec<u32>, nuc_mask: Vec<u32>
    ) -> (String, String, u32, u32) {
        let mut stage_counts: Vec<String> = Vec::new();
        let mut stage_scores: Vec<String> = Vec::new();
        let mut index_nfr_score = 0;
        let mut index_nuc_score = 0;
        for stage in self.stages.iter() {
            let stage_insert_map = self.data.get(stage).unwrap();
            // TODO: slow enough to warrant SIMD?
            let stage_count = stage_insert_map.endpoints[start0..start0 + chain_length].iter().map(|&x| x as u32).sum::<u32>();
            let stage_nfr_score = stage_insert_map.endpoints[start0..start0 + chain_length].iter()
                .zip(nfr_mask.iter())
                .map(|(&count, &mask)| count as u32 * mask )
                .sum::<u32>();
            let stage_nuc_score = stage_insert_map.nuc_centers[start0..start0 + chain_length].iter()
                .zip(nuc_mask.iter())
                .map(|(&count, &mask)| count as u32 * mask )
                .sum::<u32>();
            stage_counts.push(stage_count.to_string());
            stage_scores.push((stage_nfr_score + stage_nuc_score).to_string());
            if stage == index_stage {
                index_nfr_score = stage_nfr_score; // use the previously score calculated for the index stage
                index_nuc_score = stage_nuc_score;
            }
        }
        (stage_scores.join("\t"), stage_counts.join("\t"), index_nfr_score, index_nuc_score)
    }
}

// find the best positioned dinucleosome in a slice of best scores at the best gap lengths
// start with entire chromosome of best scores
// if a positioned dinucleosome is found above the score threshold:
//      record the position, score, and gap length of the dinucleosome
//      remove the entire dinucleosome span (but not the flanks) from further consideration
//      recursively search the remaining chromosome spans to the left and right of the dinucleosome
// use rayon join to parallelize the search, falling back to serial search when the slice falls below a threshold length
fn find_positioned_dinuc(
    best_scores: &[(i16, i16)], // the same reference to best_scores in all calls; offset and length define the working slice
    offset_i0: usize, // indices are relative to best_scores, which does not include all coordinate positions
    length: usize,    // number of best scores in the slice, again, not the number of chromosome positions
    dinuc_centers: &mut Vec<(u32, i16, i16)>,
    min_chrom_i0: usize,
    min_score: i16,
) {
    if length == 0 { return; } // terminate recursion when slice is empty ...

    // find the best dinucleosome center that does not conflict with prior high-scoring nucleosome placements
    // a first step toward dinuc collision avoidance was established by the left/right half slicing below
    // but a gap center outside of a blocked dinucleosome span may still place a second nucleosome too far into a prior committed central gap
    // this situation cannot be fully avoided prior to maximum finding here as it depends on gap length of the new left/right candidate dinuc
    // thus, the max_by_key() closure returns a zero score when it would place a nucleosome too far off the end of the working slice
    let min_nuc_start_bp = (offset_i0 * CHROM_POS_STEP).saturating_sub((BLOCK_HALF_LEN + COLLISION_TOLERANCE) as usize);
    let max_nuc_end_bp   = (offset_i0 + length) * CHROM_POS_STEP + (BLOCK_HALF_LEN + COLLISION_TOLERANCE) as usize;
    let (slice_i0, (score, gap_len)) = best_scores[offset_i0..offset_i0 + length]
        .iter().enumerate()
        .max_by_key(|(slice_i0, (score, gap_len))| {
            let gap_center_bp = (offset_i0 + slice_i0) * CHROM_POS_STEP;
            let half_width_bp = (NUC_LEN + (gap_len - 1) / 2) as usize;
            let left_nuc_start_bp = gap_center_bp.saturating_sub(half_width_bp);
            let right_nuc_end_bp  = gap_center_bp + half_width_bp;
            if left_nuc_start_bp > min_nuc_start_bp &&
               right_nuc_end_bp  < max_nuc_end_bp { 
                *score 
            } else { 
                0
            }
        }) // thus, the best score in the working slice
        .unwrap(); // cannot be empty, length > 0
    if *score < min_score { return; } // ... or when the best score is below the minimum score threshold
    let best_scores_i0 = offset_i0 + slice_i0; // index of the best score in the best_scores vector (not a chrom position)

    // commit the best dinucleosome center found in this slice
    dinuc_centers.push(((min_chrom_i0 + best_scores_i0 * CHROM_POS_STEP) as u32, *score, *gap_len));

    // establish the boundaries of the left and right flanks of the dinucleosome for further searching
    // ensures that subsequent gap centers cannot fall within the nucleosomes or central gap of this committed and now blocked dinucleosome
    // this begins but is insufficient for avoiding dinucleosome collisions as desribed above
    let blocked_half_width_bp = (BLOCK_HALF_LEN + (gap_len - 1) / 2) as usize; // number of blocked coordinate positions on each side of the dinucleosome center
    let blocked_half_width_idx= blocked_half_width_bp / CHROM_POS_STEP; // blocked half-width in best_scores indices
    let dinuc_start_i0 = best_scores_i0.saturating_sub(blocked_half_width_idx);
    let dinuc_end_i1   = (best_scores_i0 + blocked_half_width_idx + 1).min(best_scores.len());
    let left_len  = dinuc_start_i0.saturating_sub(offset_i0);
    let right_len = length.saturating_sub(dinuc_end_i1 - offset_i0);

    // when the slice is small enough, continue the search for positioned dinucleosomes in series over the left and right halves
    if length <= MIN_DINUC_SEARCH_LEN {
        find_positioned_dinuc(
            best_scores,
            offset_i0,
            left_len,
            dinuc_centers,
            min_chrom_i0,
            min_score,
        );
        find_positioned_dinuc(
            best_scores,
            dinuc_end_i1, // right side offset_i0 is first position after the dinucleosome end
            right_len,
            dinuc_centers,
            min_chrom_i0,
            min_score,
        );

    // when the slice is large enough, continue the search for positioned dinucleosomes in parallel by join left and right halves
    } else {
        let mut left_dinuc_centers:  Vec<(u32, i16, i16)> = Vec::new();
        let mut right_dinuc_centers: Vec<(u32, i16, i16)> = Vec::new();
        join(
            || find_positioned_dinuc(
                best_scores,
                offset_i0,
                left_len,
                &mut left_dinuc_centers,
                min_chrom_i0,
                min_score,
            ),
            || find_positioned_dinuc(
                best_scores,
                dinuc_end_i1,
                right_len,
                &mut right_dinuc_centers,
                min_chrom_i0,
                min_score,
            )
        );
        dinuc_centers.extend(left_dinuc_centers);
        dinuc_centers.extend(right_dinuc_centers);
    }
}

// // extract the component parts of the final dinuc scores
// // these are useful for distinguishing true positioned dinucleosomes from simpler high accessibility regions
// fn parse_score_part(data: &Vec<i16>, offset_i0: usize, length: usize, weight: i16) -> i16 {
//     data[offset_i0..offset_i0 + length].iter()
//         .map(|&n| n * weight)
//         .sum::<i16>()
// }
