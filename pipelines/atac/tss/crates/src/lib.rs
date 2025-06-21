// dependencies
use std::io::Write;
use rayon::prelude::*;
use rayon::join;
use wide::i16x8;

// set constants established from biological principles and observations
const NUC_LEN:      u16   = 147; // nucleosome span in bp, u16 to match the size of the masks
const MAX_CHROM_LEN:usize = 250_000_000; // 250 Mbp to hold human chromosome 1
const MIN_GAP_LEN:  u8    = 51;  // minimum allowed central NFR length
const MAX_GAP_LEN:  u8    = 255; // maximum allowed central NFR length, fits into u8 so gap_len ~ vector index
const TOLERANCE:    u16   = 5;   // tolerance for imprecise placement of flanking nucleosomes across a population of cells
const FLANK_LEN:    u16   = 100; // the length of the NFRs flanking the positioned nucleosomes
const MIN_SCORE:    i16   = 25;  // minimum score to report a positioned dinucleosome, i16 to match the size of the masks

// set constants that define parallelization implementation
const MIN_PAR_ITER_LEN: usize = 1_000_000; // minimum length of a parallel iterator to use parallel processing
const N_LANES: usize = 8; // number of lanes in the SIMD vector, i.e., 8 for i16x8
const MIN_DINUC_SEARCH_LEN: usize = 1_000_000; // minimum length of a slice to use parallel search for positioned dinucleosomes
const GAP_LEN_STEP: usize = 2; // step size for the central NFR gap lengths, only odd numbers are used
const CHROM_POS_STEP: usize = 2; // step size for the chromosome positions, only every other position is processed to reduce iterations

// configure a struct to hold and process insert endpoint counts
pub struct InsertMap {

    // vectors of insert endpoint counts, held one chromosome at a time
    // one vector for each endpoint type, i.e., start and end, to support distinct masks
    // each vector element represents one coordinate position on the working chromosome
    // each count is held as i16 to support negative scores and a max count >128
    pub starts: Vec<i16>,
    pub ends:   Vec<i16>,

    // central NFR, i.e., gap lengths to assess, limited to odd numbers only
    pub gap_lens: Vec<u8>,   // allowed central NFR lengths, i.e., the sizes of the gaps between positioned nucleosomes
    pub odd_j0s: Vec<usize>, // same as gap_lens as usize for indexing

    // pre-calculated masks for scoring insert endpoints
    // one mask for each NFR gap length from 0 to 255 in the outer vector (only some of these are used)
    // each inner mask is a variable-length vector of roughly 600-700 bases
    // masks of +1/0/-1 i16 values are multiplied by windows moved along the starts and ends vectors
    pub start_masks: Vec<Vec<i16>>, // pre-calculated masks for start counting in scores, by gap length
    pub end_masks:   Vec<Vec<i16>>, // pre-calculated masks for end   counting in scores, by gap length

    // the half-widths of the masks, i.e., the total number of positions on each side of the central base
    pub half_mask_widths: Vec<usize>,

    // the maximum half-width of all masks, used to determine the working size of the chromosome score matrix
    pub max_half_mask_width: usize,

    // combined mask lengths over all scored spans
    pub mask_lens: Vec<usize>,

    // the first 0-referenced index past the last full SIMD lane
    pub mask_tail_i0s: Vec<usize>,

    // the current working chromosome
    pub working_chrom: String, 

    // the range of endpoint positions encountered on the current working chromosome
    pub min_start0: usize,
    pub max_end1:   usize, 
}
impl InsertMap {

    // instantiate a new InsertMap
    pub fn new() -> Self {
        let half_mask_widths: Vec<usize> = (0..=MAX_GAP_LEN).into_iter()
            .map(|gap_len| {
                (FLANK_LEN + NUC_LEN + (gap_len as u16 - 1) / 2) as usize
            })
            .collect();
        InsertMap {
            starts: vec![0; MAX_CHROM_LEN], // allocated once, reset on each chromosome
            ends:   vec![0; MAX_CHROM_LEN],
            gap_lens: (MIN_GAP_LEN..=MAX_GAP_LEN).step_by(GAP_LEN_STEP).collect(),
            odd_j0s:  (MIN_GAP_LEN..=MAX_GAP_LEN).step_by(GAP_LEN_STEP).map(|j| j as usize).collect(),
            start_masks:  Self::get_start_masks(),
            end_masks:    Self::get_end_masks(),
            half_mask_widths:    half_mask_widths.clone(),
            max_half_mask_width: half_mask_widths[MAX_GAP_LEN as usize],
            mask_lens: (0..=MAX_GAP_LEN).into_iter().map(|gap_len| {
                (FLANK_LEN as usize + NUC_LEN as usize) * 2 + gap_len as usize
            }).collect(),
            mask_tail_i0s: (0..=MAX_GAP_LEN).into_iter().map(|gap_len| {
                let mask_len = (FLANK_LEN as usize + NUC_LEN as usize) * 2 + gap_len as usize;
                mask_len / N_LANES * N_LANES
            }).collect(),
            working_chrom: "".to_string(),
            min_start0: 0,
            max_end1:   0, 
        }
    }

    // start and end masks are multiplied by a query span of the chromosome insert endpoints counts
    // they enforce:
    //      positive scores for endpoints matching the positioned nucleosome pattern
    //      negative penalties for endpoints that do not match the pattern
    //      zero scores for endpoints that are neither awarded nor penalized
    fn get_start_masks() -> Vec<Vec<i16>> {
        let mut masks: Vec<Vec<i16>> = Vec::new();
        (0..=MAX_GAP_LEN).into_iter().for_each(|gap_len| {
            let mask: Vec<i16> = [
                vec![ 1; FLANK_LEN as usize],
                vec![ 0; TOLERANCE as usize],
                vec![-1; NUC_LEN   as usize - TOLERANCE as usize * 2],
                vec![ 0; TOLERANCE as usize],
                vec![ 1; gap_len   as usize],
                vec![ 0; TOLERANCE as usize],
                vec![-1; NUC_LEN   as usize - TOLERANCE as usize * 2],
                vec![ 0; TOLERANCE as usize],
                vec![ 0; FLANK_LEN as usize],
            ].concat();
            masks.push(mask);
        });
        masks
    }
    fn get_end_masks() -> Vec<Vec<i16>> {
        let mut masks: Vec<Vec<i16>> = Vec::new();
        (0..=MAX_GAP_LEN).into_iter().for_each(|gap_len| {
            let mask: Vec<i16> = [
                vec![ 0; FLANK_LEN as usize],
                vec![ 0; TOLERANCE as usize],
                vec![-1; NUC_LEN   as usize - TOLERANCE as usize * 2],
                vec![ 0; TOLERANCE as usize],
                vec![ 1; gap_len   as usize],
                vec![ 0; TOLERANCE as usize],
                vec![-1; NUC_LEN   as usize - TOLERANCE as usize * 2],
                vec![ 0; TOLERANCE as usize],
                vec![ 1; FLANK_LEN as usize],
            ].concat();
            masks.push(mask);
        });
        masks
    }

    // add an insert to the map, processing the previous chromosome as necessary
    pub fn add_insert (&mut self, chrom: String, start0: usize, end1: usize) {

        // if the chromosome has changed, process the previous chromosome
        // reset to the new chromosome before adding the new insert
        if chrom != self.working_chrom {
            if self.working_chrom != "" {
                self.process_chrom();
                self.starts.fill(0);
                self.ends.fill(0);
            }
            self.working_chrom = chrom;
            self.min_start0 = start0;
            self.max_end1 = 0;
        }

        // increment the endpoint counts for the current insert
        self.starts[start0] += 1;
        self.ends[end1 - 1] += 1;
        if end1 > self.max_end1 { self.max_end1 = end1; }
    }

    // move a scoring mask over the mapped counts for all central NFR gap lengths in use
    // establish a score for the evidence that each position is a nucleosome-free region at each gap length
    // choose the best score and its corresponding gap length at each chromosome position
    // extract non-overlapping local maxima above a threshold score as positioned dinucleosomes
    // i is the index for the chromosome position, j is the index for the central NFR gap length
    pub fn process_chrom(&self) {
        let min_chrom_i0 = self.min_start0 + self.max_half_mask_width;
        let max_chrom_i1 = self.max_end1   - self.max_half_mask_width;
        let i_par_iter = (min_chrom_i0..max_chrom_i1).into_par_iter()
            .step_by(CHROM_POS_STEP) // every other coordinate position to reduce the number of iterations
            .with_min_len(MIN_PAR_ITER_LEN);
        let best_scores: Vec<(i16, u8)> = i_par_iter.map(|i0| {
            let gap_scores = self.odd_j0s.iter().map(|j0| {
                self.dot_product(i0, *j0)
            }).collect::<Vec<i16>>();
            gap_scores.into_iter()
                .zip(self.gap_lens.iter().cloned())
                .max_by_key(|(score, _)| *score)
                .unwrap()
        }).collect();

        // recursively find positioned dinucleosomes in the best scores
        // the net algorithm optimizes local combinations of NFR center position and mid length
        let mut dinuc_centers: Vec<(u32, i16, u8)> = Vec::new();
        find_positioned_dinuc(
            &best_scores, 
            0, 
            best_scores.len(), 
            &mut dinuc_centers
        );

        // print the positioned dinucleosomes found in the current chromosome to STDOUT
        let stdout = std::io::stdout();
        let mut writer = std::io::BufWriter::new(stdout.lock());
        for (best_scores_i0, score, gap_len) in dinuc_centers {
            let chrom_coord = min_chrom_i0 as u32 + best_scores_i0 * CHROM_POS_STEP as u32;
            write!(
                &mut writer,
                "{}\t{}\t{}\t{}\n", 
                self.working_chrom, chrom_coord, score, gap_len
            ).unwrap();
        }
    }

    // use SIMD  for efficient dot product calculation of counts * mask
    fn dot_product(&self, i0: usize, j0: usize) -> i16 {
        let min_chrom_i0 = i0 - self.half_mask_widths[j0];

        // calculate the SIMD dot product on the full chunks of N_LANES
        let mut sum = i16x8::splat(0);
        for mask_i0 in (0..self.mask_tail_i0s[j0]).step_by(N_LANES) {
            let chrom_mask_i0 = min_chrom_i0 + mask_i0;
            sum += i16x8::from_slice_unaligned(&self.starts[   chrom_mask_i0..chrom_mask_i0 + N_LANES]) * 
                   i16x8::from_slice_unaligned(&self.start_masks[j0][mask_i0..mask_i0       + N_LANES]) +
                   i16x8::from_slice_unaligned(&self.ends[     chrom_mask_i0..chrom_mask_i0 + N_LANES]) * 
                   i16x8::from_slice_unaligned(&self.end_masks[j0][  mask_i0..mask_i0       + N_LANES]);
        }

        // calculate the horizontal sum of the SIMD lanes
        // the nature of the data ensures that the sum fits in i16
        let mut dot_product: i16 = sum.as_array_ref().iter().sum();

        // handle any remaining elements that did not fit into a full chunk
        if self.mask_tail_i0s[j0] < self.mask_lens[j0] {
            for mask_i0 in self.mask_tail_i0s[j0]..self.mask_lens[j0] {
                let chrom_mask_i0 = min_chrom_i0 + mask_i0;
                dot_product += self.starts[chrom_mask_i0] * self.start_masks[j0][mask_i0] +
                               self.ends[  chrom_mask_i0] * self.end_masks[  j0][mask_i0];
            }
        }

        // return the final dot product result
        dot_product
    }

}


// find the best positioned dinucleosome in a slice of best scores at the best gap lengths
// start with entire chromosome of best scores
// if a positioned dinucleosome is found above the score threshold:
//      record the position, score, and gap length of the dinucleosome
//      remove the entire dinucleosome span from further consideration
//      recursively search the remaining chromosome spans to the left and right of the dinucleosome
// use rayon join to parallelize the search, but fall back to serial search when the slice falls below a threshold length
fn find_positioned_dinuc(
    best_scores: &[(i16, u8)], 
    offset_i0: usize, // indices are relative to best_scores, which does not include all coordiante positions
    length: usize,    // number of best scores, not number of chromosome positions
    dinuc_centers: &mut Vec<(u32, i16, u8)>
) {
    if length <= 0 { return; }
    let (best_scores_i0, (score, gap_len)) = best_scores[offset_i0..offset_i0 + length].iter().enumerate()
        .max_by_key(|(_, (score, _))| *score)
        .unwrap();
    if *score < MIN_SCORE { return; }
    dinuc_centers.push((best_scores_i0 as u32, *score, *gap_len));
    let included_half_width_bp = (*gap_len as usize - 1) / 2 + NUC_LEN as usize;
    let included_half_width_idx = included_half_width_bp / CHROM_POS_STEP; // half-width in indices, i.e., number of positions on each side of the dinucleosome center
    let left_i0  = best_scores_i0 - included_half_width_idx; // indices of the ends of the dinucleosome span
    let right_i1 = best_scores_i0 + included_half_width_idx + 1;
    if length <= MIN_DINUC_SEARCH_LEN {
        find_positioned_dinuc(
            best_scores,
            offset_i0,
            left_i0 - offset_i0,
            dinuc_centers
        );
        find_positioned_dinuc(
            best_scores,
            right_i1,
            length - (right_i1 - offset_i0),
            dinuc_centers
        );
    } else {
        let mut left_dinuc_centers:  Vec<(u32, i16, u8)> = Vec::new();
        let mut right_dinuc_centers: Vec<(u32, i16, u8)> = Vec::new();
        join(
            || find_positioned_dinuc(
                best_scores,
                offset_i0,
                left_i0 - offset_i0,
                &mut left_dinuc_centers
            ),
            || find_positioned_dinuc(
                best_scores,
                right_i1,
                length - (right_i1 - offset_i0),
                &mut right_dinuc_centers
            )
        );
        dinuc_centers.extend(left_dinuc_centers);
        dinuc_centers.extend(right_dinuc_centers);
    }
}
