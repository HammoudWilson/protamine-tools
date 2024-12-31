use strict;
use warnings;

# pull just the genes from an annotation GTF and print as simple BED
# called as perl -n gtf_to_genes_bed.pl

$_ =~ m/^#/ and next;
chomp;
my ($chrom, $source_, $feature, $start, $end, $score, $strand, $frame, $attr) = split("\t");
$feature eq "gene" or next;
$chrom =~ m/chr(\d+|X|Y)/ or next; # standard chromosomes only
$attr =~ m/gene_id "(.+?)";/ or next;
my $geneId = $1;
$attr =~ m/gene_name "(.+?)";/ or next;
my $geneNme = $1;     
print join("\t", $chrom, $start-1, $end, "$geneId/$geneNme", 0, $strand), "\n";

# chr13   HAVANA  gene    108207439       108218368       .       -       .
# gene_id "ENSG00000174405.13"; gene_type "protein_coding";
# gene_status "KNOWN"; gene_name "LIG4"; level 2; havana_gene "OTTHUMG00000017328.5";

#=========================================================================
# GTF, used for annotations (GENCODE below)
#-------------------------------------------------------------------------
#  1	chromosome_name         chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M}
#  2	annotation_source	    {ENSEMBL,HAVANA}
#  3	feature_type	        {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
#  4	genomic_start_location	integer-value (1-based)
#  5	genomic_end_location	integer-value
#  6	score_(not used) 	    .
#  7	genomic_strand	        {+,-}
#  8	genomic_phase_(for CDS) {0,1,2,.}
#  9	additional-information as key-value pairs, see below
#-------------------------------------------------------------------------
#  gene_id              ENSGXXXXXXXXXXX *
#  transcript_id        ENSTXXXXXXXXXXX *
#  gene_type            list of biotypes
#  gene_status          {KNOWN, NOVEL, PUTATIVE}
#  gene_name            string
#  transcript_type      list of biotypes
#  transcript_status    {KNOWN, NOVEL, PUTATIVE}
#  transcript_name      string
#  exon_number          indicates the biological position of the exon in the transcript
#  exon_id              ENSEXXXXXXXXXXX *
#  level                1 (verified loci),
#                       2 (manually annotated loci),
#                       3 (automatically annotated loci)
#=========================================================================
