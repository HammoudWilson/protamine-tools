# action:
#     create genome regions BED files to be used for ATAC-seq enrichment analysis
# organization:
#     both this folder (pipelines/enrichment/bed) and the output folder ($TASK_DIR/regions_bed)
#     have subfolders each defining one _class_ of related BED region types, e.g., transcription, histone, etc.
#       - each subfolder of pipelines/enrichment/bed contains one or more scripts that each create one BED file _type_
#       - each subfolder of $TASK_DIR/regions_bed contains output BED files created by the scripts 
#            in the corresponding subfolder of pipelines/enrichment/bed
#     !!! BED region class/types must be specified in region_types.txt file in this folder to be created !!!
# example:
#     BED creation script:
#       pipelines/enrichment/bed/transcription/tss_active.sh
#     creates the content written to output BED file:
#       $TASK_DIR/regions_bed/transcription/tss_active.bed
#     where 'transcription' is the regions _class_ and 'tss_active' is the regions _type_
#     and the entry in region_types.txt would be 'transcription/tss_active'
# template:
#     use the template.sh script in this folder for constructing a new BED creation script
# inputs:
#     input files are defined in the BED creation scripts
#     follow the template.sh script to name and excerpt the input content to the job log
# outputs:
#     BED creation scripts must write their output to STDOUT
#     this script capture, sorts, and writes that content to:
#        $REGION_TYPE_BED == $REGIONS_BED_DIR/<region_class>/<region_type>.bed
#     where:
#        $REGIONS_BED_DIR was set upstream to $TASK_DIR/regions_bed
#        $REGION_TYPE_BED is set here prior to calling the BED creation script
#     along the way, it ensures that the output BED file only retains included regions,
#     merges overlapping regions, and reports genome and region bp counts
#     thus, the fraction of the genome in the regions = region bp / included bp (not total genome bp)

# create the primary genome inclusions file as needed
if [ ! -f ${PRIMARY_GENOME_INCLUSIONS_BED} ]; then
    echo "extracting primary genome inclusion regions file"
    cat ${GENOME_INCLUSIONS_BED} |
    awk '$1 ~ /-'${PRIMARY_GENOME}'/' |
    sed 's/-'${PRIMARY_GENOME}'//' > ${PRIMARY_GENOME_INCLUSIONS_BED}
    checkPipe
    head ${PRIMARY_GENOME_INCLUSIONS_BED}
fi

# report genome metrics to the job log
echo
echo "primary genome: ${PRIMARY_GENOME}"
echo "total bp: ${PRIMARY_GENOME_SIZE}"
echo "included bp: "`cat ${PRIMARY_GENOME_INCLUSIONS_BED} | awk '{sum += $3 - $2} END {print sum}'`
echo

# collect the work to do
cd ${ACTION_DIR}
REGION_TYPES=`cat region_types.txt | grep -v '^#' | grep -v '^$'`

# create regione BED files as needed
for REGION_TYPE in ${REGION_TYPES}; do
    echo
    echo "--------------------------------------------------------------------------------"
    echo "${REGION_TYPE}"
    echo "--------------------------------------------------------------------------------"

    REGION_CLASS=`echo ${REGION_TYPE} | cut -d'/' -f1`
    REGION_CLASS_DIR=${REGIONS_BED_DIR}/${REGION_CLASS}
    mkdir -p ${REGION_CLASS_DIR}
    REGIONS_TYPE_BED=${REGIONS_BED_DIR}/${REGION_TYPE}.bed # includes class in path
    if [[ ! -f ${REGIONS_TYPE_BED} || ${FORCE_BED} == ${REGION_TYPE} || ${FORCE_ALL_BED} != "0" ]]; then
        BED_SCRIPT=${REGION_TYPE}.sh
        source ${BED_SCRIPT} | 
        {
            read HEADER
            echo "${HEADER}" | 
            cut -f 1-5 | 
            awk 'BEGIN{OFS = "\t"}{ print $1, $2, $3, $4, (NF==5 ? $5 : "IGNORE")}' > ${REGIONS_TYPE_BED}
            cut -f 1-5 |
            awk 'BEGIN{OFS = "\t"}{ print $1, $2, $3, $4, (NF==5 ? $5 : 0)}' |
            bedtools intersect -a - -b ${PRIMARY_GENOME_INCLUSIONS_BED} | 
            sort -k1,1 -k2,2n -k3,3n |
            bedtools merge -c 4,5 -o distinct,max -i - >> ${REGIONS_TYPE_BED} 
        }
        checkPipe
        echo
        echo "output file: ${REGIONS_TYPE_BED}"
        echo `cat ${REGIONS_TYPE_BED} | wc -l`" lines (regions plus header)"
        echo "region bp: "`cat ${REGIONS_TYPE_BED} | awk '$1 != "chrom" {sum += $3 - $2} END {print sum}'`
        head ${REGIONS_TYPE_BED}
    else
        echo "BED file exists, skipping; use --force-bed or --force-all-bed to overwrite"
        continue
    fi
done

echo
echo "done"
