#!/bin/bash

# common UCSC URLs
UCSC_GOLDEN_PATH=https://hgdownload.soe.ucsc.edu/goldenPath
HG38_ANALYSIS_SET="${UCSC_GOLDEN_PATH}/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz"

# shared download action support
download_and_unzip_file () {
    if [ ! -f $TO_FILE ]; then
        if [ ! -f $GZ_FILE ]; then
            echo "  downloading $FILE_DESCRIPTION"
            echo "    from: $URL"
            echo "    to:   $GZ_FILE"
            wget -O $GZ_FILE $URL
        fi
        if [ "$LEAVE_ZIPPED" == "" ]; then
            echo "  unzipping $FILE_DESCRIPTION"
            gunzip $GZ_FILE
        fi
    else 
        echo "  already exists: $FILE_DESCRIPTION"
    fi    
}
download_and_unzip_UCSC () {
    UCSC_FILE=$1
    TO_FILE=$2
    FILE_DESCRIPTION=$3
    LEAVE_ZIPPED=$4
    if [ "$LEAVE_ZIPPED" == "" ]; then
        GZ_FILE=${TO_FILE}.gz
    else 
        GZ_FILE=${TO_FILE}
    fi
    URL=${UCSC_GOLDEN_PATH}/${GENOME}/${UCSC_FILE}
    if [ "$FORCE_URL" != "" ]; then
        URL=${FORCE_URL}
        FORCE_URL=""
    fi
    download_and_unzip_file
}
download_and_index_genome () {
    FORCE_URL=$1 # to support downloading an analysis set instead of just the genome
    download_and_unzip_UCSC \
        bigZips/${GENOME}.fa.gz \
        ${GENOME_FASTA} \
        "genome fasta"
    if [ ! -f $GENOME_FASTA.fai ]; then
        if [ "$APPEND_EBV" != "" ]; then
            FORCE_URL=${HG38_ANALYSIS_SET}
            WRK_DIR=`dirname ${GENOME_FASTA}`
            HG38_FILE=${WRK_DIR}/hg38.fa
            EBV_FILE=${WRK_DIR}/chrEBV.fa
            if [ ! -f ${EBV_FILE} ]; then # get the EBV genome from the hg38 analysis set
                download_and_unzip_UCSC \
                    NA \
                    ${HG38_FILE} \
                    "hg38 analysis set"
                echo "  indexing analysis set"
                samtools faidx ${HG38_FILE}
                echo "  extracting chrEBV"
                samtools faidx ${HG38_FILE} chrEBV > ${EBV_FILE}
                rm ${HG38_FILE}
                rm ${HG38_FILE}.fai # leave chrEBV.fa but remove the larger, now obsolete hg38 file
            fi
            echo "  appending chrEBV"
            cat ${EBV_FILE} >> ${GENOME_FASTA}
        fi
        echo "  indexing genome fasta"
        samtools faidx $GENOME_FASTA
    else 
        echo "  already exists: genome fasta index"
    fi
}
download_genome_gaps () {
    download_and_unzip_UCSC \
        database/gap.txt.gz \
        ${GENOME_GAPS_FILE} \
        "gaps file"
}
download_gc5Base () {
    download_and_unzip_UCSC \
        bigZips/${GENOME}.gc5Base.wigVarStep.gz \
        ${GENOME_GC5BASE_WIG} \
        "gc5Base file" \
        LEAVE_ZIPPED
}
construct_gc5Base () { # for newer UCSC genomes that do not offer gc5Base.wigVarStep.gz
    TO_FILE=${GENOME_GC5BASE_WIG}
    if [ ! -f $TO_FILE ]; then
        echo "  pending: gc5Base construction for genome $GENOME"
    fi
}
download_ENCODE_exclusions () {
    TO_FILE=${GENOME_EXCLUSIONS_BED}
    GZ_FILE=${TO_FILE}.gz
    FILE_DESCRIPTION="excluded regions"
    LEAVE_ZIPPED=""
    URL=https://github.com/Boyle-Lab/Blacklist/raw/master/lists/${GENOME}-blacklist.v2.bed.gz
    download_and_unzip_file
}
download_excluderegions () { # for newer genomes lacking Boyle lab exclusions
    REMOTE_FILE=$1
    TO_FILE=${GENOME_EXCLUSIONS_BED}
    if [ ! -f $TO_FILE ]; then
        echo
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "at present, the best choice for obtaining exclusion regions for genome ${GENOME}"
        echo "is to manually download the bed format (not bed.gz format) file:"
        echo "    $REMOTE_FILE"
        echo "from Google Drive:"
        echo "    https://drive.google.com/drive/folders/1sF9m8Y3eZouTZ3IEEywjs2kfHOWFBSJT/"
        echo "to this exact path and file name:"
        echo "    $TO_FILE"
        echo "for details on the source of these files, see:"
        echo "    https://github.com/dozmorovlab/excluderanges"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo
    fi
}
download_gtf () {
    SPECIES=$1
    RELEASE=$2
    VERSION=v${RELEASE}
    GTF_PATH=Gencode_${SPECIES}/release_${RELEASE}/gencode.${VERSION}.basic.annotation.gtf.gz
    GENCODE_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/${GTF_PATH}
    if [ ! -f ${GENOME_GTF} ]; then
        echo "  downloading GTF file"
        echo "    from: ${GENCODE_URL}"
        echo "    to:   ${GENOME_GTF}"
        wget -O ${GENOME_GTF} ${GENCODE_URL}
    else 
        echo "  already exists: genome GTF"
    fi
    if [ ! -f ${GENES_BED} ]; then
        echo "  writing genes BED file"
        echo "    from: ${GENOME_GTF}"
        echo "    to:   ${GENES_BED}"
        zcat ${GENOME_GTF} | 
        perl -n ${ACTION_DIR}/gtf_to_genes_bed.pl | 
        sort -k1,1 -k2,2n -k3,3n |
        pigz -c > ${GENES_BED}
    else 
        echo "  already exists: genes BED"
    fi
}

# execute download actions customized for each specific supported genome
echo
echo "downloading resource files for reference genome ${GENOME}"
if [ "$GENOME" == "hs1" ]; then
    APPEND_EBV="TRUE"
    download_and_index_genome ""
    echo "  touching genome gaps (hs1/CHM13 has no gaps!)"
    touch ${GENOME_GAPS_FILE} # CHM13 has no gaps! create the empty file anyway for consistency
    construct_gc5Base
    download_excluderegions T2T.excluderanges.bed

elif [ "$GENOME" == "hg38" ]; then
    APPEND_EBV=""
    download_and_index_genome ${HG38_ANALYSIS_SET}
    download_genome_gaps
    download_gc5Base
    download_ENCODE_exclusions
    download_gtf human ${GENCODE_RELEASE} # e.g., 47
elif [ "$GENOME" == "mm39" ]; then
    APPEND_EBV=""
    download_and_index_genome ""
    download_genome_gaps
    construct_gc5Base
    download_excluderegions mm39.excluderanges.bed
    download_gtf mouse ${GENCODE_RELEASE} # e.g., M36
elif [ "$GENOME" == "dm6" ]; then
    APPEND_EBV=""
    download_and_index_genome ""
    download_genome_gaps
    construct_gc5Base
    download_ENCODE_exclusions
else
    echo "unsupported genome: ${GENOME}"
    exit 1
fi

echo
echo "done"
