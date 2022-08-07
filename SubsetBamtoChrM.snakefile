import os

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
        expand("results/SubsetBamtoChrM/{tumor}/{tumor}.bam", tumor=config["pairings"]),
        expand("results/RevertSam/{tumor}/{tumor}.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}/{tumor}_mba.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}/{tumor}_md.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}/{tumor}.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}/{tumor}.metrics", tumor=config["pairings"])
        
rule SubsetBamtoChrM:
    input:
        tumor_filepath = lambda wildcards: config["samples"][wildcards.tumor],
        #normal_filepath = lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
    output:
        bam = "results/SubsetBamtoChrM/{tumor}/{tumor}.bam",
        bai = "results/SubsetBamtoChrM/{tumor}/{tumor}.bai"
    params:
        gatk = config["gatk_path"],
        contig_name = config["contig_name"]
    log:
        "logs/SubsetBamtoChrM/{tumor}.txt"
    shell:
        "({params.gatk} PrintReads \
        -L {params.contig_name} \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -I {input.tumor_filepath} \
        -O {output.bam}) 2> {log}"

rule RevertSam:
    input:
        bam = "results/SubsetBamtoChrM/{tumor}/{tumor}.bam"
    output:
        bam = "results/RevertSam/{tumor}/{tumor}.bam",
    params:
        java = config["java"],
        picard_jar = config["picard_jar"]
    log:
        "logs/RevertSam/{tumor}.txt"
    shell:
        "({params.java} -Xmx1000m -jar {params.picard_jar} \
        RevertSam \
        INPUT={input.bam} \
        OUTPUT_BY_READGROUP=false \
        OUTPUT={output.bam} \
        VALIDATION_STRINGENCY=LENIENT \
        ATTRIBUTE_TO_CLEAR=FT \
        ATTRIBUTE_TO_CLEAR=CO \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=false) 2> {log}"

rule AlignAndMarkDuplicates:
    input:
        bam = "results/RevertSam/{tumor}/{tumor}.bam",
    output:
        bwa_log = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.bwa.stderr.log",
        mba_bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}_mba.bam",
        md_bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}_md.bam",
        bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.bam",
        metrics = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.metrics"
    params:
        bwa = config["bwa"],
        java = config["java"],
        picard_jar = config["picard_jar"],
        gatk = config["gatk_path"],
        reference_genome = config["mt_ref"]
    log:
        "logs/AlignAndMarkDuplicates/{tumor}.txt"
    shell:
        """(set -o pipefail
         set -e
         
         bwa_version=$({params.bwa} 2>&1 | grep -e '^Version' | sed 's/Version: //')
         
         {params.java} -Xms5000m -jar {params.picard_jar} \
         SamToFastq \
         INPUT={input.bam} \
         FASTQ=/dev/stdout \
         INTERLEAVE=true \
         NON_PF=true | \
         {params.bwa} mem -K 100000000 -p -v 3 -t 2 -Y {params.reference_genome} /dev/stdin - 2> >(tee {output.bwa_log} >&2) | \
         {params.java} -Xms3000m -jar {params.picard_jar} \
         MergeBamAlignment \
         VALIDATION_STRINGENCY=SILENT \
         EXPECTED_ORIENTATIONS=FR \
         ATTRIBUTES_TO_RETAIN=X0 \
         ATTRIBUTES_TO_REMOVE=NM \
         ATTRIBUTES_TO_REMOVE=MD \
         ALIGNED_BAM=/dev/stdin \
         UNMAPPED_BAM={input.bam} \
         OUTPUT={output.mba_bam} \
         REFERENCE_SEQUENCE={params.reference_genome} \
         PAIRED_RUN=true \
         SORT_ORDER="unsorted" \
         IS_BISULFITE_SEQUENCE=false \
         ALIGNED_READS_ONLY=false \
         CLIP_ADAPTERS=false \
         MAX_RECORDS_IN_RAM=2000000 \
         ADD_MATE_CIGAR=true \
         MAX_INSERTIONS_OR_DELETIONS=-1 \
         PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
         PROGRAM_RECORD_ID="bwamem" \
         PROGRAM_GROUP_VERSION=$({params.bwa} 2>&1 | grep -e '^Version' | sed 's/Version: //') \
         PROGRAM_GROUP_COMMAND_LINE="{params.bwa} mem -K 100000000 -p -v 3 -t 2 -Y {params.reference_genome}" \
         PROGRAM_GROUP_NAME="bwamem" \
         UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
         ALIGNER_PROPER_PAIR_FLAGS=true \
         UNMAP_CONTAMINANT_READS=true \
         ADD_PG_TAG_TO_READS=false) 2> {log}"""
         
