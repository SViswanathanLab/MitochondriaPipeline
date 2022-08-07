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
        expand("results/AlignAndMarkDuplicates/{tumor}/{tumor}.metrics", tumor=config["pairings"]),
        expand("results/AlignShiftedMTAndMarkDuplicates/{tumor}/{tumor}_mba.bam", tumor=config["pairings"]),
        expand("results/AlignShiftedMTAndMarkDuplicates/{tumor}/{tumor}_md.bam", tumor=config["pairings"]),
        expand("results/AlignShiftedMTAndMarkDuplicates/{tumor}/{tumor}.bam", tumor=config["pairings"]),
        expand("results/AlignShiftedMTAndMarkDuplicates/{tumor}/{tumor}.metrics", tumor=config["pairings"]),
        expand("results/CollectWgsMetrics/{tumor}/metrics.txt", tumor=config["pairings"]),
        expand("results/CollectWgsMetrics/{tumor}/theoretical_sensitivity.txt", tumor=config["pairings"]),
        expand("results/CollectWgsMetrics/{tumor}/mean_coverage.txt", tumor=config["pairings"]),
        expand("results/CollectWgsMetrics/{tumor}/median_coverage.txt", tumor=config["pairings"]),
        expand("results/CallMt/{tumor}/{tumor}.vcf.gz", tumor=config["pairings"]),
        expand("results/CallMt/{tumor}/{tumor}.vcf.gz.tbi", tumor=config["pairings"]),
        expand("results/CallMt/{tumor}/{tumor}.vcf.gz.stats", tumor=config["pairings"]),
        expand("results/CallMt/{tumor}/{tumor}_bamout.bam", tumor=config["pairings"])
        
        
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
        mt_ref = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta",
        mt_ref_index = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.fai",
        mt_ref_dict = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.dict",
        mt_sa = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.sa",
        mt_amb = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.amb",
        mt_ann = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.ann",
        mt_bwt = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.bwt",
        mt_pac = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.pac"
    output:
        bwa_log = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.bwa.stderr.log",
        mba_bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}_mba.bam",
        md_bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}_md.bam",
        bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.bam",
        metrics = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.metrics",
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
         PROGRAM_GROUP_VERSION="$({params.bwa} 2>&1 | grep -e '^Version' | sed 's/Version: //')" \
         PROGRAM_GROUP_COMMAND_LINE="{params.bwa} mem -K 100000000 -p -v 3 -t 2 -Y {params.reference_genome}" \
         PROGRAM_GROUP_NAME="bwamem" \
         UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
         ALIGNER_PROPER_PAIR_FLAGS=true \
         UNMAP_CONTAMINANT_READS=true \
         ADD_PG_TAG_TO_READS=false
         
         {params.java} -Xms4000m -jar {params.picard_jar} \
         MarkDuplicates \
         INPUT={output.mba_bam} \
         OUTPUT={output.md_bam} \
         METRICS_FILE={output.metrics} \
         VALIDATION_STRINGENCY=SILENT \
         OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
         ASSUME_SORT_ORDER="queryname" \
         CLEAR_DT="false" \
         ADD_PG_TAG_TO_READS=false

         {params.java} -Xms4000m -jar {params.picard_jar} \
         SortSam \
         INPUT={output.md_bam} \
         OUTPUT={output.bam} \
         SORT_ORDER="coordinate" \
         CREATE_INDEX=true \
         MAX_RECORDS_IN_RAM=300000) 2> {log}"""
         
rule AlignShiftedMTAndMarkDuplicates:
    input:
        bam = "results/RevertSam/{tumor}/{tumor}.bam",
        mt_ref = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta",
        mt_ref_index = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai",
        mt_ref_dict = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict",
        mt_sa = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa",
        mt_amb = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb",
        mt_ann = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann",
        mt_bwt = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt",
        mt_pac = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac"
    output:
        bwa_log = "results/AlignShiftedMTAndMarkDuplicates/{tumor}/{tumor}.bwa.stderr.log",
        mba_bam = "results/AlignShiftedMTAndMarkDuplicates/{tumor}/{tumor}_mba.bam",
        md_bam = "results/AlignShiftedMTAndMarkDuplicates/{tumor}/{tumor}_md.bam",
        bam = "results/AlignShiftedMTAndMarkDuplicates/{tumor}/{tumor}.bam",
        metrics = "results/AlignShiftedMTAndMarkDuplicates/{tumor}/{tumor}.metrics",
    params:
        bwa = config["bwa"],
        java = config["java"],
        picard_jar = config["picard_jar"],
        gatk = config["gatk_path"],
        reference_genome = config["mt_shifted_ref"]
    log:
        "logs/AlignShiftedMTAndMarkDuplicates/{tumor}.txt"
    shell:
        """(set -o pipefail
         set -e
         
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
         PROGRAM_GROUP_VERSION="$({params.bwa} 2>&1 | grep -e '^Version' | sed 's/Version: //')" \
         PROGRAM_GROUP_COMMAND_LINE="{params.bwa} mem -K 100000000 -p -v 3 -t 2 -Y {params.reference_genome}" \
         PROGRAM_GROUP_NAME="bwamem" \
         UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
         ALIGNER_PROPER_PAIR_FLAGS=true \
         UNMAP_CONTAMINANT_READS=true \
         ADD_PG_TAG_TO_READS=false
         
         {params.java} -Xms4000m -jar {params.picard_jar} \
         MarkDuplicates \
         INPUT={output.mba_bam} \
         OUTPUT={output.md_bam} \
         METRICS_FILE={output.metrics} \
         VALIDATION_STRINGENCY=SILENT \
         OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
         ASSUME_SORT_ORDER="queryname" \
         CLEAR_DT="false" \
         ADD_PG_TAG_TO_READS=false

         {params.java} -Xms4000m -jar {params.picard_jar} \
         SortSam \
         INPUT={output.md_bam} \
         OUTPUT={output.bam} \
         SORT_ORDER="coordinate" \
         CREATE_INDEX=true \
         MAX_RECORDS_IN_RAM=300000) 2> {log}"""

rule CollectWgsMetrics:
    input:
        bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.bam",
        bai = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.bai"
    output:
        metrics = "results/CollectWgsMetrics/{tumor}/metrics.txt",
        theoretical_sensitivity = "results/CollectWgsMetrics/{tumor}/theoretical_sensitivity.txt",
        mean_coverage = "results/CollectWgsMetrics/{tumor}/mean_coverage.txt",
        median_coverage = "results/CollectWgsMetrics/{tumor}/median_coverage.txt"
    params:
        coverage_cap = config["coverage_cap"],
        mt_ref = config["mt_ref"],
        mt_ref_index = config["mt_ref_index"],
        java = config["java"],
        picard_jar = config["picard_jar"],
        read_length_for_optimization = config["read_length_for_optimization"],
        coverage_script = config["coverage_script"]
    log:
        "logs/CollectWgsMetrics/{tumor}.txt"
    shell:
        """(set -e

        {params.java} -Xms2000m -jar {params.picard_jar} \
        CollectWgsMetrics \
        INPUT={input.bam} \
        VALIDATION_STRINGENCY=SILENT \
        REFERENCE_SEQUENCE={params.mt_ref} \
        OUTPUT={output.metrics} \
        USE_FAST_ALGORITHM=true \
        READ_LENGTH={params.read_length_for_optimization} \
        COVERAGE_CAP= {params.coverage_cap}\
        INCLUDE_BQ_HISTOGRAM=true \
        THEORETICAL_SENSITIVITY_OUTPUT={output.theoretical_sensitivity}
        
        Rscript {params.coverage_script} --input {output.metrics} --mean_output {output.mean_coverage} --median_output {output.median_coverage}) 2> {log}"""

rule CallMt:
    input:
        bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.bam"
    output:
        vcf = protected("results/CallMt/{tumor}/{tumor}.vcf.gz"),
        tbi = protected("results/CallMt/{tumor}/{tumor}.vcf.gz.tbi"),
        stats = protected("results/CallMt/{tumor}/{tumor}.vcf.gz.stats"),
        bamout = protected("results/CallMt/{tumor}/{tumor}_bamout.bam")
    params:
        gatk = config["gatk_path"],
        max_reads_per_alignment_start = config["max_reads_per_alignment_start"],
        mt_ref = config["mt_ref"]
    log:
    shell:
        """(set -e

        touch {output.bamout}

        {params.gatk} --java-options "-Xmx2000m" Mutect2 \
        -R {params.mt_ref} \
        -I {input.bam} \
        -L chrM:576-16024 \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O {output.vcf} \
        --annotation StrandBiasBySample \
        --bam-output {output.bamout} \
        --mitochondria-mode \
        --max-reads-per-alignment-start {params.max_reads_per_alignment_start} \
        --max-mnp-distance 0) 2> {log}"""

