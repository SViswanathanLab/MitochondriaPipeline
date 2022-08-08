import os

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
        expand("results/SubsetBamtoChrM/{tumor}.bam", tumor=config["pairings"]),
        expand("results/RevertSam/{tumor}.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}_mba.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}_md.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}.metrics", tumor=config["pairings"]),
        expand("results/AlignShiftedMTAndMarkDuplicates/{tumor}_mba.bam", tumor=config["pairings"]),
        expand("results/AlignShiftedMTAndMarkDuplicates/{tumor}_md.bam", tumor=config["pairings"]),
        expand("results/AlignShiftedMTAndMarkDuplicates/{tumor}.bam", tumor=config["pairings"]),
        expand("results/AlignShiftedMTAndMarkDuplicates/{tumor}.metrics", tumor=config["pairings"]),
        expand("results/CollectWgsMetrics/{tumor}_metrics.txt", tumor=config["pairings"]),
        expand("results/CollectWgsMetrics/{tumor}_theoretical_sensitivity.txt", tumor=config["pairings"]),
        expand("results/CollectWgsMetrics/{tumor}_mean_coverage.txt", tumor=config["pairings"]),
        expand("results/CollectWgsMetrics/{tumor}_median_coverage.txt", tumor=config["pairings"]),
        expand("results/CallMt/{tumor}.vcf.gz", tumor=config["pairings"]),
        expand("results/CallMt/{tumor}.vcf.gz.tbi", tumor=config["pairings"]),
        expand("results/CallMt/{tumor}.vcf.gz.stats", tumor=config["pairings"]),
        expand("results/CallMt/{tumor}_bamout.bam", tumor=config["pairings"]),
        expand("results/CallShiftedMt/{tumor}.vcf.gz", tumor=config["pairings"]),
        expand("results/CallShiftedMt/{tumor}.vcf.gz.tbi", tumor=config["pairings"]),
        expand("results/CallShiftedMt/{tumor}.vcf.gz.stats", tumor=config["pairings"]),
        expand("results/CallShiftedMt/{tumor}_bamout.bam", tumor=config["pairings"]),
        expand("results/LiftoverAndCombineVcfs/{tumor}.shifted_back.vcf", tumor=config["pairings"]),
        expand("results/LiftoverAndCombineVcfs/{tumor}.rejected.vcf", tumor=config["pairings"]),
        expand("results/LiftoverAndCombineVcfs/{tumor}.merged.vcf", tumor=config["pairings"]),
        expand("results/MergeStats/{tumor}.raw.combined.stats", tumor=config["pairings"]),
        expand("results/InitialFilter/{tumor}.filtered.vcf", tumor=config["pairings"]),
        expand("results/InitialFilter/{tumor}.vcf", tumor=config["pairings"]),
        expand("results/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_splitAndPassOnly.vcf", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_headers.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_output_data.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_contamination.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_major_hg.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_minor_hg.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_mean_het_major.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_mean_het_minor.txt", tumor=config["pairings"]),
        expand("results/FilterContamination/{tumor}.vcf", tumor=config["pairings"]),
        expand("results/FilterContamination/{tumor}.filtered.vcf", tumor=config["pairings"]),
        expand("results/CoverageAtEveryBase/{tumor}_per_base_coverage.tsv", tumor=config["pairings"]),
        expand("results/CoverageAtEveryBase/{tumor}_non_control_region.metrics", tumor=config["pairings"]),
        expand("results/CoverageAtEveryBase/{tumor}_control_region_shifted.metrics", tumor=config["pairings"]),
        expand("results/CoverageAtEveryBase/{tumor}_non_control_region.tsv", tumor=config["pairings"]),
        expand("results/CoverageAtEveryBase/{tumor}_control_region_shifted.tsv", tumor=config["pairings"]),
        expand("results/SplitMultiAllelicSites/{tumor}.vcf", tumor=config["pairings"])
        
               
rule SubsetBamtoChrM:
    input:
        tumor_filepath = lambda wildcards: config["samples"][wildcards.tumor],
        #normal_filepath = lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
    output:
        bam = "results/SubsetBamtoChrM/{tumor}.bam",
        bai = "results/SubsetBamtoChrM/{tumor}.bai"
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
        bam = "results/SubsetBamtoChrM/{tumor}.bam"
    output:
        bam = "results/RevertSam/{tumor}.bam",
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
        bam = "results/RevertSam/{tumor}.bam",
        mt_ref = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta",
        mt_ref_index = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.fai",
        mt_ref_dict = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.dict",
        mt_sa = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.sa",
        mt_amb = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.amb",
        mt_ann = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.ann",
        mt_bwt = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.bwt",
        mt_pac = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.fasta.pac"
    output:
        bwa_log = "results/AlignAndMarkDuplicates/{tumor}.bwa.stderr.log",
        mba_bam = "results/AlignAndMarkDuplicates/{tumor}_mba.bam",
        md_bam = "results/AlignAndMarkDuplicates/{tumor}_md.bam",
        bam = "results/AlignAndMarkDuplicates/{tumor}.bam",
        bai = "results/AlignAndMarkDuplicates/{tumor}.bai",
        metrics = "results/AlignAndMarkDuplicates/{tumor}.metrics",
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
        bam = "results/RevertSam/{tumor}.bam",
        mt_ref = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta",
        mt_ref_index = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai",
        mt_ref_dict = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict",
        mt_sa = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa",
        mt_amb = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb",
        mt_ann = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann",
        mt_bwt = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt",
        mt_pac = "/mnt/storage/labs/sviswanathan/mt_gs_files/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac"
    output:
        bwa_log = "results/AlignShiftedMTAndMarkDuplicates/{tumor}.bwa.stderr.log",
        mba_bam = "results/AlignShiftedMTAndMarkDuplicates/{tumor}_mba.bam",
        md_bam = "results/AlignShiftedMTAndMarkDuplicates/{tumor}_md.bam",
        bam = "results/AlignShiftedMTAndMarkDuplicates/{tumor}.bam",
        bai = "results/AlignShiftedMTAndMarkDuplicates/{tumor}.bai",
        metrics = "results/AlignShiftedMTAndMarkDuplicates/{tumor}.metrics",
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
        bam = "results/AlignAndMarkDuplicates/{tumor}.bam",
        bai = "results/AlignAndMarkDuplicates/{tumor}.bai"
    output:
        metrics = "results/CollectWgsMetrics/{tumor}_metrics.txt",
        theoretical_sensitivity = "results/CollectWgsMetrics/{tumor}_theoretical_sensitivity.txt",
        mean_coverage = "results/CollectWgsMetrics/{tumor}_mean_coverage.txt",
        median_coverage = "results/CollectWgsMetrics/{tumor}_median_coverage.txt"
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
        bam = "results/AlignAndMarkDuplicates/{tumor}.bam",
        bai = "results/AlignAndMarkDuplicates/{tumor}.bai"
    output:
        vcf = protected("results/CallMt/{tumor}.vcf.gz"),
        tbi = protected("results/CallMt/{tumor}.vcf.gz.tbi"),
        stats = protected("results/CallMt/{tumor}.vcf.gz.stats"),
        bamout = protected("results/CallMt/{tumor}_bamout.bam")
    params:
        gatk = config["gatk_path"],
        max_reads_per_alignment_start = config["max_reads_per_alignment_start"],
        mt_ref = config["mt_ref"]
    log:
        "logs/CallMt/{tumor}.txt"
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

rule CallShiftedMt:
    input:
        bam = "results/AlignShiftedMTAndMarkDuplicates/{tumor}.bam",
        bai = "results/AlignShiftedMTAndMarkDuplicates/{tumor}.bai"
    output:
        vcf = protected("results/CallShiftedMt/{tumor}.vcf.gz"),
        tbi = protected("results/CallShiftedMt/{tumor}.vcf.gz.tbi"),
        stats = protected("results/CallShiftedMt/{tumor}.vcf.gz.stats"),
        bamout = protected("results/CallShiftedMt/{tumor}_bamout.bam")
    params:
        gatk = config["gatk_path"],
        max_reads_per_alignment_start = config["max_reads_per_alignment_start"],
        mt_ref = config["mt_shifted_ref"]
    log:
        "logs/CallShiftedMt/{tumor}.txt"
    shell:
        """(set -e

        touch {output.bamout}

        {params.gatk} --java-options "-Xmx2000m" Mutect2 \
        -R {params.mt_ref} \
        -I {input.bam} \
        -L chrM:8025-9144 \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O {output.vcf} \
        --annotation StrandBiasBySample \
        --bam-output {output.bamout} \
        --mitochondria-mode \
        --max-reads-per-alignment-start {params.max_reads_per_alignment_start} \
        --max-mnp-distance 0) 2> {log}"""

rule LiftoverAndCombineVcfs:
    input:
        vcf = "results/CallMt/{tumor}.vcf.gz",
        shifted_vcf = "results/CallShiftedMt/{tumor}.vcf.gz"
    output:
        shifted_back = "results/LiftoverAndCombineVcfs/{tumor}.shifted_back.vcf",
        rejected = "results/LiftoverAndCombineVcfs/{tumor}.rejected.vcf",
        merged = "results/LiftoverAndCombineVcfs/{tumor}.merged.vcf"
    params:
        mt_ref = config["mt_ref"],
        java = config["java"],
        picard_jar = config["picard_jar"],
        shift_back_chain = config["shift_back_chain"]
    log:
        "logs/LiftoverAndCombineVcfs/{tumor}.txt"
    shell:
        """(set -e

        {params.java} -jar {params.picard_jar} LiftoverVcf \
        I={input.shifted_vcf} \
        O={output.shifted_back} \
        R={params.mt_ref} \
        CHAIN={params.shift_back_chain} \
        REJECT={output.rejected}

        {params.java} -jar {params.picard_jar} MergeVcfs \
        I={input.shifted_vcf} \
        I={input.vcf} \
        O={output.merged}) 2> {log}"""
 
rule MergeStats:
    input:
        shifted_stats = "results/CallShiftedMt/{tumor}.vcf.gz.stats",
        non_shifted_stats = "results/CallMt/{tumor}.vcf.gz.stats"
    output:
        raw_combined_stats = "results/MergeStats/{tumor}.raw.combined.stats"
    params:
        gatk = config["gatk_path"]
        
    log:
        "logs/MergeStats/{tumor}.txt"
    shell:
        """(set -e

        {params.gatk} MergeMutectStats --stats {input.shifted_stats} --stats {input.non_shifted_stats} -O {output.raw_combined_stats}) 2> {log}"""

rule InitialFilter:
    input:
        raw_vcf = "results/LiftoverAndCombineVcfs/{tumor}.merged.vcf",
        raw_vcf_stats = "results/MergeStats/{tumor}.raw.combined.stats"
    output:
        filtered_vcf = "results/InitialFilter/{tumor}.filtered.vcf",
        output_vcf = "results/InitialFilter/{tumor}.vcf",
        bamout = ""
    params:
        gatk = config["gatk_path"],
        mt_ref = config["mt_ref"],
        max_alt_allele_count = config["max_alt_allele_count"],
        vaf_filter_threshold = config["vaf_filter_threshold"],
        blacklisted_sites = config["blacklisted_sites"]
    log:
        "logs/InitialFilter/{tumor}.txt"
    shell:
        """(set -e
        
        # We need to create these files regardless, even if they stay empty
        touch {output.bamout}

        {params.gatk} --java-options "-Xmx2500m" FilterMutectCalls \
        -V {input.raw_vcf} \
        -R {params.mt_ref} \
        -O {output.filtered_vcf} \
        --stats {input.raw_vcf_stats} \
        --max-alt-allele-count {params.max_alt_allele_count} \
        --mitochondria-mode \
        --min-allele-fraction {params.vaf_filter_threshold} 

        {params.gatk} VariantFiltration \
        -V {output.filtered_vcf}\
        -O {output.output_vcf} \
        --apply-allele-specific-filters \
        --mask {params.blacklisted_sites} \
        --mask-name "blacklisted_site") 2> {log}"""

rule SplitMultiAllelicsAndRemoveNonPassSites:
    input:
        filtered_vcf = "results/InitialFilter/{tumor}.vcf"
    output:
        split_vcf = temp("results/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_split.vcf"),
        split_vcf_idx = temp("results/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_split.vcf.idx"),
        splitAndPassOnly_vcf = "results/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_splitAndPassOnly.vcf"
    log:
        "logs/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}.txt"
    params:
        gatk = config["gatk_path"],
        mt_ref = config["mt_ref"]
    shell:
        """(set -e

        {params.gatk} LeftAlignAndTrimVariants \
        -R {params.mt_ref} \
        -V {input.filtered_vcf} \
        -O {output.split_vcf} \
        --split-multi-allelics \
        --dont-trim-alleles \
        --keep-original-ac

        {params.gatk} SelectVariants \
        -V {output.split_vcf} \
        -O {output.splitAndPassOnly_vcf} \
        --exclude-filtered) 2> {log}"""

rule GetContamination:
    input:
        input_vcf = "results/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_splitAndPassOnly.vcf"
    output:
        outputs = "results/GetContamination/{tumor}/output",
        output_noquotes = "results/GetContamination/{tumor}/output_noquotes",
        headers = "results/GetContamination/{tumor}/{tumor}_headers.txt",
        output_data = "results/GetContamination/{tumor}/{tumor}_output_data.txt",
        contamination = "results/GetContamination/{tumor}/{tumor}_contamination.txt",
        major_hg = "results/GetContamination/{tumor}/{tumor}_major_hg.txt",
        minor_hg = "results/GetContamination/{tumor}/{tumor}_minor_hg.txt",
        mean_het_major = "results/GetContamination/{tumor}/{tumor}_mean_het_major.txt",
        mean_het_minor = "results/GetContamination/{tumor}/{tumor}_mean_het_minor.txt"
    params:
        java = config["java"],
        picard_jar = config["picard_jar"],
        haplocheckCLI_path = config["haplocheckCLI_path"]
    log:
        "logs/GetContamination/{tumor}.txt"
    shell:
        """(set -e
        touch {output.headers}
        touch {output.output_data}
        touch {output.contamination}
        touch {output.major_hg}
        touch {output.minor_hg}
        touch {output.mean_het_major}
        touch {output.mean_het_minor}

        {params.java} -jar {params.haplocheckCLI_path} "$(dirname "{input.input_vcf}")" > {output.outputs}
        sed 's/\\\"//g' {output.outputs} > {output.output-noquotes}
        
        grep "SampleID" {output.output-noquotes} > {output.headers}
        if [ `awk '{{print $2}}' {output.headers}` != \"Contamination\" ]; then
          echo \"Bad contamination file format\"; exit 1
        fi
        if [ `awk '{{print $6}}' {output.headers}` != "HgMajor" ]; then
          echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{{print $8}}' {output.headers}` != "HgMinor" ]; then
          echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{{print $14}}' {output.headers}` != "MeanHetLevelMajor" ]; then
          echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{{print $15}}' {output.headers}` != "MeanHetLevelMinor" ]; then
          echo $FORMAT_ERROR; exit 1
        fi
        
        grep -v "SampleID" {output.output-noquotes}  > {output.output_data}
        awk -F \"\\t\" '{{print $2}}' {output.output_data} > {output.contamination}
        awk -F \"\\t\" '{{print $6}}' {output.output_data} > {output.major_hg}
        awk -F \"\\t\" '{{print $8}}' {output.output_data} > {output.minor_hg}
        awk -F \"\\t\" '{{print $14}}' {output.output_data} > {output.mean_het_major}
        awk -F \"\\t\" '{{print $15}}' {output.output_data} > {output.mean_het_minor}) 2> {log}"""
        
rule FilterContamination:
    input:
        hasContamination = "results/GetContamination/{tumor}/{tumor}_contamination.txt",
        contamination_major = "results/GetContamination/{tumor}/{tumor}_mean_het_major.txt",
        contamination_minor = "results/GetContamination/{tumor}/{tumor}_mean_het_minor.txt",
        filtered_vcf = "results/InitialFilter/{tumor}.vcf",
        raw_vcf_stats = "results/MergeStats/{tumor}.raw.combined.stats"
    output:
        output_vcf = "results/FilterContamination/{tumor}.vcf",
        filtered_vcf = "results/FilterContamination/{tumor}.filtered.vcf",
        bamout = ""
    params:
        gatk = config["gatk_path"],
        mt_ref = config["mt_ref"],
        max_alt_allele_count = config["max_alt_allele_count"],
        vaf_filter_threshold = config["vaf_filter_threshold"],
        blacklisted_sites = config["blacklisted_sites"],
        contamination_flag = config["contamination_flag"]
    log:
        "logs/FilterContamination/{tumor}.txt"
    shell:
        """
        (set -e
        
        if [ {params.contamination_flag} == "NO" ]; then
            exit 1
        fi
        
        # We need to create these files regardless, even if they stay empty
        touch {output.bamout}
        
        CONTAMINATION="$(cat "{input.hasContamination}")"
        MAJOR="$(cat "{input.contamination_major}")"
        MINOR="$(cat "{input.contamination_minor}")"
        if [ $MAJOR == 0.000 ]; then
          max_contamination=$MINOR
        else
          max_contamination=$((1.0-$MAJOR))
        fi
        
        {params.gatk} --java-options "-Xmx2500m" FilterMutectCalls \
        -V {input.filtered_vcf} \
        -R {params.mt_ref} \
        -O {output.filtered_vcf} \
        --stats {input.raw_vcf_stats} \
        --max-alt-allele-count {params.max_alt_allele_count} \
        --mitochondria-mode \
        --contamination-estimate $max_contamination\
        --min-allele-fraction {params.vaf_filter_threshold} 

        {params.gatk} VariantFiltration \
        -V {output.filtered_vcf}\
        -O {output.output_vcf} \
        --apply-allele-specific-filters \
        --mask {params.blacklisted_sites} \
        --mask-name "blacklisted_site") 2> {log}
        """

rule CoverageAtEveryBase:
    input:
        normal_bam = "results/AlignAndMarkDuplicates/{tumor}.bam",
        normal_bai = "results/AlignAndMarkDuplicates/{tumor}.bai",
        shifted_bam = "results/AlignShiftedMTAndMarkDuplicates/{tumor}.bam",
        shifted_bai = "results/AlignShiftedMTAndMarkDuplicates/{tumor}.bai"
    output:
        table = "results/CoverageAtEveryBase/{tumor}_per_base_coverage.tsv",
        non_control_regions = "results/CoverageAtEveryBase/{tumor}_non_control_region.metrics",
        control_region_shifted = "results/CoverageAtEveryBase/{tumor}_control_region_shifted.metrics",
        non_control_region_tsv = "results/CoverageAtEveryBase/{tumor}_non_control_region.tsv",
        control_region_shifted_tsv = "results/CoverageAtEveryBase/{tumor}_control_region_shifted.tsv"
    params:
        control_region_shifted_reference_interval_list = config["control_region_shifted_reference_interval_list"],
        non_control_region_interval_list = config["non_control_region_interval_list"],
        mt_ref = config["mt_ref"],
        mt_shifted_ref = config["mt_shifted_ref"],
        shift_back_chain = config["shift_back_chain"],
        java = config["java"],
        picard_jar = config["picard_jar"],
        CoverageAtEveryBase = config["CoverageAtEveryBase"]
    log:
        "logs/CoverageAtEveryBase/{tumor}.txt"
    shell:
        """(set -e

        {params.java} -jar {params.picard_jar} CollectHsMetrics \
        I={input.normal_bam} \
        R={params.mt_ref} \
        PER_BASE_COVERAGE={output.non_control_region_tsv} \
        O={output.non_control_regions}\
        TI={params.non_control_region_interval_list} \
        BI={params.non_control_region_interval_list} \
        COVMAX=20000 \
        SAMPLE_SIZE=1
        
        {params.java} -jar {params.picard_jar} CollectHsMetrics \
        I={input.shifted_bam} \
        R={params.mt_shifted_ref} \
        PER_BASE_COVERAGE={output.control_region_shifted_tsv} \
        O={output.control_region_shifted} \
        TI={params.control_region_shifted_reference_interval_list} \
        BI={params.control_region_shifted_reference_interval_list} \
        COVMAX=20000 \
        SAMPLE_SIZE=1
        
        Rscript {params.CoverageAtEveryBase} --control_region_shifted {output.control_region_shifted_tsv} --non_control_region {output.non_control_region_tsv} --output {output.table}) 2> {log}
        """

rule SplitMultiAllelicSites:
    input:
        input_vcf = "results/FilterContamination/{tumor}.vcf"
    output:
        split_vcf = "results/SplitMultiAllelicSites/{tumor}.vcf"
    params:
        gatk = config["gatk_path"],
        mt_ref = config["mt_ref"]
    log:
        "logs/SplitMultiAllelicSites/{tumor}.txt"
    shell:
        """(set -e
        
        {params.gatk} LeftAlignAndTrimVariants \
        -R {params.mt_ref} \
        -V {input.input_vcf} \
        -O {output.split_vcf} \
        --split-multi-allelics \
        --dont-trim-alleles \
        --keep-original-ac) 2> {log}
        """
