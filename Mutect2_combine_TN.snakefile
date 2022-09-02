import os

configfile: "config/config.yaml"
configfile: "config/sample_combine_TN.yaml"

rule all:
    input:
        expand("results2/CallMt/{tumor}.vcf.gz", tumor=config["pairings"]),
        expand("results2/CallMt/{tumor}.vcf.gz.tbi", tumor=config["pairings"]),
        expand("results2/CallMt/{tumor}.vcf.gz.stats", tumor=config["pairings"]),
        expand("results2/CallMt/{tumor}_bamout.bam", tumor=config["pairings"]),
        expand("results2/CallShiftedMt/{tumor}.vcf.gz", tumor=config["pairings"]),
        expand("results2/CallShiftedMt/{tumor}.vcf.gz.tbi", tumor=config["pairings"]),
        expand("results2/CallShiftedMt/{tumor}.vcf.gz.stats", tumor=config["pairings"]),
        expand("results2/CallShiftedMt/{tumor}_bamout.bam", tumor=config["pairings"]),
        expand("results2/LiftoverAndCombineVcfs/{tumor}.shifted_back.vcf", tumor=config["pairings"]),
        expand("results2/LiftoverAndCombineVcfs/{tumor}.rejected.vcf", tumor=config["pairings"]),
        expand("results2/LiftoverAndCombineVcfs/{tumor}.merged.vcf", tumor=config["pairings"]),
        expand("results2/MergeStats/{tumor}.raw.combined.stats", tumor=config["pairings"]),
        expand("results2/InitialFilter/{tumor}.filtered.vcf", tumor=config["pairings"]),
        expand("results2/InitialFilter/{tumor}.vcf", tumor=config["pairings"]),
        expand("results2/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_splitAndPassOnly.vcf", tumor=config["pairings"]),
        expand("results2/GetContamination/{tumor}/{tumor}_headers.txt", tumor=config["pairings"]),
        expand("results2/GetContamination/{tumor}/{tumor}_output_data.txt", tumor=config["pairings"]),
        expand("results2/GetContamination/{tumor}/{tumor}_contamination.txt", tumor=config["pairings"]),
        expand("results2/GetContamination/{tumor}/{tumor}_major_hg.txt", tumor=config["pairings"]),
        expand("results2/GetContamination/{tumor}/{tumor}_minor_hg.txt", tumor=config["pairings"]),
        expand("results2/GetContamination/{tumor}/{tumor}_mean_het_major.txt", tumor=config["pairings"]),
        expand("results2/GetContamination/{tumor}/{tumor}_mean_het_minor.txt", tumor=config["pairings"]),
        expand("results2/FilterContamination/{tumor}.vcf", tumor=config["pairings"]),
        expand("results2/FilterContamination/{tumor}.filtered.vcf", tumor=config["pairings"]),
        expand("results2/SplitMultiAllelicSites/{tumor}.vcf", tumor=config["pairings"]),
        expand("results2/vcf2maf/{tumor}.maf", tumor=config["pairings"])

rule CallMt:
    input:
        bam_T_filepath = lambda wildcards: config["samples"][wildcards.tumor],
        bam_N_filepath = lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
    output:
        vcf = protected("results2/CallMt/{tumor}.vcf.gz"),
        tbi = protected("results2/CallMt/{tumor}.vcf.gz.tbi"),
        stats = protected("results2/CallMt/{tumor}.vcf.gz.stats"),
        bamout = protected("results2/CallMt/{tumor}_bamout.bam")
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
        -I {input.bam_T_filepath} \
        -I {input.bam_N_filepath} \
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
        bam_T_filepath = lambda wildcards: config["samples2"][wildcards.tumor],
        bam_N_filepath = lambda wildcards: config["samples2"][config["pairings"][wildcards.tumor]]
    output:
        vcf = protected("results2/CallShiftedMt/{tumor}.vcf.gz"),
        tbi = protected("results2/CallShiftedMt/{tumor}.vcf.gz.tbi"),
        stats = protected("results2/CallShiftedMt/{tumor}.vcf.gz.stats"),
        bamout = protected("results2/CallShiftedMt/{tumor}_bamout.bam")
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
        -I {input.bam_T_filepath} \
        -I {input.bam_N_filepath} \
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
        vcf = "results2/CallMt/{tumor}.vcf.gz",
        shifted_vcf = "results2/CallShiftedMt/{tumor}.vcf.gz"
    output:
        shifted_back = "results2/LiftoverAndCombineVcfs/{tumor}.shifted_back.vcf",
        rejected = "results2/LiftoverAndCombineVcfs/{tumor}.rejected.vcf",
        merged = "results2/LiftoverAndCombineVcfs/{tumor}.merged.vcf"
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
        shifted_stats = "results2/CallShiftedMt/{tumor}.vcf.gz.stats",
        non_shifted_stats = "results2/CallMt/{tumor}.vcf.gz.stats"
    output:
        raw_combined_stats = "results2/MergeStats/{tumor}.raw.combined.stats"
    params:
        gatk = config["gatk_path"]

    log:
        "logs/MergeStats/{tumor}.txt"
    shell:
        """(set -e

        {params.gatk} MergeMutectStats --stats {input.shifted_stats} --stats {input.non_shifted_stats} -O {output.raw_combined_stats}) 2> {log}"""

rule InitialFilter:
    input:
        raw_vcf = "results2/LiftoverAndCombineVcfs/{tumor}.merged.vcf",
        raw_vcf_stats = "results2/MergeStats/{tumor}.raw.combined.stats"
    output:
        filtered_vcf = "results2/InitialFilter/{tumor}.filtered.vcf",
        output_vcf = "results2/InitialFilter/{tumor}.vcf",
        bamout = "results2/InitialFilter/{tumor}_bamout.bam"
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
        filtered_vcf = "results2/InitialFilter/{tumor}.vcf"
    output:
        split_vcf = temp("results2/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_split.vcf"),
        split_vcf_idx = temp("results2/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_split.vcf.idx"),
        splitAndPassOnly_vcf = "results2/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_splitAndPassOnly.vcf"
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
        input_vcf = "results2/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_splitAndPassOnly.vcf"
    output:
        outputs = "results2/GetContamination/{tumor}/output",
        output_noquotes = "results2/GetContamination/{tumor}/output-noquotes",
        headers = "results2/GetContamination/{tumor}/{tumor}_headers.txt",
        output_data = "results2/GetContamination/{tumor}/{tumor}_output_data.txt",
        contamination = "results2/GetContamination/{tumor}/{tumor}_contamination.txt",
        major_hg = "results2/GetContamination/{tumor}/{tumor}_major_hg.txt",
        minor_hg = "results2/GetContamination/{tumor}/{tumor}_minor_hg.txt",
        mean_het_major = "results2/GetContamination/{tumor}/{tumor}_mean_het_major.txt",
        mean_het_minor = "results2/GetContamination/{tumor}/{tumor}_mean_het_minor.txt"
    params:
        java = config["java"],
        picard_jar = config["picard_jar"],
        haplocheckCLI_path = config["haplocheckCLI_path"],
        haplocheckCLI_newpath = "results/GetContamination/{tumor}/",
        MitochondriaPipeline_path = config["MitochondriaPipeline_path"]
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

        cd {params.haplocheckCLI_newpath}
        {params.java} -jar {params.haplocheckCLI_path} "$(dirname "{params.MitochondriaPipeline_path}{input.input_vcf}")"
        #mv output {output.outputs}
        sed 's/\\\"//g' {params.MitochondriaPipeline_path}{output.outputs} > {params.MitochondriaPipeline_path}{output.output_noquotes}
        cd {params.MitochondriaPipeline_path}
        grep "SampleID" {output.output_noquotes} > {output.headers}
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

        grep -v "SampleID" {output.output_noquotes}  > {output.output_data}
        awk -F \"\\t\" '{{print $2}}' {output.output_data} > {output.contamination}
        awk -F \"\\t\" '{{print $6}}' {output.output_data} > {output.major_hg}
        awk -F \"\\t\" '{{print $8}}' {output.output_data} > {output.minor_hg}
        awk -F \"\\t\" '{{print $14}}' {output.output_data} > {output.mean_het_major}
        awk -F \"\\t\" '{{print $15}}' {output.output_data} > {output.mean_het_minor}) 2> {log}"""

rule FilterContamination:
    input:
        hasContamination = "results2/GetContamination/{tumor}/{tumor}_contamination.txt",
        contamination_major = "results2/GetContamination/{tumor}/{tumor}_mean_het_major.txt",
        contamination_minor = "results2/GetContamination/{tumor}/{tumor}_mean_het_minor.txt",
        filtered_vcf = "results2/InitialFilter/{tumor}.vcf",
        raw_vcf_stats = "results2/MergeStats/{tumor}.raw.combined.stats"
    output:
        output_vcf = "results2/FilterContamination/{tumor}.vcf",
        filtered_vcf = "results2/FilterContamination/{tumor}.filtered.vcf",
        bamout = "results2/FilterContamination/{tumor}_bamout.bam"
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
          max_contamination=$(echo "1 - $MAJOR" | bc)
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

rule SplitMultiAllelicSites:
    input:
        input_vcf = "results2/FilterContamination/{tumor}.vcf"
    output:
        split_vcf = "results2/SplitMultiAllelicSites/{tumor}.vcf"
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

rule vcf2maf:
    input:
        tumor_vcf = "results2/SplitMultiAllelicSites/{tumor}.vcf"
    output:
        maf = "results2/vcf2maf/{tumor}.maf"
    params:
        perl = config["perl"],
        vcf2maf = config["vcf2maf"],
        vep = config["vep"],
        vep_cache = config["vep_cache"],
        samtools = config["samtools"],
        reference_genome = config["mt_ref"]
    log:
        "logs/vcf2maf/{tumor}.txt"
    shell:
        """(set -e

        export PATH=$PATH:{params.samtools}
        export PERL5LIB={params.vep}
        export HTSLIB_DIR={params.vep}/htslib
        export LD_LIBRARY_PATH=/mnt/storage/apps/htslib/1.9

        {params.perl} {params.vcf2maf} \
        --ref-fasta  {params.reference_genome} \
        --vep-path {params.vep} \
        --vep-data {params.vep_cache} \
        --filter-vcf 0 \
        --ncbi-build GRCh38 \
        --input-vcf {input.tumor_vcf} \
        --output-maf {output.maf}) 2> {log}
        """

