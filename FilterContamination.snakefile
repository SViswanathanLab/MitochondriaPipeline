import os

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
        expand("results/FilterContamination/{tumor}.vcf", tumor=config["pairings"]),
        expand("results/FilterContamination/{tumor}.filtered.vcf", tumor=config["pairings"]),
        expand("results/SplitMultiAllelicSites/{tumor}.vcf", tumor=config["pairings"]),
        expand("results/vcf2maf/{tumor}.maf", tumor=config["pairings"])

     
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
        bamout = "results/FilterContamination/{tumor}_bamout.bam"
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

rule vcf2maf:
    input:
        tumor_vcf = "results/SplitMultiAllelicSites/{tumor}.vcf"
    output:
        maf = "results/vcf2maf/{tumor}.maf"
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
        --tumor-id `grep "tumor_sample" {input.tumor_vcf} | awk -F "=" '{{print $2}}'` \
        --input-vcf {input.tumor_vcf} \
        --output-maf {output.maf}) 2> {log}
        """
