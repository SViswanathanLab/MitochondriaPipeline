### Based on WDL script in: https://github.com/OpenGenomics/vcf2maf-tools/blob/master/vcf2maf.wdl


import os

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
         expand("results/SplitMultiAllelicSites/{tumor}.vcf", tumor=config["pairings"]),
         expand("results/vcf2maf/{tumor}.maf", tumor=config["pairings"])
                  
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
        --input-vcf {input.tumor_vcf} \
        --output-maf {output.maf}) 2> {log}
        """
