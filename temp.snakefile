    
rule InitialFilter:

rule SplitMultiAllelicsAndRemoveNonPassSites:
    input:
    output:
        split_vcf = 
        splitAndPassOnly_vcf
        
    log:
        "logs/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}.txt"
    params:
        gatk = config["gatk_path"]
    shell:
        """(set -e

        {params.gatk} LeftAlignAndTrimVariants \
        -R {params.ref_fasta} \
        -V ~{filtered_vcf} \
        -O {output.split_vcf} \
        --split-multi-allelics \
        --dont-trim-alleles \
        --keep-original-ac

        {params.gatk} SelectVariants \
        -V {output.split_vcf} \
        -O {output.splitAndPassOnly_vcf} \
        --exclude-filtered) 2> {log}"""

rule GetContamination:
    params:
        java = config["java"],
        picard_jar = config["picard_jar"]
    log:
        "logs/GetContamination/{tumor}.txt"
    shell:
        """(set -e
        PARENT_DIR="$(dirname "~{input_vcf}")"
        {params.java} -jar /usr/mtdnaserver/haplocheckCLI.jar "${PARENT_DIR}"

        sed 's/\"//g' output > output-noquotes

        grep "SampleID" output-noquotes > headers
        FORMAT_ERROR="Bad contamination file format"
        if [ `awk '{print $2}' headers` != "Contamination" ]; then
          echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{print $6}' headers` != "HgMajor" ]; then
          echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{print $8}' headers` != "HgMinor" ]; then
          echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{print $14}' headers` != "MeanHetLevelMajor" ]; then
          echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{print $15}' headers` != "MeanHetLevelMinor" ]; then
          echo $FORMAT_ERROR; exit 1
        fi

        grep -v "SampleID" output-noquotes > output-data
        awk -F "\t" '{print $2}' output-data > contamination.txt
        awk -F "\t" '{print $6}' output-data > major_hg.txt
        awk -F "\t" '{print $8}' output-data > minor_hg.txt
        awk -F "\t" '{print $14}' output-data > mean_het_major.txt
        awk -F "\t" '{print $15}' output-data > mean_het_minor.txt) 2> {log}"""

rule FilterContamination:
    
rule CoverageAtEveryBase:

rule SplitMultiAllelicSites:
