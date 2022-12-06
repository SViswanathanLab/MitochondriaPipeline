import os

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
        #expand("results/GetContamination/{tumor}/output", tumor=config["pairings"]),
        #expand("results/GetContamination/{tumor}/output-noquotes", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_headers.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_output_data.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_contamination.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_major_hg.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_minor_hg.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_mean_het_major.txt", tumor=config["pairings"]),
        expand("results/GetContamination/{tumor}/{tumor}_mean_het_minor.txt", tumor=config["pairings"])

rule GetContamination:
    input:
        input_vcf = "results/SplitMultiAllelicsAndRemoveNonPassSites/{tumor}/{tumor}_splitAndPassOnly.vcf"
    output:
        outputs = "results/GetContamination/{tumor}/output",
        output_noquotes = "results/GetContamination/{tumor}/output-noquotes",
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
        #mv --force output {params.MitochondriaPipeline_path}{output.outputs}
        sed 's/\\\"//g' output > {params.MitochondriaPipeline_path}{output.output_noquotes}
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
