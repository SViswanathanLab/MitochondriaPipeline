FROM java:openjdk-8-jre

WORKDIR /home/mi724/Tools/MitochondriaPipeline

RUN wget https://github.com/leklab/haplocheckCLI/raw/master/haplocheckCLI.jar
RUN jar cvfe haplocheckCLI/haplocheckCLI.jar haplocheck_contam *
