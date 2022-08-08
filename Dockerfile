FROM java:openjdk-8-jre

WORKDIR /home/mi724/Tools/MitochondriaPipeline/haplocheckCLI

RUN jar cvfe haplocheckCLI.jar haplocheck_contam *
