# GenomeAnnotation
Genome annotation in an easy way

For avians, RNA samples are usually not easy to acquire. 

This project aims to create an easy way to annotate genomes.

We provide two pipelines for avian genome annotation, one for with RNA-seq and one for without RNA-seq.


Personally suggested that input transcripts from Hiast2+TransDecoder(cufflinks_gtf_to_alignment_gff3.pl) to PASA

###########update#################

V2.sh: More simplicity, faster, and no frameshift 

RNA-based: Hisat2+TransDdecoder
Ab-initio: braker3
homology: Complete structure from Miniprot 

#################################################

1.RNA-seq+homology+ab-initio based annotation.sh
![annotation1](https://user-images.githubusercontent.com/57522086/180604868-19489cc6-d8c4-4b64-a885-8e0fca3e33b9.png)


2.Homology+Ab-initio.sh based annotation.sh

Ab initio gene prediction using BUSCO single-copy genes as training sets. 

GMAP, Exonerate and GTH were used in homology-based annotation

AUGSUTUS, SNAP(MAKER3), GeneID were used in ab initio annotation


RNA-seq+homology+ab-initio based annotation.sh have been used in the great bustard genome




