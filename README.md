# GenomeAnnotation
Genome annotation in an easy way


This project aims to create an easy way to annotate genomes.

We provide two pipelines for avian genome annotation, one for with RNA-seq and one for without RNA-seq.


Personally suggested that input transcripts from Hiast2+TransDecoder(cufflinks_gtf_to_alignment_gff3.pl) to PASA

###########update#################

V2.sh: More simplicity, faster, and no frameshift 

RNA-based: Hisat2+TransDdecoder

Ab-initio: braker3

Homology: Complete structure from Miniprot 

########### Post PASA ############
###Soft-masked genome would result in TE contained in annotation. We used OrthoFinder to filter annotation results to the retention of orthologous genes and remove non-orthologous with 1 or 2 exons.
##Target species (Gene ID ) in Orthogroups.GeneCount.tsv $3 with reference species in $2 and $4

orthofinder -f orthof -og -M msa -t 12 -S blast_gz

cd orthof/*/Orthogroups/

cat Orthogroups.GeneCount.tsv | awk '{if ($2 > 0 || $4 >0) print}' | awk '{if ($3 > 0) print}' > nny.allortho.count.tsv

awk 'FNR==NR {a[$1]=$0;next} $1 in a {print a[$1],$0}'  nny.allortho.count.tsv Orthogroups.tsv > nny.merge.tsv

grep -o "NNYC[0-9]*\.[0-9]*" nny.merge.tsv | cut -f1 -d "." > nny.orthogene.list

#nny.orthogene.list contained all orthologous genes that should kept. Non-orthologous with 1 or 2 exons can be found via TBTools GXF STAT or any other method you prefer.



#################################################

1.RNA-seq+homology+ab-initio based annotation.sh
![annotation1](https://user-images.githubusercontent.com/57522086/180604868-19489cc6-d8c4-4b64-a885-8e0fca3e33b9.png)


2.Homology+Ab-initio.sh based annotation.sh

Ab initio gene prediction using BUSCO single-copy genes as training sets. 

GMAP, Exonerate and GTH were used in homology-based annotation

AUGSUTUS, SNAP(MAKER3), GeneID were used in ab initio annotation


RNA-seq+homology+ab-initio based annotation.sh have been used in the great bustard genome



