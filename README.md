# GenomeAnnotation


This project aims to create an easy way to annotate genomes.


Personally suggested that input transcripts from Hiast2+TransDecoder(cufflinks_gtf_to_alignment_gff3.pl) to PASA


###########update#################

V2.sh: More simplicity, faster, and no frameshift 

RNA-based: Hisat2+TransDdecoder

Ab-initio: braker3

Homology: Complete structure from Miniprot


########### Post PASA ############

This [
](https://github.com/hrluo93/python4bio/blob/main/false-gene-model.py) script can check if annotation contained false gene model.

###Soft-masked genome would result in TE contained in annotation. We used OrthoFinder to filter annotation results to retain orthologous genes and remove non-orthologous with 1 or 2 exons.


##Target species (Gene ID like NNYC0000010.1 )in Orthogroups.GeneCount.tsv $3 with reference species in $2 and $4

orthofinder -f orthof -og -M msa -t 12 -S blast_gz

cd orthof/*/Orthogroups/

cat Orthogroups.GeneCount.tsv | awk '{if ($2 > 0 || $4 >0) print}' | awk '{if ($3 > 0) print}' > nny.allortho.count.tsv

awk 'FNR==NR {a[$1]=$0;next} $1 in a {print a[$1],$0}'  nny.allortho.count.tsv Orthogroups.tsv > nny.merge.tsv

grep -o ![image](https://github.com/user-attachments/assets/2c3c5925-c08a-4e5f-ad6a-85f51b9cc063) nny.merge.tsv | cut -f1 -d "." > nny.orthogene.list

#nny.orthogene.list contained all orthologous genes that should kept. Non-orthologous with 1 or 2 exons can be found via TBTools GXF STAT or any other method you prefer.



######################################################

RNA-seq+homology+ab-initio based annotation.sh used in the great bustard genome.
![annotation1](https://user-images.githubusercontent.com/57522086/180604868-19489cc6-d8c4-4b64-a885-8e0fca3e33b9.png)




