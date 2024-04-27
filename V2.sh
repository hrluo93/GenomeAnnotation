1###RNA-based
#Hisat2.sh batch maaping

ls /media/perimeter/data/WHXWZB-2023080075A/raw_data/illumina/RNA/*/*/*.R1.fastq.gz | cut -f 1-2 -d '.'| while read loc;do
spname=`ls ${loc}.R1.fastq.gz | cut -f 9 -d '/'`
/media/perimeter/r2/srcs/fastp -i ${loc}.R1.fastq.gz -I ${loc}.R2.fastq.gz -o ${loc}.R1.clean.fastq.gz -O ${loc}.R2.clean.fastq.gz -q 20
hisat2 -x eeu2 -1 ${loc}.R1.clean.fastq.gz -2 ${loc}.R2.clean.fastq.gz --dta -p 16 | samtools sort -@ 16 -O BAM -o ${spname}.sort.bam  
echo ${loc} ${spname} >> suc.samples
done

samtools merge -@ 16 merged.bam *.sort.bam
stringtie merged.bam -p 10 -o merged.gtf

#soft-masked or no mask assembly
gawn/01_scripts/TransDecoder/util/cufflinks_gtf_genome_to_cdna_fasta.pl merged.gtf sp.fasta.masked > sptranscripts.fasta
gawn/01_scripts/TransDecoder/util/cufflinks_gtf_to_alignment_gff3.pl merged.gtf > sptranscripts.gff3
TransDecoder.LongOrfs -t sptranscripts.fasta
TransDecoder.Predict -t sptranscripts.fasta
gawn/01_scripts/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl sptranscripts.fasta.transdecoder.gff3 sptranscripts.gff3 sptranscripts.fasta > sptranscripts.fasta.transdecoder.genome.gff3

2####Ab-initio
##Running braker3 singularity, soft-masked assembly required
cp -r /home/perimeter/miniconda3/envs/braker/config/ /media/perimeter/r2/eeu/
export AUGUSTUS_CONFIG_PATH=/media/perimeter/r2/eeu/config/
singularity exec -B $PWD:$PWD /media/perimeter/r2/srcs/braker3.sif braker.pl --genome=/media/perimeter/r2/eeu/eeu.final2.fasta.masked --species=eeub1 --prot_seq=/media/perimeter/r2/eeu/homofaa2-adcy.fasta --bam=/media/perimeter/r2/eeu/trans/merged.bam --gff3 --threads 32 --AUGUSTUS_CONFIG_PATH=/media/perimeter/r2/eeu/config/

3###homology
##Running in Miniprot, try to use as closest relative species homology proteins as you can 
miniprot -t8 --gff eeu.final2.fasta homofaa2-adcy.fasta > eeu.mini.gff
grep -e "stop_codon" /media/perimeter/r2/eeu/trans/eeu.mini.gff > miniprot.compelte.list
grep -o "MP[0-9]*" miniprot.compelte.list > miniprot.compelte.id
#https://github.com/jorvis/biocode/blob/9f043705e436db4a17dace4a2f70be0da5dfc3b5/gff/filter_gff3_by_id_list.py
python /media/perimeter/r2/eeu/trans/filter_gff3_by_id_list.py -l /media/perimeter/r2/eeu/trans/miniprot.compelte.id -i /media/perimeter/r2/eeu/trans/eeu.mini.gff -o eeu.mini.complete.gff

#####Then convert to EVM gff formation as described in RNA-seq+homology+ab-initio based annotation.sh
#####EVM weights.txt

PROTEIN	miniprot	5
TRANSCRIPT	transdecoder	10
ABINITIO_PREDICTION	AUGUSTUS	2
ABINITIO_PREDICTION	gmst	1
ABINITIO_PREDICTION	GeneMark.hmm3	1

#####PASA pipeline as described in RNA-seq+homology+ab-initio based annotation.sh. Personally suggested that input transcripts from Hiast2+TransDecoder(cufflinks_gtf_to_alignment_gff3.pl) to PASA would be better than trinity assembled (less TE)


