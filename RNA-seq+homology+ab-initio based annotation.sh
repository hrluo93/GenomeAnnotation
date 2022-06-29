#Homology
#1 Genomethreader
/gth-1.7.1/bin/gth -genomic sp.fasta.masked -intermediate -protein homo.faa -gff3out -species chicken -o sp_gth.gff3
#2 Exonerate 2.2.0 (We found some memery use bugs in 2.4.0)
exonerate -t sp.fasta.masked -q homo.faa --querytype protein --targettype dna --model protein2genome --bestn 1 --showtargetgff yes --showalignment no -M 256000

#RNA-seq Based 
hisat2-build sp.fasta.masked spmasked 
hisat2 -x spmasked -1 d1.fq.gz -2 d2.fq.gz --dta -p 16 | samtools sort -@ 16 -O BAM -o datuijirou.sort.bam
hisat2 -x spmasked -1 x1.fq.gz -2 x2.fq.gz --dta -p 16 | samtools sort -@ 16 -O BAM -o xin.sort.bam
hisat2 -x spmasked -1 g1.fq.gz -2 g2.fq.gz --dta -p 16 | samtools sort -@ 16 -O BAM -o gan.sort.bam
hisat2 -x spmasked -1 f1.fq.gz -2 f2.fq.gz --dta -p 16 | samtools sort -@ 16 -O BAM -o fei.sort.bam
hisat2 -x spmasked -1 n1.fq.gz -2 n2.fq.gz --dta -p 16 | samtools sort -@ 16 -O BAM -o danao.sort.bam
samtools merge -@ 16 merged.bam *.sort.bam
stringtie /home/tilapia/xuluohao/dabao/asm/anno/stringtie/merged.bam -p 10 -o stringtieout/merged.gtf

gawn/01_scripts/TransDecoder/util/cufflinks_gtf_genome_to_cdna_fasta.pl merged.gtf sp.fasta.masked > sptranscripts.fasta
gawn/01_scripts/TransDecoder/util/cufflinks_gtf_to_alignment_gff3.pl merged.gtf > sptranscripts.gff3
TransDecoder.LongOrfs -t sptranscripts.fasta
TransDecoder.Predict -t sptranscripts.fasta
gawn/01_scripts/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl sptranscripts.fasta.transdecoder.gff3 sptranscripts.gff3 sptranscripts.fasta > sptranscripts.fasta.transdecoder.genome.gff3
