#Homology
#1 Genomethreader
#genometools 
gt seqtransform -addstopaminos yes homoraw.faa > homo.faa
/gth-1.7.1/bin/gth -genomic sp.fasta.masked -intermediate -protein homo.faa -gff3out -species chicken -o sp_gth.gff3
EVidenceModeler-1.1.1/EvmUtils/misc/genomeThreader_to_evm_gff3.pl sp_gth.gff3 > sp_gth.prot.gff3
#2 Exonerate 2.2.0 (We found some memery usage bugs in 2.4.0)
exonerate -t sp.fasta.masked -q homo.faa --querytype protein --targettype dna --model protein2genome --bestn 1 --showtargetgff yes --showalignment no -M 256000 -o sp_exonerate.gff
EVidenceModeler-1.1.1/EvmUtils/misc/exonerate_gff_to_alignment_gff3.pl sp_exonerate.gff > sp_exonerate.prot.gff
#Or you can use maker output 
cat sp.masked.fasta.maker.output/sp.masked.fasta_datastore/*/*/scaffold_*/theVoid.scaffold_*/evidence_*.gff > evi-sca-all.gff
#maker match to evm input gff 
EVidenceModeler-1.1.1/EvmUtils/misc/maker_match_gff_to_gene_gff3.pl allevi.gff > allmakergene.gff
awk '{if ($2=="est2genome") print }' allmakergene.gff > sp.est2genome.est.gff
awk '{if ($2=="protein2genome") print }' allmakergene.gff > sp.prot2genome.prot.gff
awk '{if ($2=="blastx") print }' allmakergene.gff > sp.blastx2genome.blastx.gff

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

#ab-initio
#1 AUGUSTUS can directly trained in BuscoV5
busco -i sp.fasta.masked -o sp -l aves_odb10 -m genome -c 12 -f --long --augustus --augustus_parameters='--progress=true'
# The training result can be found in run_aves_odb10/augustus_output/retraining_parameters/BUSCO_sp
# CP BUSCO_sp to AUGUSTUS speciesdir
augustus --species=BUSCO_sp --gff3=on sp.masked.fa > sp_aug.gff
#Converted to evm.input
/home/tilapia/srcs/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl sp_aug.gff > sp_aug.ab.EVM.gff

#2 GeneID 
#Train From sptranscripts.fasta.transdecoder.genome.gff3
#Converted to Gff2
#Like this
#scaffold222	transdecoder	CDS	294413	294689	.	+	.	g1
#scaffold222	transdecoder	CDS	295082	295257	.	+	.	g1
#scaffold222	transdecoder	CDS	295572	295889	.	+	.	g1
#scaffold215	transdecoder	CDS	678615	678830	.	+	.	g2
#scaffold215	transdecoder	CDS	682004	682167	.	+	.	g2
#scaffold215	transdecoder	CDS	682257	682392	.	+	.	g2
#scaffold215	transdecoder	CDS	682773	682972	.	+	.	g2
#scaffold215	transdecoder	CDS	683046	683295	.	+	.	g2
#Make sure you have permission run Docker, otherwise you can choose Fgenesh 
service docker start      
docker build -t geneidtrainerdocker .
docker run -u $(id -u):$(id -g) -v /your/gff2&fa/path:/data -w /data geneidtrainerdocker -species Sp -gff sp.gff2 -fastas sp.masked.fasta -results ./output/ -reduced no
geneid -3 -P SP.geneid.param sp.masked.fasta > sp_genid.gff
EVidenceModeler-1.1.1/EvmUtils/misc/GeneID_to_gff3.pl sp_genid.gff > sp_genid.ab.EVM.gff

#3 SNAP
#In Maker3
Trinity --genome_guided_bam merged.bam --genome_guided_max_intron 10000 --max_memory 50G --CPU 20
#Maker
maker -CTL
vim maker_opts.ctl
#genome=sp.masked.fasta
#est=Trinity-GG.fasta
#protein=homo.faa
#-----Repeat Masking (leave values blank to skip repeat masking)
#model_org=
#rmlib= 
#repeat_protein= 
#rm_gff= 
#prok_rm=0
#softmask=0 
#est2genome=1 
#protein2genome=1
mpiexec -n 40 maker/bin/maker > runsnapr1.log
maker/bin/gff3_merge -d sp.masked.fasta.maker.output/sp.masked.fasta_master_datastore_index.log
maker/bin/maker2zff -c 0.8 -e 0.8 -o 0.8 -x 0.2 sp.masked.fasta.all.gff
fathom genome.ann genome.dna -gene-stats > stats.log
fathom genome.ann genome.dna -validate > validate.log
fathom genome.ann genome.dna -categorize 1000 > categorize.log
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log
mkdir paramsr1
cd paramsr1
forge ../export.ann ../export.dna > forge.log
cd ..
hmm-assembler.pl sp_snap1 /paramsr1 > sp_snap1.hmm
#2nd round training
cp sp.masked.fasta.all.gff sp.round1.masked.fasta.all.gff
vim maker_opts.ctl
#est2genome=0
#protein2genome=0
#snaphmm=sp_snap1.hmm
mpiexec -n 40 maker/bin/maker > runsnapr2.log
maker/bin/gff3_merge -d sp.masked.fasta.maker.output/sp.masked.fasta_master_datastore_index.log
maker/bin/maker2zff -c 0.8 -e 0.8 -o 0.8 -x 0.2 sp.masked.fasta.all.gff
fathom genome.ann genome.dna -gene-stats > stats.log
fathom genome.ann genome.dna -validate > validate.log
fathom genome.ann genome.dna -categorize 1000 > categorize.log
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log
mkdir paramsr2
cd paramsr2
forge ../export.ann ../export.dna > forge.log
cd ..
hmm-assembler.pl sp_snap2 /paramsr2 > sp_snap2.hmm
#3rd training
cp sp.masked.fasta.all.gff sp.round2.masked.fasta.all.gff
vim maker_opts.ctl
#snaphmm=sp_snap2.hmm
mpiexec -n 40 maker/bin/maker > runsnapr3.log
maker/bin/gff3_merge -d sp.masked.fasta.maker.output/sp.masked.fasta_master_datastore_index.log
maker/bin/maker2zff -c 0.8 -e 0.8 -o 0.8 -x 0.2 sp.masked.fasta.all.gff
fathom genome.ann genome.dna -gene-stats > stats.log
fathom genome.ann genome.dna -validate > validate.log
fathom genome.ann genome.dna -categorize 1000 > categorize.log
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log
mkdir paramsr3
cd paramsr3
forge ../export.ann ../export.dna > forge.log
cd ..
hmm-assembler.pl sp_snap3 /paramsr3 > sp_snap3.hmm
#Final of SNAP
vim maker_opts.ctl
#snaphmm=sp_snap3.hmm
mpiexec -n 40 maker/bin/maker > runsnapf.log
maker/bin/gff3_merge -d sp.masked.fasta.maker.output/sp.masked.fasta_master_datastore_index.log
awk '{if ($2=="maker") print }' sp.masked.fasta.all.gff > sp_snap.ab.gff

#EVM
cat *.ab.gff > ab.evm.gff
cat *.prot.gff > prot.evm.gff
vim weights.txt
PROTEIN	protein2genome	5
PROTEIN	genomeThreader	5
TRANSCRIPT	transdecoder	10
ABINITIO_PREDICTION	Augustus	3
ABINITIO_PREDICTION	maker	2
ABINITIO_PREDICTION	geneid_v1.4	1


EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome sp.masked.fasta --protein_alignments prot.evm.gff --gene_predictions prot.evm.gff --transcript_alignments sptranscripts.fasta.transdecoder.genome.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome sp.masked.fasta --protein_alignments prot.evm.gff --gene_predictions prot.evm.gff --transcript_alignments sptranscripts.fasta.transdecoder.genome.gff3 --weights full/path/weights.txt --output_file_name evm.out  --partitions partitions_list.out >  commands.list
EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list
EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output_file_name evm.out --genome sp.fasta.masked
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3

#PASA update annotation and add UTR
cp Trinity-GG.fasta all_transcripts.fasta
makeblastdb -in UniVec -dbtype nucl -parse_seqids -title UniVec -out UniVec
pasa-2.5.2/bin/tools/seqclean all_transcripts.fasta -v /home/PASA/UniVec
Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g sp.masked.fasta -t all_transcripts.fasta.clean -T -u all_transcripts.fasta --ALIGNERS blat,gmap --CPU 20 --TRANSDECODER
Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g genome_sample.fasta -P EVM.all.gff3
#SQLite path in alignAssembly.config and annotCompare.config
Launch_PASA_pipeline.pl -c annotCompare.config -A -g sp.masked.fasta -t all_transcripts.fasta.clean --CPU 20
#2nd 
Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g genome_sample.fasta -P pasa.round1.gff3
Launch_PASA_pipeline.pl -c annotCompare.config -A -g sp.masked.fasta -t all_transcripts.fasta.clean --CPU 20

#Done
#Bye!
