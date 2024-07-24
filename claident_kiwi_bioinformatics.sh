# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Claident script
# Kiwi gut microbiomes through time
# Priscilla San Juan
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


# BACTERIA 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# Demultiplex bacteria
clsplitseq --runname=kiwiMic --truncateN=enable 
--index1file=index1.txt ---index2file=index2.txt 
--primerfile=forward.primer.16S.txt --reverseprimerfile=reverse.primer.16S.txt --minqualtag=20 --numthreads=8 Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_I1_001.fastq.gz Undetermined_S0_L001_I2_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz Demultiplexed_16S
 
# Removing "undetermined" files.
cd ../Demultiplexed_16S
rm *undetermined*fastq.gz

# Merging of paired-end reads with PEAR
module load PEAR
for f in `ls *.forward.fastq.gz | grep -P -o '^[^\.]+'`; do temp=(`echo $f|grep -P -o '[^\__]+___L001'`); pear -p 0.001 -u 0 -j 32 -f $f.forward.fastq.gz -r $f.reverse.fastq.gz -o $f; done

# Filtering of merged reads
for f in `ls *.assembled.fastq | grep -o -P '^[^\.]+'`
do clfilterseq --minqual=30 --minlen=150 --maxplowqual=0.1 --numthreads=8 $f.assembled.fastq $f.filtered.fastq
done

# Detection/removal of noisy reads
module load VSEARCH
for f in `ls *filtered.fastq | grep -o -P '^[^\.]+'`
do clcleanseqv --derepmode=FULLLENGTH --primarymaxnmismatch=1 --secondarymaxnmismatch=3 --pnoisycluster=0.5 --numthreads=8 $f.filtered.fastq clcleanseqv.$f
done

# Clustering with VSEARCH
module load VSEARCH
clclassseqv --minident=0.97 --numthreads=8 
clcleanseqv.*/primarycluster.denoised.fasta.gz clclassseqv.16S


# Mapping of raw reads that were removed in the previous noisy-read detection process but are mapped to the clustered OTUs with a given cutoff similarity
clrecoverseqv --minident=0.97 --centroid=clclassseqv.16S/clustered.fasta --numthreads=8 clcleanseqv.*/primarycluster.fasta.gz clrecoverseqv.16S

# Change directory
cd ./clrecoverseqv.16S

# De-novo chimera removal with UCHIME. 
vsearch  --uchime_denovo clustered.fasta --nonchimeras otus_nonchimeras.fasta --chimeras otus_chimera.fasta --uchimeout otus_chimeras.uc --threads 4 --sizein --sizeout --fasta_width 0 --log chimera.log 

# Reference-based chimera removal with UCHIME
vsearch --referencedb=silva132LSUref.fasta chimera.denovo/nonchimeras.fasta chimera.ref

# Generating a sample x OTU summary table. Only the OTUs that appear in the above chimera-excluded fasta file will appear in the summary matrix.
clsumclass clustered.otu.gz summary.16S.txt

# Removing OTUs with low abundance from the summary matrix
clfiltersum --minntotalseqotu=10 summary.16S.txt summary.16S.10more.txt

# Removing rare OTUs from the fasta file.
clfilterseq --otufile=clustered.otu.gz --minnseq=10 chimera.ref/nonchimeras.fasta nonchimera.10more.fasta

# Tax assignment
module load BLAST/2.10.0-GCC-9.2.0
clidentseq ---blastdb=prokaryota_16S_genus --numthreads=32 clustered.fasta bacteria_kiwi_2020_genus

classigntax --taxdb=prokaryota_16S_species bacteria_kiwi_2020_genus bacteria_kiwi_2020_tax_table.txt

