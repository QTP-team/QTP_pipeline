# QTP pipeline
## 1. Assemble
### 1.1 bowtie2 (bowtie2 v2.5.1)
```
seqkit seq -m 1500 {sample_contig} > {sample_contig_filter}
bowtie2-build {sample_contig_filter} {sample_contig_index} 1> {sample_bowtie2_index.log}
bowtie2 -p 8 -x {sample_contig_index} -1 {sample.rmhost.1.fq.gz} -2 {sample.rmhost.2.fq.gz} > {smaple.sam}
samtools sort -@ 8 {sample.sam} -O BAM > {sample.bam}
samtools index {sample.bam} {smaple.bam.bai}
```

### 1.2 concoct (concoct v1.1.0) 
```
cut_up_fasta.py {sample_contig_filter} -c 10000 -o 0 --merge_last -b {sample_contigs_10K.bed} > {smaple_contigs_10K.fa}
concoct_coverage_table.py {sample_contigs_10K.bed} {sample.bam} > {sample_coverage_table.tsv}
concoct -t 8 --composition_file {sample_contigs_10K.fa} --coverage_file {sample_coverage_table.tsv} -b {sample_concoct_output}
merge_cutup_clustering.py {sample_clustering_gt1000.csv} > {sample_clustering_merged.csv}
extract_fasta_bins.py {sample_contig_filter} {sample_clustering_merged.csv} --output_path {fasta_bins}
```

### 1.3 maxbin2 (maxbin2 v2.2.7)
```
pileup.sh in={sample.sam} out={sample_cov.txt}
rm {sample.sam}
awk '{print $1"\t"$2}' {sample_cov.txt} | grep -v '^#' > {sample_abundance.txt}
run_MaxBin.pl -thread 8 -contig {sample_contig_filter} -out {sample_maxbin2_out} -abund {sample_abundance.txt}
```

### 1.4 metabat2 (metabat2 v2.15)
```
jgi_summarize_bam_contig_depths {sample.bam} --outputDepth {sample.depth}
metabat2 -i {sample_contig_filter} -a {sample.depth} -o {sample.binning} -m 1500 -t 8 > {sample_binning_metabat2.log}
```

### 1.5 MAGScoT (MAGScoT v1.0.0)
```
prodigal {sample_contig_filter} -p meta -a {sample.prodigal.faa} -d {sample.prodigal.ffn} -o {sample_tmpfile}
hmmsearch -o {sample.hmm.tigr.out} --tblout {sample.hmm.tigr.hit.out} --noali --notextw --cut_nc --cpu 8 gtdbtk_rel214_tigrfam.hmm {sample.prodigal.faa}
hmmsearch -o {sample.hmm.pfam.out} --tblout {sample.hmm.pfam.hit.out} --noali --notextw --cut_nc --cpu 8 gtdbtk_rel214_Pfam-A.hmm {sample.prodigal.faa}
cat {sample.hmm.tigr.hit.out} | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > {sample.tigr}
cat {sample.hmm.pfam.hit.out} | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > {sample.pfam}
cat {sample.pfam} {sample.tigr} > {sample.hmm}
sh Fasta_to_Contig2Bin.sh -i {sample_maxbin2_bins} -e fasta > {sample.maxbin2.contigs2bin.tsv}
sh Fasta_to_Contig2Bin.sh -i {sample_metabat2_bins} -e fa > {sample.metabat2.contigs2bin.tsv}
sh Fasta_to_Contig2Bin.sh -i {sample_concoct_bins} -e fa > {sample.concoct.contigs2bin.tsv}
awk '{print $2"\t"$1"\tmaxbin2"}'  {sample.maxbin2.contigs2bin.tsv} >> {sample.contigs_to_bin.tsv}
awk '{print $2"\t"$1"\tmetabat2"}'  {sample.metabat2.contigs2bin.tsv} >> {sample.contigs_to_bin.tsv}
awk '{print $2"\t"$1"\tconcoct"}'  {sample.concoct.contigs2bin.tsv} >> {sample.contigs_to_bin.tsv}
Rscript MAGScoT.R -i {sample.contigs_to_bin.tsv} --hmm {sample.hmm} -o {sample.MAGScoT}
```

## 2. Genome de-redundancy (https://github.com/QTP-team/dRep)


## 3. SGBs profile (https://github.com/QTP-team/meta_profile)


## 4. SGBs annotation
```
gtdbtk classify_wf --genome_dir SGBs/ --out_dir result/ --cpus 8 --extension fa
```


## 5. Build phylogenetic tree (https://github.com/QTP-team/build_phylogenetic_tree)


## 6. Mapping rate
### 6.1 build kraken index (https://github.com/QTP-team/build_kraken_index)
### 6.2 kraken
```
kraken2 --db kraken_database --threads 4 --report report/sample.report --output /dev/null --paired sample.1.fq.gz sample.2.fq.gz
```

## 7. Gene cataloge
### 7.1 Gene prediction
```
prodigal -c -m -p single -i genome.fa -a gene.faa > /dev/null
```

### 7.2 Gene de-redundancy
```
sh script/mmseq2_clust.sh library_ID input.fasta 0.95 output_dir 32 > log
```

## 8. Calculate pN/pS
```
bowtie2 -p 8 -x genome_index -1 sample.1.fq.gz -2 sample.2.fq.gz 2> logs/sample.log | samtools view -S -@ 8 -b > sample.bam
inStrain profile sample.bam genome.fa -o sample.IS -p 8 -l 0.95 -c 4 -f 0.01 -g genome.gene.fna -s genome.stb
```

## 9. Calculate Orthogroup Sequences
```
orthofinder -f Sequence/ -s diamond -t 20
```

## 10. Calculate Ka/Ks
```
perl script/ParaAT.pl -h Orthogroups.txt -n Orthogroup_Sequences.fna -a Orthogroup_Sequences.faa -p proc -f axt -g -k -c 11 -o Ka_Ks_results
```

## 11. Reconciles host and symbiont phylogenies
```
python script/AnGST.py Sample.input 1>sample.o 2>sample.e
```
