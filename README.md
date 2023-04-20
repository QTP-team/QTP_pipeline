# QTP pipeline
## 1. Assemble (https://github.com/QTP-team/meta_assemble)


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
sh script/mmseq2_clust.sh library_ID input.fasta 0.9 output_dir 32 > log
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
