# Trioza anthrisci 
Contains commands run by T.Heaven in assembly of Trioza anthrisci

Unless stated otherwise commands were performed from the directory /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae

## Collect data
```bash
#HiFi reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAn22_hifi_reads.fastq.gz raw_data/T_anthrisci/HiFi/anthrisci_hifi-reads.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAn22_hifi_3rdSMRTcell.fastq.gz raw_data/T_anthrisci/HiFi/anthrisci_hifi-3rdSMRTcell.fastq.gz

#Tellseq reads, from 3rd tellseq run:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_I1_T508.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_anthrisci/TellSeq/anthrisci_T508_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R1_T508.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_anthrisci/TellSeq/anthrisci_T508_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R2_T508.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_anthrisci/TellSeq/anthrisci_T508_R2.fastq.gz

#HiC reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_anthrisci/anthrisci-286154_S3HiC_R1.fastq\(1\)-002.gz raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_anthrisci/anthrisci-286154_S3HiC_R2.fastq\(1\)-003.gz raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
```
#### longQC
```bash
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi/*.fastq.gz); do
Datatype=pb-hifi
OutDir=$(dirname $Reads)/longqc/$(basename $Reads | cut -d '.' -f1)
OutFile=$(basename $Reads | cut -d '.' -f1)
ProgDir=~/git_repos/Wrappers/NBI
echo ${OutDir}/${OutFile}
mkdir $(dirname $Reads)/longqc
sbatch $ProgDir/run_longqc.sh $Reads $OutDir $OutFile $Datatype
done #55614592,3
```

## Hifiasm

### Default settings
Tom Mathers has performed hifiasm assembly using the HiFi reads with default settings.
```bash
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/Trant_default.p_ctg.fa assembly/genome/T_anthrisci/hifiasm/default/Trant_default.p_ctg.fa
```
#### Abyss
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
8639    8639    1111    5968    66931    196278    413917    269446    2788957    787.3e6    Trant_default.p_ctg.fa
#### KAT
Versus HiFi reads:
```bash
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/kat_comp-main.mx assembly/genome/T_anthrisci/hifiasm/default/kat/.
cp /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/kat_comp-main.mx.spectra-cn.png assembly/genome/T_anthrisci/hifiasm/default/kat/.
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/kat_comp.stats assembly/genome/T_anthrisci/hifiasm/default/kat/.
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/kat_comp.dist_analysis.json assembly/genome/T_anthrisci/hifiasm/default/kat/.
```
Versus TellSeq reads:
```bash
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/kat_comp_vs_tellseq-main.mx assembly/genome/T_anthrisci/hifiasm/default/kat/.
cp /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/kat_comp_vs_tellseq-main.mx.spectra-cn.png assembly/genome/T_anthrisci/hifiasm/default/kat/.
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/kat_comp_vs_tellseq.stats assembly/genome/T_anthrisci/hifiasm/default/kat/.
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/kat_comp_vs_tellseq.dist_analysis.json assembly/genome/T_anthrisci/hifiasm/default/kat/.
```
#### Jellyfish/kmc
```bash
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/jellyfish
Outfile=default_HiFi
Reads=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi/anthrisci_hifi-reads.fastq.gz
Reads2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi/anthrisci_hifi-3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_jellyfish.sh $OutDir $Outfile $Reads $Reads2
#55426797
```
Giving default_HiFi_19mer_out.histo to GenomeScope a size of 796,778,905 is estimated (Av. read length 9,355, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_21mer_out.histo to GenomeScope a size of 820,517,224 is estimated (Av. read length 9,355, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_25mer_out.histo to GenomeScope a size of 810,450,769 is estimated (Av. read length 9,355, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_31mer_out.histo to GenomeScope a size of 804,201,406 is estimated (Av. read length 9,355, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_39mer_out.histo to GenomeScope a size of 430,421,090 is estimated (Av. read length 9,355, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_49mer_out.histo to GenomeScope a size of 439,361,795 is estimated (Av. read length 9,355, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_61mer_out.histo to GenomeScope a size of 449,427,878 is estimated (Av. read length 9,355, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_75mer_out.histo to GenomeScope a size of 457,185,238 is estimated (Av. read length 9,355, 10,000,000 max. k-mer coverage from jellyfish settings).
```R
setwd("C:/Users/did23faz/OneDrive - Norwich Bioscience Institutes/Desktop/R")
dataframe19 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/jellyfish/default_HiFi_19mer_out.histo") 
dataframe21 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/jellyfish/default_HiFi_21mer_out.histo") 
dataframe25 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/jellyfish/default_HiFi_25mer_out.histo") 
dataframe31 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/jellyfish/default_HiFi_31mer_out.histo") 

#Plot kmer distribution:
plot(dataframe19[1:200,], type="l")
plot(dataframe21[1:200,], type="l")
plot(dataframe25[1:200,], type="l")
plot(dataframe31[1:200,], type="l")

#Ignore low frequency kmers:
plot(dataframe19[3:100,], type="l")
points(dataframe19[3:100,])
#Single copy region looks to be 5 - 45, peak = 12
sum(as.numeric(dataframe19[3:36566,1]*dataframe19[3:36566,2]))/12
#1,547,263,306
sum(as.numeric(dataframe19[3:45,1]*dataframe19[3:45,2]))/12
#425,033,858

plot(dataframe21[3:100,], type="l")
points(dataframe21[3:100,])
#Single copy region looks to be 5 - 45, peak = 11
sum(as.numeric(dataframe21[3:35197,1]*dataframe21[3:35197,2]))/11
#1,663,330,734
sum(as.numeric(dataframe21[3:45,1]*dataframe21[3:45,2]))/11
#484,823,960

plot(dataframe25[3:100,], type="l")
points(dataframe25[3:100,])
#Single copy region looks to be 5 - 45, peak = 10
sum(as.numeric(dataframe25[3:32814,1]*dataframe25[3:32814,2]))/10
#1,781,581,866
sum(as.numeric(dataframe25[3:45,1]*dataframe25[3:45,2]))/10
#556,806,278

plot(dataframe31[3:100,], type="l")
points(dataframe31[3:100,])
#Single copy region looks to be 5 - 45, peak = 9
sum(as.numeric(dataframe31[3:30041,1]*dataframe31[3:30041,2]))/9
#1,911,690,449
sum(as.numeric(dataframe31[3:45,1]*dataframe31[3:45,2]))/9
#641,797,281

dataframe39 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/jellyfish/default_HiFi_39mer_out.histo")
dataframe49 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/jellyfish/default_HiFi_49mer_out.histo")
dataframe61 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/jellyfish/default_HiFi_61mer_out.histo")
dataframe75 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/jellyfish/default_HiFi_75mer_out.histo")

plot(dataframe39[1:200,], type="l")
plot(dataframe49[1:200,], type="l")
plot(dataframe61[1:200,], type="l")
plot(dataframe75[1:200,], type="l")

plot(dataframe39[6:100,], type="l")
points(dataframe39[6:100,])
#Single copy region looks to be 6 - 50, peak = 12
sum(as.numeric(dataframe39[6:100000,1]*dataframe39[6:100000,2]))/12
#965,416,764
sum(as.numeric(dataframe39[6:50,1]*dataframe39[6:50,2]))/12
#666,341,284

plot(dataframe49[6:100,], type="l")
points(dataframe49[6:100,])
#Single copy region looks to be 6 - 50, peak = 11
sum(as.numeric(dataframe49[6:100000,1]*dataframe49[6:100000,2]))/11
#1,077,106,864
sum(as.numeric(dataframe49[6:50,1]*dataframe49[6:50,2]))/11
#762,257,618

plot(dataframe61[6:100,], type="l")
points(dataframe61[6:100,])
#Single copy region looks to be 6 - 50, peak = 10
sum(as.numeric(dataframe61[6:100000,1]*dataframe61[6:100000,2]))/10
#1,202,449,818
sum(as.numeric(dataframe61[6:50,1]*dataframe61[6:50,2]))/10
#870,638,819

plot(dataframe75[6:100,], type="l")
points(dataframe75[6:100,])
#Single copy region looks to be 6 - 50, peak = 9
sum(as.numeric(dataframe75[6:100000,1]*dataframe75[6:100000,2]))/9
#1,346,154,274
sum(as.numeric(dataframe75[6:50,1]*dataframe75[6:50,2]))/9
#995,246,302
```
K-mer estimates of genome size range from 440-960, genomescope plots with smaller and larger sizes both look like they capture all of the observed kmers.

Genomescope and 21-mer distribution plots give ~ the same genome size.
```bash
  for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/anthrisci_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/anthrisci_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm@g')/831m
    OutFile=T_anthrisci_1
    Haploid_Genomesize=831m #based on 21mers
    Homozygous_coverage=30 #based upon KAT spectra (0x approaches 0)
    Min_contig=2 #default
    Purge_haplotigs_level=3 #default (strict)
    Kmer_cuttoff=5.0 #default
    Overlap_iterations=200 #increasing can improve assembly quality
    Kmer_size=51 #default
    Similarity_threshold=0.75 #default
    mkdir -p $OutDir
    sbatch $ProgDir/run_hifiasm.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2
  done #55339299, 55339621

#n		n:500	L50		min		N80		N50		N20		max		sum
#8,649	8,649	1,126	7,787	66,214	193,738	412,133	2,801,313	1,082,000,000

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/831m/T_anthrisci_1.bp.p_ctg.fa); do
  Jobs=$(squeue -u did23faz| grep 'kat_comp'  | wc -l)
  echo x
  while [ $Jobs -gt 0 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'kat_comp'  | wc -l)
  done
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_HiFi_reads
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi/anthrisci_hifi-reads.fastq.gz
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi/anthrisci_hifi-3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1
done
```

```bash
  for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/anthrisci_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/anthrisci_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm@g')/773m
    OutFile=T_anthrisci_773m
    Haploid_Genomesize=773m #based on 19mers
    Homozygous_coverage=30 #based upon KAT spectra (0x approaches 0)
    Min_contig=2 #default
    Purge_haplotigs_level=3 #default (strict)
    Kmer_cuttoff=5.0 #default
    Overlap_iterations=200 #increasing can improve assembly quality
    Kmer_size=51 #default
    Similarity_threshold=0.75 #default
    mkdir -p $OutDir
    sbatch $ProgDir/run_hifiasm.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2
  done #55352511
```

