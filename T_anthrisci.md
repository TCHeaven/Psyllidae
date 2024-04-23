# Trioza anthrisci 
Contains commands run by T.Heaven in assembly of Trioza anthrisci

Unless stated otherwise commands were performed from the directory /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae

## Collect data
```bash
#HiFi reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAn22_hifi_reads.fastq.gz raw_data/T_anthrisci/HiFi/anthrisci_hifi-reads.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAn22_hifi_3rdSMRTcell.fastq.gz raw_data/T_anthrisci/HiFi/anthrisci_hifi-3rdSMRTcell.fastq.gz

#HiC reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_anthrisci/anthrisci-286154_S3HiC_R1.fastq\(1\)-002.gz raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_anthrisci/anthrisci-286154_S3HiC_R2.fastq\(1\)-003.gz raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz

#Tellseq reads, from 3rd tellseq run:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_I1_T508.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_anthrisci/TellSeq/anthrisci_T508_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R1_T508.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_anthrisci/TellSeq/anthrisci_T508_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R2_T508.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_anthrisci/TellSeq/anthrisci_T508_R2.fastq.gz

#RNASeq
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/N11/N11_1.fq.gz raw_data/T_anthrisci/RNASeq/T_anthrisci_N11_1.fq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/N11/N11_2.fq.gz raw_data/T_anthrisci/RNASeq/T_anthrisci_N11_2.fq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/N7/N7_1.fq.gz raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_1.fq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/N7/N7_2.fq.gz raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_2.fq.gz
```
Remove HiFi reads containing adapters:
```bash
for InFile in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi/*.fastq.gz); do
OutDir=$(dirname $InFile)/filtered
OutFile=$(basename $InFile | sed 's@.fastq.gz@@g')_filtered
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_HiFiAdapterFilt.sh $InFile $OutDir $OutFile
done 
#57334949,50
```
TellSeq reads converted to 10X format:
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run3_T_anthrisci_T_apicales/10x_conversion/T508/T508_S1_L001_R1_001.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x/T508_S1_L001_R1_001.fastq.gz
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run3_T_anthrisci_T_apicales/10x_conversion/T508/T508_S1_L001_R2_001.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x/T508_S1_L001_R2_001.fastq.gz
```
TellSeq reads - remove the internal barcodes and adapters:
```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger

mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/longranger/longranger-2.2.2/longranger-cs/2.2.2/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/longranger/longranger-2.2.2/longranger-cs/2.2.2/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt.backup
cp /nbi/software/testing/supernova/2.1.1_TCM/x86_64/bin/supernova-cs/2.1.1/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/longranger/longranger-2.2.2/longranger-cs/2.2.2/tenkit/lib/python/tenkit/barcodes/.

for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x/T508_S1_L001_R1_001.fastq.gz); do
Fread=$Reads
Rread=$(echo $Reads | sed 's@_S1_L001_R1_001.fastq.gz@_S1_L001_R2_001.fastq.gz@g')
Run=$(basename $Reads | cut -d '_' -f1)
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger
OutFile=Tant_${Run}
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_longranger_basic.sh $Fread $Rread $Run $OutDir $OutFile
done 
#57161309
```
#### Fastqc
```bash
for ReadFile in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/*/*.fastq.gz); do
OutDir=$(dirname $ReadFile)/fastqc
OutFile=$(basename $ReadFile | sed 's@.fastq.gz@@g')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_fastqc.sh $ReadFile $OutDir $OutFile
done
#57205238-57205244
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
done 
#57206546-7
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
#### BUSCO
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/Trant_default.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 

rm temp.txt
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default); do
echo "" >> temp.txt
echo "" >> temp.txt
echo $assembly >> temp.txt
cat ${assembly}/*abyss_report.txt >> temp.txt
echo Arthropoda: >> temp.txt
cat ${assembly}/BUSCO/*arthropoda_odb10_short_summary.txt | grep 'C:' >> temp.txt
echo Insecta: >> temp.txt
cat ${assembly}/BUSCO/*insecta_odb10_short_summary.txt | grep 'C:' >> temp.txt
echo Hemiptera: >> temp.txt
cat ${assembly}/BUSCO/*hemiptera_odb10_short_summary.txt | grep 'C:' >> temp.txt
done

cat temp.txt 
#n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
#8639    8639    1111    5968    66931    196278    413917    269446    2788957    787.3e6    Trant_default.p_ctg.fa
#Arthropoda:
#        C:92.4%[S:72.9%,D:19.5%],F:3.9%,M:3.7%,n:1013
#Insecta:
#        C:91.9%[S:73.1%,D:18.8%],F:4.5%,M:3.6%,n:1367
#Hemiptera:
#        C:92.7%[S:74.2%,D:18.5%],F:3.5%,M:3.8%,n:2510
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

#n    n:500 L50   min   N80   N50   N20   max   sum
#8,649  8,649 1,126 7,787 66,214  193,738 412,133 2,801,313 1,082,000,000

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

Assembly was tried with a range of genome sizes based upon the k-mer distribution and genomescope estimates:

```bash
Haploid_Genomesize_values=("440m" "485m" "530m" "600m" "675m" "770m" "820m" "890m" "955m")
Homozygous_coverage_values=(40)
Purge_haplotigs_level_values=(3)
Kmer_cuttoff_values=("5.0")
Similarity_threshold_values=("0.75")

for Haploid_Genomesize in "${Haploid_Genomesize_values[@]}"
do
for Homozygous_coverage in "${Homozygous_coverage_values[@]}"
do
for Purge_haplotigs_level in "${Purge_haplotigs_level_values[@]}"
do
for Kmer_cuttoff in "${Kmer_cuttoff_values[@]}"
do
for Similarity_threshold in "${Similarity_threshold_values[@]}"
do
ReadDir=$(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi)
ProgDir=~/git_repos/Wrappers/NBI
Run1=$(ls $ReadDir/anthrisci_hifi-reads.fastq.gz)
Run2=$(ls $ReadDir/anthrisci_hifi-3rdSMRTcell.fastq.gz)
OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm_19.5@g')/${Haploid_Genomesize}/${Homozygous_coverage}/${Purge_haplotigs_level}/${Kmer_cuttoff}/${Similarity_threshold}
OutFile=T_anthrisci_${Haploid_Genomesize}_${Homozygous_coverage}_${Purge_haplotigs_level}_${Kmer_cuttoff}_${Similarity_threshold}
Min_contig=2 #default
Kmer_size=51 #default
Overlap_iterations=200 #increasing can improve assembly quality 
Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
echo x
while [ $Jobs -gt 11 ]; do
sleep 900s
printf "."
Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
done
echo ${OutDir}/$OutFile >> logs/anthrisci_hifilog.txt
if [ -s "${OutDir}/${OutFile}.bp.p_ctg.fa" ]; then
echo Already done for: $OutFile
else 
echo Running for: $OutFile
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/anthrisci_hifilog.txt
fi
done
done
done
done
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/*/*/*/*/*/*.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 

rm temp.txt
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/*/*/*/*/*); do
echo "" >> temp.txt
echo "" >> temp.txt
echo $assembly >> temp.txt
cat ${assembly}/*abyss_report.txt >> temp.txt
echo Arthropoda: >> temp.txt
cat ${assembly}/BUSCO/*arthropoda_odb10_short_summary.txt | grep 'C:' >> temp.txt
echo Insecta: >> temp.txt
cat ${assembly}/BUSCO/*insecta_odb10_short_summary.txt | grep 'C:' >> temp.txt
echo Hemiptera: >> temp.txt
cat ${assembly}/BUSCO/*hemiptera_odb10_short_summary.txt | grep 'C:' >> temp.txt
done

cat temp.txt 
```
HiFiasm gives different estimates for homozygous coverage:
```bash
grep -m 1 'peak_hom:' slurm.56260036.err slurm.56260037.err slurm.56260038.err slurm.56260039.err slurm.56262037.err slurm.56262425.err slurm.56262426.err slurm.56266456.err slurm.56266457.err

#440m: 
#slurm.56260036.err:[M::ha_ft_gen] peak_hom: 53; peak_het: 12
#485m: 
#slurm.56260037.err:[M::ha_ft_gen] peak_hom: 48; peak_het: 12
#530m: 
#slurm.56260038.err:[M::ha_ft_gen] peak_hom: 44; peak_het: 12
#600m: 
#slurm.56260039.err:[M::ha_ft_gen] peak_hom: 39; peak_het: 12
#675m: 
#slurm.56262037.err:[M::ha_ft_gen] peak_hom: 34; peak_het: 12
#770m: 
#slurm.56262425.err:[M::ha_ft_gen] peak_hom: 30; peak_het: 12
#820m: 
#slurm.56262426.err:[M::ha_ft_gen] peak_hom: 28; peak_het: 12
#890m: 
#slurm.56266456.err:[M::ha_ft_gen] peak_hom: 26; peak_het: 12
#955m:
#slurm.56266457.err:[M::ha_ft_gen] peak_hom: 24; peak_het: 12

#If the heterozygous peak is 12 the homozygous peak should be 24?
```
The assembly results from this range are very similar, however 820m has both the highest BUSCO completeness and lowest contigs, whilst 890m has highest N50.
```bash
Haploid_Genomesize_values=("820m" "890m")
Homozygous_coverage_values=(24 28 39 48)
Purge_haplotigs_level_values=(0 1 2 3)
Kmer_cuttoff_values=("3.0" "4.0" "5.0" "10.0")
Similarity_threshold_values=("0.75" "0.50" "0.25")

for Haploid_Genomesize in "${Haploid_Genomesize_values[@]}"
do
for Homozygous_coverage in "${Homozygous_coverage_values[@]}"
do
for Purge_haplotigs_level in "${Purge_haplotigs_level_values[@]}"
do
for Kmer_cuttoff in "${Kmer_cuttoff_values[@]}"
do
for Similarity_threshold in "${Similarity_threshold_values[@]}"
do
ReadDir=$(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi)
ProgDir=~/git_repos/Wrappers/NBI
Run1=$(ls $ReadDir/anthrisci_hifi-reads.fastq.gz)
Run2=$(ls $ReadDir/anthrisci_hifi-3rdSMRTcell.fastq.gz)
OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm_19.5@g')/${Haploid_Genomesize}/${Homozygous_coverage}/${Purge_haplotigs_level}/${Kmer_cuttoff}/${Similarity_threshold}
OutFile=T_anthrisci_${Haploid_Genomesize}_${Homozygous_coverage}_${Purge_haplotigs_level}_${Kmer_cuttoff}_${Similarity_threshold}
Min_contig=2 #default
Kmer_size=51 #default
Overlap_iterations=200 #increasing can improve assembly quality 
Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
echo x
while [ $Jobs -gt 11 ]; do
sleep 900s
printf "."
Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
done
echo ${OutDir}/$OutFile >> logs/anthrisci_hifilog.txt
if [ -s "${OutDir}/${OutFile}.bp.p_ctg.fa" ]; then
echo Already done for: $OutFile
else 
echo Running for: $OutFile
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/anthrisci_hifilog.txt
fi
done
done
done
done
done

ls: cannot access /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/890m/48/2/10.0/0.25/T_anthrisci_890m_48_2_10.0_0.25.bp.p_ctg.fa: No such file or directory

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/*/*/*/*/*/*.bp.p_ctg.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/*/*/*/*/*/*.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 

rm temp.txt
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/*/*/*/*/*); do
echo "" >> temp.txt
echo "" >> temp.txt
echo $assembly >> temp.txt
cat ${assembly}/*abyss_report.txt >> temp.txt
echo Arthropoda: >> temp.txt
cat ${assembly}/BUSCO/*arthropoda_odb10_short_summary.txt | grep 'C:' >> temp.txt
echo Insecta: >> temp.txt
cat ${assembly}/BUSCO/*insecta_odb10_short_summary.txt | grep 'C:' >> temp.txt
echo Hemiptera: >> temp.txt
cat ${assembly}/BUSCO/*hemiptera_odb10_short_summary.txt | grep 'C:' >> temp.txt
done

cp temp.txt Reports/anthrisci_assembly_report3.txt
```
Amongst hifiasm_19.5 assemblies the highist BUSCO score is 93.3% of Hemiptera BUSCOs, of 30 assemblies with this score the most contiguous is T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa with 9,575 contigs, N50=185,165. The most contiguous hifiasm_19.5 assembly was T_anthrisci_820m_48_3_10.0_0.25.bp.p_ctg.fa with 7,714 contigs, N50=240,906 and 89.2% of Hemiptera BUSCOs. The best N50 is 244,067 for T_anthrisci_890m_48_3_4.0_0.25.bp.p_ctg.fa (90.1% of Hemiptera BUSCOs, 7,739 contigs). 

Based on Merqury/KAT plot the true homozygous coverage is ~24

```bash
#These assemblies were kept:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/default/Trant_default.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/3/10.0/0.25/T_anthrisci_820m_48_3_10.0_0.25.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/890m/48/3/4.0/0.25/T_anthrisci_890m_48_3_4.0_0.25.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/890m/24/0/5.0/0.75/T_anthrisci_890m_24_0_5.0_0.75.bp.p_ctg.fa

#All other assemblies were deleted:
for assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/*/*/*/*/*/*.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm/*/*.fa | grep -v 'Trant_default.p_ctg.fa\|T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa\|T_anthrisci_820m_48_3_10.0_0.25.bp.p_ctg.fa\|T_anthrisci_890m_48_3_4.0_0.25.bp.p_ctg.fa\|T_anthrisci_890m_24_0_5.0_0.75.bp.p_ctg.fa'); do
rm $assembly  
done
```
#### KAT
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/*/*/*/*/*/*.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_hifi_reads
F1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAn22_hifi_reads.fastq.gz
R1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAn22_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp.sh $OutDir $Outfile $Genome $F1 $R1
done #56881242, 56881243, 56881244

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/*/*/*/*/*/*.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_tellseq_reads
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/anthrisci_T508_R1.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/anthrisci_T508_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1
done #56881247, 56881248, 56881249
```
#### Merqury
Kmer plots with all reads; Hifi, HiC, Tellseq
```bash
#Prepare meryl kmer counts
source /nbi/software/staging/RCSUPPORT-2452/stagingloader
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/*/*.fastq.gz); do
OutDir=$(dirname $Reads)/meryl
mkdir $OutDir
meryl k=21 count output ${OutDir}/$(basename $Reads | sed 's@.fastq.gz@@g').meryl $Reads 
done

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
mkdir $(dirname $Assembly)/meryl
meryl union-sum output $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/*/meryl/*.meryl
mkdir $(dirname $Assembly)/meryl/HiFi 
meryl union-sum output $(dirname $Assembly)/meryl/HiFi/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi/meryl/*.meryl
mkdir $(dirname $Assembly)/meryl/HiC
meryl union-sum output $(dirname $Assembly)/meryl/HiC/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/meryl/*.meryl
mkdir $(dirname $Assembly)/meryl/Tellseq
meryl union-sum output $(dirname $Assembly)/meryl/Tellseq/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/meryl/*.meryl

#Plot with merqury
mkdir $(dirname $Assembly)/merqury/
cd $(dirname $Assembly)/merqury/
ln -s $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl
ln -s $Assembly
meryl histogram $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl.hist
merqury.sh $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl/ $(basename $Assembly) $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g') #NOTE: needs >25GB memory to complete, errors will output to log/ file not to slurm.err

mkdir $(dirname $Assembly)/merqury/HiFi
cd $(dirname $Assembly)/merqury/HiFi
ln -s $(dirname $Assembly)/meryl/HiFi/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl
ln -s $Assembly
meryl histogram $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl.hist
merqury.sh $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl/ $(basename $Assembly) $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')

mkdir $(dirname $Assembly)/merqury/HiC
cd $(dirname $Assembly)/merqury/HiC
ln -s $(dirname $Assembly)/meryl/HiC/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl
ln -s $Assembly
meryl histogram $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl.hist
merqury.sh $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl/ $(basename $Assembly) $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')

mkdir $(dirname $Assembly)/merqury/TellSeq
cd $(dirname $Assembly)/merqury/TellSeq
ln -s $(dirname $Assembly)/meryl/Tellseq/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl
ln -s $Assembly
meryl histogram $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl.hist
merqury.sh $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl/ $(basename $Assembly) $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
#57159096

source /nbi/software/staging/RCSUPPORT-2452/stagingloader
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
cd $(dirname $Assembly)/merqury/
spectra_cn=$(ls *.spectra-cn.hist)
cn_only_hist=$(ls *.only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_1000 -z $cn_only_hist -m 1000 -n 4000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_300 -z $cn_only_hist -m 300 -n 4000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_100 -z $cn_only_hist -m 100 -n 4000000 
spectra_asm=$(ls *.spectra-asm.hist)
asm_only_host=$(ls *.dist_only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_1000 -z $asm_only_host -m 1000 -n 6000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_300 -z $asm_only_host -m 300 -n 6000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_100 -z $asm_only_host -m 100 -n 6000000 

cd $(dirname $Assembly)/merqury/HiFi
spectra_cn=$(ls *.spectra-cn.hist)
cn_only_hist=$(ls *.only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_1000 -z $cn_only_hist -m 1000 -n 12500000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_300 -z $cn_only_hist -m 300 -n 12500000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_100 -z $cn_only_hist -m 100 -n 12500000 
spectra_asm=$(ls *.spectra-asm.hist)
asm_only_host=$(ls *.dist_only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_1000 -z $asm_only_host -m 1000 -n 16000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_300 -z $asm_only_host -m 300 -n 16000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_100 -z $asm_only_host -m 100 -n 16000000 

cd $(dirname $Assembly)/merqury/HiC
spectra_cn=$(ls *.spectra-cn.hist)
cn_only_hist=$(ls *.only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_1000 -z $cn_only_hist -m 100 -n 30000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_300 -z $cn_only_hist -m 30 -n 30000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_100 -z $cn_only_hist -m 10 -n 30000000 
spectra_asm=$(ls *.spectra-asm.hist)
asm_only_host=$(ls *.dist_only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_1000 -z $asm_only_host -m 100 -n 20000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_300 -z $asm_only_host -m 30 -n 20000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_100 -z $asm_only_host -m 10 -n 20000000 

cd $(dirname $Assembly)/merqury/TellSeq
spectra_cn=$(ls *.spectra-cn.hist)
cn_only_hist=$(ls *.only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_1000 -z $cn_only_hist -m 1000 -n 6000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_300 -z $cn_only_hist -m 300 -n 6000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_100 -z $cn_only_hist -m 100 -n 6000000 
spectra_asm=$(ls *.spectra-asm.hist)
asm_only_host=$(ls *.dist_only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_1000 -z $asm_only_host -m 1000 -n 10000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_300 -z $asm_only_host -m 300 -n 10000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_100 -z $asm_only_host -m 100 -n 10000000 
#57192972

#NOTE: meryl databases take up a lot of space so remove after using them
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/meryl/HiC/T_anthrisci_820m_48_1_10.0_0.25.meryl
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/meryl/HiFi/T_anthrisci_820m_48_1_10.0_0.25.meryl
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/meryl/Tellseq/T_anthrisci_820m_48_1_10.0_0.25.meryl
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/meryl/T_anthrisci_820m_48_1_10.0_0.25.meryl

######################################################################################################################################################################
source /nbi/software/staging/RCSUPPORT-2452/stagingloader
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_*/TellSeq/longranger/*.fastq.gz); do
OutDir=$(dirname $Reads)/meryl
mkdir $OutDir
meryl k=21 count output ${OutDir}/$(basename $Reads | sed 's@.fastq.gz@@g').meryl $Reads 
done
#57186506

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
mkdir $(dirname $Assembly)/meryl/Tellseq_trimmed
meryl union-sum output $(dirname $Assembly)/meryl/Tellseq_trimmed/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger/meryl/*.meryl

mkdir -p $(dirname $Assembly)/merqury/Tellseq_trimmed
cd $(dirname $Assembly)/merqury/Tellseq_trimmed
ln -s $(dirname $Assembly)/meryl/Tellseq_trimmed/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl
ln -s $Assembly
meryl histogram $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl.hist
merqury.sh $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl/ $(basename $Assembly) $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
rm -r $(dirname $Assembly)/meryl/Tellseq_trimmed/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl

spectra_cn=$(ls *.spectra-cn.hist)
cn_only_hist=$(ls *.only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_1000 -z $cn_only_hist -m 1000 -n 6000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_300 -z $cn_only_hist -m 300 -n 6000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_100 -z $cn_only_hist -m 100 -n 6000000 
spectra_asm=$(ls *.spectra-asm.hist)
asm_only_host=$(ls *.dist_only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_1000 -z $asm_only_host -m 1000 -n 10000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_300 -z $asm_only_host -m 300 -n 10000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_100 -z $asm_only_host -m 100 -n 10000000 

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae
mkdir temp
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi temp/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC temp/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger temp/.
mkdir $(dirname $Assembly)/meryl/All_Tellseq_trimmed
meryl union-sum output $(dirname $Assembly)/meryl/All_Tellseq_trimmed/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl temp/*/meryl/*.meryl
rm -r temp

mkdir $(dirname $Assembly)/merqury/All_Tellseq_trimmed
cd $(dirname $Assembly)/merqury/All_Tellseq_trimmed
ln -s $(dirname $Assembly)/meryl/All_Tellseq_trimmed/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl
ln -s $Assembly
meryl histogram $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl.hist
merqury.sh $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl/ $(basename $Assembly) $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
rm -r $(dirname $Assembly)/meryl/All_Tellseq_trimmed/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl

spectra_cn=$(ls *.spectra-cn.hist)
cn_only_hist=$(ls *.only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_1000 -z $cn_only_hist -m 1000 -n 4000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_300 -z $cn_only_hist -m 300 -n 4000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_cn -o $(echo $spectra_cn | sed 's@.hist@@g')_100 -z $cn_only_hist -m 100 -n 4000000 
spectra_asm=$(ls *.spectra-asm.hist)
asm_only_host=$(ls *.dist_only.hist)
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_1000 -z $asm_only_host -m 1000 -n 6000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_300 -z $asm_only_host -m 300 -n 6000000 
Rscript /opt/software/merqury/plot/plot_spectra_cn.R -f $spectra_asm -o $(echo $spectra_asm | sed 's@.hist@@g')_100 -z $asm_only_host -m 100 -n 6000000 
#57195977

#NOTE: meryl databases take up a lot of space so remove after using them
rm /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_*/*/meryl/*.meryl/*
rm /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_*/TellSeq/longranger/meryl/*.meryl/*
```
#### Blobtools
Blobtools with Hifi reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAn22_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAn22_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2 #56939403
sbatch $ProgDir/run_qualimap.sh $(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam $Assembly $OutDir #57111024

#Blast
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/blast2.12.0/3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blastn.sh $Assembly $Database $OutDir $OutPrefix 
#57182845,57308572, 57398064,57044991
#parse_seqids is required for -taxid_map however when running databases created with -parse_seqids blast runs fail with c++ errors

#BUSCO - keeping output files
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done 
#57207205, 57307069, 57307080, 57307097

#Diamond blast - BUSCO regions
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
BUSCOgff=$(dirname $Assembly)/BUSCO/hemiptera_odb10/run_hemiptera_odb10/metaeuk_output/rerun_results/T_anthrisci_820m_48_1_10_hemiptera_odb10.fa.gff
awk '$3 == "gene" {print $0}' $BUSCOgff > $(echo $BUSCOgff | sed 's@.fa.gff@_genes_only.fa.gff@g')
bedtools getfasta -fi $Assembly -bed $(echo $BUSCOgff | sed 's@.fa.gff@_genes_only.fa.gff@g') -fo $(dirname $BUSCOgff)/busco_regions.fasta
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/uniprot_01102023/Uniprot_01102023_reference_proteomes.dmnd
ProgDir=~/git_repos/Wrappers/NBI
mkdir -p $OutDir
sbatch $ProgDir/run_diamond_blastx.sh $(dirname $BUSCOgff)/busco_regions.fasta $Database $OutDir $OutPrefix 
#57307107, 57315795 - no deletion of WorkDir

#Diamond blast - whole genome
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/fasta_to_bed.py $Assembly > $(dirname $Assembly)/$(basename $Assembly | sed 's@.fa@.bed@g')
bedtools makewindows -g $(dirname $Assembly)/$(basename $Assembly | sed 's@.fa@.bed@g') -w 25000 > $(dirname $Assembly)/25000regions.bed
bedtools intersect -a $(dirname $Assembly)/25000regions.bed -b $(echo $BUSCOgff | sed 's@.fa.gff@_genes_only.fa.gff@g') -wa -u > $(dirname $BUSCOgff)/overlapping_windows.bed
bedtools subtract -a $(dirname $Assembly)/25000regions.bed -b $(dirname $BUSCOgff)/overlapping_windows.bed > $(dirname $Assembly)/buscofiltered25000regions.bed
bedtools getfasta -fi $Assembly -bed $(dirname $Assembly)/buscofiltered25000regions.bed -fo $(dirname $Assembly)/buscofiltered25000regions.fasta
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/diamond0.9.29_blastx/nonBUSCO_regions
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/uniprot_01102023/Uniprot_01102023_reference_proteomes.dmnd
ProgDir=~/git_repos/Wrappers/NBI
mkdir -p $OutDir
sbatch $ProgDir/run_diamond_blastx.sh $(dirname $Assembly)/buscofiltered25000regions.fasta $Database $OutDir $OutPrefix 
#57307108, 57315796

BUSCODiamond=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions/*.diamondblastx.out
ElseDiamond=$(dirname $Assembly)/diamond0.9.29_blastx/nonBUSCO_regions/*.diamondblastx.out
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/normalise_blast.py $BUSCODiamond $(echo $BUSCODiamond | sed 's@x.out@x_2.out@g')
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/normalise_blast.py $ElseDiamond $(echo $ElseDiamond | sed 's@x.out@x_2.out@g')

#Tiara
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/tiara
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_tiara.sh $Assembly $OutDir $OutPrefix
#57401244

#Blobtools
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
MappingFile=$(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
BlastFile=$(dirname $Assembly)/blast2.12.0/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').vs.nt.mts1.hsp1.1e25.megablast.out
OutDir=$(dirname $Assembly)/blobtools1.1.1
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
ColourFile=NA
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobplot.sh $Assembly $MappingFile $BlastFile $OutDir $OutPrefix $ColourFile
#57189468, 57191145

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
Record_type=contig
MappingFile=$(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
BlastFile=$(dirname $Assembly)/blast2.12.0/3/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').vs.nt.mts1.hsp1.1e25.megablast.out
BUSCOFile=$(dirname $Assembly)/BUSCO/hemiptera_odb10/run_hemiptera_odb10/full_table.tsv
BUSCODiamond=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions/*.diamondblastx_2.out
ElseDiamond=$(dirname $Assembly)/diamond0.9.29_blastx/nonBUSCO_regions/*.diamondblastx_2.out
Tiara=$(dirname $Assembly)/tiara/*.tiara
OutDir=$(dirname $Assembly)/blobtoolkit4.2.1
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Genus=Trioza
Species=anthrisci
TaxID=2023874
alias=Tant820481100025
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57080570, 57092797

cp -r $OutDir/Tant820481100025_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
#The contents of /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs were subsequently copied to local machine WSL-Ubuntu: \\wsl.localhost\Ubuntu\home\did23faz
```
```bash
ubuntu
conda activate btk
#From \\wsl.localhost\Ubuntu\home\did23faz
apptainer exec blobtoolkit.sif blobtools host BlobDirs
```
Blobtools with HiC reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_HiC
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#57207038
MappingFile=$(ls ${OutDir}/*_HiC.bam)
sbatch $ProgDir/run_qualimap.sh $MappingFile $Assembly $OutDir #57111425

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
Record_type=contig
MappingFile=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_HiC.bam
BlastFile=$(dirname $Assembly)/blast2.12.0/3/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').vs.nt.mts1.hsp1.1e25.megablast.out
BUSCOFile=$(dirname $Assembly)/BUSCO/hemiptera_odb10/run_hemiptera_odb10/full_table.tsv
BUSCODiamond=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions/*.diamondblastx_2.out
ElseDiamond=$(dirname $Assembly)/diamond0.9.29_blastx/nonBUSCO_regions/*.diamondblastx_2.out
Tiara=$(dirname $Assembly)/tiara/*.tiara
OutDir=$(dirname $Assembly)/blobtoolkit4.2.1
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Genus=Trioza
Species=anthrisci
TaxID=2023874
alias=Tant820481100025_HiC
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57217715

cp -r $OutDir/Tant820481100025_HiC_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
```
Blobtools with TellSeq reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/anthrisci_T508_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/anthrisci_T508_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#57206460
MappingFile=$(ls ${OutDir}/*_Tellseq.bam)
sbatch $ProgDir/run_qualimap.sh $MappingFile $Assembly $OutDir #57111432

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
Record_type=contig
MappingFile=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq.bam
BlastFile=$(dirname $Assembly)/blast2.12.0/3/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').vs.nt.mts1.hsp1.1e25.megablast.out
BUSCOFile=$(dirname $Assembly)/BUSCO/hemiptera_odb10/run_hemiptera_odb10/full_table.tsv
BUSCODiamond=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions/*.diamondblastx_2.out
ElseDiamond=$(dirname $Assembly)/diamond0.9.29_blastx/nonBUSCO_regions/*.diamondblastx_2.out
Tiara=$(dirname $Assembly)/tiara/*.tiara
OutDir=$(dirname $Assembly)/blobtoolkit4.2.1
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Genus=Trioza
Species=anthrisci
TaxID=2023874
alias=Tant820481100025__Tellseq
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57218784

#########################################################################################################################
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger/Tant_T508_barcoded.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $Read1
#57186563
MappingFile=$(ls ${OutDir}/*_Tellseq_trimmed.bam)
sbatch $ProgDir/run_qualimap.sh $MappingFile $Assembly $OutDir #57317736

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
Record_type=contig
MappingFile=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed.bam
BlastFile=$(dirname $Assembly)/blast2.12.0/3/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').vs.nt.mts1.hsp1.1e25.megablast.out
BUSCOFile=$(dirname $Assembly)/BUSCO/hemiptera_odb10/run_hemiptera_odb10/full_table.tsv
BUSCODiamond=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions/*.diamondblastx_2.out
ElseDiamond=$(dirname $Assembly)/diamond0.9.29_blastx/nonBUSCO_regions/*.diamondblastx_2.out
Tiara=$(dirname $Assembly)/tiara/*.tiara
OutDir=$(dirname $Assembly)/blobtoolkit4.2.1
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Genus=Trioza
Species=anthrisci
TaxID=2023874
alias=Tant820481100025__Tellseq_trimmed
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57331163

cp -r $OutDir/Tant820481100025__Tellseq_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
cp -r $OutDir/Tant820481100025__Tellseq_trimmed_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
```
#### Kraken 
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_kraken2nt
OutDir=$(dirname $Assembly)/kraken2.1.3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/kraken/nt_14092023
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_kraken2.sh $Assembly $Database $OutDir $OutPrefix 2>&1 >> ${OutDir}/log.txt
```
Contigs classified to non-Arthropoda phylum level were considered contamininants. Contigs classified to Arthropoda taxa or above (eg. Eukaryota) were considered Psyllid contigs.
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
cd $(dirname $Assembly)/kraken2.1.3/
nano contaminantlist.txt #Edit with contaminant names
grep -f contaminantlist.txt $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_kraken2nt_output.txt > contaminantcontigs.txt
awk -F'\t' '{print $2}' contaminantcontigs.txt > contaminantcontignames.txt
rm contaminantcontigs.txt
seqtk subseq $Assembly contaminantcontignames.txt | gzip > kraken2_contaminants.fa.gz
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/filter.py  $Assembly contaminantcontignames.txt > filtered_$(basename $Assembly)
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' filtered_$(basename $Assembly) > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_krakenfiltered.fa
rm contaminantcontignames.txt
rm filtered_$(basename $Assembly)

#BUSCO
for Genome in $(ls $(dirname $Assembly)/kraken2.1.3/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_krakenfiltered.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 
```
##### Tiara
Contigs identified as bacteria and archaea by tiara were highlighted, those unknown to tiara that did not blast to arthropoda were also highlighted for removal
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
cd $(dirname $Assembly)/tiara/
nano contaminantcontignames.txt #edit with bacteria, archaea and unknown contig names from tiara
seqtk subseq $Assembly contaminantcontignames.txt | gzip > tiara_contaminants.fa.gz
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/filter.py  $Assembly contaminantcontignames.txt > filtered_$(basename $Assembly)
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' filtered_$(basename $Assembly) > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_tiarafiltered.fa
rm filtered_$(basename $Assembly)

#BUSCO
for Genome in $(ls $(dirname $Assembly)/tiara/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')*filtered.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 
```
##### BTK
Contigs not within the main arthropoda cluster in blobtoolkit were highlighted for removal
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
cd $(dirname $Assembly)/blobtoolkit4.2.1/
nano contaminantcontignames.txt #edit with contig names from blobtoolkit 
seqtk subseq $Assembly contaminantcontignames.txt | gzip > btk_contaminants.fa.gz
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/filter.py  $Assembly contaminantcontignames.txt > filtered_$(basename $Assembly)
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' filtered_$(basename $Assembly) > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_btkfiltered.fa
rm filtered_$(basename $Assembly)

#BUSCO
for Genome in $(ls $(dirname $Assembly)/blobtoolkit4.2.1/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')*filtered.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 
```
Sliders:
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
cd $(dirname $Assembly)/blobtoolkit4.2.1/
nano contignames2.txt #edit with contig names from blobtoolkit 
grep -A 1 -F -f contignames2.txt $Assembly > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_btkfiltered2.fa

#BUSCO
for Genome in $(ls $(dirname $Assembly)/blobtoolkit4.2.1/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')*filtered2.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 

grep '>' $Assembly | sed 's@>@@g' > $(dirname $Assembly)/allcontigs.txt
awk 'NR==FNR {busco[$0]=1; next} !($0 in busco)' "$(dirname $Assembly)/blobtoolkit4.2.1/contignames2.txt" "$(dirname $Assembly)/allcontigs.txt" > "$(dirname $Assembly)/blobtoolkit4.2.1/contaminantcontignames.txt"
```
##### blastn
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
cd $(dirname $Assembly)/blast2.12.0/
nano contaminantcontignames.txt #edit with contig names which blast to non-ascomycota phyla 
seqtk subseq $Assembly contaminantcontignames.txt | gzip > blast_contaminants.fa.gz
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/filter.py $Assembly contaminantcontignames.txt > filtered_$(basename $Assembly)
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' filtered_$(basename $Assembly) > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_blastfiltered.fa
rm filtered_$(basename $Assembly)

#BUSCO
for Genome in $(ls $(dirname $Assembly)/blast2.12.0/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')*filtered.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 
```
##### All
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
cd $(dirname $Assembly)
cat kraken2.1.3/contaminantcontignames.txt > contaminantcontignames.txt
cat tiara/contaminantcontignames.txt >> contaminantcontignames.txt
cat blobtoolkit4.2.1/contaminantcontignames.txt >> contaminantcontignames.txt
cat blast2.12.0/contaminantcontignames.txt >> contaminantcontignames.txt

seqtk subseq $Assembly contaminantcontignames.txt | gzip > contaminants.fa.gz
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/filter.py $Assembly contaminantcontignames.txt > filtered_$(basename $Assembly)
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' filtered_$(basename $Assembly) > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_filtered.fa
rm filtered_$(basename $Assembly)

#BUSCO
for Genome in $(ls $(dirname $Assembly)/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')*filtered.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 
```
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/blobtoolkit4.2.1/T_anthrisci_820m_48_1_10.0_0.25_btkfiltered2.fa
cd $(dirname $Assembly)
cd ..
cat kraken2.1.3/contaminantcontignames.txt > contaminantcontignames.txt
cat tiara/contaminantcontignames.txt >> contaminantcontignames.txt
cat blast2.12.0/contaminantcontignames.txt >> contaminantcontignames.txt

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/filter.py $Assembly contaminantcontignames.txt > filtered_$(basename $Assembly)
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' filtered_$(basename $Assembly) > $(basename $Assembly | sed 's@_btkfiltered2.fa@@g')_allfiltered.fa
rm filtered_$(basename $Assembly)

#BUSCO
for Genome in $(ls $(dirname $Assembly)/../$(basename $Assembly | sed 's@_btkfiltered2.fa@@g')_allfiltered.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 
```
T_anthrisci_820m_48_1_10.0_0.25_hemiptera_odb10_short_summary.txt
        C:92.3%[S:72.6%,D:19.7%],F:3.6%,M:4.1%,n:2510

T_anthrisci_820m_48_1_10.0_0.25_krakenfiltered_hemiptera_odb10_short_summary.txt
        C:91.8%[S:72.6%,D:19.2%],F:3.5%,M:4.7%,n:2510
T_anthrisci_820m_48_1_10.0_0.25_blastfiltered_hemiptera_odb10_short_summary.txt
        C:91.6%[S:72.6%,D:19.0%],F:3.7%,M:4.7%,n:2510
T_anthrisci_820m_48_1_10.0_0.25_tiarafiltered_hemiptera_odb10_short_summary.txt
        C:92.1%[S:72.6%,D:19.5%],F:3.7%,M:4.2%,n:2510

T_anthrisci_820m_48_1_10.0_0.25_btkfiltered_hemiptera_odb10_short_summary.txt
        C:92.1%[S:72.6%,D:19.5%],F:3.6%,M:4.3%,n:2510
T_anthrisci_820m_48_1_10.0_0.25_filtered_hemiptera_odb10_short_summary.txt
        C:2.1%[S:1.6%,D:0.5%],F:0.0%,M:97.9%,n:2510

T_anthrisci_820m_48_1_10.0_0.25_btkfiltered2_hemiptera_odb10_short_summary.txt
        C:92.3%[S:72.6%,D:19.7%],F:3.6%,M:4.1%,n:2510
#T_anthrisci_820m_48_1_10.0_0.25_allfiltered_hemiptera_odb10_short_summary.txt
        C:91.2%[S:72.6%,D:18.6%],F:3.6%,M:5.2%,n:2510
```bash
#BUSCO - keeping output files
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done 
#57044702

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done 
#57044706

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
for file in $(ls $(dirname $Assembly)/BUSCO/*_odb10/run_*_odb10/full_table.tsv); do
grep 'Complete\|Duplicated' $file | awk '!/^#/ {split($3, arr, ":"); print arr[1]}' >> $(dirname $Assembly)/BUSCO/buscocontignames.txt
done

cd $(dirname $Assembly)
cat kraken2.1.3/contaminantcontignames.txt > contaminantcontignames.txt
cat tiara/contaminantcontignames.txt >> contaminantcontignames.txt
cat blast2.12.0/contaminantcontignames.txt >> contaminantcontignames.txt
cat blobtoolkit4.2.1/contaminantcontignames.txt >> contaminantcontignames.txt
mkdir filtered
awk 'NR==FNR {busco[$0]=1; next} !($0 in busco)' "$(dirname $Assembly)/BUSCO/buscocontignames.txt" "$(dirname $Assembly)/contaminantcontignames.txt" | sort | uniq > "$(dirname $Assembly)/filtered/contaminantcontignames.txt"
cat $(dirname $Assembly)/filtered/contaminantcontignames.txt | wc -l #4020

seqtk subseq $Assembly $(dirname $Assembly)/filtered/contaminantcontignames.txt | gzip > $(dirname $Assembly)/contaminants.fa.gz
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/filter.py $Assembly $(dirname $Assembly)/filtered/contaminantcontignames.txt > $(dirname $Assembly)/filtered/filtered_$(basename $Assembly)
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $(dirname $Assembly)/filtered/filtered_$(basename $Assembly) > $(dirname $Assembly)/filtered/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_filtered.fa
rm $(dirname $Assembly)/filtered/filtered_$(basename $Assembly)

#BUSCO
for Genome in $(ls $(dirname $Assembly)/filtered/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_filtered.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 
```
T_anthrisci_820m_48_1_10.0_0.25_filtered_hemiptera_odb10_short_summary
  C:92.3%[S:72.6%,D:19.7%],F:3.5%,M:4.2%,n:2510

### Purge Dups
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
#Purge with Tellseq (Illumina) Reads:
MappingFile=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57347690

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4384    4384    766     6088    104935  230299  446600  2788957 595.6e6 T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa

#Purge with HiFi reads:
MappingFile=$(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
Type=long
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57344889

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4372    4372    808     8481    105996  228568  430982  1522448 614.8e6 T_anthrisci_820m_48_1_10.0_0.25_HiFiPurged.fa
```
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_*Purged.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_*/hifiasm_19.5/*/*/*/*/*/purge_dups/BUSCO/*_short_summary.txt); do
echo $(basename $file) >> temp.txt
grep 'C:' $file >> temp.txt
done
```
T_anthrisci_820m_48_1_10.0_0.25_filtered_hemiptera_odb10_short_summary
  C:92.3%[S:72.6%,D:19.7%],F:3.5%,M:4.2%,n:2510

T_anthrisci_820m_48_1_10.0_0.25_HiFiPurged_arthropoda_odb10_short_summary.txt
        C:91.2%[S:85.0%,D:6.2%],F:4.5%,M:4.3%,n:1013
T_anthrisci_820m_48_1_10.0_0.25_HiFiPurged_hemiptera_odb10_short_summary.txt
        C:91.7%[S:85.3%,D:6.4%],F:3.8%,M:4.5%,n:2510
T_anthrisci_820m_48_1_10.0_0.25_HiFiPurged_insecta_odb10_short_summary.txt
        C:91.1%[S:85.2%,D:5.9%],F:4.7%,M:4.2%,n:1367
        
T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_arthropoda_odb10_short_summary.txt
        C:91.0%[S:89.1%,D:1.9%],F:4.7%,M:4.3%,n:1013
T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_hemiptera_odb10_short_summary.txt
        C:91.6%[S:89.3%,D:2.3%],F:3.7%,M:4.7%,n:2510
T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_insecta_odb10_short_summary.txt
        C:90.7%[S:88.5%,D:2.2%],F:4.8%,M:4.5%,n:1367

Purging with TellSeq reads seems to be more effective than using HiFi reads
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/tellseq
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/tellseq/.
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/hifi
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_HiFiPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/hifi/.

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/*/*Purged.fa); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done #57291752,3
```
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/filtered/T_anthrisci_820m_48_1_10.0_0.25_filtered.fa
#Purge with Tellseq (Illumina) Reads:
MappingFile=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57051904

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4001    4001    725     8068    109745  235314  448219  2788957 571.2e6 T_anthrisci_820m_48_1_10.0_0.25_filtered_TellSeqPurged.fa

#Purge with HiFi reads:
MappingFile=$(dirname $Assembly)/../minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
Type=long
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57051905

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4001    4001    725     8068    109745  235314  448219  2788957 571.2e6 T_anthrisci_820m_48_1_10.0_0.25_filtered_HiFiPurged.fa
```
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/filtered/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_*Purged.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    echo $OutFile >> logs/buscolog.txt
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_*/hifiasm_19.5/*/*/*/*/*/filtered/purge_dups/BUSCO/*_short_summary.txt); do
echo $(basename $file) >> temp.txt
grep 'C:' $file >> temp.txt
done

/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_*/hifiasm_19.5/*/*/*/*/*/filtered/purge_dups/T_*Purged.fa
```
T_anthrisci_820m_48_1_10.0_0.25_filtered_hemiptera_odb10_short_summary
  C:92.3%[S:72.6%,D:19.7%],F:3.5%,M:4.2%,n:2510

T_anthrisci_820m_48_1_10.0_0.25_filtered_HiFiPurged_arthropoda_odb10_short_summary.txt
        C:91.1%[S:89.2%,D:1.9%],F:4.7%,M:4.2%,n:1013
T_anthrisci_820m_48_1_10.0_0.25_filtered_HiFiPurged_insecta_odb10_short_summary.txt
        C:90.8%[S:88.6%,D:2.2%],F:4.8%,M:4.4%,n:1367
T_anthrisci_820m_48_1_10.0_0.25_filtered_TellSeqPurged_arthropoda_odb10_short_summary.txt
        C:91.1%[S:89.2%,D:1.9%],F:4.7%,M:4.2%,n:1013
T_anthrisci_820m_48_1_10.0_0.25_filtered_TellSeqPurged_hemiptera_odb10_short_summary.txt
        C:91.7%[S:89.4%,D:2.3%],F:3.7%,M:4.6%,n:2510
T_anthrisci_820m_48_1_10.0_0.25_filtered_TellSeqPurged_insecta_odb10_short_summary.txt
        C:90.8%[S:88.6%,D:2.2%],F:4.8%,M:4.4%,n:1367

```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/filtered/purge_dups/tellseq
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/filtered/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_filtered_TellSeqPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/filtered/purge_dups/tellseq/.
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/filtered/purge_dups/hifi
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/filtered/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_filtered_HiFiPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/filtered/purge_dups/hifi/.

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/filtered/purge_dups/*/*Purged.fa); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done #57291737,8
```
```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_*/hifiasm_19.5/*/*/*/*/*/purge_dups/BUSCO/*); do
basename $file >> temp.txt
grep 'C:' $file >> temp.txt
done
```
Purging with/without contaminants seems to make no difference to the final scores, neither does using Tellseq/Hifi reads.

### Scaffolding
Tellseq reads have are linked via 18bp barcodes, Tom Mathers has already used the conversion software prodived by Universal sequencing and the 4Mwith-alts-february-2016.txt barcode whitelist file to convert to 16bp barcoded versions of the reads that are compatible with 10x genomics linked read software. The files also have to be in a standardised naming format https://cdn.shopify.com/s/files/1/0654/0378/1341/files/100027-USG_TELL-Seq_Software_Roadmap_User_Guide_v1.0_4d2abbf8-3d19-4899-84c2-9a79594eb7f3.pdf?v=1684422029

Scaff10X can now be used to scaffold the HiFi assembly contigs using the TellSeq reads.

Longranger is run in order to remove the barcodes from the reads so that they can simly be treated as large insert illumina reads and used for assembly polishing.

Purging should be done before scaffolding as presence of haplotigs will confuse scaffolder, however the order of pilon,scaff10x,break10x,YAHS or whether using at all will improve the assembly is unclear.

#### 3DDNA
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa
OutDir=$(dirname $Assembly)/3ddna
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_3dDNA.sh $Assembly $OutDir $OutFile $Read1 $Read2
#57364468, 57615781, 57625368
#NOTE: 3ddna output is very large ~600GB therefore only final files kept

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1612    1612    6       977     27.51e6 37.42e6 46.85e6 58.62e6 554.4e6 T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.FINAL.fasta
```

#### 0
Scaffolding without purging

##### Scaff10X
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234151

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#9147    9147    1071    6088    60145   214611  457160  2788957 829.7e6 output_scaffolds.fasta
```
##### YAHS
Generate mapping file
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57245207

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57278990

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#8819    8819    773     1000    63856   270247  700275  2914693 829.7e6 T_anthrisci_820m_48_1_10.0_0.25_scaffolds_final.fa
```
#### 1
Purge
```bash
#Already performed above
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234129

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4353    4353    762     6088    106632  232287  447197  2788957 595.6e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57288677

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57297212

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1304    1304    5       1000    20.61e6 52.29e6 62.67e6 77e6    595.6e6 T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57245203

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57279053, 57302560

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1305    1305    3       1000    20.16e6 74.73e6 120.4e6 120.4e6 595.6e6 T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57285519

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1297    1297    3       1000    20.16e6 74.73e6 120.4e6 120.4e6 595.6e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57298234

Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds_mapped.PT.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.PT.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57301306, 57303597

#No dedup:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap2.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57301397

Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@.bam@@')_nodedup
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57333204

#Sanger
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')_222
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap2.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57561855
```
Remove 13 and 156:
```python
from Bio import SeqIO

# Define the input and output filenames
input_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta"
output_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds_edit.fasta"

# Define the list of contigs to remove
contigs_to_remove = ["tarseq_13", "tarseq_156"]

# Open the input and output files
with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
    # Iterate over the records in the input FASTA file
    for record in SeqIO.parse(input_handle, "fasta"):
        # Check if the contig is not in the list of contigs to remove
        if record.id not in contigs_to_remove:
            # Write the record to the output file
            SeqIO.write(record, output_handle, "fasta")
```
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds_edit.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57308023

Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final_scaff10xscaffolds_edit_mapped.PT.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57333147
```
#### 2
Purge
```bash
#Already performed above
```
Pilon
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa
Alignment=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@_TellSeqPurged.fa@@g')_Tellseq_trimmed.bam
source package 638df626-d658-40aa-80e5-14a275b7464b
samtools sort -o $(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_sorted_Tellseq_trimmed.bam@g') $Alignment
rm $Alignment

source switch-institute ei
source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
source package /tgac/software/testing/bin/picardtools-2.1.1
chmod 777 $(dirname $Alignment)
java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar MarkDuplicates I=$(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_sorted_Tellseq_trimmed.bam@g') O=$(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_markdups_Tellseq_trimmed.bam@g') M=$(dirname $Assembly)/marked_dup_metrics.txt
rm $(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_sorted_Tellseq_trimmed.bam@g')
cd $(dirname $Alignment)
samtools index $(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_markdups_Tellseq_trimmed.bam@g')

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa
Alignment=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@_TellSeqPurged.fa@@g')_sorted_markdups_Tellseq_trimmed.bam
Index=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@_TellSeqPurged.fa@@g')_sorted_markdups_Tellseq_trimmed.bam.bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57234146
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa_pilon.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa_pilon.fasta@_pilon@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234831

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4325    4325    752     6088    107835  235400  452506  2788759 595.3e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/pilon/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_pilon_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa_pilon.fasta@_pilon@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57297142

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/pilon/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_pilon_scaff10xscaffolds.fasta
Alignment=
Alignment_Index=
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa_pilon.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa_pilon.fasta@_pilon@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57245209

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa_pilon.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_pilon_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_pilon_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57279782

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1278    1278    4       1000    21.37e6 56.26e6 140e6   140e6   595.3e6 T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa_pilon_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/pilon/yahs/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa_pilon_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa_pilon@_pilon@' | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57297220
```
#### 3
Break10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/break10x
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_break10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234077

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#9790    9790    1366    6088    59062   173895  353115  2788957 829.7e6 T_anthrisci_820m_48_1_10.0_0.25_break.fa

#BUSCO
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/T_anthrisci_820m_48_1_10.0_0.25_break.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done #
```
Purge
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/T_anthrisci_820m_48_1_10.0_0.25_break.fa
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger/Tant_T508_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1
#57245190

#Purge with Tellseq (Illumina) Reads:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/T_anthrisci_820m_48_1_10.0_0.25_break.fa
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/bwa/T_anthrisci_820m_48_1_10.0_0.25_break.fa_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57252456

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4531    4531    824     6088    100360  215484  406833  2788957 593.8e6 T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa

sbatch ~/git_repos/Pipelines/Trioza_merqury.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa  
#57291755

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done #57381425

#T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_hemiptera_odb10_short_summary.txt
#        C:91.7%[S:89.3%,D:2.4%],F:3.7%,M:4.6%,n:2510
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57258652

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4499    4499    817     6088    101965  218217  410989  2788957 593.8e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57286927

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57291004

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1323    1323    6       1000    12.37e6 34.83e6 68.88e6 88.56e6 593.8e6 T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/scaff10x/yahs/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57295617
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57276685, 57276779, 57276823, 57279785, 57279829

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57285414

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1299    1299    6       1000    20.08e6 36.36e6 64.41e6 92.34e6 593.8e6 T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/yahs/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57291010

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1292    1292    6       1000    20.08e6 36.36e6 64.41e6 92.65e6 593.8e6 output_scaffolds.fasta
```
#### 3.1
Final curated version from T.Mathers
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa.gz
```
#### 4
Break10x
```bash
#Already performed above
```
Purge
```bash
#Already performed above
```
Pilon
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa
Alignment=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@_TellSeqPurged.fa@@g')_Tellseq_trimmed.bam
source package 638df626-d658-40aa-80e5-14a275b7464b
samtools sort -o $(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_sorted_Tellseq_trimmed.bam@g') $Alignment
rm $Alignment

source switch-institute ei
source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
source package /tgac/software/testing/bin/picardtools-2.1.1
chmod 777 $(dirname $Alignment)
java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar MarkDuplicates I=$(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_sorted_Tellseq_trimmed.bam@g') O=$(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_markdups_Tellseq_trimmed.bam@g') M=$(dirname $Assembly)/marked_dup_metrics.txt
rm $(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_sorted_Tellseq_trimmed.bam@g')
cd $(dirname $Alignment)
samtools index $(echo $Alignment | sed 's@_Tellseq_trimmed.bam@_markdups_Tellseq_trimmed.bam@g')

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/bwa/T_anthrisci_820m_48_1_10.0_0.25_break_sorted_markdups_Tellseq_trimmed.bam
Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/bwa/T_anthrisci_820m_48_1_10.0_0.25_break_sorted_markdups_Tellseq_trimmed.bam.bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57271544

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon.fasta); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done #57381431
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57276691

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4471    4471    811     6088    103508  220397  414291  2788759 593.5e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57285462

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57290955

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1273    1273    5       1000    23.1e6  42.79e6 76.03e6 98.71e6 593.5e6 T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/scaff10x/yahs/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291027

Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/scaff10x/yahs/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_mapped.PT.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.PT.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57297211 - 57306167
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57279828

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57285415

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1289    1289    3       1000    17.79e6 63.65e6 177.7e6 177.7e6 593.5e6 T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/yahs/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57285479

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1282    1282    3       1000    17.79e6 63.65e6 177.7e6 177.7e6 593.5e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/yahs/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57297158
```
#### 5
Break10x
```bash
#Already performed above
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/T_anthrisci_820m_48_1_10.0_0.25_break.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57252795

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#9337    9337    1149    6088    59062   202992  420426  2788957 829.7e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291765

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57297196
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/T_anthrisci_820m_48_1_10.0_0.25_break.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57252799

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/T_anthrisci_820m_48_1_10.0_0.25_break.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/T_anthrisci_820m_48_1_10.0_0.25_break_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/T_anthrisci_820m_48_1_10.0_0.25_break_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57279792

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#8806    8806    762     1000    64026   270007  707093  3120170 829.7e6 T_anthrisci_820m_48_1_10.0_0.25_break_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/yahs/T_anthrisci_820m_48_1_10.0_0.25_break_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57297221
```
#### 6
Pilon
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/T_anthrisci_820m_48_1_10.0_0.25.bp.p_ctg.fa
Alignment=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_sorted_markdups_Tellseq_trimmed.bam
Index=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_sorted_markdups_Tellseq_trimmed.bam.bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57227862
```
Break10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/T_anthrisci_820m_48_1_10.0_0.25.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/break10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_break10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234076

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#9789    9789    1367    6088    59062   173762  353037  2788759 829.3e6 T_anthrisci_820m_48_1_10.0_0.25_break.fa
```
Purge
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/T_anthrisci_820m_48_1_10.0_0.25_break.fa
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger/Tant_T508_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1
#57245192

#Purge with Tellseq (Illumina) Reads:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/T_anthrisci_820m_48_1_10.0_0.25_break.fa
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/bwa/T_anthrisci_820m_48_1_10.0_0.25_break_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57252544

sbatch ~/git_repos/Pipelines/Trioza_merqury.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa
#57291761
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57276686

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4502    4502    818     6088    101605  218822  407978  2788759 593.8e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291804

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/scaff10x/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57297199
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@_break_TellSeqPurged.fa@_pilon_break_TellSeqPurged@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57285416

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_pilon_break_TellSeqPurged_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_pilon_break_TellSeqPurged_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57285416

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1297    1297    4       1000    25.6e6  62.38e6 81.19e6 99.91e6 593.8e6 T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/pilon/break10x/purge_dups/scaff10x/yahs/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57297228
```
### Filtering
#### MitoHiFi
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa
ReferenceFasta=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/NCBI/Trant_mito_NC_038141.1.fasta
ReferenceGenebank=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/NCBI/Trant_mito_NC_038141.1.gb
PercentOverlap=50
Code=5
Kingdom=animal
Good_Reference=Y
OutDir=$(dirname $Assembly)/MitoHifi
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_mitohifi_contigs.sh $Assembly $ReferenceFasta $ReferenceGenebank $PercentOverlap $Code $Kingdom $Good_Reference $OutDir 
#58068711

OutFile=T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated
```
scaffold_397 identified as the mitochondrial genome.
```bash
echo scaffold_397 > id.txt
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_rm.py  --id_file id.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito.fa

#Mitochondrial genome:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/final_mitogenome.fasta
```

#### Filter contaminants from curated assembly

Investigate the large contaminant contig/scaffold 17.
```bash
#Extract scaffold 17:
echo scaffold_17_3 > temp_search.txt
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/git_repos/Scripts/NBI/seq_get.py --id_file temp_search.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa --output ant_17.fa

#Align scaffold 17 to the liberibacter genome from the apicales assembly, and the reference published liberibacter genome:
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
Query=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/ant_17.fa
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae
OutFile=scaff17
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58078409
source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/scaff17.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/scaff17.delta

Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna
Query=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/ant_17.fa
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae
OutFile=scaff17-2
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58089601
source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/scaff17-2.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/scaff17-2.delta
#There was no alignment to liberibacter


#Blast scaffold 17 to identify its origin:
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/ant_17.fa
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_blastn
Max=9999999
OutFile=ant_17
sbatch ~/git_repos/Wrappers/NBI/run_blastn.sh $InFile $Database $OutDir $OutFile $Max
#58092051
awk '{print $5}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_blastn/ant_17.vs.nt.mts1.hsp1.1e25.megablast.out | sort | uniq -c | sort -nr | head -n 1
#blast hits suggest scaffodl 17 is a Staphylococcus xylosus genome

#Align scaffold 17 to the reference Staphylococcus xylosus genome:
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/ant_17.fa
Query=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Staphylococcus/xylosus/GCF_000709415.1/GCA_000709415.1_ASM70941v1_genomic.fna
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae
OutFile=scaff17-3
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58098311
source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/scaff17-3.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/scaff17-3.delta
#Alignment is very strong
```
Scaffold 17 looks like it is Staphylococcus xylosus.

#### Kraken and Blobtools
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito.fa
Species=$(echo $Assembly | cut -d '/' -f10)
Genus=Trioza
TaxID=2023874
alias=$(basename $Assembly | cut -d '_' -f1,2 | sed 's@_@@g' | cut -c1-4)
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger/Tant_T508_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI

#Kraken 
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')_kraken2nt
OutDir=$(dirname $Assembly)/kraken2.1.3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/kraken/nt_14092023
mkdir $OutDir
sbatch $ProgDir/run_kraken2.sh $Assembly $Database $OutDir $OutPrefix
#58072788

#Tiara
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')
OutDir=$(dirname $Assembly)/tiara
mkdir $OutDir
sbatch $ProgDir/run_tiara.sh $Assembly $OutDir $OutPrefix
#58072789

#Alignment
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1 
#58072790,58089920

#Blast
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')
OutDir=$(dirname $Assembly)/blast2.12.0
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
mkdir $OutDir
sbatch $ProgDir/run_blastn.sh $Assembly $Database $OutDir $OutPrefix 
#58072804

#BUSCO - keeping output files
OutDir=$(dirname $Assembly)/BUSCO
Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
OutFile=$(basename $Assembly | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
mkdir $OutDir 
sbatch $ProgDir/run_busco_keep.sh $Assembly $Database $OutDir $OutFile 
#58098588

#Blobtoolkit
Record_type=contig
MappingFile=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed.bam
BlastFile=$(dirname $Assembly)/blast2.12.0/$(basename $Assembly | sed 's@.fa@@g').vs.nt.mts1.hsp1.1e25.megablast.out
BUSCOFile=$(dirname $Assembly)/BUSCO/hemiptera_odb10/run_hemiptera_odb10/full_table.tsv
Tiara=$(ls $(dirname $Assembly)/tiara/${Species}*.tiara)
OutDir=$(dirname $Assembly)/blobtoolkit4.2.1
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')
Alias=$(echo $alias)_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/run_btk4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $Alias
#58102771

cp -r ${OutDir}/Tant_Tellseq_trimmed_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
#The contents of /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs were subsequently copied to local machine WSL-Ubuntu: \\wsl.localhost\Ubuntu\home\did23faz
```
```bash
ubuntu
conda activate btk
#From \\wsl.localhost\Ubuntu\home\did23faz
apptainer exec blobtoolkit.sif blobtools host BlobDirs
```
Kraken and blobtools were used to screen the assembly for contaminants resulting in the removal of 122 scaffolds, all scaffolds without kraken classification to an arthropoda taxa outside of the range 0.3 - 0.4 GC, 40 - 70x coverage were removed, within this range scaffolds without classification to arthropoda taxa were kept if there classification was more general (eg. eukaryota), unclassified, or to a taxa it is implausible the sample to be contaminatied with (eg. tuna). Additionally, 24 scaffolds scaffolds were removed that were classified to arthropoda taxa by kraken, but which were identified by either tiara or blast as potential contaminants, these all fell outside of the 0.3 - 0.4 GC, 40 - 70x coverage range.

Blast and kraken both identify contamination from Rickettsia bacteria at high coverage, these are likely symbionts, as well as Ignavibacterium at low coverage. They also confirm Staphylococcus xylosus contamination.

Suspected contaminants were removed to a seperate fasta file:
```bash
nano /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/suspected_contaminant_names.txt

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/suspected_contaminant_names.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito.fa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_suspected_contaminants.fa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_rm.py --id_file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/suspected_contaminant_names.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito.fa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered.fa
```

##### Carsonella

```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_HiFiPurged_curated_break_scaffolds_final_nomito_filtered_corrected.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_suspected_contaminants.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_suspected_contaminants.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_HiFiPurged_curated_break_scaffolds_final_nomito_suspected_contaminants.fa); do
name=$(basename $file | rev | cut -d '/' -f1 | rev | cut -d '_' -f1,2)
echo $name
cp $file ./x.fa
sed -i "s/>/>${name}_/g" x.fa
cat x.fa >> db3.fa
done

makeblastdb -in db3.fa -input_type fasta -dbtype nucl -title psyllid3  -parse_seqids -out psyllid3

blastn -query ../Symbionts/Candidatus/Carsonella/ruddii/GCA_000287275.1/GCA_000287275.1_ASM28727v1_genomic.fna -db psyllid3 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/Carsonella_ruddii_results -evalue 1e-5 -outfmt 6 -num_threads 1

#T_anthrisci_scaffold_54 tens of thousands of bp alignment
#T_anthrisci_scaffold_858 6701bp alignment
#T_anthrisci_scaffold_75 <1220bp alignment
#T_anthrisci_scaffold_285 <795bp alignment
#T_anthrisci_scaffold_286 <891bp alignment
#T_anthrisci_scaffold_17_3 <590bp alignment
#T_anthrisci_scaffold_17_2 <590bp alignment
#T_anthrisci_scaffold_254 <574bp alignment
#T_anthrisci_scaffold_250 64bp alignment

awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_54' > temp_scaffold_54.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_858' > temp_scaffold_858.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_75' > temp_scaffold_75.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_285' > temp_scaffold_285.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_17_3' > temp_scaffold_17_3.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_254' > temp_scaffold_254.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_17_2' > temp_scaffold_17_2.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_250' > temp_scaffold_250.fasta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_54.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_54
mummerplot -color scaffold_54.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_54.fasta -p scaffold_54-2
mummerplot -color scaffold_54-2.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_858.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_858
mummerplot -color scaffold_858.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_858.fasta -p scaffold_858-2
mummerplot -color scaffold_858-2.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_75.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_75
mummerplot -color scaffold_75.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_75.fasta -p scaffold_75-2
mummerplot -color scaffold_75-2.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_285.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_285
mummerplot -color scaffold_285.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_285.fasta -p scaffold_285-2
mummerplot -color scaffold_285-2.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_17_2.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_17_2
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_17_3.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_17_3
mummerplot -color scaffold_17_3.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_17_3.fasta -p scaffold_17_3-2
mummerplot -color scaffold_17_3-2.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_250.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_250
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_254.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_254
mummerplot -color scaffold_254.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_254.fasta -p scaffold_254-2
mummerplot -color scaffold_254-2.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_858.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_54.fasta -p scaffold_54_858
mummerplot -color scaffold_54_858.delta

cat temp_scaffold_54.fasta temp_scaffold_858.fasta temp_scaffold_75.fasta temp_scaffold_285.fasta temp_scaffold_17_3.fasta temp_scaffold_17_2.fasta temp_scaffold_250.fasta temp_scaffold_254.fasta > temp_all.fasta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_all.fasta -p temp_all_ant
mummerplot -color temp_all_ant.delta

nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna ../Symbionts/Candidatus/Carsonella/ruddii/GCA_000010365.1/GCA_000010365.1_ASM1036v1_genomic.fna -p temp
mummerplot -color temp.delta
```
Alignments show that scaffold 54 consists of ~140,000bp of carsonella sequence, some sections of the genome are covered twice. Scaffold 858 is also ~50% carsonella sequence, all of which is within scaffold 54, the other scaffolds have minimal alignment. 
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/carsonella_names.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/carsonella/carsonella_contigs.fa

ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/carsonella/carsonella_contigs.fa
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAn22_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAn22_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2
#58772573

source package c92263ec-95e5-43eb-a527-8f1496d56f1a
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/carsonella/minimap2
samtools view -h carsonella_contigs.bam -o carsonella_contigs.sam
samtools fastq -@32 carsonella_contigs.sam > reads.fastq
#58776517

ProgDir=~/git_repos/Wrappers/NBI
Reads=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/carsonella/minimap2/reads.fastq
OutDir=$(dirname $Reads)/hifiasm_19.5.2
OutFile=carsonella_ant
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_default.sh $OutDir $OutFile $Reads
#58778305

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#9452    9452    1152    6429    60185   198829  430445  2788957 831.3e6 carsonella_api.bp.p_ctg.fa

source package d6092385-3a81-49d9-b044-8ffb85d0c446
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/carsonella/minimap2/hifiasm_19.5.2/carsonella_api.bp.p_ctg.fa -input_type fasta -dbtype nucl -title carant  -parse_seqids -out carant
blastn -query ../Symbionts/Candidatus/Carsonella/ruddii/GCA_000287275.1/GCA_000287275.1_ASM28727v1_genomic.fna -db carant -out carant -evalue 1e-5 -outfmt 6 -num_threads 1
awk '{print $2}' carant | sort | uniq | wc -l #222
```
#### FCS
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered.fa
TAXID=2023874
OutDir=$(dirname $Genome)/fcs
OutFile=$(basename $Genome | sed 's@.fa@@g')
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_fcs.sh $Genome $TAXID $OutDir $OutFile
#59164346
```

### Final assembly assessment

Fresh blobplots were prepared for the filtered assembly:
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered.fa
Species=$(echo $Assembly | cut -d '/' -f10)
Genus=Trioza
TaxID=872318
alias=$(basename $Assembly | cut -d '_' -f1,2 | sed 's@_@@g' | cut -c1-4)
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger/Tant_T508_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI

#Tiara
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')
OutDir=$(dirname $Assembly)/tiara
mkdir $OutDir
sbatch $ProgDir/run_tiara.sh $Assembly $OutDir $OutPrefix
#58113113

#Alignment
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1 
#58113116

#Blast
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')
OutDir=$(dirname $Assembly)/blast2.12.0
Max_target=10
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
mkdir $OutDir
sbatch $ProgDir/run_blastn.sh $Assembly $Database $OutDir $OutPrefix $Max_target
#58746494

#BUSCO - keeping output files
OutDir=$(dirname $Assembly)/BUSCO
Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
OutFile=$(basename $Assembly | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
mkdir $OutDir 
sbatch $ProgDir/run_busco_keep.sh $Assembly $Database $OutDir $OutFile 
#58113126

#Blobtoolkit
Record_type=contig
MappingFile=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed.bam
BlastFile=$(dirname $Assembly)/blast2.12.0/$(basename $Assembly | sed 's@.fa@@g').vs.nt.mts1.hsp1.1e25.megablast.out
BUSCOFile=$(dirname $Assembly)/BUSCO/hemiptera_odb10/run_hemiptera_odb10/full_table.tsv
Tiara=$(ls $(dirname $Assembly)/tiara/${Species}*.tiara)
OutDir=$(dirname $Assembly)/blobtoolkit4.2.1
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')
Alias=$(echo $alias)_Tellseq_trimmed_filtered
mkdir $OutDir
sbatch $ProgDir/run_btk4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $Alias
#58760372

cp -r ${OutDir}/Tant_Tellseq_trimmed_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
#The contents of /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs were subsequently copied to local machine WSL-Ubuntu: \\wsl.localhost\Ubuntu\home\did23faz
```
#### Inspector

The filtered assembly was assessed and polished with inspector:
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered.fa
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)
Datatype=hifi
Correct_Datatype=pacbio-hifi
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAn22_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAn22_hifi_3rdSMRTcell.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_inspector.sh $OutFile $OutDir $Genome $Datatype $Correct_Datatype $Read1 $Read2
#58153902
```
Inspector reported the following summary statistics for the assembly:

Statics of contigs:
Number of contigs       361
Number of contigs > 1000 bp     361
Number of contigs >1000000 bp   14
Total length    587809615
Total length of contigs > 1000 bp       587809615
Total length of contigs >1000000bp      571400068
Longest contig  67459246
Second longest contig length    66647579
N50     48836234
N50 of contigs >1Mbp    48836234


Read to Contig alignment:
Mapping rate /% 79.6
Split-read rate /%      20.19
Depth   31.8739
Mapping rate in large contigs /%        78.07
Split-read rate in large contigs /%     20.22
Depth in large conigs   32.1703


Structural error        462
Expansion       301
Collapse        64
Haplotype switch        79
Inversion       18


Small-scale assembly error /per Mbp     98.53360428614288
Total small-scale assembly error        57919
Base substitution       45288
Small-scale expansion   6959
Small-scale collapse    5672

QV      29.932235272211926

#### Juicebox
The final assembly was visualised using juicebox:
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
OutDir=$(dirname $Assembly)/juicebox
OutFile=$(basename $Assembly | sed 's@.fa@@g')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_juicebox-prep.sh $Assembly $Read1 $Read2 $OutDir $OutFile
#58971621

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
OutDir=$(dirname $Assembly)/3ddna
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_3dDNA.sh $Assembly $OutDir $OutFile $Read1 $Read2
#58971633

source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
java17 -Xmx48000m -Djava.awt.headless=true -jar ~/git_repos/Scripts/NBI/juicer_tools_1.22.01.jar pre --threads 32 mapped.pairs hic.hic genome.genome


#From omni-c guide
source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
java17 -Xmx48000m -Djava.awt.headless=true -jar ~/git_repos/Scripts/NBI/juicer_tools_1.22.01.jar pre --threads 16 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicebox/*.pairs /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicebox/*.hic /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicebox/*.genome
sbatch ~/git_repos/Wrappers/NBI/temp3.sh
#57554716, 57816287, 57818847
 ~/git_repos/Scripts/NBI/generate-assembly-file-from-fasta.awk /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta

source lastz-1.03.73;source gnu_parallel-20180322;source jdk-1.7.0_25
~/3ddna-master/visualize/run-assembly-visualizer.sh draft.assembly aligned/merged_nodups.txt

#From T.Mathers
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicer
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicer
mkdir restriction_sites;mkdir references;mkdir fastq
cd fastq
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz .
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz .
cd ../references
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa .
source bwa-0.7.17
bwa index T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
cd ../restriction_sites
Enzyme=DpnII
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/generate_site_positions.py $Enzyme OutFile ../references/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' OutFile_DpnII.txt > OutFile_DpnII.chrom.sizes
cd ..
mkdir -p scripts/common
cp ~/git_repos/Scripts/NBI/chimeric_blacklist.awk scripts/common/.
cp ~/git_repos/Scripts/NBI/countligations.sh scripts/common/. 
srun -p jic-long  -c 32 --mem 50G --pty bash
source jdk-1.7.0_25; source bwa-0.7.17;source samtools-1.6;source gnu_parallel-20180322
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/juicer.sif juicer.sh -D /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicer -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicer -g OutFile -z references/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa -y restriction_sites/OutFile_DpnII.txt -s DpnII -t 32 -p restriction_sites/OutFile_DpnII.chrom.sizes 
#58959622 pwd=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicer/references

source lastz-1.03.73;source gnu_parallel-20180322;source jdk-1.7.0_25


source lastz-1.03.73;source gnu_parallel-20180322;source jdk-1.7.0_25

awk -f ~/git_repos/Scripts/NBI/generate-assembly-file-from-fasta.awk 
~/git_repos/Scripts/NBI/run-assembly-visualizer.sh


~/3ddna-master/visualize/run-assembly-visualizer.sh draft.assembly aligned/merged_nodups.txt


bwa mem -SP5M

submit-slurm_v1.1.pl -q ei-medium -m 100000 -c 2 -t 2-
00:00 -e -j Smis_HiC_file -i "source lastz-1.03.73;source gnu_parallel-20180322;source jdk-1.7.0_25;~/3ddna-
master/visualize/run-assembly-visualizer.sh draft.assembly aligned/merged_nodups.txt"

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/juicer.sif juicer.sh -t 16 -d $WorkDir/juicer -g saundersiae -z genome_wrapped.fa -y ${OutFile}_Phase.txt -p genome_wrapped.fa.fai -D /opt/juicer-1.6.2/CPU


ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa ###
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_HiC
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#59554224

awk -f ~/git_repos/Scripts/NBI/generate-assembly-file-from-fasta.awk /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa > tant.assembly
srun -p jic-long  -c 32 --mem 128G --pty bash
source lastz-1.03.73;source gnu_parallel-20180322;source jdk-1.7.0_25
~/git_repos/Scripts/NBI/run-assembly-visualizer.sh tant.assembly /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_59638370/juicer/aligned/merged_nodups.txt
#NOTE: this destroys the input merged_nodups.txt file

mv tant* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/3ddna-vis/.
```
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
OutDir=$(dirname $Assembly)/3ddna-vis
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiC/anthrisci_286154-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_3dDNA.sh $Assembly $OutDir $OutFile $Read1 $Read2
#59638370,59644856

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_59644856/juicer/aligned/merged_nodups.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/3ddna-vis/*
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_59644856/juicer/aligned/merged_nodups.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/3ddna-vis/.
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_59638370/juicer/aligned/genome_wrapped.rawchrom* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/3ddna-vis/.
```
## Annotation
### Repeatmasking

#### Repeatmodeler
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/repeatmodeler
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmodeler.sh $Genome $OutFile $OutDir
#58170701
```
#### Repeatmasker
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
Species=hemiptera
Repeat_library=$(ls $(dirname $Genome)/repeatmodeler/*/consensi.fa.classified)
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/repeatmasker
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmasker4.1.5.sh $Genome $OutFile $OutDir $Species $Repeat_library
#58186045
```
#### EarlGreyTE
```bash
source package 14fbfadb-9fe7-419a-9f20-cd5f458c0fff

source package /tgac/software/testing/bin/transposonPSI-08222010

bedtools maskfasta -fi [genome.fasta] -bed [teannotations.bed] -fo [genome.softmasked.fa] -soft

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/earlgrey4.0.6.sif earlGrey 

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/earlgrey
RMsearch=arthropoda
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_earlgrey.sh $Genome $OutFile $OutDir $RMsearch
#58797251, 58852735 ran out of memory with 184GB, 58929433 time out at 2 days, 58940974, 59071710

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/earlgrey4.0.6.sif earlGrey -g $Genome -s $OutFile -o $OutDir -t 1 -d yes -r arthropoda
```

### RNASeq 
#### Trimmomatic
```bash
for Rawdata in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_*.fq.gz); do
OutDir=$(dirname $Rawdata)/fastqc
OutFile=$(basename $Rawdata | sed 's@.fq.gz@@g')
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_fastqc.sh $Rawdata $OutDir $OutFile
done

mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_3.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_4.fq.gz
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f3 | rev)
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    Fread2=$(ls ${ReadDir}*_3.fq.gz)
    Rread2=$(ls ${ReadDir}*_4.fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@dna_qc@g')trim_galore
    OutFile=${sample}_trimmed
    Quality=20
    Length=50
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3 $Fread4 $Rread4 $Fread5 $Rread5
done 
#58156405

for QCdata in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/RNASeq/trim_galore/T_anthrisci_*.fq.gz); do
OutDir=$(dirname $QCdata)/fastqc
OutFile=$(basename $QCdata | sed 's@.fq.gz@@g')
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_fastqc.sh $QCdata $OutDir $OutFile
done
#58171201,2
```
#### HiSat2
Hisat2 is the aligner used internally by braker, it is slower than star but requires less memory.
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/RNASeq/trim_galore/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    InGenome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    OutDir=$(dirname $InGenome)/hisat2
    OutFile=$(basename $InGenome | sed 's@.fa@@g')
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_HiSat2.sh $InGenome $Fread $Rread $OutDir $OutFile
done
#58233931
```
100164402 reads; of these:
  100164402 (100.00%) were paired; of these:
    37485981 (37.42%) aligned concordantly 0 times
    59758386 (59.66%) aligned concordantly exactly 1 time
    2920035 (2.92%) aligned concordantly >1 times
    ----
    37485981 pairs aligned concordantly 0 times; of these:
      197173 (0.53%) aligned discordantly 1 time
    ----
    37288808 pairs aligned 0 times concordantly or discordantly; of these:
      74577616 mates make up the pairs; of these:
        57021940 (76.46%) aligned 0 times
        16951489 (22.73%) aligned exactly 1 time
        604187 (0.81%) aligned >1 times
71.54% overall alignment rate

#### STAR
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/RNASeq/trim_galore/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    InGenome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    OutDir=$(dirname $InGenome)/star
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_star_align.sh $InGenome $Fread $Rread $OutDir 
done
#59434622
```
#### Braker
```bash
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda.fa > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
for file in $(ls ../Genomes/*/*/*/protein.faa | grep -v 'Frankliniella\|Thrips\|Megalurothrips\|Tribolium\|Drosophila\|Bombyx\|Apis\|Anopheles'); do
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $file | sed 's@*@@g'| sed 's@|@@g'| sed 's@ @@g'>> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
done

mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa.masked /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa


Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
Species=D_anthrisci
OutDir=$(dirname $Masked_Genome)/braker3
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker3.0.7.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#58236703

samtools sort -n /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
#58379020
Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=NA
Species=D_anthrisci_rna3
OutDir=$(dirname $Masked_Genome)/braker3-rna
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker3.0.7.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#58527346

grep '>' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3-soft/braker.aa | wc -l
#13980
grep '>' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.aa | wc -l
#13765

Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
Species=D_anthrisci
OutDir=$(dirname $Masked_Genome)/braker2
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker2.1.6.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59430194

Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/star/star_aligmentAligned.sortedByCoord.out.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
Species=D_anthrisci-rna
OutDir=$(dirname $Masked_Genome)/braker1
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker2.1.6-rna.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59438025


conda activate braker
Assembly=/home/theaven/scratch/uncompressed/hogenhout/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
OutDir=/home/theaven/scratch/uncompressed/hogenhout/T_anthrisci_braker1
AcceptedHits=/home/theaven/scratch/uncompressed/hogenhout/T_anthrisci_star.bam
GeneModelName=tant
ProgDir=/home/theaven/scratch/apps/braker
sbatch $ProgDir/braker1.sh $Assembly $OutDir $GeneModelName $AcceptedHits 
conda deactivate
#19524348

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker1/braker.gff3
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff3@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein

grep '>' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker1/braker.faa | wc -l
#26,406

conda activate braker
Assembly=/home/theaven/scratch/uncompressed/hogenhout/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
OutDir=/home/theaven/scratch/uncompressed/hogenhout/T_anthrisci_braker2
AcceptedHits=/home/theaven/scratch/uncompressed/hogenhout/Arthropoda+hemiptera.fa
GeneModelName=tant2
ProgDir=/home/theaven/scratch/apps/braker
sbatch $ProgDir/braker2.sh $Assembly $OutDir $GeneModelName $AcceptedHits 
conda deactivate
#19521886

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker2/braker.gff3
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff3@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein

grep '>' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker2/braker.faa | wc -l
#27,975

source package c1ac247d-c4d1-4747-9817-bf03617f979b
Dir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask
kallisto index -i ${Dir}/braker3/transcriptome_index.idx ${Dir}/braker3/braker.codingseq
kallisto quant -t 8 -i ${Dir}/braker3/transcriptome_index.idx -o ${Dir}/braker3/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N11_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N11_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_3.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_4.fq.gz
#Mapping Rate = (Mapped Reads / Total Reads) * 100

kallisto index -i ${Dir}/T_anthrisci_braker1/transcriptome_index.idx ${Dir}/T_anthrisci_braker1/augustus.hints.codingseq
kallisto quant -t 8 -i ${Dir}/T_anthrisci_braker1/transcriptome_index.idx -o ${Dir}/T_anthrisci_braker1/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N11_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N11_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_3.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_4.fq.gz

kallisto index -i ${Dir}/T_anthrisci_braker2/transcriptome_index.idx ${Dir}/T_anthrisci_braker2/augustus.hints.codingseq
kallisto quant -t 8 -i ${Dir}/T_anthrisci_braker2/transcriptome_index.idx -o ${Dir}/T_anthrisci_braker2/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N11_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N11_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_3.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/RNASeq/T_anthrisci_N7_4.fq.gz
#59554876, 59609885

cat ${Dir}/T_anthrisci_braker1/run_info.json
cat ${Dir}/T_anthrisci_braker2/run_info.json
cat ${Dir}/braker3/run_info.json

#Braker1: 38.6% pseudoaligned
#Braker2: 38.3% pseudoaligned
#Braker3: 35.9% pseudoaligned

seqtk seq -U /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.faa > temp.faa && mv temp.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.faa 
seqtk seq -U /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker1/braker.faa > temp.faa && mv temp.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker1/braker.faa
seqtk seq -U /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/T_urticae_braker1/braker.faa > temp.faa && mv temp.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/T_urticae_braker1/braker.faa

awk '{print $1}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.faa > temp.faa && mv temp.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.faa 
awk '{print $1}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker1/braker.faa > temp.faa && mv temp.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker1/braker.faa
awk '{print $1}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/T_urticae_braker1/braker.faa > temp.faa && mv temp.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/T_urticae_braker1/braker.faa

awk '/^>/ { print $0 "_t_apicales" } !/^>/ { print $0 }' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.faa > temp.faa && mv temp.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.faa 
awk '/^>/ { print $0 "_t_anthrisci" } !/^>/ { print $0 }' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker1/braker.faa > temp.faa && mv temp.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker1/braker.faa
awk '/^>/ { print $0 "_t_urticae" } !/^>/ { print $0 }' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/T_urticae_braker1/braker.faa > temp.faa && mv temp.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/T_urticae_braker1/braker.faa

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_braker1/braker.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/T_urticae_braker1/braker.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+braker1.fa
seqkit rmdup /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+braker1.fa -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+braker1.faa
seqtk seq -U /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+braker1.faa > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+braker1.fa

Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+braker1.fa
Species=D_anthrisci
OutDir=$(dirname $Masked_Genome)/braker3+1
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker3.0.7.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59727992,59729476,59729651
```
#### Helixer
```bash
singularity exec ~/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --model-filepath ../databases/helixer/invertebrate_v0.3_a_0500/invertebrate_v0.3_a_0500.h5 --subsequence-length 213840 --overlap-offset 106920 --overlap-core-length 160380 --fasta-path /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa  \
  --species Dyspera_anthrisci --gff-output-path Dyspera_anthrisci_helixer_invertebrate_05.gff3

fasta=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
model_filepath=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/helixer/invertebrate_v0.3_m_0100/invertebrate_v0.3_m_0100.h5
lineage=invertebrate
species=T_anthrisci
outfile=T_anthrisci
outdir=$(dirname $fasta)/helixer
ProgDir=~/git_repos/Wrappers/NBI
mkdir $outdir
sbatch $ProgDir/run_helixer.sh $fasta $model_filepath $lineage $species $outfile $outdir
#58575598

grep 'gene' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.gff | wc -l #17,077
```
#### Swissprot
```bash
Proteome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.aa  
OutDir=$(dirname $Proteome)/swissprot
SwissDbDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/Uniprot/swissprot_2024_March_10
SwissDbName=uniprot_sprot
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName 
#58956715

Proteome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.aa  
OutDir=$(dirname $Proteome)/swissprot
SwissDbDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/Uniprot/swissprot_2024_March_10
SwissDbName=uniprot_sprot
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName 
#59328939
```
#### Interproscan
```bash
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.aa
OutDir=$(dirname $InFile)/interproscan
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_interproscan.sh $InFile $OutDir
#58949316

InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.aa  
OutDir=$(dirname $InFile)/interproscan
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_interproscan.sh $InFile $OutDir
#59328996
```
#### Pretextgraph
```bash
#NOTE: mapped scaffolds were produced by T.mathers from /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged.fa, therefore we do not have the scaffolds themselves
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/scaff10x/yahs/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_anthrisci/TellSeq/longranger/Tant_T508_barcoded.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $Read1
#57822383, 57831161, 

ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/pilon/scaff10x/yahs/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAn22_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAn22_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2
#57822435, 57831162, 

source package 6daf0c37-1c5e-4cd6-9884-2d0b4d5f9d8f
bamCoverage -b reads.bam -o coverage.bw --outFileFormat bedgraph

source /jic/software/staging/RCSUPPORT-2245/stagingloader
PretextGraph /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/anthrisci2.pretext

zcat bedgraph.file.gz | PretextGraph -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/anthrisci2.pretext 
bigWigToBedGraph bigwig.file /dev/stdout | PretextGraph -i input.pretext -n "graph name"
```
```bash
echo scaffold_17_3 > temp_search.txt
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/git_repos/Scripts/NBI/seq_get.py --id_file temp_search.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/anthrisci2.curated_primary.no_mt.unscrubbed.fa --output ant_17.fa

Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/ant_17.fa
Query=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Staphylococcus/xylosus/GCF_000709415.1/GCA_000709415.1_ASM70941v1_genomic.fna
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae
OutFile=scaff17-3
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58098311
source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/scaff17-3.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/scaff17-3.delta

#The contaminant contig scaffold_17 is Staphylococcus xylosus, almost the full genome aligns.
```
#### Curation
Tom Mathers says anthrisci2.pretext.savestate_4 is best














































```bash
source package /nbi/software/testing/bin/seqtk-1.2
seqtk seq -a "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/yahs/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged_scaffolds_final.fa" | awk '/^>/ {if (seq) print length(seq); seq=""; print $0; next} { seq = seq $0 } END { if (seq) print length(seq) }'
```

#### Canu
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/anthrisci_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/anthrisci_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@canu@g')/750m
    OutFile=T_anthrisci_750m
    Genomesize=750m
    DataType=pacbio-hifi
    mkdir -p $OutDir
    sbatch $ProgDir/run_canu.sh $OutDir $OutFile $Genomesize $DataType $Run1 $Run2
done #57221669
``` 