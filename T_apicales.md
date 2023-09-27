# Trioza apicales
Contains run by T.Heaven in assembly of Trioza apicales

Unless stated otherwise commands were performed from the directory /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae

## Collect data
```bash
#HiFi reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz raw_data/T_apicales/HiFi/apicales_hifi-reads.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz raw_data/T_apicales/HiFi/apicales_hifi-3rdSMRTcell.fastq.gz

#Tellseq reads,  from 3rd tellseq run:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_I1_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T505_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R1_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T505_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R2_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T505_R2.fastq.gz

ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_I1_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T507_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R1_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T507_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R2_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T507_R2.fastq.gz

#HiC reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_apicales/apicales-286172_S3HiC_R1.fastq-002.gz raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_apicales/apicales-286172_S3HiC_R2.fastq-001.gz raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
```
#### Fastqc
```bash
for ReadFile in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/*/*.fastq.gz); do
OutDir=$(dirname $ReadFile)/fastqc
OutFile=$(basename $ReadFile | sed 's@.fastq.gz@@g')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_fastqc.sh $ReadFile $OutDir $OutFile
done
#57205247-57205256
```
#### longQC
```bash
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/*.fastq.gz); do
Datatype=pb-hifi
OutDir=$(dirname $Reads)/longqc/$(basename $Reads | cut -d '.' -f1)
OutFile=$(basename $Reads | cut -d '.' -f1)
ProgDir=~/git_repos/Wrappers/NBI
echo ${OutDir}/${OutFile}
mkdir $(dirname $Reads)/longqc
sbatch $ProgDir/run_longqc.sh $Reads $OutDir $OutFile $Datatype
done #57206555-6
```
## Hifiasm

### Default settings
Tom Mathers has performed hifiasm assembly using the HiFi reads with default settings.
```bash
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trapi/Trapi_default.p_ctg.fa assembly/genome/T_apicales/hifiasm/default/Trapi_default.p_ctg.fa
```
#### Abyss
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
5789    5789    1018    6062    114841    253936    484828    324707    2804342    868.6e6    Trapi_default.p_ctg.fa

#### KAT
Versus HiFi reads:
```bash
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trapi/kat_comp-main.mx assembly/genome/T_apicales/hifiasm/default/kat/.
cp /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trapi/kat_comp-main.mx.spectra-cn.png assembly/genome/T_apicales/hifiasm/default/kat/.
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trapi/kat_comp.stats assembly/genome/T_apicales/hifiasm/default/kat/.
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trapi/kat_comp.dist_analysis.json assembly/genome/T_apicales/hifiasm/default/kat/.
```
Versus TellSeq reads:
```bash
source anaconda2-5.1.0;source activate kat
cd /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trapi #kat cannot find the assembly file when not in the same directory???
~/git_repos/Wrappers/NBI/submit-slurm_v1.1.pl -q nbi-long -m 200000 -c 32 -t 5-00:00:00 -e -j kat_comp -i "source package 7f4fb852-f5c2-4e4b-b9a6-e648176d5543;kat comp -t 32 -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/kat/kat_comp_vs_tellseq_reads -m 31 -H 100000000 -I 100000000 '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R2.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R2.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_I1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_I1.fastq.gz' 'Trapi_default.p_ctg.fa'"
#55340356
```
#### BUSCO
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/Trapi_default.p_ctg.fa); do
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
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default); do
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
#5789    5789    1018    6062    114841    253936    484828    324707    2804342    868.6e6    Trapi_default.p_ctg.fa
#Arthropoda:
#        C:94.4%[S:45.2%,D:49.2%],F:3.0%,M:2.6%,n:1013
#Insecta:
#        C:93.6%[S:45.8%,D:47.8%],F:3.3%,M:3.1%,n:1367
#Hemiptera:
#        C:94.0%[S:47.0%,D:47.0%],F:3.1%,M:2.9%,n:2510
```
#### Jellyfish/kmc
```bash
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/jellyfish
Outfile=default_HiFi
Reads=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/apicales_hifi-reads.fastq.gz
Reads2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/apicales_hifi-3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_jellyfish.sh $OutDir $Outfile $Reads $Reads2
#55426796
```
Giving default_HiFi_19mer_out.histo to GenomeScope a size of 450,475,659 is estimated (Av. read length 8,549, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_21mer_out.histo to GenomeScope a size of 448,208,575 is estimated (Av. read length 8,549, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_25mer_out.histo to GenomeScope a size of 443,089,611 is estimated (Av. read length 8,549, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_31mer_out.histo to GenomeScope a size of 878,468,665 is estimated (Av. read length 8,549, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_39mer_out.histo to GenomeScope a size of 363,563,900 is estimated (Av. read length 8,549, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_49mer_out.histo to GenomeScope a size of 381,631,324 is estimated (Av. read length 8,549, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_61mer_out.histo to GenomeScope a size of 392,026,162 is estimated (Av. read length 8,549, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_75mer_out.histo to GenomeScope a size of 821,182,569 is estimated (Av. read length 8,549, 10,000,000 max. k-mer coverage from jellyfish settings).

Full genomescope model does not capture the peak of the observed kmers except with k = 31.
```R
setwd("C:/Users/did23faz/OneDrive - Norwich Bioscience Institutes/Desktop/R")
dataframe19 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/jellyfish/default_HiFi_19mer_out.histo") 
dataframe21 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/jellyfish/default_HiFi_21mer_out.histo") 
dataframe25 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/jellyfish/default_HiFi_25mer_out.histo") 
dataframe31 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/jellyfish/default_HiFi_31mer_out.histo") 

#Plot kmer distribution:
plot(dataframe19[1:200,], type="l")
plot(dataframe21[1:200,], type="l")
plot(dataframe25[1:200,], type="l")
plot(dataframe31[1:200,], type="l")

#Ignore low frequency kmers:
plot(dataframe19[3:100,], type="l")
points(dataframe19[3:100,])
#Single copy region looks to be 3 - 80, peak = 16
sum(as.numeric(dataframe19[3:41343,1]*dataframe19[3:41343,2]))/16
#1,224,553,110
sum(as.numeric(dataframe19[3:80,1]*dataframe19[3:80,2]))/16
#594,229,709

plot(dataframe21[3:100,], type="l")
points(dataframe21[3:100,])
#Single copy region looks to be 3 - 80, peak = 15
sum(as.numeric(dataframe21[3:39451,1]*dataframe21[3:39451,2]))/15
#1,287,959,928
sum(as.numeric(dataframe21[3:80,1]*dataframe21[3:80,2]))/15
#662,954,879

plot(dataframe25[3:100,], type="l")
points(dataframe25[3:100,])
#Single copy region looks to be 3 - 80, peak = 14
sum(as.numeric(dataframe25[3:36326,1]*dataframe25[3:36326,2]))/14
#1,344,716,739
sum(as.numeric(dataframe25[3:80,1]*dataframe25[3:80,2]))/14
#743,172,579

plot(dataframe31[3:100,], type="l")
points(dataframe31[3:100,])
#Single copy region looks to be 3 - 80, peak = 13
sum(as.numeric(dataframe31[3:32265,1]*dataframe31[3:32265,2]))/14
#1,299,595,401
sum(as.numeric(dataframe31[3:80,1]*dataframe31[3:80,2]))/14
#773,044,449

dataframe39 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/jellyfish/default_HiFi_39mer_out.histo")
dataframe49 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/jellyfish/default_HiFi_49mer_out.histo")
dataframe61 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/jellyfish/default_HiFi_61mer_out.histo")
dataframe75 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/jellyfish/default_HiFi_75mer_out.histo")

plot(dataframe39[1:200,], type="l")
plot(dataframe49[1:200,], type="l")
plot(dataframe61[1:200,], type="l")
plot(dataframe75[1:200,], type="l")

plot(dataframe39[7:100,], type="l")
points(dataframe39[7:100,])
#Single copy region looks to be 6 - 80, peak = 18
sum(as.numeric(dataframe39[7:100000,1]*dataframe39[7:100000,2]))/18
#1,026,948,894
sum(as.numeric(dataframe39[7:80,1]*dataframe39[7:80,2]))/18
#784,809,759

plot(dataframe49[7:100,], type="l")
points(dataframe49[7:100,])
#Single copy region looks to be 6 - 80, peak = 17
sum(as.numeric(dataframe49[7:100000,1]*dataframe49[7:100000,2]))/17
#1,112,206,740
sum(as.numeric(dataframe49[7:80,1]*dataframe49[7:80,2]))/17
#871,583,801

plot(dataframe61[7:100,], type="l")
points(dataframe61[7:100,])
#Single copy region looks to be 6 - 80, peak = 16
sum(as.numeric(dataframe61[7:100000,1]*dataframe61[7:100000,2]))/16
#1,197,961,568
sum(as.numeric(dataframe61[7:80,1]*dataframe61[7:80,2]))/16
#960,803,812

plot(dataframe75[7:100,], type="l")
points(dataframe75[7:100,])
#Single copy region looks to be 6 - 80, peak = 16
sum(as.numeric(dataframe75[7:100000,1]*dataframe75[7:100000,2]))/16
#1,204,614,045
sum(as.numeric(dataframe75[7:80,1]*dataframe75[7:80,2]))/16
#986,780,093
```
K-mer estimates of genome size range from 607-672-880, genomescope plots with smaller sizes do not look to capture all of the observed kmers.

Genomescope and kmer distribution plots give different genome sizes.


```bash
  for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/apicales_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/apicales_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm@g')/645m
    OutFile=T_apicales_645m
    Haploid_Genomesize=645m #based on 21mers
    Homozygous_coverage=30 #based upon KAT spectra (0x approaches 0)
    Min_contig=2 #default
    Purge_haplotigs_level=3 #default (strict)
    Kmer_cuttoff=5.0 #default
    Overlap_iterations=200 #increasing can improve assembly quality 
    Kmer_size=51 #default
    Similarity_threshold=0.75 #default
    mkdir -p $OutDir
    sbatch $ProgDir/run_hifiasm.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2
  done #55339301, 55339623

#n    n:500 n:N50 min   N80   N50   N20   max   sum
#7735 7735  1455  4192  98560 208448  399581  2804342 1.023e9 assembly/genome/T_apicales/hifiasm/645m/T_apicales_645m.bp.p_ctg.fa

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/645m/T_apicales_645m.bp.p_ctg.fa); do
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
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/apicales_hifi-reads.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/apicales_hifi-3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1
done
```

```bash
  for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/apicales_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/apicales_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm@g')/878m
    OutFile=T_apicales_878m
    Haploid_Genomesize=878m #based on genomescope
    Homozygous_coverage=30 #based upon KAT spectra (0x approaches 0)
    Min_contig=2 #default
    Purge_haplotigs_level=3 #default (strict)
    Kmer_cuttoff=5.0 #default
    Overlap_iterations=200 #increasing can improve assembly quality 
    Kmer_size=51 #default
    Similarity_threshold=0.75 #default
    mkdir -p $OutDir
    sbatch $ProgDir/run_hifiasm.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2
  done #55341159

#n  n:500 L50 min N80 N50 N20 max sum
#7,745  7,745 1,457 7,321 98,475  209,373 397,104 2,804,342 1,022,000,000

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/878m/T_apicales_878m.bp.p_ctg.fa); do
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
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/apicales_hifi-reads.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/apicales_hifi-3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1
done
```     

Assembly was tried with a range of genome sizes based upon the k-mer distribution and genomescope estimates:

```bash
Haploid_Genomesize_values=("360m" "510m" "555m" "600m" "645m" "675m" "820m" "880m")
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
ReadDir=$(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi)
ProgDir=~/git_repos/Wrappers/NBI
Run1=$(ls $ReadDir/apicales_hifi-reads.fastq.gz)
Run2=$(ls $ReadDir/apicales_hifi-3rdSMRTcell.fastq.gz)
OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm_19.5@g')/${Haploid_Genomesize}/${Homozygous_coverage}/${Purge_haplotigs_level}/${Kmer_cuttoff}/${Similarity_threshold}
OutFile=T_apicales_${Haploid_Genomesize}_${Homozygous_coverage}_${Purge_haplotigs_level}_${Kmer_cuttoff}_${Similarity_threshold}
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
echo ${OutDir}/$OutFile >> logs/apicales_hifilog.txt
if [ -s "${OutDir}/${OutFile}.bp.p_ctg.fa" ]; then
echo Already done for: $OutFile
else 
echo Running for: $OutFile
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/apicales_hifilog.txt
fi
done
done
done
done
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/*/*/*/*/*/*.bp.p_ctg.fa); do
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
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/*/*/*/*/*); do
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
grep -m 1 'peak_hom:' slurm.56260026.err slurm.56260027.err slurm.56260028.err slurm.56260029.err slurm.56260030.err slurm.56260031.err slurm.56260032.err slurm.56260033.err

#360m:
#slurm.56260026.err:[M::ha_ft_gen] peak_hom: 68; peak_het: 19
#510m:
#slurm.56260027.err:[M::ha_ft_gen] peak_hom: 48; peak_het: 19
#555m: 
#slurm.56260028.err:[M::ha_ft_gen] peak_hom: 44; peak_het: 19
#600m: 
#slurm.56260029.err:[M::ha_ft_gen] peak_hom: 40; peak_het: 19
#645m: 
#slurm.56260030.err:[M::ha_ft_gen] peak_hom: 38; peak_het: 19
#675m: 
#slurm.56260031.err:[M::ha_ft_gen] peak_hom: 36; peak_het: 19
#820m: 
#slurm.56260032.err:[M::ha_ft_gen] peak_hom: 29; peak_het: 19
#880m:
#slurm.56260033.err:[M::ha_ft_gen] peak_hom: 19; peak_het: -1
```
The results from this range are very similar, 880m gives slightly better BUSCO completeness and 510m gives slightly better N50. 
```bash
Haploid_Genomesize_values=("510m" "880m")
Homozygous_coverage_values=(19 29 38 48)
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
ReadDir=$(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi)
ProgDir=~/git_repos/Wrappers/NBI
Run1=$(ls $ReadDir/apicales_hifi-reads.fastq.gz)
Run2=$(ls $ReadDir/apicales_hifi-3rdSMRTcell.fastq.gz)
OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm_19.5@g')/${Haploid_Genomesize}/${Homozygous_coverage}/${Purge_haplotigs_level}/${Kmer_cuttoff}/${Similarity_threshold}
OutFile=T_apicales_${Haploid_Genomesize}_${Homozygous_coverage}_${Purge_haplotigs_level}_${Kmer_cuttoff}_${Similarity_threshold}
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
echo ${OutDir}/$OutFile >> logs/apicales_hifilog.txt
if [ -s "${OutDir}/${OutFile}.bp.p_ctg.fa" ]; then
echo Already done for: $OutFile
else 
echo Running for: $OutFile
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/apicales_hifilog.txt
fi
done
done
done
done 
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/*/*/*/*/*/*.bp.p_ctg.fa); do
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
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/*/*/*/*/*); do
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

cp temp.txt Reports/apicales_assembly_report3.txt
```
Amongst hifiasm_19.5 assemblies the highist BUSCO score is 94.7% of Hemiptera BUSCOs, 105 assemblies have this score, amongst these T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa has highest contiguity with 7,411 contigs and N50=223,272. The most contiguous hifiasm_19.5 assembly was T_apicales_510m_48_3_10.0_0.25.bp.p_ctg.fa with 3,811 contigs, N50=330,047, 93.4% of Hemiptera BUSCOs, this is also the highest N50. T_apicales_880m_48_2_5.0_0.25.bp.p_ctg.fa has 4,505 contigs, N50=298,439, 94% of Hemiptera BUSCOs.

```bash
#These assemblies were kept:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/Trapi_default.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/510m/48/3/10.0/0.25/T_apicales_510m_48_3_10.0_0.25.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/48/2/5.0/0.25/T_apicales_880m_48_2_5.0_0.25.bp.p_ctg.fa

#All other assemblies were deleted:
for assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/*/*/*/*/*/*.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/*/*.fa | grep -v 'Trapi_default.p_ctg.fa\|T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa\|T_apicales_510m_48_3_10.0_0.25.bp.p_ctg.fa\|T_apicales_880m_48_2_5.0_0.25.bp.p_ctg.fa'); do
rm $assembly  
done
```
#### KAT
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/*/*/*/*/*/*.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_hifi_reads
F1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
R1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp.sh $OutDir $Outfile $Genome $F1 $R1
done #56881253, 56881254, 56881255

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/*/*/*/*/*/*.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_tellseq_reads
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R1.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R2.fastq.gz
F2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R1.fastq.gz
R2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1
done #56881258, 56881259, 56881260
```
#### Blobtools
Blobtools with Hifi reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutDir=$(dirname $Reference)/minimap2
Outfile=$(basename $Reference | sed 's@.bp.p_ctg.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Reference $Read1 $Read2 #56939444

#Blast
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/blast2.12.0
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_23092023/5/nt
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blastn.sh $Assembly $Database $OutDir $OutPrefix 
#57175082, 57182846

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/blast2.7.1
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_23092023/4/nt
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blastn_2.7.1.sh $Assembly $Database $OutDir $OutPrefix 
#57193492

#BUSCO - keeping output files
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done 
#57207206

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
MappingFile=$(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
BlastFile=$(dirname $Assembly)/blast2.12.0/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').vs.nt.mts1.hsp1.1e25.megablast.out
OutDir=$(dirname $Assembly)/blobtools1.1.1
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
ColourFile=NA
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobplot.sh $Assembly $MappingFile $BlastFile $OutDir $OutPrefix $ColourFile
#57160576, 57193462
```
Blobtools with HiC reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutDir=$(dirname $Reference)/bwa
Outfile=$(basename $Reference | sed 's@.bp.p_ctg.fa@@g')_HiC
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Reference $Read1 $Read2 
#57207042
```
Blobtools with TellSeq reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutDir=$(dirname $Reference)/bwa
Outfile=$(basename $Reference | sed 's@.bp.p_ctg.fa@@g')_Tellseq
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R2.fastq.gz
Read3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R1.fastq.gz
Read4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Reference $Read1 $Read2 $Read3 $Read4 
#57207047
```
#### Kraken 
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_kraken2nt
OutDir=$(dirname $Assembly)/kraken2.1.3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/kraken/nt_14092023
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_kraken2.sh $Assembly $Database $OutDir $OutPrefix 2>&1 >> ${OutDir}/log.txt
```
#### Purge assembly or remove haplotigs in pretext?
```bash
source purge_dups-git05062020
split_fa Trapi_default.p_ctg.fa > Trapi_default.p_ctg.fa.split
```
#### Scaffold with tellseq reads

#### 3dDNA



















































































1:
/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/apicales_hifi-reads.fastq.gz
2:
/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/apicales_hifi-3rdSMRTcell.fastq.gz

### Notes on T.Mathers Trioza anthrisci
default:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
8639    8639    1111    5968    66931    196278    413917    269446    2788957    787.3e6    Trant_default.p_ctg.fa

    C:92.6%[S:78.5%,D:14.1%],F:3.0%,M:4.4%,n:1066

    987    Complete BUSCOs (C)
    837    Complete and single-copy BUSCOs (S)
    150    Complete and duplicated BUSCOs (D)
    32    Fragmented BUSCOs (F)
    47    Missing BUSCOs (M)
    1066    Total BUSCO groups searched

without purging:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
12813    12813    2433    4826    44086    99656    207050    142852    2799219    854.2e6    Trant_l0.p_ctg.fa

Flye with scaffolding:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
13530    13527    2613    514    40623    86028    179472    118970    690345    794.8e6    assembly.fasta

Manually selecting homo coverage from kmers:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
8648    8648    1116    5968    66820    196059    411004    268626    2788957    787.5e6    Trant_hom_cov_25.p_ctg.fa

hg-size 750m:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
7907    7907    964    8704    67507    206341    463316    287377    2801313    729.5e6    Trant_hg-size_set.p_ctg.fa

    C:94.9%[S:85.1%,D:9.8%],F:1.7%,M:3.4%,n:1066

    1011    Complete BUSCOs (C)
    907    Complete and single-copy BUSCOs (S)
    104    Complete and duplicated BUSCOs (D)
    18    Fragmented BUSCOs (F)
    37    Missing BUSCOs (M)
    1066    Total BUSCO groups searched

3dDNA:
juicer -> 3dDNA - without HiC data???
fragmented assembly - try with different ERC (expected run count) values

different -s values:
  default:
  n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
  8004    8004    958    5968    66414    207299    456948    292192    2788957    730.8e6    Trant_default.p_ctg.fa

  0.5:
  n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
  7979    7979    951    5968    66383    208789    457683    291866    2788957    727.9e6    Trant_default_hom_cov_30_s_0.5.p_ctg.fa

  0.25
  n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
  7497    7497    883    5968    71044    221944    466650    303088    2788957    704.8e6    Trant_default_hom_cov_30_s_0.25.p_ctg.fa

default, -D 10 --hg-size 750m:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
7782    7782    946    8819    71213    211226    460235    292099    2801313    728.2e6    Trant_hg-size_750m_D10.bp.p_ctg.fa

default, -D 20 -N200 --hg-size 750m:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
7577    7577    904    8650    72684    222221    479180    307570    2801313    727.8e6    Trant_hg-size_750m_D20_N200.bp.p_ctg.fa

mashmap to D. citri

Flye:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
8911    8910    1667    500    61309    130262    267373    178648    1216246    758.1e6    assembly.fasta

### T. Mathers Trioza apicales cont.
default:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
5804    5804    985    5105    102461    230316    469459    307658    2815354    791.1e6    Trapi_default.p_ctg.fa


default after purging duplicates:

split assembly

get mapping stats and auto cutoffs

minimap self align

purge dups

n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
4189    4189    694    7139    111199    259480    506701    332135    2815354    612.2e6    purged.fa
5804    5804    985    5105    102461    230316    469459    307658    2815354    791.1e6    Trapi_default.p_ctg.fa

 C:95.5%[S:88.5%,D:7.0%],F:1.6%,M:2.9%,n:1066

    1018    Complete BUSCOs (C)
    943    Complete and single-copy BUSCOs (S)
    75    Complete and duplicated BUSCOs (D)
    17    Fragmented BUSCOs (F)
    31    Missing BUSCOs (M)
    1066    Total BUSCO groups searched

Scaffold with tell-seq reads

scaffold with default settings - scaff10x gives segmentation fault
using tigmint longranger basic - no scaffolds made:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
6295    5819    739    501    102722    240538    481378    313202    2815354    612.1e6    draft.tigmint.arcs.fa
4189    4189    694    7139    111199    259480    506701    332135    2815354    612.2e6    draft.fa

With more data, default:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
5789    5789    1018    6062    114841    253936    484828    324707    2804342    868.6e6    Trapi_default.p_ctg.fa

--hg-size to 650m
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
6378    6378    1180    4192    110990    235242    444264    308940    2804342    932.9e6    Trapi_hg-size_650m.bp.p_ctg.fa
5789    5789    1018    6062    114841    253936    484828    324707    2804342    868.6e6    ../Trapi_default.p_ctg.fa

--hg-size to 650m -D 10