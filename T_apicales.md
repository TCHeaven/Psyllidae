# Trioza apicales
Contains run by T.Heaven in assembly of Trioza apicales

Unless stated otherwise commands were performed from the directory /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae

## Collect data
```bash
#HiFi reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz raw_data/T_apicales/HiFi/apicales_hifi-reads.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz raw_data/T_apicales/HiFi/apicales_hifi-3rdSMRTcell.fastq.gz

#HiC reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_apicales/apicales-286172_S3HiC_R1.fastq-002.gz raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_apicales/apicales-286172_S3HiC_R2.fastq-001.gz raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
for file in raw_data/T_apicales/HiC/apicales_286172-S3HiC*.fastq.gz; do echo -n "$file: "; zcat "$file" | wc -l | awk '{print $1/4}'; done

#Tellseq reads,  from 3rd tellseq run:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_I1_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T505_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R1_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T505_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R2_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T505_R2.fastq.gz

ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_I1_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T507_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R1_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T507_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_T_anthrisci_T_apicales_May_2022/220505_NB501793_0306_AHHMK5BGXK/Caliber_tellseq_run3_T_anthrisci_T_apicales_R2_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_apicales/TellSeq/apicales_T507_R2.fastq.gz

for file in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T505_barcoded.fastq.gz; do echo -n "$file: "; zcat "$file" | wc -l | awk '{print $1/4}'; done
for file in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T507_barcoded.fastq.gz; do echo -n "$file: "; zcat "$file" | wc -l | awk '{print $1/4}'; done

#RNA
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/Pf/Pf_1.fq.gz raw_data/T_apicales/RNASeq/T_apicales_Pf_1.fq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/Pf/Pf_2.fq.gz raw_data/T_apicales/RNASeq/T_apicales_Pf_2.fq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/Pm/Pm_1.fq.gz raw_data/T_apicales/RNASeq/T_apicales_Pm_1.fq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/Pm/Pm_2.fq.gz raw_data/T_apicales/RNASeq/T_apicales_Pm_2.fq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/P9/P9_1.fq.gz raw_data/T_apicales/RNASeq/T_apicales_P9_1.fq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/RNASeq/Trioza_psyllids_2022/RNAseqMay2022/X204SC22051079-Z01-F001/raw_data/P9/P9_2.fq.gz raw_data/T_apicales/RNASeq/T_apicales_P9_2.fq.gz
```
Remove HiFi reads containing adapters:
```bash
for InFile in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/*.fastq.gz); do
OutDir=$(dirname $InFile)/filtered
OutFile=$(basename $InFile | sed 's@.fastq.gz@@g')_filtered
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_HiFiAdapterFilt.sh $InFile $OutDir $OutFile
done 
#57334905,6
```
TellSeq reads converted to 10X format:
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run3_T_anthrisci_T_apicales/10x_conversion/T505/*fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x/.

ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run3_T_anthrisci_T_apicales/10x_conversion/T507/*fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x/.
```
TellSeq reads - remove the internal barcodes and adapters:
```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger

for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x/*_S1_L001_R1_001.fastq.gz); do
Fread=$Reads
Rread=$(echo $Reads | sed 's@_S1_L001_R1_001.fastq.gz@_S1_L001_R2_001.fastq.gz@g')
Run=$(basename $Reads | cut -d '_' -f1)
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger
OutFile=Trapi_${Run}
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_longranger_basic.sh $Fread $Rread $Run $OutDir $OutFile
done 
#57160255-6
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

Based on Merqury/KAT plot the true homozygous coverage is ~38

```bash
#These assemblies were kept:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/default/Trapi_default.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/510m/48/3/10.0/0.25/T_apicales_510m_48_3_10.0_0.25.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/48/2/5.0/0.25/T_apicales_880m_48_2_5.0_0.25.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/38/0/5.0/0.75/T_apicales_880m_38_0_5.0_0.75.bp.p_ctg.fa

#All other assemblies were deleted:
for assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/*/*/*/*/*/*.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm/*/*.fa | grep -v 'Trapi_default.p_ctg.fa\|T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa\|T_apicales_510m_48_3_10.0_0.25.bp.p_ctg.fa\|T_apicales_880m_48_2_5.0_0.25.bp.p_ctg.fa\|T_apicales_880m_38_0_5.0_0.75.bp.p_ctg.fa'); do
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
#### Merqury
Kmer plots with all reads; Hifi, HiC, Tellseq
```bash
#Prepare meryl kmer counts
source /nbi/software/staging/RCSUPPORT-2452/stagingloader
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/*/*.fastq.gz); do
OutDir=$(dirname $Reads)/meryl
mkdir $OutDir
meryl k=21 count output ${OutDir}/$(basename $Reads | sed 's@.fastq.gz@@g').meryl $Reads 
done

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
mkdir $(dirname $Assembly)/meryl
meryl union-sum output $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/*/meryl/*.meryl
mkdir $(dirname $Assembly)/meryl/HiFi 
meryl union-sum output $(dirname $Assembly)/meryl/HiFi/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi/meryl/*.meryl
mkdir $(dirname $Assembly)/meryl/HiC
meryl union-sum output $(dirname $Assembly)/meryl/HiC/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/meryl/*.meryl
mkdir $(dirname $Assembly)/meryl/Tellseq
meryl union-sum output $(dirname $Assembly)/meryl/Tellseq/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/meryl/*.meryl

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
#57159097

source /nbi/software/staging/RCSUPPORT-2452/stagingloader
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/meryl/HiC/T_apicales_880m_29_3_3.0_0.75.meryl
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/meryl/HiFi/T_apicales_880m_29_3_3.0_0.75.meryl
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/meryl/Tellseq/T_apicales_880m_29_3_3.0_0.75.meryl
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/meryl/T_apicales_880m_29_3_3.0_0.75.meryl

######################################################################################################################################################################
source /nbi/software/staging/RCSUPPORT-2452/stagingloader
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_*/TellSeq/longranger/*.fastq.gz); do
OutDir=$(dirname $Reads)/meryl
mkdir $OutDir
meryl k=21 count output ${OutDir}/$(basename $Reads | sed 's@.fastq.gz@@g').meryl $Reads 
done
#57186506

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
mkdir -p $(dirname $Assembly)/meryl/Tellseq_trimmed
meryl union-sum output $(dirname $Assembly)/meryl/Tellseq_trimmed/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/meryl/*.meryl

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
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi temp/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC temp/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger temp/.
mkdir -p $(dirname $Assembly)/meryl/All_Tellseq_trimmed
meryl union-sum output $(dirname $Assembly)/meryl/All_Tellseq_trimmed/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl temp/*/meryl/*.meryl
rm -r temp

mkdir -p $(dirname $Assembly)/merqury/All_Tellseq_trimmed
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2 #56939444
sbatch $ProgDir/run_qualimap.sh $(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam $Assembly $OutDir #57111022

#Blast
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/blast2.12.0/3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blastn.sh $Assembly $Database $OutDir $OutPrefix 
#57175082, 57182846,57308571,57398067,57044987
#parse_seqids is required for -taxid_map however when running databases created with -parse_seqids blast runs fail with c++ errors

#BUSCO - keeping output files
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done 
#57207206, 57307068, 57307079, 57307096

#Diamond blast - BUSCO regions
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
BUSCOgff=$(dirname $Assembly)/BUSCO/hemiptera_odb10/run_hemiptera_odb10/metaeuk_output/rerun_results/T_apicales_880m_29_3_3_hemiptera_odb10.fa.gff
awk '$3 == "gene" {print $0}' $BUSCOgff > $(echo $BUSCOgff | sed 's@.fa.gff@_genes_only.fa.gff@g')
bedtools getfasta -fi $Assembly -bed $(echo $BUSCOgff | sed 's@.fa.gff@_genes_only.fa.gff@g') -fo $(dirname $BUSCOgff)/busco_regions.fasta
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/uniprot_01102023/Uniprot_01102023_reference_proteomes.dmnd
ProgDir=~/git_repos/Wrappers/NBI
mkdir -p $OutDir
sbatch $ProgDir/run_diamond_blastx.sh $(dirname $BUSCOgff)/busco_regions.fasta $Database $OutDir $OutPrefix 
#57307105, 57315793

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
#57307106, 57315794

BUSCODiamond=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions/*.diamondblastx.out
ElseDiamond=$(dirname $Assembly)/diamond0.9.29_blastx/nonBUSCO_regions/*.diamondblastx.out
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/normalise_blast.py $BUSCODiamond $(echo $BUSCODiamond | sed 's@x.out@x_2.out@g')
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/normalise_blast.py $ElseDiamond $(echo $ElseDiamond | sed 's@x.out@x_2.out@g')

#Tiara
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/tiara
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_tiara.sh $Assembly $OutDir $OutPrefix
#57401241

#Blobtools
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

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
Species=apicalis
TaxID=872318
alias=Tapi88029330075
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57127861

cp -r $OutDir/Tapi88029330075_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
#The contents of /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs were subsequently copied to local machine WSL-Ubuntu: \\wsl.localhost\Ubuntu\home\did23faz
```
```bash
conda activate btk
#From \\wsl.localhost\Ubuntu\home\did23faz
apptainer exec blobtoolkit.sif blobtools host BlobDirs
```
Blobtools with HiC reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_HiC
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#57207042
MappingFile=$(ls ${OutDir}/*_HiC.bam)
sbatch $ProgDir/run_qualimap.sh $MappingFile $Assembly $OutDir #57111435

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
Species=apicalis
TaxID=872318
alias=Tapi88029330075_HiC
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57217721

cp -r $OutDir/Tapi88029330075_HiC_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
```
Blobtools with TellSeq reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R2.fastq.gz
Read3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R1.fastq.gz
Read4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Read1 $Read2 $Read3 $Read4 
#57207047
MappingFile=$(ls ${OutDir}/*_Tellseq.bam)
sbatch $ProgDir/run_qualimap.sh $MappingFile $Assembly $OutDir #57111437

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
Species=apicalis
TaxID=872318
alias=Tapi88029330075_Tellseq
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57218783
#########################################################################################################################
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T505_barcoded.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T507_barcoded.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#57186561
MappingFile=$(ls ${OutDir}/*_Tellseq_trimmed.bam)
sbatch $ProgDir/run_qualimap.sh $MappingFile $Assembly $OutDir #57317734

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
Species=apicalis
TaxID=872318
alias=Tapi88029330075_Tellseq_trimmed
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57331164

cp -r $OutDir/Tapi88029330075_Tellseq_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
cp -r $OutDir/Tapi88029330075_Tellseq_trimmed_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
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
Contigs classified to non-Arthropoda phylum level were considered contamininants. Contigs classified to Arthropoda taxa or above (eg. Eukaryota) were considered Psyllid contigs.
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
cd $(dirname $Assembly)/kraken2.1.3/
nano contaminantlist.txt #Edit with contaminant names
grep -f contaminantlist.txt $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_kraken2nt_output.txt > contaminantcontigs.txt
awk -F'\t' '{print $2}' contaminantcontigs.txt > contaminantcontignames.txt
rm contaminantcontigs.txt
seqtk subseq $Assembly contaminantcontignames.txt | gzip > kraken2_contaminants.fa.gz
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/filter.py  $Assembly contaminantcontignames.txt > filtered_$(basename $Assembly)
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' filtered_$(basename $Assembly) > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_krakenfiltered.fa
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
Contigs identified as bacteria and archaea do not form a cluster, however several are very large, there is only one contig that is identified as bacteria by tiara and arthropoda by blast
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
cd $(dirname $Assembly)/tiara/
nano contaminantcontignames.txt #edit with bacteria and archaea contig names from tiara
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/blobtoolkit4.2.1/T_apicales_880m_29_3_3.0_0.75_btkfiltered2.fa
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
T_apicales_880m_29_3_3.0_0.75_hemiptera_odb10_short_summary.txt
        C:94.7%[S:26.3%,D:68.4%],F:2.7%,M:2.6%,n:2510

T_apicales_880m_29_3_3.0_0.75_krakenfiltered_hemiptera_odb10_short_summary.txt
        C:94.5%[S:26.5%,D:68.0%],F:2.7%,M:2.8%,n:2510
T_apicales_880m_29_3_3.0_0.75_blastfiltered_hemiptera_odb10_short_summary.txt
        C:94.4%[S:26.8%,D:67.6%],F:2.7%,M:2.9%,n:2510
T_apicales_880m_29_3_3.0_0.75_tiarafiltered_hemiptera_odb10_short_summary.txt
        C:94.7%[S:26.3%,D:68.4%],F:2.7%,M:2.6%,n:2510

T_apicales_880m_29_3_3.0_0.75_btkfiltered_hemiptera_odb10_short_summary.txt
        C:94.6%[S:26.2%,D:68.4%],F:2.7%,M:2.7%,n:2510
T_apicales_880m_29_3_3.0_0.75_filtered_hemiptera_odb10_short_summary.txt
        C:94.3%[S:26.7%,D:67.6%],F:2.7%,M:3.0%,n:2510

T_apicales_880m_29_3_3.0_0.75_btkfiltered2_hemiptera_odb10_short_summary.txt
        C:94.7%[S:26.3%,D:68.4%],F:2.7%,M:2.6%,n:2510
T_apicales_880m_29_3_3.0_0.75_allfiltered_hemiptera_odb10_short_summary.txt
        C:94.4%[S:26.8%,D:67.6%],F:2.7%,M:2.9%,n:2510
```bash
#BUSCO - keeping output files
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done 
#57044709

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done 
#57044749

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
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
cat $(dirname $Assembly)/filtered/contaminantcontignames.txt | wc -l #206

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
T_apicales_880m_29_3_3.0_0.75_filtered_hemiptera_odb10_short_summary
  C:94.7%[S:26.3%,D:68.4%],F:2.7%,M:2.6%,n:2510

### Purge Dups
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
#Purge with Tellseq (Illumina) Reads:
MappingFile=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57347671

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4278    4278    676     6289    118431  271325  498877  2845242 603.9e6 T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.fa

#Purge with HiFi reads:
MappingFile=$(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
Type=long
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57344866

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4463    4463    789     7321    119249  266404  488568  2845242 683.2e6 T_apicales_880m_29_3_3.0_0.75_HiFiPurged.fa
```
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_*Purged.fa); do
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
```
T_apicales_880m_29_3_3.0_0.75_filtered_hemiptera_odb10_short_summary
  C:94.7%[S:26.3%,D:68.4%],F:2.7%,M:2.6%,n:2510

T_apicales_880m_29_3_3.0_0.75_HiFiPurged_arthropoda_odb10_short_summary.txt
        C:93.8%[S:74.6%,D:19.2%],F:3.5%,M:2.7%,n:1013
T_apicales_880m_29_3_3.0_0.75_HiFiPurged_hemiptera_odb10_short_summary.txt
        C:93.8%[S:76.2%,D:17.6%],F:3.2%,M:3.0%,n:2510
T_apicales_880m_29_3_3.0_0.75_HiFiPurged_insecta_odb10_short_summary.txt
        C:93.2%[S:75.6%,D:17.6%],F:3.4%,M:3.4%,n:1367
        
T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_arthropoda_odb10_short_summary.txt
        C:92.8%[S:89.1%,D:3.7%],F:4.0%,M:3.2%,n:1013
T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_hemiptera_odb10_short_summary.txt
        C:93.3%[S:89.4%,D:3.9%],F:3.5%,M:3.2%,n:2510
T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_insecta_odb10_short_summary.txt
        C:92.1%[S:88.5%,D:3.6%],F:4.2%,M:3.7%,n:1367

Purging with TellSeq reads seems to be more effective than using HiFi reads
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/tellseq
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/tellseq/.
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/hifi
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_HiFiPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/hifi/.

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/*/*Purged.fa); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done #57291744,45
```
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/filtered/T_apicales_880m_29_3_3.0_0.75_filtered.fa
#Purge with Tellseq (Illumina) Reads:
MappingFile=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57051896

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4172    4172    663     6289    118442  273181  504047  2845242 595.8e6 T_apicales_880m_29_3_3.0_0.75_filtered_TellSeqPurged.fa

#Purge with HiFi reads:
MappingFile=$(dirname $Assembly)/../minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
Type=long
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57051897

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4172    4172    663     6289    118442  273181  504047  2845242 595.8e6 T_apicales_880m_29_3_3.0_0.75_filtered_HiFiPurged.fa
```
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/filtered/purge_dups/T_apicales_880m_29_3_3.0_0.75_*Purged.fa); do
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
```
T_apicales_880m_29_3_3.0_0.75_filtered_hemiptera_odb10_short_summary
  C:94.7%[S:26.3%,D:68.4%],F:2.7%,M:2.6%,n:2510

T_apicales_880m_29_3_3.0_0.75_filtered_HiFiPurged_arthropoda_odb10_short_summary.txt
        C:92.8%[S:89.1%,D:3.7%],F:4.0%,M:3.2%,n:1013
T_apicales_880m_29_3_3.0_0.75_filtered_HiFiPurged_hemiptera_odb10_short_summary.txt
        C:93.2%[S:89.4%,D:3.8%],F:3.5%,M:3.3%,n:2510
T_apicales_880m_29_3_3.0_0.75_filtered_HiFiPurged_insecta_odb10_short_summary.txt
        C:92.2%[S:88.5%,D:3.7%],F:4.1%,M:3.7%,n:1367
T_apicales_880m_29_3_3.0_0.75_filtered_TellSeqPurged_arthropoda_odb10_short_summary.txt
        C:92.8%[S:89.1%,D:3.7%],F:4.0%,M:3.2%,n:1013
T_apicales_880m_29_3_3.0_0.75_filtered_TellSeqPurged_hemiptera_odb10_short_summary.txt
        C:93.2%[S:89.4%,D:3.8%],F:3.5%,M:3.3%,n:2510
T_apicales_880m_29_3_3.0_0.75_filtered_TellSeqPurged_insecta_odb10_short_summary.txt
        C:92.2%[S:88.5%,D:3.7%],F:4.1%,M:3.7%,n:1367

```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/filtered/purge_dups/tellseq
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/filtered/purge_dups/T_apicales_880m_29_3_3.0_0.75_filtered_TellSeqPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/filtered/purge_dups/tellseq/.
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/filtered/purge_dups/hifi
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/filtered/purge_dups/T_apicales_880m_29_3_3.0_0.75_filtered_HiFiPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/filtered/purge_dups/hifi/.

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/filtered/purge_dups/*/*Purged.fa); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done #57291750,1
```
Purging with/without contaminants seems to make no difference to the final scores, neither does using Tellseq/Hifi reads.

### Scaffolding
Tellseq reads have are linked via 18bp barcodes, Tom Mathers has already used the conversion software prodived by Universal sequencing and the 4Mwith-alts-february-2016.txt barcode whitelist file to convert to 16bp barcoded versions of the reads that are compatible with 10x genomics linked read software. The files also have to be in a standardised naming format https://cdn.shopify.com/s/files/1/0654/0378/1341/files/100027-USG_TELL-Seq_Software_Roadmap_User_Guide_v1.0_4d2abbf8-3d19-4899-84c2-9a79594eb7f3.pdf?v=1684422029

Scaff10X can now be used to scaffold the HiFi assembly contigs using the TellSeq reads.

Longranger is run in order to remove the barcodes from the reads so that they can simply be treated as large insert illumina reads and used for assembly polishing.

Purging should be done before scaffolding as presence of haplotigs will confuse scaffolder, however the order of pilon,scaff10x,break10x,YAHS or whether using at all will improve the assembly is unclear.

#### 3DDNA
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.fa
OutDir=$(dirname $Assembly)/3ddna
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_3dDNA.sh $Assembly $OutDir $OutFile $Read1 $Read2
#57364471, 57615784
#NOTE: 3ddna output is very large ~600GB therefore only final files kept

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1989    1989    1       977     515.1e6 515.1e6 515.1e6 515.1e6 569.6e6 T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.FINAL.fasta
#One enourmous chromosome...

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/3ddna/juicer/aligned/genome_wrapped.FINAL.hic /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/.

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/
samtools faidx T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
cut -f1,2 T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta.fai > T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds.genome

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/juicer.sif juicer.sh pre --threads 1 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds_mapped.pairs /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds_contact_map.hic /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds.genome

submit-slurm_v1.1.pl -q ei-medium -m 200000 -c 32 -t 5-
00:00 -e -j Smis_v1_juicer -i "source jdk-1.7.0_25; source bwa-0.7.17;source samtools-1.6;~/juicer/juicer-
1.6.2/CPU/juicer.sh -D /hpchome/
tmathers/JIC_TM_scratch_DIR/Sitobion_miscanthi_pub_QC_and_assemble/HIC_map_published_
asm -d /hpchome/
tmathers/JIC_TM_scratch_DIR/Sitobion_miscanthi_pub_QC_and_assemble/HIC_map_published_
asm -g Smis_v1 -z references/Lachesis_assembly.fasta -y restriction_sites/Smis_v1_DpnII.txt -s DpnII -t
32 -p restriction_sites/Smis_v1_DpnII.chrom.sizes"

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/3ddna/juicer/aligned
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/juicer.sif run-assembly-visualizer.sh -q 10 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/3ddna/juicer/aligned/genome_wrapped.FINAL.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/3ddna/juicer/aligned/merged_nodups.txt

#From omni-c guide
source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
java17 -Xmx48000m -Djava.awt.headless=true -jar ~/git_repos/Scripts/NBI/juicer_tools_1.22.01.jar pre --threads 16 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds_mapped.pairs /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds_contact_map.hic /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds.genome
sbatch ~/git_repos/Wrappers/NBI/temp3.sh
#57554716, 57816287, 57818847
 ~/git_repos/Scripts/NBI/generate-assembly-file-from-fasta.awk /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta

source lastz-1.03.73;source gnu_parallel-20180322;source jdk-1.7.0_25
~/3ddna-master/visualize/run-assembly-visualizer.sh draft.assembly aligned/merged_nodups.txt
```

#### 0
Scaffolding without purging

##### Scaff10X
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234150

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#7411    7411    1390    6289    103691  223272  420487  2845242 1.032e9 output_scaffolds.fasta
```
##### YAHS
Generate mapping file
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57245206

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57278946

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#8832    8832    1568    1000    84335   191590  381738  2735242 1.032e9 T_apicales_880m_29_3_3.0_0.75_scaffolds_final.fa
```

#### 1
Purge
```bash
#Already performed above
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234130

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4278    4278    676     6289    118431  271325  498877  2845242 603.9e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57287409

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57297213

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1423    1423    5       1000    13.36e6 48.6e6  62.8e6  75.53e6 603.9e6 T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap3.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57677372, 57677375

Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaff10xscaffolds_scaffolds_final_mapped.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57681461

Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaff10xscaffolds_scaffolds_final_cat_mapped.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57681465
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57245204

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57279048

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1423    1423    5       1000    13.36e6 48.6e6  62.8e6  75.53e6 603.9e6 T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57285518

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1423    1423    5       1000    13.36e6 48.6e6  62.8e6  75.53e6 603.9e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57297936

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap2.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57301395



Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged_scaffolds_final_scaff10xscaffolds_mapped.PT.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.PT.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57301308, 57303608, 57307616
```
#### 2
Purge
```bash
#Already performed above
```
Pilon
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.fa
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

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.fa
Alignment=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@_TellSeqPurged.fa@@g')_sorted_markdups_Tellseq_trimmed.bam
Index=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@_TellSeqPurged.fa@@g')_sorted_markdups_Tellseq_trimmed.bam.bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@_TellSeqPurged.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57234145
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_pilon.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234829

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4278    4278    676     6289    118346  271326  498579  2845242 603.8e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/pilon/scaff10x/T_apicales_880m_29_3_3.0_0.75_pilon_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57297188

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/pilon/scaff10x/T_apicales_880m_29_3_3.0_0.75_pilon_scaff10xscaffolds.fasta
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_pilon.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57245211

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_pilon.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_pilon_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_pilon_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57279783

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1446    1446    7       1000    9620398 29.66e6 46.19e6 63.16e6 603.8e6 T_apicales_880m_29_3_3.0_0.75_pilon_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/pilon/yahs/T_apicales_880m_29_3_3.0_0.75_pilon_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57297218
```
#### 3
Break10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/break10x
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_break10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234087

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#7668    7668    1488    6289    99358   209609  392936  2845242 1.032e9 T_apicales_880m_29_3_3.0_0.75_break.fa

#BUSCO
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/T_apicales_880m_29_3_3.0_0.75_break.fa); do
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/T_apicales_880m_29_3_3.0_0.75_break.fa
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T505_barcoded.fastq.gz
T2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T507_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1 $T2
#57245191

#Purge with Tellseq (Illumina) Reads:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/T_apicales_880m_29_3_3.0_0.75_break.fa
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/bwa/T_apicales_880m_29_3_3.0_0.75_break_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57252575

sbatch ~/git_repos/Pipelines/Trioza_merqury.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa 
#57291756
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57258657

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4387    4387    714     6289    115705  259273  468815  2845242 603.4e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57287105

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57291005

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1448    1448    6       1000    14.65e6 44.7e6  63.55e6 64.22e6 603.4e6 T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap3.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57682437

Alignment=
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57276684

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57279790

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1448    1448    6       1000    14.65e6 44.7e6  63.55e6 64.22e6 603.4e6 T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa

awk '/^>/ { if (seq) { print id, length(seq) }; id = $0; seq = "" } /^[^>]/ { seq = seq $0 } END { if (seq) { print id, length(seq) } }' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa


Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap3.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57873985

Alignment=
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.PT.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57291011

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1448    1448    6       1000    14.65e6 44.7e6  63.55e6 64.22e6 603.4e6 output_scaffolds.fasta
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa
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

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/bwa/T_apicales_880m_29_3_3.0_0.75_break_sorted_markdups_Tellseq_trimmed.bam
Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/bwa/T_apicales_880m_29_3_3.0_0.75_break_sorted_markdups_Tellseq_trimmed.bam.bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57271539
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57276692

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4387    4387    714     6289    115705  259269  468768  2845242 603.2e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57285460

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57290956

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1481    1481    7       1000    12.18e6 32.16e6 63.47e6 65.84e6 603.2e6 T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291022

Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_mapped.PT.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.PT.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57297210 - 57306200

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_mapped.pairs | PretextMap -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_mapped.pairs.pretext --sortby length --sortorder descend --mapq 10 --highRes

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')_13
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap3.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57396583, 57398896, 57402359, 57437956, 57624391, 57624411, 57624413, 57624510, 57624515, 57624531, 57624561, 57624571, 57636606, 57673246
#2 57673274
#3 57673329
sbatch $ProgDir/temp4.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57673786,57674361,57676539

samtools view -h -b -q 40 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_mapped.bam > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_mapped_f.bam

samtools markdup --write-index -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_mapped.bam /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_mapped_f2.bam


for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_57673329/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_13_mapped.bam); do
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.bam@@')_v21
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57398903, 57406312, 57437893, 57450185, 57491068, 57491855, 57491907, 57493008, 57493026, 57493062, 57493069, 57493083,57561638,57561677,57624436,57624464,57624465,57625365,57625366,57625367,57677242,57677376
done
#57493189-57493192, 57493228
```
visualiser
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/juicer
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/juicer
mkdir restriction_sites;mkdir references;mkdir fastq
cd fastq
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
cd ../references
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa .
source bwa-0.7.17
bwa index T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa
#57571149
cd ../restriction_sites
Enzyme=DpnII
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/generate_site_positions.py $Enzyme OutFile ../references/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' OutFile_DpnII.txt > OutFile_DpnII.chrom.sizes
cd ..
mkdir -p scripts/common
cp ~/git_repos/Scripts/NBI/chimeric_blacklist.awk scripts/common/.
source jdk-1.7.0_25; source bwa-0.7.17;source samtools-1.6;source gnu_parallel-20180322
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/juicer.sif juicer.sh -D /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/juicer -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/scaff10x/yahs/juicer -g OutFile -z references/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa -y restriction_sites/OutFile_DpnII.txt -s DpnII -t 64 -p restriction_sites/OutFile_DpnII.chrom.sizes 
#57615361
source lastz-1.03.73;source gnu_parallel-20180322;source jdk-1.7.0_25
~/git_repos/Scripts/NBI/run-assembly-visualizer.sh


~/3ddna-master/visualize/run-assembly-visualizer.sh draft.assembly aligned/merged_nodups.txt


bwa mem -SP5M


singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/juicer.sif juicer.sh -t 16 -d $WorkDir/juicer -g saundersiae -z genome_wrapped.fa -y ${OutFile}_Phase.txt -p genome_wrapped.fa.fai -D /opt/juicer-1.6.2/CPU



```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57279810

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57285406

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1481    1481    7       1000    12.18e6 32.16e6 63.47e6 65.84e6 603.2e6 T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57285446

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1481    1481    7       1000    12.18e6 32.16e6 63.47e6 65.84e6 603.2e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/pilon/yahs/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_pilon_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291847
```
#### 5
Break10x
```bash
#Already performed above
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/T_apicales_880m_29_3_3.0_0.75_break.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57252782

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#7668    7668    1488    6289    99358   209609  392936  2845242 1.032e9 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291826

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57297197
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/T_apicales_880m_29_3_3.0_0.75_break.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57252786

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/T_apicales_880m_29_3_3.0_0.75_break.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/T_apicales_880m_29_3_3.0_0.75_break_mapped.PT.bam 
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/T_apicales_880m_29_3_3.0_0.75_break_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57285407

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#8851    8851    1588    1000    84279   190530  374636  2735242 1.032e9 T_apicales_880m_29_3_3.0_0.75_break_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/yahs/T_apicales_880m_29_3_3.0_0.75_break_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57297222
```
#### 6
Pilon
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/T_apicales_880m_29_3_3.0_0.75.bp.p_ctg.fa
Alignment=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_sorted_markdups_Tellseq_trimmed.bam
Index=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_sorted_markdups_Tellseq_trimmed.bam.bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57227855
```
Break10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/T_apicales_880m_29_3_3.0_0.75.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/break10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_break10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234086

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#7667    7667    1488    6289    99375   209586  392907  2845242 1.032e9 T_apicales_880m_29_3_3.0_0.75_break.fa
```
Purge
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/T_apicales_880m_29_3_3.0_0.75_break.fa
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T505_barcoded.fastq.gz
T2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T507_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1 $T2
#57245193

#Purge with Tellseq (Illumina) Reads:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/T_apicales_880m_29_3_3.0_0.75_break.fa
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/bwa/T_apicales_880m_29_3_3.0_0.75_break_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57252624

sbatch ~/git_repos/Pipelines/Trioza_merqury.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa
#57291759
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57276688

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#4391    4391    714     6289    115362  259269  468768  2845242 603.2e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291768

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/scaff10x/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57297198
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@_break_TellSeqPurged.fa@_pilon_break_TellSeqPurged@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57285409

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_pilon_break_TellSeqPurged_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_pilon_break_TellSeqPurged_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57285409

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#1477    1477    5       1000    20.46e6 32.93e6 72.57e6 89.14e6 603.2e6 T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57297227
```
### Filtering
#### MitoHiFi
```bash
#No reference T.apicales mitochondial genome is available, therefore both T.urticae and T.anthrisci was tried with Good_Reference=N:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/Trapi.curated_primary.no_mt.unscrubbed.fa
ReferenceFasta=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/NCBI/Trant_mito_NC_038141.1.fasta
ReferenceGenebank=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/NCBI/Trant_mito_NC_038141.1.gb
PercentOverlap=50
Code=5
Kingdom=animal
Good_Reference=N
OutDir=$(dirname $Assembly)/MitoHifi_tant
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_mitohifi_contigs.sh $Assembly $ReferenceFasta $ReferenceGenebank $PercentOverlap $Code $Kingdom $Good_Reference $OutDir 
#58510768

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/Trapi.curated_primary.no_mt.unscrubbed.fa
ReferenceFasta=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/NCBI/Turt_mito_NC_038113.1.fasta
ReferenceGenebank=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/NCBI/Turt_mito_NC_038113.1.gb
PercentOverlap=50
Code=5
Kingdom=animal
Good_Reference=N
OutDir=$(dirname $Assembly)/MitoHifi_turt
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_mitohifi_contigs.sh $Assembly $ReferenceFasta $ReferenceGenebank $PercentOverlap $Code $Kingdom $Good_Reference $OutDir 
#58510771

```
scaffold_198 was identified as the mitochondrial genome with both urticae and anthrisci references.
```bash
echo scaffold_198 > id.txt
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_rm.py  --id_file id.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/Trapi.curated_primary.no_mt.unscrubbed.fa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito.fa

#Mitochondrial genome:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi_tant/final_mitogenome.fasta
```

#### Kraken and Blobtools
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/Trapi.curated_primary.no_mt.unscrubbed.fa
Species=$(echo $Assembly | cut -d '/' -f10)
Genus=Trioza
TaxID=872318
alias=$(echo $Species | sed 's@_@@g' | cut -c1-4)
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T505_barcoded.fastq.gz
T2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T507_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI

#Kraken 
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')_kraken2nt
OutDir=$(dirname $Assembly)/kraken2.1.3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/kraken/nt_14092023
mkdir $OutDir
sbatch $ProgDir/run_kraken2.sh $Assembly $Database $OutDir $OutPrefix
#58511029

#Tiara
OutPrefix=${Species}
OutDir=$(dirname $Assembly)/tiara
mkdir $OutDir
sbatch $ProgDir/run_tiara.sh $Assembly $OutDir $OutPrefix
#58511030

#Alignment
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1 $T2
#58511031, 60032850

#Blast
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')
OutDir=$(dirname $Assembly)/blast2.12.0
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
Max_target=10
mkdir $OutDir
sbatch $ProgDir/run_blastn.sh $Assembly $Database $OutDir $OutPrefix $Max_target
#58511033

#BUSCO - keeping output files
OutDir=$(dirname $Assembly)/BUSCO
Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
OutFile=$(basename $Assembly | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
mkdir $OutDir 
sbatch $ProgDir/run_busco_keep.sh $Assembly $Database $OutDir $OutFile 
#58511036

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
#58524025

cp -r ${OutDir}/Tapi_Tellseq_trimmed_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
#The contents of /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs were subsequently copied to local machine WSL-Ubuntu: \\wsl.localhost\Ubuntu\home\did23faz
```
```bash
ubuntu
conda activate btk
#From \\wsl.localhost\Ubuntu\home\did23faz
apptainer exec blobtoolkit.sif blobtools host BlobDirs
```
Kraken and blobtools were used to screen the assembly for contaminants resulting in the removal of 72 scaffolds. All scaffolds without kraken classification to an arthropoda taxa outside of the range 0.3 - 0.4 GC, 40 - 110x coverage were removed (43), within this range scaffolds without classification to arthropoda taxa were kept if there classification was more general (eg. eukaryota), unclassified, or to a taxa it is implausible the sample to be contaminatied with (eg. tuna) (14 removed). Additionally, 15 scaffolds scaffolds were removed that were classified to arthropoda taxa by kraken, but which were identified by either tiara or blast as potential contaminants and fell outside of the 0.3 - 0.4 GC, 40 - 110x coverage range.

temp-cont4.txt

Suspected contaminants were removed to a seperate fasta file:
```bash
nano /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/suspected_contaminant_names.txt

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/suspected_contaminant_names.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito.fa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_suspected_contaminants.fa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_rm.py --id_file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/suspected_contaminant_names.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito.fa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered.fa
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

#T_apicales_SUPER_11 tens of thousands of bp alignment
#T_apicales_scaffold_49 <1220 bp alignments
#T_apicales_scaffold_18 <1009 bp alignments
#T_apicales_scaffold_15 <587 bp alignments

awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/Trapi.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>SUPER_11' > temp_SUPER_11.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/Trapi.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_49' > temp_scaffold_49.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/Trapi.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_18' > temp_scaffold_18.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/Trapi.curated_primary.no_mt.unscrubbed.fa | grep -A 1 -w '>scaffold_15' > temp_scaffold_15.fasta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_000010365.1/GCA_000010365.1_ASM1036v1_genomic.fna -p scaffold_SUPER_11
mummerplot -color scaffold_SUPER_11.delta
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_001274515.1/GCA_001274515.1_ASM127451v1_genomic.fna -p scaffold_SUPER_11
mummerplot -color scaffold_SUPER_11.delta
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_029909935.1/GCA_029909935.1_ASM2990993v1_genomic.fna -p scaffold_SUPER_11
mummerplot -color scaffold_SUPER_11.delta
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_000287275.1/GCA_000287275.1_ASM28727v1_genomic.fna -p scaffold_SUPER_11
mummerplot -color scaffold_SUPER_11.delta
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_SUPER_11
mummerplot -l -c scaffold_SUPER_11.delta
mummerplot -color scaffold_SUPER_11.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta -p scaffold_SUPER_11-2
mummerplot -l -c scaffold_SUPER_11-2.delta
mummerplot -color scaffold_SUPER_11-2.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_49.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_49
mummerplot -l -c scaffold_49.delta
mummerplot -color scaffold_49.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_49.fasta -p scaffold_49-2
mummerplot -l -c scaffold_49-2.delta
mummerplot -color scaffold_49-2.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_18
mummerplot -l -c scaffold_18.delta
mummerplot -color scaffold_18.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta -p scaffold_18-2
mummerplot -l -c scaffold_18-2.delta
mummerplot -color scaffold_18-2.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_15.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_15
mummerplot -l -c scaffold_15.delta
mummerplot -color scaffold_15.delta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_15.fasta -p scaffold_15-2
mummerplot -l -c scaffold_15-2.delta
mummerplot -color scaffold_15-2.delta

cat temp_SUPER_11.fasta temp_scaffold_49.fasta temp_scaffold_18.fasta temp_scaffold_15.fasta > temp_all.fasta
nucmer --maxmatch --nosimplify ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_all.fasta -p temp_all_api
mummerplot -color temp_all_api.delta


#unscaffolded contigs:
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa -input_type fasta -dbtype nucl -title psyllid1  -parse_seqids -out psyllid1
blastn -query ../Symbionts/Candidatus/Carsonella/ruddii/GCA_000287275.1/GCA_000287275.1_ASM28727v1_genomic.fna -db psyllid1 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/Carsonella_ruddii_results_2 -evalue 1e-5 -outfmt 6 -num_threads 1

awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa | grep -A 1 -w '>tarseq_44' > temp_SUPER_44.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa | grep -A 1 -w '>tarseq_2188' > temp_scaffold_2188.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa | grep -A 1 -w '>tarseq_107' > temp_scaffold_107.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa | grep -A 1 -w '>tarseq_152' > temp_scaffold_152.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/pilon/break10x/purge_dups/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged.fa | grep -A 1 -w '>tarseq_1132' > temp_scaffold_1132.fasta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_1132.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_1132
mummerplot -color scaffold_1132.delta #~70,000/25000bp
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_44.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_44
mummerplot -color scaffold_44.delta #50%
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_2188.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_2188
mummerplot -color scaffold_2188.delta
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_107.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_107
mummerplot -color scaffold_107.delta
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_152.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p scaffold_152
mummerplot -color scaffold_152.delta
```
Alignments show that a full Candidatus Carsonella ruddii has been integrated in SUPER_11 scaffold, alignments to other scaffolds are short.
```python
from Bio import SeqIO

def extract_sequence_to_file(fasta_file, start_pos, end_pos, output_file):
  for record in SeqIO.parse(fasta_file, "fasta"):
    sequence = record.seq[start_pos-1:end_pos]
    with open(output_file, "w") as f:
      f.write(">extracted_sequence\n")
      f.write(str(sequence))

fasta_file = "temp_SUPER_11.fasta"

start_pos = 26930000
end_pos = 27118000

output_file = "extracted_sequence.fasta"

extract_sequence_to_file(fasta_file, start_pos, end_pos, output_file)
```
```bash
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/extracted_sequence.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p extracted_sequence
mummerplot -color extracted_sequence.delta

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/extracted_sequence.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/carsonella/carsonella_contigs.fa

cat extracted_sequence.fasta >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/carsonella/carsonella_contigs.fa


ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/carsonella/carsonella_contigs.fa
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2
#58772557

source package c92263ec-95e5-43eb-a527-8f1496d56f1a
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/carsonella/minimap2
samtools view -h carsonella_contigs.bam -o carsonella_contigs.sam
samtools fastq -@32 carsonella_contigs.sam > reads.fastq
#58776516


ProgDir=~/git_repos/Wrappers/NBI
Reads=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/carsonella/minimap2/reads.fastq
OutDir=$(dirname $Reads)/hifiasm_19.5.2
OutFile=carsonella_api
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_default.sh $OutDir $OutFile $Reads
#58778299
```
Coverage of integrated Carsonella genome regions:
```python
from Bio import SeqIO

def extract_sequence_to_file(fasta_file, start_pos, end_pos, output_file):
  for record in SeqIO.parse(fasta_file, "fasta"):
    sequence = record.seq[start_pos-1:end_pos]
    with open(output_file, "w") as f:
      f.write(">extracted_sequence\n")
      f.write(str(sequence))

fasta_file = "temp_SUPER_11.fasta"

start_pos = 26900000
end_pos = 27150000

output_file = "extracted_sequence2.fasta"

extract_sequence_to_file(fasta_file, start_pos, end_pos, output_file)
```
```bash
for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/extracted_sequence2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_44.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_2188.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_107.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_152.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_1132.fasta); do
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/carsonella/minimap2/$(basename $Assembly | sed 's@.fasta@@g')
Outfile=$(basename $Assembly | sed 's@.fasta@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2
done
#58971549-54

source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/extracted_sequence2.fasta ../Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -p extracted_sequence2
mummerplot -color extracted_sequence2.delta 
```
Chop out sequence
```bash
for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta); do
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/carsonella/minimap2/$(basename $Assembly | sed 's@.fasta@@g')
Outfile=$(basename $Assembly | sed 's@.fasta@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2
done
#58984158

samtools sort -o sorted.bam temp_SUPER_11.bam
samtools index sorted.bam
#Input to IGV to get breakpoints

echo -e "SUPER_11\t0\t26938864" >> temp_chop.bed
echo -e "SUPER_11\t27108941\t31362224" >> temp_chop.bed

bedtools getfasta -fi /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta -fo /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11_1.fasta -bed temp_chop.bed
#This makes seperate entries for the 2 flanking regions, edited manually.

echo -e "SUPER_11\t26938865\t27108940" > temp_chop.bed
bedtools getfasta -fi /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta -fo /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/carsonella/Carsonella_1.fasta -bed temp_chop.bed

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_rm.py --id_file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/temp.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered.fa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered2.fa

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered2.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11_1.fasta > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered.fa
```
#### FCS
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered.fa
TAXID=872318
OutDir=$(dirname $Genome)/fcs
OutFile=$(basename $Genome | sed 's@.fa@@g')
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_fcs.sh $Genome $TAXID $OutDir $OutFile
#59162147
```
### Final assembly assessment

Fresh blobplots were prepared for the filtered assembly:
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered.fa
Species=$(echo $Assembly | cut -d '/' -f10)
Genus=Trioza
TaxID=872318
alias=$(echo $Species | sed 's@_@@g' | cut -c1-4)
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T505_barcoded.fastq.gz
T2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/TellSeq/longranger/Trapi_T507_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI

#Tiara
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')
OutDir=$(dirname $Assembly)/tiara
mkdir $OutDir
sbatch $ProgDir/run_tiara.sh $Assembly $OutDir $OutPrefix
#58746508

#Alignment
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1 $T2
#58746513

#Blast
OutPrefix=$(basename $Assembly | sed 's@.fa@@g')
OutDir=$(dirname $Assembly)/blast2.12.0
Max_target=10
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
mkdir $OutDir
sbatch $ProgDir/run_blastn.sh $Assembly $Database $OutDir $OutPrefix $Max_target
#58746518

#BUSCO - keeping output files
OutDir=$(dirname $Assembly)/BUSCO
Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
OutFile=$(basename $Assembly | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
mkdir $OutDir 
sbatch $ProgDir/run_busco_keep.sh $Assembly $Database $OutDir $OutFile 
#58746520

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
#58760374

cp -r ${OutDir}/Tapi_Tellseq_trimmed_filtered_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
#The contents of /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs were subsequently copied to local machine WSL-Ubuntu: \\wsl.localhost\Ubuntu\home\did23faz
```
#### Inspector

The filtered assembly was assessed and polished with inspector:
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered.fa
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/test
Datatype=hifi
Correct_Datatype=pacbio-hifi
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_inspector.sh $OutFile $OutDir $Genome $Datatype $Correct_Datatype $Read1 $Read2
#58539171, 59062467 (repeat with carsonella removed), 1079915, 1081928

seqtk seq -A /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered.fa | awk '/^>/ {if (seq) print length(seq), id; id=$0; seq=""} {seq=seq$0} END {print length(seq), id}' | sort -n | awk '{print $2, $1}'

```
Statics of contigs:
Number of contigs       559
Number of contigs > 1000 bp     559
Number of contigs >1000000 bp   15
Total length    594055691
Total length of contigs > 1000 bp       594055691
Total length of contigs >1000000bp      580396919
Longest contig  70138210
Second longest contig length    66895101
N50     50117826
N50 of contigs >1Mbp    50117826


Read to Contig alignment:
Mapping rate /% 96.21
Split-read rate /%      23.77
Depth   39.7706
Mapping rate in large contigs /%        94.51
Split-read rate in large contigs /%     23.68
Depth in large conigs   40.001


Structural error        1340
Expansion       892
Collapse        231
Haplotype switch        180
Inversion       37


Small-scale assembly error /per Mbp     335.9634509418411
Total small-scale assembly error        199581
Base substitution       156749
Small-scale expansion   23043
Small-scale collapse    19789

QV      24.282248832154515


#### Juicebox
The final assembly was visualised using juicebox:
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
OutDir=$(dirname $Assembly)/juicebox
OutFile=$(basename $Assembly | sed 's@.fa@@g')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_juicebox-prep.sh $Assembly $Read1 $Read2 $OutDir $OutFile
#58755223, 58957934

source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
java17 -Xmx48000m -Djava.awt.headless=true -jar ~/git_repos/Scripts/NBI/juicer_tools_1.22.01.jar pre --threads 16 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicebox/*.pairs /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicebox/*.hic /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/juicebox/*.genome

python ~/git_repos/Scripts/NBI/makeAgpFromFasta.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.agp

python ~/git_repos/Scripts/NBI/agp2assembly.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.agp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.assembly

ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_HiC
Gff=NA
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_apicales/apicales-286172_S3HiC_R2.fastq-001.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_apicales/apicales-286172_S3HiC_R1.fastq-002.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Gff $Read1 $Read2 
#59554223, 59618884

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/matlock_20181227--h47b34e0_7 matlock bam2 juicer /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/bwa/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_HiC.bam out.links.txt
sort -k2,2 -k6,6 out.links.txt > out.sorted.links.txt
bash 3d-dna/visualize/run-assembly-visualizer.sh -p false in.assembly out.sorted.links.txt # creates a .hic file


singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/juicer.sif




source lastz-1.03.73;source gnu_parallel-20180322;source jdk-1.7.0_25
awk -f ~/git_repos/Scripts/NBI/generate-assembly-file-from-fasta.awk /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected-reorder.fa > tapi.assembly
~/git_repos/Scripts/NBI/run-assembly-visualizer.sh tapi.assembly /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_59644228/juicer/aligned/merged_nodups.txt

mv tapi* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/3ddna-vis/.
```
The .hic format is another type of binarized, indexed and highly-compressed file (Durand et al. (2016)). It can store virtually the same information than a .cool file. However, parsing .hic files is not as straightforward as .cool files, as it does not rely on a generic file standard. Still, the straw library has been implemented in several computing languages to facilitate parsing of .hic files (Durand et al. (2016)).
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
OutDir=$(dirname $Assembly)/3ddna-vis
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_apicales/apicales-286172_S3HiC_R2.fastq-001.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_apicales/apicales-286172_S3HiC_R1.fastq-002.gz
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_3dDNA.sh $Assembly $OutDir $OutFile $Read1 $Read2
#59644228

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_59644228/juicer/aligned/merged_nodups.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/3ddna-vis/.
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/tmp_59644228/juicer/aligned/genome_wrapped.rawchrom* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/3ddna-vis/.
```
```python
def reorder_fasta(input_file, output_file):
    sequences = {}
    current_header = None
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_header = line.strip()
                sequences[current_header] = ''
            else:
                sequences[current_header] += line.strip()
    # Find the index positions of the headers
    header_order = [">SUPER_1",">SUPER_2",">SUPER_3",">SUPER_4",">SUPER_5",">SUPER_6",">SUPER_7",">SUPER_8",">SUPER_9",">SUPER_10", ">SUPER_11_1", ">SUPER_12"]
    new_order = []
    for header in header_order:
        new_order.append((header, sequences.pop(header)))
    # Write the reordered sequences to the output file
    with open(output_file, 'w') as f:
        for header, sequence in new_order:
            f.write(header + '\n')
            f.write(sequence + '\n')
        # Write remaining sequences
        for header, sequence in sequences.items():
            f.write(header + '\n')
            f.write(sequence + '\n')

# Usage example:
input_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa"
output_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected-reorder.fa"
reorder_fasta(input_file, output_file)

```

## Annotation
### Repeatmasking

#### Repeatmodeler
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/repeatmodeler
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmodeler.sh $Genome $OutFile $OutDir
#58565331
```
#### Repeatmasker
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
Species=hemiptera
Repeat_library=$(ls $(dirname $Genome)/repeatmodeler/*/consensi.fa.classified)
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/repeatmasker
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmasker4.1.5.sh $Genome $OutFile $OutDir $Species $Repeat_library
#58575055, 59062478, 59070036 (with Carsonella removed)
```
#### Earlgrey

```bash
#Test with one scaffold:
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_SUPER_11.fasta
OutFile=$(basename $Genome | sed 's@.fasta@@g')
OutDir=$(dirname $Genome)/earlgreytest
RMsearch=arthropoda
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_earlgrey.sh $Genome $OutFile $OutDir $RMsearch
#58949705, 59062468, 59092823 (rerun to see if error related to dirs existing are generated), 59093405 (rerun)

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/earlgrey
RMsearch=arthropoda
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_earlgrey.sh $Genome $OutFile $OutDir $RMsearch
#59088751 (output missing, checkpointed but cannot complete properly), 59093551 rerun in new directory
```
```bash
pwd #/home/theaven/scratch/uncompressed/hogenhout

conda activate earlgrey
for Genome in $(ls /home/theaven/scratch/uncompressed/hogenhout/psyllid-fin/*.fa); do
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/earlgrey
RMsearch=arthropoda
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch ../../apps/earlgrey/run_earlgrey.sh $Genome $OutFile $OutDir $RMsearch
done
#22751278-80


for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/download20082024/*.fa); do
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/earlgrey
RMsearch=arthropoda
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_earlgrey.sh $Genome $OutFile $OutDir $RMsearch
done
#3215923-5
```
### RNASeq 
#### Trimmomatic
```bash
for Rawdata in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_*.fq.gz); do
OutDir=$(dirname $Rawdata)/fastqc
OutFile=$(basename $Rawdata | sed 's@.fq.gz@@g')
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_fastqc.sh $Rawdata $OutDir $OutFile
done
#58575078-83

mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_3.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_4.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_5.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_6.fq.gz
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f3 | rev)
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    Fread2=$(ls ${ReadDir}*_3.fq.gz)
    Rread2=$(ls ${ReadDir}*_4.fq.gz)
    Fread3=$(ls ${ReadDir}*_5.fq.gz)
    Rread3=$(ls ${ReadDir}*_6.fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@dna_qc@g')trim_galore
    OutFile=${sample}_trimmed
    Quality=20
    Length=50
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3 $Fread4 $Rread4 $Fread5 $Rread5
done 
#58575543,59797428

for QCdata in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/RNASeq/trim_galore/T_apicales_*.fq.gz); do
OutDir=$(dirname $QCdata)/fastqc
OutFile=$(basename $QCdata | sed 's@.fq.gz@@g')
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_fastqc.sh $QCdata $OutDir $OutFile
done
#58615839,40
```
#### HiSat2
Hisat2 is the aligner used internally by braker, it is slower than star but requires less memory.
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/RNASeq/trim_galore/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    InGenome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    OutDir=$(dirname $InGenome)/hisat2
    OutFile=$(basename $InGenome | sed 's@.fa@@g')
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_HiSat2.sh $InGenome $Fread $Rread $OutDir $OutFile
done
#58615903, 59093381 (with Carsonella removed) 
```
132952361 reads; of these:
  132952361 (100.00%) were paired; of these:
    50677098 (38.12%) aligned concordantly 0 times
    76835386 (57.79%) aligned concordantly exactly 1 time
    5439877 (4.09%) aligned concordantly >1 times
    ----
    50677098 pairs aligned concordantly 0 times; of these:
      443469 (0.88%) aligned discordantly 1 time
    ----
    50233629 pairs aligned 0 times concordantly or discordantly; of these:
      100467258 mates make up the pairs; of these:
        74589020 (74.24%) aligned 0 times
        23998644 (23.89%) aligned exactly 1 time
        1879594 (1.87%) aligned >1 times
71.95% overall alignment rate

```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/RNASeq/trim_galore/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    InGenome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/carsonella/Carsonella_1.fasta
    OutDir=$(dirname $InGenome)/hisat2
    OutFile=$(basename $InGenome | sed 's@.fasta@@g')
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_HiSat2.sh $InGenome $Fread $Rread $OutDir $OutFile
done
#59824443

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/RNASeq/trim_galore/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    InGenome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_suspected_contaminants.fa
    OutDir=$(dirname $InGenome)/hisat2
    OutFile=$(basename $InGenome | sed 's@.fa@@g')
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_HiSat2.sh $InGenome $Fread $Rread $OutDir $OutFile
done
#59824446

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/RNASeq/trim_galore/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    InGenome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
    OutDir=$(dirname $InGenome)/hisat2
    OutFile=$(basename $InGenome | sed 's@.fa@@g')
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_HiSat2.sh $InGenome $Fread $Rread $OutDir $OutFile
done
#58797870
```
***Carsonella***
132952361 reads; of these:
  132952361 (100.00%) were paired; of these:
    132824740 (99.90%) aligned concordantly 0 times
    127606 (0.10%) aligned concordantly exactly 1 time
    15 (0.00%) aligned concordantly >1 times
    ----
    132824740 pairs aligned concordantly 0 times; of these:
      128 (0.00%) aligned discordantly 1 time
    ----
    132824612 pairs aligned 0 times concordantly or discordantly; of these:
      265649224 mates make up the pairs; of these:
        265642973 (100.00%) aligned 0 times
        6008 (0.00%) aligned exactly 1 time
        243 (0.00%) aligned >1 times
0.10% overall alignment rate

***Liberibacter***
132952361 reads; of these:
  132952361 (100.00%) were paired; of these:
    132871496 (99.94%) aligned concordantly 0 times
    66897 (0.05%) aligned concordantly exactly 1 time
    13968 (0.01%) aligned concordantly >1 times
    ----
    132871496 pairs aligned concordantly 0 times; of these:
      533 (0.00%) aligned discordantly 1 time
    ----
    132870963 pairs aligned 0 times concordantly or discordantly; of these:
      265741926 mates make up the pairs; of these:
        265539165 (99.92%) aligned 0 times
        196407 (0.07%) aligned exactly 1 time
        6354 (0.00%) aligned >1 times
0.14% overall alignment rate

***Contaminants***
132952361 reads; of these:
  132952361 (100.00%) were paired; of these:
    132425079 (99.60%) aligned concordantly 0 times
    487508 (0.37%) aligned concordantly exactly 1 time
    39774 (0.03%) aligned concordantly >1 times
    ----
    132425079 pairs aligned concordantly 0 times; of these:
      2482 (0.00%) aligned discordantly 1 time
    ----
    132422597 pairs aligned 0 times concordantly or discordantly; of these:
      264845194 mates make up the pairs; of these:
        264321957 (99.80%) aligned 0 times
        473631 (0.18%) aligned exactly 1 time
        49606 (0.02%) aligned >1 times
0.60% overall alignment rate

#### STAR
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_apicales/RNASeq/trim_galore/); do
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    InGenome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    OutDir=$(dirname $InGenome)/star
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_star_align.sh $InGenome $Fread $Rread $OutDir 
done
#59434553
```

#### Braker
```bash
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda.fa > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
for file in $(ls ../Genomes/*/*/*/protein.faa | grep -v 'Frankiniella\|Thrips\|Megalurothrips\|Tribolium\|Drosophila\|Bombyx\|Apis\|Anopheles'); do
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $file | sed 's@*@@g'| sed 's@|@@g'| sed 's@ @@g'>> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
done

mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa.masked /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa


Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
Species=T_apicales
OutDir=$(dirname $Masked_Genome)/braker3
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker3.0.7.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#58761739, 59160058 (with Carsonella removed)

grep '>' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.aa | wc -l #15,345

Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
Species=T_apicales
OutDir=$(dirname $Masked_Genome)/braker2
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker2.1.6.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59430270

Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/star/star_aligmentAligned.sortedByCoord.out.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera.fa
Species=T_apicales-rna
OutDir=$(dirname $Masked_Genome)/braker1
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker2.1.6-rna.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59438021

conda activate braker
Assembly=/home/theaven/scratch/uncompressed/hogenhout/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
OutDir=/home/theaven/scratch/uncompressed/hogenhout/T_apicales_braker1
AcceptedHits=/home/theaven/scratch/uncompressed/hogenhout/T_apicales_star.bam
GeneModelName=tapi
ProgDir=/home/theaven/scratch/apps/braker
sbatch $ProgDir/braker1.sh $Assembly $OutDir $GeneModelName $AcceptedHits 
conda deactivate
#19524303

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.gff3
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff3@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein

grep '>' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.faa | wc -l
#27,875

conda activate braker
Assembly=/home/theaven/scratch/uncompressed/hogenhout/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
OutDir=/home/theaven/scratch/uncompressed/hogenhout/T_apicales_braker2
AcceptedHits=/home/theaven/scratch/uncompressed/hogenhout/Arthropoda+hemiptera.fa
GeneModelName=tapi2
ProgDir=/home/theaven/scratch/apps/braker
sbatch $ProgDir/braker2.sh $Assembly $OutDir $GeneModelName $AcceptedHits 
conda deactivate
#19521898

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/braker.gff3
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff3@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein

grep '>' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/braker.faa | wc -l
#23,256

source package c1ac247d-c4d1-4747-9817-bf03617f979b
kallisto index -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/transcriptome_index.idx /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.codingseq
kallisto quant -t 8 -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/transcriptome_index.idx -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_P9_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_P9_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_3.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_4.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_5.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_6.fq.gz
#Mapping Rate = (Mapped Reads / Total Reads) * 100

kallisto index -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/transcriptome_index.idx /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/augustus.hints.codingseq
kallisto quant -t 8 -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/transcriptome_index.idx -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_P9_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_P9_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_3.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_4.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_5.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_6.fq.gz

kallisto index -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/transcriptome_index.idx /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/augustus.hints.codingseq
kallisto quant -t 8 -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/transcriptome_index.idx -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_P9_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_P9_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_3.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_4.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_5.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_6.fq.gz
#59554876

grep -c "^@" /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_P9_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_P9_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_3.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pf_4.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_5.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/RNASeq/T_apicales_Pm_6.fq.gz #

awk 'NR > 1 {sum += $4} END {print sum}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/abundance.tsv #6.44333e+07
awk 'NR > 1 {sum += $4} END {print sum}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/abundance.tsv #6.17169e+07
awk 'NR > 1 {sum += $4} END {print sum}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/abundance.tsv #6.1324e+07

#Braker1: 48.4% pseudoaligned
#Braker2: 46.4% pseudoaligned
#Braker3: 46.1% pseudoaligned

Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+braker1.fa
Species=T_apicales
OutDir=$(dirname $Masked_Genome)/braker3+1
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker3.0.7.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59727994,59729518,59729652


Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+helixer.fa
Species=D_apicales_h
OutDir=$(dirname $Masked_Genome)/braker3+helixer
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker3.0.7.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59910583

Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+braker1+helixer.fa
Species=D_apicales_h1
OutDir=$(dirname $Masked_Genome)/braker3+helixer+1
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker3.0.7.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59910585

Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+tsebra_abinitio+helixer.fa
Species=T_apicales_h2
OutDir=$(dirname $Masked_Genome)/braker3+helixer+tsebra_abinitio
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker3.0.7.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59855998

Masked_Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
RNA_alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/hisat2/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam
Protein_database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/orthodb/11/arthropoda/Arthropoda+hemiptera+tsebra_loose+helixer.fa
Species=T_apicales_h3
OutDir=$(dirname $Masked_Genome)/braker3+helixer+tsebra_loose
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_braker3.0.7.sh $Masked_Genome $RNA_alignment $Protein_database $Species $OutDir
#59856004
```
#### TSEBRA 
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/tsebra
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/tsebra
gtf1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/augustus.hints.gtf
gff1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/hintsfile.gff
gtf2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/augustus.hints.gtf
gff2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/hintsfile.gff
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/braker3.0.7.sif tsebra.py -g ${gtf1},${gtf2} \
    -e ${gff1},${gff2} \
    -o braker1+2_combined.gtf

agat_sp_extract_sequences.pl -g braker1+2_combined.gtf -f $Genome -t cds --output ${OutDir}/T_apicales_braker1+2_combined.faa --clean_final_stop --protein

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/braker3.0.7.sif tsebra.py -g ${gtf1},${gtf2} \
    -c ~/git_repos/Scripts/NBI/tsebra_default.cfg \
    -e ${gff1},${gff2} \
    -o braker1+2_combined.gtf

agat_sp_extract_sequences.pl -g braker1+2_combined.gtf -f $Genome -t cds --output ${OutDir}/T_apicales_braker1+2_combined_default.faa --clean_final_stop --protein

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/braker3.0.7.sif tsebra.py -g ${gtf1},${gtf2} \
    -c ~/git_repos/Scripts/NBI/tsebra_braker3.cfg \
    -e ${gff1},${gff2} \
    -o braker1+2_combined.gtf

agat_sp_extract_sequences.pl -g braker1+2_combined.gtf -f $Genome -t cds --output ${OutDir}/T_apicales_braker1+2_combined_braker3.faa --clean_final_stop --protein

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/braker3.0.7.sif tsebra.py -g ${gtf1},${gtf2} \
    -c ~/git_repos/Scripts/NBI/tsebra_pref_braker1.cfg \
    -e ${gff1},${gff2} \
    -o braker1+2_combined.gtf

agat_sp_extract_sequences.pl -g braker1+2_combined.gtf -f $Genome -t cds --output ${OutDir}/T_apicales_braker1+2_combined_pref1.faa --clean_final_stop --protein

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/braker3.0.7.sif tsebra.py -g ${gtf1},${gtf2} \
    -c ~/git_repos/Scripts/NBI/tsebra_keep_ab_initio.cfg \
    -e ${gff1},${gff2} \
    -o braker1+2_combined.gtf

agat_sp_extract_sequences.pl -g braker1+2_combined.gtf -f $Genome -t cds --output ${OutDir}/T_apicales_braker1+2_combined_abinitio.faa --clean_final_stop --protein

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/braker3.0.7.sif tsebra.py -g ${gtf1},${gtf2} \
    -c temp.cfg \
    -e ${gff1},${gff2} \
    -o braker1+2_combined.gtf

agat_sp_extract_sequences.pl -g braker1+2_combined.gtf -f $Genome -t cds --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/tsebra/T_apicales_braker1+2_1-10-5-1_0-1_0.5-2-0.1-0.36.faa --clean_final_stop --protein

nano temp.cfg
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/braker3.0.7.sif tsebra.py -g ${gtf1},${gtf2} \
    -c temp.cfg \
    -e ${gff1},${gff2} \
    -o braker1+2_combined.gtf

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/braker3.0.7.sif tsebra.py -g ${gtf1},${gtf2} \
    -c temp2.cfg \
    -e ${gff1},${gff2} \
    -o braker1+2_combined2.gtf

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/braker3.0.7.sif tsebra.py -g ${gtf1},${gtf2} \
    -c temp3.cfg \
    -e ${gff1},${gff2} \
    -o braker1+2_combined3.gtf

agat_sp_extract_sequences.pl -g braker1+2_combined.gtf -f $Genome -t cds --output ${OutDir}/T_apicales_braker1+2_1-20-5-1_0-1_0-0-0.01-0.01.faa --clean_final_stop --protein

agat_sp_extract_sequences.pl -g braker1+2_combined2.gtf -f $Genome -t cds --output ${OutDir}/T_apicales_braker1+2_1-10-5-1_0.25-2_0.25-0.25-0.05-0.18.faa --clean_final_stop --protein

agat_sp_extract_sequences.pl -g braker1+2_combined3.gtf -f $Genome -t cds --output ${OutDir}/T_apicales_braker1+2_1-10-5-1_0-1_0-0.5-0.01-0.1.faa --clean_final_stop --protein


```
```bash
for file in $(ls ${OutDir}/*.faa | grep -v 'longest.faa'); do
Out=$(echo $file | sed 's@.faa@_longest.faa@g')
if [ ! -e ${Out} ]; then
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py $file $Out
fi
done

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/braker.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker1/T_apicales_braker1_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/braker.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_braker2/T_apicales_braker2_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.aa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/T_apicales_braker3_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales_helixer_longest.faa


singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py ${OutDir}/T_apicales_braker1+2_combined.faa ${OutDir}/T_apicales_braker1+2_combined_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py ${OutDir}/T_apicales_braker1+2_combined_default.faa ${OutDir}/T_apicales_braker1+2_combined_default_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py ${OutDir}/T_apicales_braker1+2_combined_braker3.faa ${OutDir}/T_apicales_braker1+2_combined_braker3_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py ${OutDir}/T_apicales_braker1+2_combined_pref1.faa ${OutDir}/T_apicales_braker1+2_combined_pref1_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py ${OutDir}/T_apicales_braker1+2_combined_abinitio.faa ${OutDir}/T_apicales_braker1+2_combined_abinitio_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+1/braker.aa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+1/T_apicales_braker3+1_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py ${OutDir}/T_apicales_braker1+2_1-10-5-1_0-1_0.5-2-0.1-0.36.faa ${OutDir}/T_apicales_braker1+2_1-10-5-1_0-1_0.5-2-0.1-0.36_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.aa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/T_apicalis_braker3+helixer+tsebra_loose_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_abinitio/braker.aa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_abinitio/T_apicalis_braker3+helixer+tsebra_abinitio_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer/braker.aa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer/T_apicalis_braker3+helixer_longest.faa

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/longest_variant.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+1/braker.aa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+1/T_apicalis_braker3+helixer+1_longest.faa
872318
```
```bash
for Proteome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/*/T_apicalis*longest.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/*/T_apicales*longest.faa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Proteome)/BUSCO
    mkdir $OutDir 
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    while [ $Jobs -gt 5 ]; do
      sleep 300s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Proteome | sed 's@.faa@@g')_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco-prot.sh $Proteome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done 
#59794569,59794582,59794586,59794592,59794728,59795108,59795919,59798082,59798157,59798181
```

#### Helixer
```bash
fasta=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
model_filepath=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/helixer/invertebrate_v0.3_m_0100/invertebrate_v0.3_m_0100.h5
lineage=invertebrate
species=T_apicales
outfile=T_apicales
outdir=$(dirname $fasta)/helixer
ProgDir=~/git_repos/Wrappers/NBI
mkdir $outdir
sbatch $ProgDir/run_helixer.sh $fasta $model_filepath $lineage $species $outfile $outdir
#58575605, 59099434 (with Carsonella removed)

grep 'gene' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.gff | wc -l #17,996
```
#### Swissprot
```bash
Proteome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.aa 
OutDir=$(dirname $Proteome)/swissprot
SwissDbDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/Uniprot/swissprot_2024_March_10
SwissDbName=uniprot_sprot
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName 
#58956716, 59168519 (with Carsonella removed)

Proteome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa
OutDir=$(dirname $Proteome)/swissprot
SwissDbDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/Uniprot/swissprot_2024_March_10
SwissDbName=uniprot_sprot
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName 
#59328912

Proteome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.aa
OutDir=$(dirname $Proteome)/swissprot
SwissDbDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/Uniprot/swissprot_2024_March_10
SwissDbName=uniprot_sprot
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName 
#63243636
```
#### Interproscan
```bash
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.aa
OutDir=$(dirname $InFile)/interproscan
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_interproscan.sh $InFile $OutDir
#58949512, 59168521 (with Carsonella removed) 

InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa
OutDir=$(dirname $InFile)/interproscan
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_interproscan.sh $InFile $OutDir
#59328913

/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.aa
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.aa
OutDir=$(dirname $InFile)/interproscan
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_interproscan.sh $InFile $OutDir
#63243648

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.aa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/tant.aa
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/braker3/braker.aa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/turt.aa
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.aa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/tapi.aa

grep 'GO:' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/interproscan/braker.aa.gff3 > GO_tapi.gff3
awk '{print $1}' GO_tapi.gff3 | sort | uniq | wc -l #11,059
```
#### Circos
```bash
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/circos
mkdir $OutDir

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/fasta-to-karyotype.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa > $OutDir/karyotype.txt

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/gc_skew_for_circos.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa 200000 200000 blue orange > $OutDir/gc_skew.txt

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/gff_to_circos.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.gff3 200000 gene $OutDir/karyotype.txt $OutDir gene

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/gff_to_circos.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/download20082024/earlgrey/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked_EarlGrey/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked_summaryFiles/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.filteredRepeats.gff 200000 all $OutDir/karyotype.txt $OutDir all_TE

ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#4177353

source package b0ed0698-358b-4c9b-9d21-603ea8d6e478
source package c92263ec-95e5-43eb-a527-8f1496d56f1a

samtools sort -o temp.bam $(basename $Assembly | sed 's@.fa@@g').bam && mv temp.bam $(basename $Assembly | sed 's@.fa@@g').bam
samtools index $OutDir/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam 
samtools depth -a $OutDir/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam > ${OutDir}/coverage.txt
samtools faidx $Assembly
cut -f1,2 ${Assembly}.fai > ${OutDir}/genome.txt
bedtools makewindows -g ${Assembly}.fai -w 200000 -s 200000 > $OutDir/genome.windows
bedtools multicov -bams $OutDir/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.bam -bed $OutDir/genome.windows > $OutDir/genome.cov.histogram
cp $OutDir/genome.cov.histogram /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/circos/.

OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/circos
grep 'SUPER' $OutDir/karyotype.txt > $OutDir/karyotype2.txt
awk '{if(NR>1) max=($4>max?$4:max)} END {print max}' $OutDir/gc_skew.txt #0.8989
awk 'NR==1 {next} {if (min=="") min=$4; if ($4<min) min=$4} END {print min}' $OutDir/gc_skew.txt #-0.9422
awk '{if(NR>1) max=($4>max?$4:max)} END {print max}' $OutDir/gene_count.tsv #56
awk '{if(NR>1) max=($4>max?$4:max)} END {print max}' $OutDir/all_TE_density.tsv #0.938
cp /hpc-home/did23faz/git_repos/temp/bands.conf $OutDir/.
cp /hpc-home/did23faz/git_repos/temp/ideogram.conf $OutDir/.
cp /hpc-home/did23faz/git_repos/temp/ideogram.label.conf $OutDir/.
cp /hpc-home/did23faz/git_repos/temp/ideogram.position.conf $OutDir/.
cp /hpc-home/did23faz/git_repos/temp/ticks.conf $OutDir/.

cd $OutDir
circos -conf /hpc-home/did23faz/git_repos/temp/apicalis_circos.conf
```
#### HGT

```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/blast
blastp -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.aa -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nr_08042024/nr -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/blast/HGT_results.out -evalue 1e-10 -outfmt 6 -num_threads 32
#944000
```














#### EarlGreyTE
```bash
source package 14fbfadb-9fe7-419a-9f20-cd5f458c0fff

source package /tgac/software/testing/bin/transposonPSI-08222010

bedtools maskfasta -fi [genome.fasta] -bed [teannotations.bed] -fo [genome.softmasked.fa] -soft

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/earlgrey4.0.6.sif earlGrey 

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
OutFile=$(basename $Genome | sed 's@.fasta@@g')
OutDir=$(dirname $Genome)/temp_18_te
RMsearch=NA
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_earlgrey.sh $Genome $OutFile $OutDir $RMsearch
#58135217,58135228

srun -p short -c 8 --mem 50G --pty bash
echo "y" | singularity exec earlgrey4.0.6.sif earlGrey -g /home/theaven/scratch/uncompressed/temp/temp_scaffold_18.fasta -s temp_scaffold_18 -o /home/theaven/scratch/uncompressed/temp/ -t 8 -d yes > output.log 2>&1
sbatch run_earlgreyte.sh /home/theaven/scratch/uncompressed/temp/temp_scaffold_18.fasta temp_scaffold_18 /home/theaven/scratch/uncompressed/temp/ NA
#15840390


Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
OutFile=$(basename $Genome | sed 's@.fasta@@g')
OutDir=$(dirname $Genome)/temp_18_RM
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmasker.sh $Genome $OutFile $OutDir
#58135895

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
OutFile=$(basename $Genome | sed 's@.fasta@@g')
OutDir=$(dirname $Genome)/temp_18_RM415
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmasker4.1.5.sh $Genome $OutFile $OutDir
#58149444



source /jic/software/staging/RCSUPPORT-2681/stagingloader 
earlGrey -g /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta -s temp_scaffold_18-2 -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_18_te -t 1 > output.log 2>&1
```





























#### Polish
Improve cv from Merqury before running MitoHiFi
Polish ideally with good coverage of illumina reads

#### MitoHiFi
```bash
Assembly=
ReferenceFasta=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/NCBI/Trant_mito_NC_038141.1.fasta
ReferenceGenebank=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_anthrisci/NCBI/Trant_mito_NC_038141.1.gb
PercentOverlap=50
Code=5
Kingdom=animal
Good_Reference=
OutDir=
OutFile=
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_mitohifi_contigs.sh $Assembly $ReferenceFasta $ReferenceGenebank $PercentOverlap $Code $Kingdom $Good_Reference $OutDir $OutFile

Assembly=
ReferenceFasta=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/NCBI/Turt_mito_NC_038113.1.fasta
ReferenceGenebank=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/NCBI/Turt_mito_NC_038113.1.gb
PercentOverlap=50
Code=5
Kingdom=animal
Good_Reference=
OutDir=
OutFile=
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_mitohifi_contigs.sh $Assembly $ReferenceFasta $ReferenceGenebank $PercentOverlap $Code $Kingdom $Good_Reference $OutDir $OutFile
```
#### Purge assembly or remove haplotigs in pretext?
```bash
source purge_dups-git05062020
split_fa Trapi_default.p_ctg.fa > Trapi_default.p_ctg.fa.split
```
#### Scaffold with tellseq reads

#### 3dDNA
#### Canu
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/apicales_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/apicales_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@canu@g')/650m
    OutFile=T_apicales_650m
    Genomesize=650m
    DataType=pacbio-hifi
    mkdir -p $OutDir
    sbatch $ProgDir/run_canu.sh $OutDir $OutFile $Genomesize $DataType $Run1 $Run2
done #57221674
``` 


















































































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