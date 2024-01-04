# Trioza urticae assembly
Contains commands run by T.Heaven in assembly of Trioza urticae

Unless stated otherwise commands were performed from the directory /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae

## Collect data
```bash
#HiFi reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz raw_data/T_urticae/HiFi/urticae_hifi-reads.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz raw_data/T_urticae/HiFi/urticae_hifi-3rdSMRTcell.fastq.gz

#HiC reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_urticae/urticae-286170_S3HiC_R1.fastq-001.gz raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_urticae/urticae-286170_S3HiC_R2.fastq-002.gz raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz

#Tellseq reads, from 2nd tellseq run:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_I1_T502.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T502_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_R1_T502.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T502_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_R2_T502.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T502_R2.fastq.gz

ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_I1_T504.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T504_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_R1_T504.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T504_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_R2_T504.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T504_R2.fastq.gz

#No reads here?:
ls /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_Dec_2021/211210_NB501793_0280_AHWCHWBGXH/

#Tellseq reads from 1st tellseq run -NOTE: these are from a different sample to the other data
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_I1_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T505_I1.fastq.gz
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_R1_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T505_R1.fastq.gz
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_R2_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T505_R2.fastq.gz

ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_I1_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T507_I1.fastq.gz
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_R1_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T507_R1.fastq.gz
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_R2_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T507_R2.fastq.gz
```
Remove HiFi reads containing adapters:
```bash
for InFile in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/*.fastq.gz); do
OutDir=$(dirname $InFile)/filtered
OutFile=$(basename $InFile | sed 's@.fastq.gz@@g')_filtered
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_HiFiAdapterFilt.sh $InFile $OutDir $OutFile
done 
#57334903,4
```
TellSeq reads converted to 10X format:
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x

for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_*_R1.fastq.gz); do
Run=$(basename $Reads | cut -d '_' -f2)
Fread=$Reads
Rread=$(echo $Reads | sed 's@_R1.fastq.gz@_R2.fastq.gz@g')
Sread=$(echo $Reads | sed 's@_R1.fastq.gz@_I1.fastq.gz@g')
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_TellSeq10XConversion.sh $Fread $Rread $Sread $OutDir $Run
done 
#57160151-4
```
TellSeq reads - remove the internal barcodes and adapters:
```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger

for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x/T502_S1_L001_R1_001.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x/T507_S1_L001_R1_001.fastq.gz); do
Fread=$Reads
Rread=$(echo $Reads | sed 's@_S1_L001_R1_001.fastq.gz@_S1_L001_R2_001.fastq.gz@g')
Run=$(basename $Reads | cut -d '_' -f1)
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger
OutFile=Turt_${Run}
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_longranger_basic.sh $Fread $Rread $Run $OutDir $OutFile
done 
#57161625-8, 57165657-8
```
#### Fastqc
```bash
for ReadFile in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/*/*.fastq.gz); do
OutDir=$(dirname $ReadFile)/fastqc
OutFile=$(basename $ReadFile | sed 's@.fastq.gz@@g')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_fastqc.sh $ReadFile $OutDir $OutFile
done
#57205220-57205235
```
#### longQC
```bash
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/*.fastq.gz); do
Datatype=pb-hifi
OutDir=$(dirname $Reads)/longqc/$(basename $Reads | cut -d '.' -f1)
OutFile=$(basename $Reads | cut -d '.' -f1)
ProgDir=~/git_repos/Wrappers/NBI
echo ${OutDir}/${OutFile}
mkdir $(dirname $Reads)/longqc
sbatch $ProgDir/run_longqc.sh $Reads $OutDir $OutFile $Datatype
done #57206533-4
```
## Hifiasm

### Default settings
Tom Mathers has performed hifiasm assembly using the HiFi reads with default settings.
```bash
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/Trurt_default.bp.p_ctg.fa assembly/genome/T_urticae/hifiasm/default/Trurt_default.bp.p_ctg.fa
```
#### Abyss
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
20013    20013    4360    6253    34674    69264    130990    89516    877047    1.009e9    Trurt_default.bp.p_ctg.fa
#### KAT
Versus HiFi reads:
```bash
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/kat_comp-main.mx assembly/genome/T_urticae/hifiasm/default/kat/.
cp /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/kat_comp-main.mx.spectra-cn.png assembly/genome/T_urticae/hifiasm/default/kat/.
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/kat_comp.stats assembly/genome/T_urticae/hifiasm/default/kat/.
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/kat_comp.dist_analysis.json assembly/genome/T_urticae/hifiasm/default/kat/.


#Repeat to compare:
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/Trurt_default.bp.p_ctg.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_hifi_reads_repeat
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-reads.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1
done #55476213

for Genome in $(ls /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/Trurt_default.bp.p_ctg.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/kat
Outfile=kat_comp_vs_hifi_reads_repeated
F1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz
R1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp.sh $OutDir $Outfile $Genome $F1 $R1
done #55497585

for Genome in $(ls /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/Trurt_default.bp.p_ctg.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/kat
Outfile=kat_comp_vs_hifi_reads_repeated2
F1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz
R1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp.sh $OutDir $Outfile $Genome $F1 $R1
done #55507928

for Genome in $(ls /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/Trurt_default.bp.p_ctg.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/kat
Outfile=kat_comp_vs_hifi_reads_repeated3
F1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz
R1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp.sh $OutDir $Outfile $Genome $F1 $R1
done #55508185

for Genome in $(ls /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/Trurt_default.bp.p_ctg.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/kat
Outfile=kat_comp_vs_hifi_reads_repeated4
F1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz
R1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp.sh $OutDir $Outfile $Genome $F1 $R1
done #55509882, 55510209

for Genome in $(ls /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/Trurt_default.bp.p_ctg.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/kat
Outfile=kat_comp_vs_hifi_reads_repeated5
F1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz
#R1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp.sh $OutDir $Outfile $Genome $F1 $R1
done #55509883, 55510199, 55516398, 55516470, 55516726, 55517981 The only way to get KAT to run properly is to only give it one read file, otherwise it keep comparing between the 2 read files rather than to the assembly...
```
Versus TellSeq reads:
Versus TellSeq reads:
```bash
cd /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt #kat cannot find the assembly file when not in the same directory???
~/git_repos/Wrappers/NBI/submit-slurm_v1.1.pl -q nbi-long -m 200000 -c 32 -t 5-00:00:00 -e -j kat_comp -i "source package 7f4fb852-f5c2-4e4b-b9a6-e648176d5543;kat comp -t 32 -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/kat/kat_comp_vs_tellseq_reads -m 31 -H 100000000 -I 100000000 '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T505_R1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T507_R1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R2.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R2.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T505_R2.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T507_R2.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_I1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_I1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T505_I1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T507_I1.fastq.gz' 'Trurt_default.bp.p_ctg.fa'"
#55350142 (13th)

~/git_repos/Wrappers/NBI/submit-slurm_v1.1.pl -q nbi-long -m 200000 -c 32 -t 5-00:00:00 -e -j kat_comp -i "source package 7f4fb852-f5c2-4e4b-b9a6-e648176d5543;kat comp -t 32 -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/kat/kat_comp_vs_tellseq_reads -m 31 -H 100000000 -I 100000000 '       ' 'Trurt_default.bp.p_ctg.fa'"
#55298674, 55340599

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/Trurt_default.bp.p_ctg.fa); do
  Jobs=$(squeue -u did23faz| grep 'kat_comp'  | wc -l)
  echo x
  while [ $Jobs -gt 0 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'kat_comp'  | wc -l)
  done
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_tellseq_reads
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R1.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R2.fastq.gz
F2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R1.fastq.gz
R2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R2.fastq.gz
F3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T505_R1.fastq.gz
R3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T505_R2.fastq.gz
F4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T507_R1.fastq.gz
R4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T507_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1 $F2 $R2 $F3 $R3 $F4 $R4
done #55352335, 55367815 - runs out of memory even with 500GB

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/Trurt_default.bp.p_ctg.fa); do
  Jobs=$(squeue -u did23faz| grep 'kat_comp'  | wc -l)
  echo x
  while [ $Jobs -gt 0 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'kat_comp'  | wc -l)
  done
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_tellseq_reads
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R1.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R2.fastq.gz
F2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R1.fastq.gz
R2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1 $F2 $R2
done #55414265 - 51mers still not run - try again 55476047 - not enough memory even with 500GB for 51mers

#Assuming that homozygous peak is the largest in the spectra with frequency of: 24x
#Homozygous peak index: 2
#CAUTION: the following estimates are based on having a clean spectra and having identified the correct homozygous peak!
#Estimated genome size: 1037.06 Mbp
#Estimated heterozygous rate: 0.16%
#xEstimated assembly completeness: 62.11%


for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/Trurt_default.bp.p_ctg.fa); do
  Jobs=$(squeue -u did23faz| grep 'kat_comp'  | wc -l)
  echo x
  while [ $Jobs -gt 0 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'kat_comp'  | wc -l)
  done
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_tellseq_reads
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R1.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R2.fastq.gz
F2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R1.fastq.gz
R2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/temp.sh $OutDir $Outfile $Genome $F1 $R1 $F2 $R2
done #55426786 - just make preparatory files
```
#### Jellyfish/kmc
```bash
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/jellyfish
Outfile=default_HiFi
Reads=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-reads.fastq.gz
Reads2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_jellyfish.sh $OutDir $Outfile $Reads $Reads2
#55426793

sed -i 's/\t/ /g' input.tsv

```
Giving default_HiFi_19mer_out.histo to GenomeScope a size of 436,626,497 is estimated (Av. read length 7,573, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_21mer_out.histo to GenomeScope a size of 423,480,356 is estimated (Av. read length 7,573, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_25mer_out.histo to GenomeScope a failed to converge.
Giving default_HiFi_31mer_out.histo to GenomeScope a failed to converge.
Giving default_HiFi_39mer_out.histo to GenomeScope a size of 345,144,661 is estimated (Av. read length 7,573, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_49mer_out.histo to GenomeScope a size of 343,452,083 is estimated (Av. read length 7,573, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_61mer_out.histo to GenomeScope a size of 775,686,720 is estimated (Av. read length 7,573, 10,000,000 max. k-mer coverage from jellyfish settings).
Giving default_HiFi_75mer_out.histo to GenomeScope a size of 798,956,485 is estimated (Av. read length 7,573, 10,000,000 max. k-mer coverage from jellyfish settings).

Full genomescope model does not capture the peak of the obseved kmers.
```R
setwd("C:/Users/did23faz/OneDrive - Norwich Bioscience Institutes/Desktop/R")
dataframe19 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/jellyfish/default_HiFi_19mer_out.histo") 
dataframe21 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/jellyfish/default_HiFi_21mer_out.histo") 
dataframe25 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/jellyfish/default_HiFi_25mer_out.histo") 
dataframe31 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/jellyfish/default_HiFi_31mer_out.histo") 

#Plot kmer distribution:
plot(dataframe19[1:200,], type="l")
plot(dataframe21[1:200,], type="l")
plot(dataframe25[1:200,], type="l")
plot(dataframe31[1:200,], type="l")

#Ignore low frequency kmers:
plot(dataframe19[4:100,], type="l")
points(dataframe19[4:100,])
#Single copy region looks to be 5 - 70, peak = 14
sum(as.numeric(dataframe19[5:42813,1]*dataframe19[5:42813,2]))/14
#1,348,240,984
sum(as.numeric(dataframe19[5:70,1]*dataframe19[5:70,2]))/14
#637,224,429

plot(dataframe21[4:100,], type="l")
points(dataframe21[4:100,])
#Single copy region looks to be 5 - 70, peak = 13
sum(as.numeric(dataframe21[5:40553,1]*dataframe21[5:40553,2]))/13
#1,431,496,506
sum(as.numeric(dataframe21[5:70,1]*dataframe21[5:70,2]))/13
#718,852,241

plot(dataframe25[4:100,], type="l")
points(dataframe25[4:100,])
#Single copy region looks to be 5 - 70, peak = 12
sum(as.numeric(dataframe25[5:36904,1]*dataframe25[5:36904,2]))/12
#1,511,706,433
sum(as.numeric(dataframe25[5:70,1]*dataframe25[5:70,2]))/12
#819,821,022

plot(dataframe31[4:100,], type="l")
points(dataframe31[4:100,])
#Single copy region looks to be 5 - 70, peak = 10
sum(as.numeric(dataframe31[5:32521,1]*dataframe31[5:32521,2]))/10
#1,754,293,609
sum(as.numeric(dataframe31[5:70,1]*dataframe31[5:70,2]))/10
#1,031,918,205

dataframe39 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/jellyfish/default_HiFi_39mer_out.histo")
dataframe49 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/jellyfish/default_HiFi_49mer_out.histo")
dataframe61 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/jellyfish/default_HiFi_61mer_out.histo")
dataframe75 <- read.table("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/jellyfish/default_HiFi_75mer_out.histo")

plot(dataframe39[1:200,], type="l")
plot(dataframe49[1:200,], type="l")
plot(dataframe61[1:200,], type="l")
plot(dataframe75[1:200,], type="l")

plot(dataframe39[9:100,], type="l")
points(dataframe39[9:100,])
#Single copy region looks to be 6 - 80, peak = 14
sum(as.numeric(dataframe39[9:100000,1]*dataframe39[9:100000,2]))/14
#1,205,868,112
sum(as.numeric(dataframe39[9:80,1]*dataframe39[9:80,2]))/14
#960,390,417

plot(dataframe49[9:100,], type="l")
points(dataframe49[9:100,])
#Single copy region looks to be 6 - 80, peak = 14
sum(as.numeric(dataframe49[9:100000,1]*dataframe49[9:100000,2]))/14
#1,235,327,529
sum(as.numeric(dataframe49[9:80,1]*dataframe49[9:80,2]))/14
#1,006,742,819

plot(dataframe61[10:100,], type="l")
points(dataframe61[10:100,])
#Single copy region looks to be 6 - 70, peak = 12
sum(as.numeric(dataframe61[10:100000,1]*dataframe61[10:100000,2]))/12
#1,442,311,829
sum(as.numeric(dataframe61[10:70,1]*dataframe61[10:70,2]))/12
#1,167,993,952

plot(dataframe75[10:100,], type="l")
points(dataframe75[10:100,])
#Single copy region looks to be 6 - 70, peak = 12
sum(as.numeric(dataframe75[10:100000,1]*dataframe75[10:100000,2]))/12
#1,446,407,637
sum(as.numeric(dataframe75[10:70,1]*dataframe75[10:70,2]))/12
#1,197,764,826
```
Setting genome size estimate based upon 21-mer estimate:
```bash
  for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm@g')/715m
    OutFile=T_urticae_715m
    Haploid_Genomesize=715m #based on 21mers
    Homozygous_coverage=40 #based upon KAT spectra (0x approaches 0)
    Min_contig=2 #default
    Purge_haplotigs_level=3 #default (strict)
    Kmer_cuttoff=5.0 #default, #increasing can improve assembly quality 
    Overlap_iterations=200 #increasing can improve assembly quality 
    Kmer_size=51 #default
    Similarity_threshold=0.75 #default
    mkdir -p $OutDir
    sbatch $ProgDir/run_hifiasm.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2
  done #55338939, 55339625

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#20704   20704   4589    4842    36154   71172   134425  832727  1.082e9 assembly/genome/T_urticae/hifiasm/715m/T_urticae_715m.bp.p_ctg.fa

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/715m/T_urticae_715m.bp.p_ctg.fa); do
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
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-reads.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1
done #55359939
```
```bash
awk '/^S/{print ">"$2;print $3}' assembly/genome/T_urticae/hifiasm/715m/T_urticae_715m.bp.p_ctg.gfa > assembly/genome/T_urticae/hifiasm/715m/T_urticae_715m.bp.p_ctg.fa
source package /tgac/software/production/bin/abyss-1.3.5
abyss-fac assembly/genome/T_urticae/hifiasm/715m/T_urticae_715m.bp.p_ctg.fa > assembly/genome/T_urticae/hifiasm/715m/abyss_report.txt
```
Setting genome size estimate based upon 31-mers:
```bash
  for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm@g')/877m
    OutFile=T_urticae_877m
    Haploid_Genomesize=877m #based on 31mers
    Homozygous_coverage=40 #based upon KAT spectra (0x approaches 0)
    Min_contig=2 #default
    Purge_haplotigs_level=3 #default (strict)
    Kmer_cuttoff=5.0 #default, #increasing can improve assembly quality 
    Overlap_iterations=200 #increasing can improve assembly quality 
    Kmer_size=51 #default
    Similarity_threshold=0.75 #default
    mkdir -p $OutDir
    sbatch $ProgDir/run_hifiasm.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2
  done #55341232

#n  n:500 L50 min N80 N50 N20 max sum
#21,066 21,066  4,683 6,253 34,557  68,405  126,932 1,000,732 1,054,000,000

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/877m/T_urticae_877m.bp.p_ctg.fa); do
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
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-reads.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1
done #55359923
```   
K-mer estimates of genome size range from 600-770-800, genomescope plots with smaller sizes do not look to capture all of the observed kmers. Running assmelby with a braod range of estimated genme sizes:
```bash
  for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi); do
  	Count=900
  	for ((i=Count; i>290; i-=10)); do
    echo "Count: $i"
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
    Haploid_Genomesize=${i}m 
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm@g')/${Haploid_Genomesize}
    OutFile=T_urticae_${Haploid_Genomesize}
    Homozygous_coverage=40 #based upon KAT spectra (0x approaches 0)
    Min_contig=2 #default
    Purge_haplotigs_level=3 #default (strict)
    Kmer_cuttoff=5.0 #default, #increasing can improve assembly quality 
    Overlap_iterations=200 #increasing can improve assembly quality 
    Kmer_size=51 #default
    Similarity_threshold=0.75 #default
    Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
  	echo x
  	while [ $Jobs -gt 3 ]; do
    	sleep 900s
    	printf "."
    	Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
  	done
  	echo ${OutDir}/$OutFile >> logs/hifilog.txt
  	mkdir -p $OutDir
    sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/hifilog.txt
    done
  done

cat assembly/genome/T_urticae/hifiasm/*/abyss_report.txt >> temp.txt

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/*/*.bp.p_ctg.fa); do
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    echo $OutFile >> logs/buscolog.txt
    sleep 30s
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
done 

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/*/*.bp.p_ctg.fa | grep -v '300m'); do
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    echo $OutFile >> logs/buscolog.txt
    sleep 30s
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/*/*.bp.p_ctg.fa); do
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    echo $OutFile >> logs/buscolog.txt
    sleep 30s
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
done

rm temp.txt
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/*);do
echo "" >> temp.txt
echo "" >> temp.txt
cat ${assembly}/*abyss_report.txt >> temp.txt
echo Arthropoda: >> temp.txt
cat ${assembly}/BUSCO/*arthropoda_odb10_short_summary.txt | grep 'C:' >> temp.txt
echo Insecta: >> temp.txt
cat ${assembly}/BUSCO/*insecta_odb10_short_summary.txt | grep 'C:' >> temp.txt
echo Hemiptera: >> temp.txt
cat ${assembly}/BUSCO/*hemiptera_odb10_short_summary.txt | grep 'C:' >> temp.txt
done

cp temp.txt Reports/urticae_assembly_report.txt
```
A genome size of 715 seems to give the most complete assembly.
```bash
  for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi); do
    Count=52
    for ((i=Count; i>8; i-=2)); do
    echo "Count: $i"
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
    Haploid_Genomesize=715m 
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm@g')/${Haploid_Genomesize}/${Homozygous_coverage}
    OutFile=T_urticae_${Haploid_Genomesize}
    Homozygous_coverage=${i} #based upon KAT spectra (0x approaches 0)
    Min_contig=2 #default
    Purge_haplotigs_level=3 #default (strict)
    Kmer_cuttoff=5.0 #default, #increasing can improve assembly quality 
    Overlap_iterations=200 #increasing can improve assembly quality 
    Kmer_size=51 #default
    Similarity_threshold=0.75 #default
    Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
    done
    echo ${OutDir}/$OutFile >> logs/hifilog.txt
    mkdir -p $OutDir
    sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/hifilog.txt
    done
  done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/715m/*/*.bp.p_ctg.fa); do
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Genome | cut -d '/' -f13)_$(echo $Database | cut -d '/' -f7)
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    echo $OutFile >> logs/buscolog.txt
    sleep 30s
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
done 

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/715m/*/*.bp.p_ctg.fa | grep -v '300m'); do
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Genome | cut -d '/' -f13)_$(echo $Database | cut -d '/' -f7)
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    echo $OutFile >> logs/buscolog.txt
    sleep 30s
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/715m/*/*.bp.p_ctg.fa); do
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Genome | cut -d '/' -f13)_$(echo $Database | cut -d '/' -f7)
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    echo $OutFile >> logs/buscolog.txt
    sleep 30s
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
done

rm temp.txt
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/715m/*/);do
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

cp temp.txt Reports/urticae_assembly_report2.txt

#24 appears to be the best assembly, however they are all very similar
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/715m/*/ | grep -v BUSCO | grep -v kat | grep -v 24);do
rm -r $assembly
done
```
All of these assemblies are worse, even the one with the same coverage setting... >:(

Have updated hifiasm from v16 to v19.5, will try with this version at different haplotig purge levels:
```bash
Haploid_Genomesize_values=("345m" "430m" "470m" "570m" "610m" "670m" "715m" "775m" "880m")
Homozygous_coverage_values=(12 30 40 50)
Purge_haplotigs_level_values=(0 1 2 3)

for Haploid_Genomesize in "${Haploid_Genomesize_values[@]}"
do
for Homozygous_coverage in "${Homozygous_coverage_values[@]}"
do
for Purge_haplotigs_level in "${Purge_haplotigs_level_values[@]}"
do
ReadDir=$(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi)
ProgDir=~/git_repos/Wrappers/NBI
Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm_19.5@g')/${Haploid_Genomesize}/${Homozygous_coverage}/${Purge_haplotigs_level}
OutFile=T_urticae_${Haploid_Genomesize}_${Homozygous_coverage}_${Purge_haplotigs_level}
Min_contig=2 #default
Kmer_cuttoff=5.0
Kmer_size=51 #default
Overlap_iterations=200 #increasing can improve assembly quality 
Similarity_threshold=0.75 #default
Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
echo x
while [ $Jobs -gt 11 ]; do
sleep 900s
printf "."
Jobs=$(squeue -u did23faz| grep 'hifiasm'  | wc -l)
done
echo ${OutDir}/$OutFile >> logs/hifilog.txt
if [ -s "${OutDir}/${OutFile}.bp.p_ctg.fa" ]; then
echo Already done for: $OutFile
else 
echo Running for: $OutFile
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/hifilog.txt
fi
done
done
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/*/*/*/*.bp.p_ctg.fa); do
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
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/*/*/*);do
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

cp temp.txt Reports/urticae_assembly_report3.txt

#Once again results are very similar accross assemblies, T_urticae_570m_12_3.bp.p_ctg.fa, is probably best overall with highest N50 of those with highest buscos, although 25992 contigs - peak_hom: 32; peak_het: 14
```
```bash
Haploid_Genomesize_values=("570m")
Homozygous_coverage_values=(12 32)
Purge_haplotigs_level_values=(3)
Kmer_cuttoff_values=("3.0" "4.0" "5.0" "10.0")
Similarity_threshold_values=("0.75" "0.5" "0.25")

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
ReadDir=$(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi)
ProgDir=~/git_repos/Wrappers/NBI
Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm_19.5@g')/${Haploid_Genomesize}/${Homozygous_coverage}/${Purge_haplotigs_level}/${Kmer_cuttoff}/${Similarity_threshold}
OutFile=T_urticae_${Haploid_Genomesize}_${Homozygous_coverage}_${Purge_haplotigs_level}_${Kmer_cuttoff}_${Similarity_threshold}
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
echo ${OutDir}/$OutFile >> logs/hifilog.txt
if [ -s "${OutDir}/${OutFile}.bp.p_ctg.fa" ]; then
echo Already done for: $OutFile
else 
echo Running for: $OutFile
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/hifilog.txt
fi
done
done
done
done
done

Haploid_Genomesize_values=("610m")
Homozygous_coverage_values=(12 32)
Purge_haplotigs_level_values=(1)
Kmer_cuttoff_values=("3.0" "4.0" "5.0" "10.0")
Similarity_threshold_values=("0.75" "0.5" "0.25")

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
ReadDir=$(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi)
ProgDir=~/git_repos/Wrappers/NBI
Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm_19.5@g')/${Haploid_Genomesize}/${Homozygous_coverage}/${Purge_haplotigs_level}/${Kmer_cuttoff}/${Similarity_threshold}
OutFile=T_urticae_${Haploid_Genomesize}_${Homozygous_coverage}_${Purge_haplotigs_level}_${Kmer_cuttoff}_${Similarity_threshold}
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
echo ${OutDir}/$OutFile >> logs/hifilog.txt
if [ -s "${OutDir}/${OutFile}.bp.p_ctg.fa" ]; then
echo Already done for: $OutFile
else 
echo Running for: $OutFile
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/hifilog.txt
fi
done
done
done
done
done

Haploid_Genomesize_values=("715m")
Homozygous_coverage_values=(12 32)
Purge_haplotigs_level_values=(2)
Kmer_cuttoff_values=("3.0" "4.0" "5.0" "10.0")
Similarity_threshold_values=("0.75" "0.5" "0.25")

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
ReadDir=$(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi)
ProgDir=~/git_repos/Wrappers/NBI
Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@hifiasm_19.5@g')/${Haploid_Genomesize}/${Homozygous_coverage}/${Purge_haplotigs_level}/${Kmer_cuttoff}/${Similarity_threshold}
OutFile=T_urticae_${Haploid_Genomesize}_${Homozygous_coverage}_${Purge_haplotigs_level}_${Kmer_cuttoff}_${Similarity_threshold}
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
echo ${OutDir}/$OutFile >> logs/hifilog.txt
if [ -s "${OutDir}/${OutFile}.bp.p_ctg.fa" ]; then
echo Already done for: $OutFile
else 
echo Running for: $OutFile
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_fa_only.sh $OutDir $OutFile $Haploid_Genomesize $Homozygous_coverage $Min_contig $Purge_haplotigs_level $Kmer_cuttoff $Overlap_iterations $Kmer_size $Similarity_threshold $Run1 $Run2 2>&1 >> logs/hifilog.txt
fi
done
done
done
done
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/4.0/0.75/*.bp.p_ctg.fa); do
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
for assembly in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/*/*/*/*/0.*); do
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

cp temp.txt Reports/urticae_assembly_report3.txt
```
The most complete assembly is T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa with 87.6% of Hemiptera BUSCOs across 26,089 contigs, whilst T_urticae_570m_32_3_10.0_0.25.bp.p_ctg.fa was most contiguous with 85.4% of BUSCOs across 18,726 contigs. Using hifiasm with default settings produced an assembly with 84.7% of BUSCOs across 20,013 contigs (Trurt_default.bp.p_ctg.fa).

n       n:500   n:N50   min     N80     N50     N20     max     sum
26089   26089   5889    3103    31608   62783   117571  753257  1.221e9 T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa

Based on Merqury/KAT plot the true homozygous coverage is ~32
```bash
#These assemblies were kept:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/Trurt_default.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/570m/32/3/10.0/0.25/T_urticae_570m_32_3_10.0_0.25.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/32/2/3.0/0.75/T_urticae_715m_32_2_3.0_0.75.bp.p_ctg.fa

#All other assemblies were deleted - 278 in total:
for assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifi*/*m/*.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifi*/*m/*/*.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifi*/*m/*/*/*.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifi*/*m/*/*/*/*/*.fa | grep -v 'Trurt_default.bp.p_ctg.fa\|T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa\|T_urticae_570m_32_3_10.0_0.25.bp.p_ctg.fa\|T_urticae_715m_32_2_3.0_0.75.bp.p_ctg.fa'); do
rm $assembly  
done
```
#### KAT
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/570m/32/3/10.0/0.25/T_urticae_570m_32_3_10.0_0.25.bp.p_ctg.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_hifi_reads
F1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz
R1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp.sh $OutDir $Outfile $Genome $F1 $R1
done #56145393, 56145394

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/570m/32/3/10.0/0.25/T_urticae_570m_32_3_10.0_0.25.bp.p_ctg.fa); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Genome)/kat
Outfile=kat_comp_vs_tellseq_reads
F1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R1.fastq.gz
R1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_kat_comp_paired.sh $OutDir $Outfile $Genome $F1 $R1
done #56145395, 56145396
```
#### Merqury
Kmer plots with all reads; Hifi, HiC, Tellseq
```bash
#Prepare meryl kmer counts
source /nbi/software/staging/RCSUPPORT-2452/stagingloader
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/*/*.fastq.gz); do
OutDir=$(dirname $Reads)/meryl
mkdir $OutDir
meryl k=21 count output ${OutDir}/$(basename $Reads | sed 's@.fastq.gz@@g').meryl $Reads 
done
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
mkdir $(dirname $Assembly)/meryl
meryl union-sum output $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/*/meryl/*.meryl
mkdir $(dirname $Assembly)/meryl/HiFi 
meryl union-sum output $(dirname $Assembly)/meryl/HiFi/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/meryl/*.meryl
mkdir $(dirname $Assembly)/meryl/HiC
meryl union-sum output $(dirname $Assembly)/meryl/HiC/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/meryl/*.meryl
mkdir $(dirname $Assembly)/meryl/Tellseq
meryl union-sum output $(dirname $Assembly)/meryl/Tellseq/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/meryl/*.meryl
#57139543

#Plot with merqury
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
mkdir $(dirname $Assembly)/merqury/
cd $(dirname $Assembly)/merqury/
ln -s $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl
ln -s $Assembly
meryl histogram $(dirname $Assembly)/meryl/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl > $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl.hist
merqury.sh $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl/ $(basename $Assembly) $(basename $Assembly | sed 's@.bp.p_ctg.fa@@g') #NOTE: needs >25GB memory to complete, errors will output to log/ file not to slurm.err
#57136861

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
#57140338

source /nbi/software/staging/RCSUPPORT-2452/stagingloader
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/meryl/HiC/T_urticae_715m_12_2_3.0_0.5.meryl
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/meryl/HiFi/T_urticae_715m_12_2_3.0_0.5.meryl
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/meryl/Tellseq/T_urticae_715m_12_2_3.0_0.5.meryl
rm -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/meryl/T_urticae_715m_12_2_3.0_0.5.meryl

######################################################################################################################################################################
source /nbi/software/staging/RCSUPPORT-2452/stagingloader
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_*/TellSeq/longranger/Turt_T502*.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_*/TellSeq/longranger/Turt_T504*.fastq.gz); do
OutDir=$(dirname $Reads)/meryl
mkdir $OutDir
meryl k=21 count output ${OutDir}/$(basename $Reads | sed 's@.fastq.gz@@g').meryl $Reads 
done
#57186506

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
mkdir $(dirname $Assembly)/meryl/Tellseq_trimmed
meryl union-sum output $(dirname $Assembly)/meryl/Tellseq_trimmed/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').meryl /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/meryl/*.meryl

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
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi temp/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC temp/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger temp/.
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
#>:( 57216689
```
```bash
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.fa
for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.fa); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done
#57305004
rm /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.fa
```
#### Blobtools
Blobtools with Hifi reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2 #56939385
sbatch $ProgDir/run_qualimap.sh $(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam $Assembly $OutDir #57111016

#Blast
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/blast2.12.0/3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blastn.sh $Assembly $Database $OutDir $OutPrefix 
#57175072, 57182844, 57308570, 57398068, 57044985
#parse_seqids is required for -taxid_map however when running databases created with -parse_seqids blast runs fail with c++ errors

#Blobtools
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf taxdump.tar.gz -C . nodes.dmp names.dmp
./blobtools nodesdb --nodes nodes.dmp --names names.dmp

#BUSCO - keeping output files
for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Assembly)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Assembly | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Assembly $Database $OutDir $OutFile 
done 
#57207203, 57307067, 57307078, 57307098

#Diamond blast - BUSCO regions
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
BUSCOgff=$(dirname $Assembly)/BUSCO/hemiptera_odb10/run_hemiptera_odb10/metaeuk_output/rerun_results/T_urticae_715m_12_2_3_hemiptera_odb10.fa.gff
awk '$3 == "gene" {print $0}' $BUSCOgff > $(echo $BUSCOgff | sed 's@.fa.gff@_genes_only.fa.gff@g')
bedtools getfasta -fi $Assembly -bed $(echo $BUSCOgff | sed 's@.fa.gff@_genes_only.fa.gff@g') -fo $(dirname $BUSCOgff)/busco_regions.fasta
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/uniprot_01102023/Uniprot_01102023_reference_proteomes.dmnd
ProgDir=~/git_repos/Wrappers/NBI
mkdir -p $OutDir
sbatch $ProgDir/run_diamond_blastx.sh $(dirname $BUSCOgff)/busco_regions.fasta $Database $OutDir $OutPrefix 
#57307103, 57315733

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
#57307104, 57315723

BUSCODiamond=$(dirname $Assembly)/diamond0.9.29_blastx/BUSCO_regions/*.diamondblastx.out
ElseDiamond=$(dirname $Assembly)/diamond0.9.29_blastx/nonBUSCO_regions/*.diamondblastx.out
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/normalise_blast.py $BUSCODiamond $(echo $BUSCODiamond | sed 's@x.out@x_2.out@g')
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/normalise_blast.py $ElseDiamond $(echo $ElseDiamond | sed 's@x.out@x_2.out@g')

#Tiara
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
OutDir=$(dirname $Assembly)/tiara
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_tiara.sh $Assembly $OutDir $OutPrefix
#57401224

#Blobtools
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
MappingFile=$(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
BlastFile=$(dirname $Assembly)/blast2.12.0/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').vs.nt.mts1.hsp1.1e25.megablast.out
OutDir=$(dirname $Assembly)/blobtools1.1.1
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
ColourFile=NA
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobplot.sh $Assembly $MappingFile $BlastFile $OutDir $OutPrefix $ColourFile
#57194209

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
Species=urticae
TaxID=121826
alias=Turt1223005
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57217567

cp -r $OutDir/Turt1223005_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_HiC
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#57207041
MappingFile=$(ls ${OutDir}/*_HiC.bam)
sbatch $ProgDir/run_qualimap.sh $MappingFile $Assembly $OutDir #57111447

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
Species=urticae
TaxID=121826
alias=Turt1223005_HiC
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57218296

cp -r $OutDir/Turt1223005_HiC_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
```
Blobtools with TellSeq reads
```bash
#alignment of reads to unfiltered assembly
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T505_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T505_R2.fastq.gz
Read3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T507_R1.fastq.gz
Read4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T507_R2.fastq.gz
Read5=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R1.fastq.gz
Read6=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R2.fastq.gz
Read7=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R1.fastq.gz
Read8=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Read1 $Read2 $Read3 $Read4 $Read5 $Read6 $Read7 $Read8 
#57207046, 57110411, 57160169
MappingFile=$(ls ${OutDir}/*_Tellseq.bam)
sbatch $ProgDir/run_qualimap.sh $MappingFile $Assembly $OutDir #57186438

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
Species=urticae
TaxID=121826
alias=Turt1223005_Tellseq
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57218782

#########################################################################################################################
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T502_barcoded.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T504_barcoded.fastq.gz
Read3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T505_barcoded.fastq.gz
Read4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T507_barcoded.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $Read1 $Read2 $Read3 $Read4 
#57186559
MappingFile=$(ls ${OutDir}/*_Tellseq_trimmed.bam)
sbatch $ProgDir/run_qualimap.sh $MappingFile $Assembly $OutDir #57317732

#Blobtoolkit
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
Species=urticae
TaxID=121826
alias=Turt1223005_Tellseq_trimmed
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_blobtoolkit4.2.1.sh $Assembly $Record_type $MappingFile $BlastFile $BUSCOFile $BUSCODiamond $ElseDiamond $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#57334957

cp -r $OutDir/Turt1223005_Tellseq_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
cp -r $OutDir/Turt1223005_Tellseq_trimmed_blobdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blobtools/BlobDirs/.
```
#### Kraken 
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_kraken2nt
OutDir=$(dirname $Assembly)/kraken2.1.3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/kraken/nt_14092023
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_kraken2.sh $Assembly $Database $OutDir $OutPrefix 2>&1 >> ${OutDir}/log.txt
#57115815, 57115857
```
#### Filtering
##### Kraken
Contigs classified to non-Arthropoda phylum level were considered contamininants. Contigs classified to Arthropoda taxa or above (eg. Eukaryota) were considered Psyllid contigs.
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
There is a set of contigs that are identified as bacteria by Tiara and group together in blobtoolkit, these do not have bacterial blast identidy
```bash
source package /nbi/software/testing/bin/seqtk-1.2
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
cd $(dirname $Assembly)/tiara/
nano contaminantcontignames.txt #edit with contig names from tiara cluster
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/blobtoolkit4.2.1/T_urticae_715m_12_2_3.0_0.5_btkfiltered2.fa
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
T_urticae_715m_12_2_3_hemiptera_odb10_short_summary.txt
        C:87.6%[S:32.7%,D:54.9%],F:4.6%,M:7.8%,n:2510

T_urticae_715m_12_2_3.0_0.5_krakenfiltered_hemiptera_odb10_short_summary.txt
        C:87.1%[S:33.4%,D:53.7%],F:4.6%,M:8.3%,n:2510
T_urticae_715m_12_2_3.0_0.5_blastfiltered_hemiptera_odb10_short_summary.txt
        C:87.6%[S:33.1%,D:54.5%],F:4.6%,M:7.8%,n:2510
T_urticae_715m_12_2_3.0_0.5_tiarafiltered_hemiptera_odb10_short_summary.txt
        C:87.6%[S:33.0%,D:54.6%],F:4.6%,M:7.8%,n:2510

T_urticae_715m_12_2_3.0_0.5_btkfiltered_hemiptera_odb10_short_summary.txt
        C:87.6%[S:32.7%,D:54.9%],F:4.6%,M:7.8%,n:2510
T_urticae_715m_12_2_3.0_0.5_filtered_hemiptera_odb10_short_summary.txt
        C:87.1%[S:33.3%,D:54.2%],F:4.6%,M:7.9%,n:2510

T_urticae_715m_12_2_3.0_0.5_btkfiltered2_hemiptera_odb10_short_summary.txt
        C:87.6%[S:32.7%,D:54.9%],F:4.6%,M:7.8%,n:2510
T_urticae_715m_12_2_3.0_0.5_allfiltered_hemiptera_odb10_short_summary.txt
        C:87.0%[S:33.7%,D:53.3%],F:4.6%,M:8.4%,n:2510
```bash
#BUSCO - keeping output files
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done 
#57044988

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done 
#57045017

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
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
cat $(dirname $Assembly)/filtered/contaminantcontignames.txt | wc -l #759

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
T_urticae_715m_12_2_3.0_0.5_filtered_hemiptera_odb10_short_summary
  C:87.6%[S:32.7%,D:54.9%],F:4.5%,M:7.9%,n:2510

### Purge Dups
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
#Purge with Tellseq (Illumina) Reads:
MappingFile=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57138576

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#18864   18864   3796    3103    28042   63188   128064  753257  831.4e6 T_urticae_715m_12_2_3.0_0.5_TellSeqPurged.fa

#Purge with HiFi reads:
MappingFile=$(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
Type=long
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57344837

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#16842   16842   3812    5898    35223   69271   129722  753257  867.5e6 T_urticae_715m_12_2_3.0_0.5_HiFiPurged.fa

######################################################################################################################

#Try Purge_Haplotigs
#TellSeq
interactive
source /nbi/software/staging/RCSUPPORT-2652/stagingloader
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
AssemblyIndex=$(echo $Assembly).fai
MappingFile=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_sorted_markdups_Tellseq_trimmed.bam
purge_haplotigs hist -b $MappingFile -g $Assembly -t 1
#Assess plot for low,mid,high cuttoffs
MappingIndex=$(echo $MappingFile).bai
low=10
mid=95
high=220
diploid_cuttoff=80
junk_cuttoff=80
hap_percent_cuttoff=70
rep_percent_cuttoff=250
OutDir=$(dirname $Assembly)/purge_haplotigs/
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_haplotigs.sh $Assembly $AssemblyIndex $MappingFile $MappingIndex $low $mid $high $diploid_cuttoff $junk_cuttoff $hap_percent_cuttoff $rep_percent_cuttoff $OutDir $OutPrefix
#57373614

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#12750   12750   2960    5898    41312   78619   138007  753257  733.3e6 T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_curated.fasta

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_haplotigs/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_curated.fasta); do
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
#57381722

#T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_curated_hemiptera_odb10_short_summary.txt
#        C:84.9%[S:76.4%,D:8.5%],F:5.6%,M:9.5%,n:2510

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_haplotigs/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_curated.fasta); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done
#

#HiFi
interactive
source /nbi/software/staging/RCSUPPORT-2652/stagingloader
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
AssemblyIndex=$(echo $Assembly).fai
MappingFile=$(dirname $Assembly)/minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_sorted_markdups_HiFi.bam
purge_haplotigs hist -b $MappingFile -g $Assembly -t 1
#Assess plot for low,mid,high cuttoffs
MappingIndex=$(echo $MappingFile).bai
low=3
mid=33
high=85
diploid_cuttoff=80
junk_cuttoff=80
hap_percent_cuttoff=70
rep_percent_cuttoff=250
OutDir=$(dirname $Assembly)/purge_haplotigs
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_haplotigs.sh $Assembly $AssemblyIndex $MappingFile $MappingIndex $low $mid $high $diploid_cuttoff $junk_cuttoff $hap_percent_cuttoff $rep_percent_cuttoff $OutDir $OutPrefix
#57364802, 57368497

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#12753   12753   2960    5898    41286   78608   138007  753257  733.2e6 T_urticae_715m_12_2_3.0_0.5_HiFiPurged_curated.fasta

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_haplotigs/T_urticae_715m_12_2_3.0_0.5_HiFiPurged_curated.fasta); do
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
#57372000

#T_urticae_715m_12_2_3.0_0.5_HiFiPurged_curated_hemiptera_odb10_short_summary.txt
#        C:85.0%[S:76.5%,D:8.5%],F:5.6%,M:9.4%,n:2510

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_haplotigs/T_urticae_715m_12_2_3.0_0.5_HiFiPurged_curated.fasta); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done
#57372815

minimap2 -t 4 -ax map-hifi /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-3rdSMRTcell.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/urticae_hifi-reads.fastq.gz --secondary=no \
    | samtools sort -m 1G -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/minimap2/T_urticae_715m_12_2_3.0_0.5_aligned.bam -T tmp.ali

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/minimap2/
samtools faidx T_urticae_715m_12_2_3.0_0.5_aligned.bam
#57362837,57372005

interactive
source /nbi/software/staging/RCSUPPORT-2652/stagingloader
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
AssemblyIndex=$(echo $Assembly).fai
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/minimap2/T_urticae_715m_12_2_3.0_0.5_aligned.bam
purge_haplotigs hist -b $MappingFile -g $Assembly -t 1
#Assess plot for low,mid,high cuttoffs
MappingIndex=$(echo $MappingFile).bai
low=3
mid=33
high=85
diploid_cuttoff=80
junk_cuttoff=80
hap_percent_cuttoff=70
rep_percent_cuttoff=250
OutDir=$(dirname $Assembly)/purge_haplotigs
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_HiFiPurged_2
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_haplotigs.sh $Assembly $AssemblyIndex $MappingFile $MappingIndex $low $mid $high $diploid_cuttoff $junk_cuttoff $hap_percent_cuttoff $rep_percent_cuttoff $OutDir $OutPrefix
#57373613

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#12753   12753   2960    5898    41304   78611   138007  753257  733.2e6 T_urticae_715m_12_2_3.0_0.5_HiFiPurged_2_curated.fasta

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_haplotigs/T_urticae_715m_12_2_3.0_0.5_HiFiPurged_2_curated.fasta); do
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
#57381723

#T_urticae_715m_12_2_3.0_0.5_HiFiPurged_2_curated_hemiptera_odb10_short_summary.txt
#        C:85.0%[S:76.5%,D:8.5%],F:5.6%,M:9.4%,n:2510

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_haplotigs/T_urticae_715m_12_2_3.0_0.5_HiFiPurged_2_curated.fasta); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done
#57382079
```
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_*Purged.fa); do
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
T_urticae_715m_12_2_3.0_0.5_filtered_hemiptera_odb10_short_summary
  C:87.6%[S:32.7%,D:54.9%],F:4.5%,M:7.9%,n:2510

T_urticae_715m_12_2_3.0_0.5_HiFiPurged_arthropoda_odb10_short_summary.txt
        C:84.1%[S:64.3%,D:19.8%],F:9.3%,M:6.6%,n:1013
T_urticae_715m_12_2_3.0_0.5_HiFiPurged_hemiptera_odb10_short_summary.txt
        C:85.9%[S:67.5%,D:18.4%],F:5.3%,M:8.8%,n:2510
T_urticae_715m_12_2_3.0_0.5_HiFiPurged_insecta_odb10_short_summary.txt
        C:83.5%[S:64.4%,D:19.1%],F:9.4%,M:7.1%,n:1367
T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_arthropoda_odb10_short_summary.txt
        C:83.3%[S:71.8%,D:11.5%],F:9.4%,M:7.3%,n:1013
T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_hemiptera_odb10_short_summary.txt
        C:85.5%[S:74.3%,D:11.2%],F:5.6%,M:8.9%,n:2510
T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_insecta_odb10_short_summary.txt
        C:82.4%[S:70.3%,D:12.1%],F:10.1%,M:7.5%,n:1367

Purging with TellSeq reads seems to be more effective than using HiFi reads, duplication is still high
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/tellseq
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/tellseq/.
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/hifi
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_HiFiPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/hifi/.

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/*/*Purged.fa); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done #57291746,7
```
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/T_urticae_715m_12_2_3.0_0.5_filtered.fa
#Purge with Tellseq (Illumina) Reads:
MappingFile=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57138700

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#18364   18364   3715    3103    28332   63809   128086  753257  816.3e6 T_urticae_715m_12_2_3.0_0.5_filtered_TellSeqPurged.fa

#Purge with HiFi reads:
MappingFile=$(dirname $Assembly)/../minimap2/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g').bam
Type=long
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57051889

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#18364   18364   3715    3103    28332   63809   128086  753257  816.3e6 T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa
```
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/T_urticae_715m_12_2_3.0_0.5_*Purged.fa); do
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
T_urticae_715m_12_2_3.0_0.5_filtered_hemiptera_odb10_short_summary
  C:87.6%[S:32.7%,D:54.9%],F:4.5%,M:7.9%,n:2510

T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_arthropoda_odb10_short_summary.txt
        C:83.4%[S:71.9%,D:11.5%],F:9.2%,M:7.4%,n:1013
T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_hemiptera_odb10_short_summary.txt
        C:85.5%[S:74.5%,D:11.0%],F:5.5%,M:9.0%,n:2510
T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_insecta_odb10_short_summary.txt
        C:82.4%[S:70.4%,D:12.0%],F:9.9%,M:7.7%,n:1367
T_urticae_715m_12_2_3.0_0.5_filtered_TellSeqPurged_arthropoda_odb10_short_summary.txt
        C:83.4%[S:71.9%,D:11.5%],F:9.2%,M:7.4%,n:1013
T_urticae_715m_12_2_3.0_0.5_filtered_TellSeqPurged_hemiptera_odb10_short_summary.txt
        C:85.5%[S:74.5%,D:11.0%],F:5.5%,M:9.0%,n:2510
T_urticae_715m_12_2_3.0_0.5_filtered_TellSeqPurged_insecta_odb10_short_summary.txt
        C:82.4%[S:70.4%,D:12.0%],F:9.9%,M:7.7%,n:1367
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/tellseq
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/T_urticae_715m_12_2_3.0_0.5_filtered_TellSeqPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/tellseq/.
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/hifi
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/hifi/.

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/*/*Purged.fa); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done #57291748,9
```

Purging with/without contaminants seems to make no difference to the final scores, neither does using Tellseq/Hifi reads.

Purge_dups and purge_haplotigs:
```bash
interactive
source /nbi/software/staging/RCSUPPORT-2652/stagingloader
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa
AssemblyIndex=$(echo $Assembly).fai
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/minimap2/T_urticae_715m_12_2_3.0_0.5_aligned.bam
purge_haplotigs hist -b $MappingFile -g $Assembly -t 1
#Assess plot for low,mid,high cuttoffs
MappingIndex=$(echo $MappingFile).bai
low=3
mid=33
high=85
diploid_cuttoff=80
junk_cuttoff=80
hap_percent_cuttoff=70
rep_percent_cuttoff=250
OutDir=$(dirname $Assembly)/purge_haplotigs
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_haplotigs.sh $Assembly $AssemblyIndex $MappingFile $MappingIndex $low $mid $high $diploid_cuttoff $junk_cuttoff $hap_percent_cuttoff $rep_percent_cuttoff $OutDir $OutPrefix
#57449666

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#12489   12489   2820    5898    37959   74436   136409  753257  682.4e6 T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged_curated.fasta

low=6
mid=33
high=85
diploid_cuttoff=80
junk_cuttoff=80
hap_percent_cuttoff=70
rep_percent_cuttoff=250
OutDir=$(dirname $Assembly)/purge_haplotigs
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_HiFiPurged2
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_haplotigs.sh $Assembly $AssemblyIndex $MappingFile $MappingIndex $low $mid $high $diploid_cuttoff $junk_cuttoff $hap_percent_cuttoff $rep_percent_cuttoff $OutDir $OutPrefix
#57449669

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#12430   12430   2814    5898    38083   74522   136409  753257  681.5e6 T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged2_curated.fasta

low=6
mid=33
high=75
diploid_cuttoff=50
junk_cuttoff=80
hap_percent_cuttoff=70
rep_percent_cuttoff=250
OutDir=$(dirname $Assembly)/purge_haplotigs
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_HiFiPurged3
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_haplotigs.sh $Assembly $AssemblyIndex $MappingFile $MappingIndex $low $mid $high $diploid_cuttoff $junk_cuttoff $hap_percent_cuttoff $rep_percent_cuttoff $OutDir $OutPrefix
#57449683

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#12491   12491   2814    5898    38123   74836   137726  753257  687.1e6 T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged3_curated.fasta

low=6
mid=33
high=75
diploid_cuttoff=50
junk_cuttoff=80
hap_percent_cuttoff=50
rep_percent_cuttoff=250
OutDir=$(dirname $Assembly)/purge_haplotigs
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')_HiFiPurged4
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_haplotigs.sh $Assembly $AssemblyIndex $MappingFile $MappingIndex $low $mid $high $diploid_cuttoff $junk_cuttoff $hap_percent_cuttoff $rep_percent_cuttoff $OutDir $OutPrefix
#57449694

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#11200   11200   2565    5898    40481   78255   140412  753257  642.7e6 T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged4_curated.fasta

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged*curated.fasta); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    OutFile=$(basename $Genome | cut -d '.' -f1,2,3,4)_$(echo $Database | cut -d '/' -f7)
    if [ ! -e ${OutDir}/${OutFile}_short_summary.txt ]; then
    echo Running BUSCO for: $OutFile
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 
    sleep 30s
    else 
    echo Already done for: $OutFile
    fi
done
#57456278, 57456310, 57456322, 57456351

#T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged_curated_hemiptera_odb10_short_summary.txt
#        C:84.6%[S:80.9%,D:3.7%],F:5.8%,M:9.6%,n:2510

#T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged2_curated_hemiptera_odb10_short_summary.txt
#        C:84.5%[S:80.8%,D:3.7%],F:5.8%,M:9.7%,n:2510

#T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged3_curated_hemiptera_odb10_short_summary.txt
#        C:84.5%[S:80.8%,D:3.7%],F:5.8%,M:9.7%,n:2510

#T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged4_curated_hemiptera_odb10_short_summary.txt
#        C:82.9%[S:80.3%,D:2.6%],F:5.8%,M:11.3%,n:2510
```
#### Check HiC match assemblies
```bash
for HiCreads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_*/HiC/*.fastq.gz); do
OutDir=$(dirname $HiCreads)
Assembly1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged.fa
Assembly2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/purge_dups/T_apicales_880m_29_3_3.0_0.75_TellSeqPurged.fa
Assembly3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/purge_dups/T_anthrisci_820m_48_1_10.0_0.25_TellSeqPurged.fa
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_mash.sh $HiCreads $OutDir $Assembly1 $Assembly2 $Assembly3
done
#57308291-6

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_*/HiC/*_mashscreen.txt); do
echo $(basename $file | sed 's@_mashscreen.txt@@g').fastq.gz : >> tempmash.txt
cat $file >> tempmash.txt
done
```

### Scaffolding
Tellseq reads have are linked via 18bp barcodes, Tom Mathers has already used the conversion software prodived by Universal sequencing and the 4Mwith-alts-february-2016.txt barcode whitelist file to convert to 16bp barcoded versions of the reads that are compatible with 10x genomics linked read software. The files also have to be in a standardised naming format https://cdn.shopify.com/s/files/1/0654/0378/1341/files/100027-USG_TELL-Seq_Software_Roadmap_User_Guide_v1.0_4d2abbf8-3d19-4899-84c2-9a79594eb7f3.pdf?v=1684422029

Scaff10X can now be used to scaffold the HiFi assembly contigs using the TellSeq reads.

Longranger is run in order to remove the barcodes from the reads so that they can simly be treated as large insert illumina reads and used for assembly polishing.

Purging should be done before scaffolding as presence of haplotigs will confuse scaffolder, however the order of pilon,scaff10x,break10x,YAHS or whether using at all will improve the assembly is unclear.

#### 3DDNA
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
OutDir=$(dirname $Assembly)/3ddna
OutFile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_3dDNA.sh $Assembly $OutDir $OutFile $Read1 $Read2
#57291096
#57301681

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged.fa
OutDir=$(dirname $Assembly)/3ddna
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_3dDNA.sh $Assembly $OutDir $OutFile $Read1 $Read2
#57548996, 57615776
#NOTE: 3ddna output is very large ~600GB therefore only final files kept

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#8659    8659    1       977     730168  612.6e6 612.6e6 612.6e6 780.5e6 genome_wrapped.FINAL.fasta
#One enourmous chromosome...

#######################################################################################################################################################################
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged_curated.fasta
OutDir=$(dirname $Assembly)/3ddna
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_3dDNA.sh $Assembly $OutDir $OutFile $Read1 $Read2
#57624373

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#5266    5266    17      977     2872641 12.03e6 23.58e6 37.04e6 659.1e6 T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged_curated.FINAL.fasta
```

#### 0
Scaffolding without purging

##### Scaff10X
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234149

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#24239   24239   4047    3103    31608   65011   212852  1144978 1.221e9 T_urticae_715m_12_2_3.0_0.5_scaff10xscaffolds.fasta
```
##### YAHS
Generate mapping file
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57119211

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57276768, 57278888, 57278931

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#28653   28653   6421    1000    29777   57531   106980  753257  1.221e9 T_urticae_715m_12_2_3.0_0.5_scaffolds_final.fa
```

#### 1
Purge
```bash
#Already performed above
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234131

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17369   17369   2310    3103    28042   66834   271119  1429785 831.4e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57288023

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57297214

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17635   17635   1958    1000    27236   85729   296983  1429785 831.4e6 T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57245205

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57278999

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#18570   18570   2699    1000    27169   70150   203022  894820  831.4e6 T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57285516

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17658   17658   1817    1000    27169   79480   352049  1398806 831.4e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57297383

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57303753

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17909   17909   1631    1000    26826   98900   359011  1503673 831.4e6 T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds_scaffolds_final.fa

Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds_mapped.PT.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.PT.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57301312

#no contig correction, no lower contig limit:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs_no
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57308051

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#16165   16165   1459    1000    27961   122976  392097  1398806 831.4e6 T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds_scaffolds_final.fa

seqtk seq -a "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/yahs_no/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds_scaffolds_final.fa" | awk '/^>/ {if (seq) print length(seq); seq=""; print $0; next} { seq = seq $0 } END { if (seq) print length(seq) }' | head -n 20

#NO dedup
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap2.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57301391
```
#### 2
Purge
```bash
#Already performed above
```
Pilon
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged.fa
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

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged.fa
Alignment=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@_TellSeqPurged.fa@@g')_sorted_markdups_Tellseq_trimmed.bam
Index=$(dirname $Assembly)/../bwa/$(basename $Assembly | sed 's@_TellSeqPurged.fa@@g')_sorted_markdups_Tellseq_trimmed.bam.bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@_TellSeqPurged.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57234144
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_pilon.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57234826

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17348   17348   2288    3103    28028   66896   275511  1429565 831e6   output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/scaff10x/T_urticae_715m_12_2_3.0_0.5_pilon_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291067

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/scaff10x/T_urticae_715m_12_2_3.0_0.5_pilon_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/scaff10x/T_urticae_715m_12_2_3.0_0.5_pilon_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/scaff10x/T_urticae_715m_12_2_3.0_0.5_pilon_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57297215
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_pilon.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57245212

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_pilon.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_pilon_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_pilon_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57279132

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#18599   18599   2711    1000    27149   69430   199205  804660  831e6   T_urticae_715m_12_2_3.0_0.5_pilon_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/pilon/yahs/T_urticae_715m_12_2_3.0_0.5_pilon_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57297217
```
#### 3
Break10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/break10x
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_break10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57103744

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#26257   26257   6044    3103    31608   62171   113346  753257  1.221e9 T_urticae_715m_12_2_3.0_0.5_break.fa

#BUSCO
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/T_urticae_715m_12_2_3.0_0.5_break.fa); do
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
done #57228868
```
Purge
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/T_urticae_715m_12_2_3.0_0.5_break.fa
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T502_barcoded.fastq.gz
T2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T504_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1 $T2
#57231173, 57231332

#Purge with Tellseq (Illumina) Reads:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/T_urticae_715m_12_2_3.0_0.5_break.fa
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/bwa/T_urticae_715m_12_2_3.0_0.5_break_sorted_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57233947

sbatch ~/git_repos/Pipelines/Trioza_merqury.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
#57291757
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57258663

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17453   17453   2389    3103    28006   65999   256603  1467321 830.6e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57287410

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57291036

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17620   17620   1933    1000    27203   84539   304039  1429785 830.6e6 T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/scaff10x/yahs/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaff10xscaffolds_scaffolds_final.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57298252
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57276683

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57279788

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#18552   18552   2697    1000    27144   70273   201772  894820  830.6e6 T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/yahs/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57291012

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17626   17626   1803    1000    27144   80146   346630  1308351 830.6e6 output_scaffolds.fasta
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
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

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/bwa/T_urticae_715m_12_2_3.0_0.5_break_sorted_markdups_Tellseq_trimmed.bam
Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/bwa/T_urticae_715m_12_2_3.0_0.5_break_sorted_markdups_Tellseq_trimmed.bam.bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57276696
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57279793

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17443   17443   2379    3103    27998   65874   258926  1466934 830.2e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57285451

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57290957

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17606   17606   1908    1000    27197   84000   309042  1429565 830.2e6 T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/scaff10x/yahs/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291017

Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/scaff10x/yahs/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaff10xscaffolds_scaffolds_final_mapped.PT.bam
OutDir=$(dirname $Alignment)
OutFile=$(basename $Alignment | sed 's@_mapped.PT.bam@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pretextmap.sh $Alignment $OutDir $OutFile
#57297209 - 57306220
```
YAHS -> Scaff10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57279794

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57285405

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#18578   18578   2707    1000    27136   69724   200944  804660  830.2e6 T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/yahs/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57285437

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17661   17661   1816    1000    27136   78817   342002  1168469 830.2e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/purge_dups/pilon/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_pilon_scaffolds_final_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291810
```
#### 5
Break10x
```bash
#Already performed above
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/T_urticae_715m_12_2_3.0_0.5_break.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57252744

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#24385   24385   4177    3103    31608   64079   205296  1378212 1.221e9 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57291870

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_scaff10xscaffolds.fasta
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/T_urticae_715m_12_2_3.0_0.5_break.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57252750

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/T_urticae_715m_12_2_3.0_0.5_break.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/T_urticae_715m_12_2_3.0_0.5_break_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/T_urticae_715m_12_2_3.0_0.5_break_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57285404

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#28698   28698   6469    1000    29774   57418   105867  753257  1.221e9 T_urticae_715m_12_2_3.0_0.5_break_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/break10x/yahs/T_urticae_715m_12_2_3.0_0.5_break_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57297225
```
#### 6
Pilon
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
Alignment=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_sorted_markdups_Tellseq_trimmed.bam
Index=$(dirname $Assembly)/bwa/$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_sorted_markdups_Tellseq_trimmed.bam.bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57227844
```
Break10x
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/T_urticae_715m_12_2_3.0_0.5.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/break10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_break10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57233953

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#26256   26256   6042    3103    31606   62143   113314  753257  1.22e9  T_urticae_715m_12_2_3.0_0.5_break.fa
```
Purge
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/T_urticae_715m_12_2_3.0_0.5_break.fa
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T502_barcoded.fastq.gz
T2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T504_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fa@@g')_Tellseq_trimmed
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1 $T2
#57245194

#Purge with Tellseq (Illumina) Reads:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/T_urticae_715m_12_2_3.0_0.5_break.fa
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/bwa/T_urticae_715m_12_2_3.0_0.5_break_Tellseq_trimmed.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57252713

sbatch ~/git_repos/Pipelines/Trioza_merqury.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
#57291758
```
Scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57276689

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#17448   17448   2387    3103    27991   65744   256559  1466891 829.8e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57293720

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/purge_dups/scaff10x/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaff10xscaffolds.fasta
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@_break_TellSeqPurged.fa@_pilon_break_TellSeqPurged@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57285401

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_pilon_break_TellSeqPurged_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/purge_dups/T_urticae_715m_12_2_3.0_0.5_pilon_break_TellSeqPurged_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57285403

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#18528   18528   2678    1000    27115   71000   202177  804660  829.8e6 T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaffolds_final.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/pilon/break10x/purge_dups/yahs/T_urticae_715m_12_2_3.0_0.5_break_TellSeqPurged_scaffolds_final.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57297226
```
#### 7
Filtered -> Purge_dups -> pilon -> break10x -> scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/bwa/T_urticae_715m_12_2_3.0_0.5_sorted_markdups_Tellseq_trimmed.bam
Index=$(echo $Alignment).bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57491054

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/pilon/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_pilon.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/break10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_break10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57495010

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/pilon/break10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_pilon_break.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57495930

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/pilon/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_pilon_break_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57508776

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/pilon/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_pilon_break_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/pilon/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_pilon_break_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/pilon/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_pilon_break_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57508814

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#16954   16954   2302    3103    28327   66600   262598  1429565 816e6   T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_pilon_break_scaff10xscaffolds_scaffolds_final.fa
```
Filtered -> Purge_dups -> purge_haplotigs -> pilon -> break10x -> scaff10x -> YAHS
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged_curated.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/bwa/T_urticae_715m_12_2_3.0_0.5_sorted_markdups_Tellseq_trimmed.bam
Index=$(echo $Alignment).bai
OutDir=$(dirname $Assembly)/pilon
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_pilon.sh $Assembly $Alignment $Index $OutDir $OutPrefix
#57491048

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/pilon/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged_curated_pilon.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/break10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_break10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57495009

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/pilon/break10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged_curated_pilon_break.fa
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57495931

sleep 36000s
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/pilon/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_HiFiPurged_curated_pilon_break.fa_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57508795

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/pilon/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_HiFiPurged_curated_pilon_break.fa_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/pilon/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged_curated_pilon_break_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/pilon/break10x/scaff10x/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged.fa_HiFiPurged_curated_pilon_break_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57517791
```
#### TreeVal
```bash
switch-institute EI
export PATH=/ei/software/cb/bin:$PATH
source eihic-0.2.1_CBG

echo /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz > samples.csv #HiC reads forward 
echo /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz >> samples.csv #HiC reverse reads
echo /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/purge_dups/yahs/scaff10x/T_urticae_715m_12_2_3.0_0.5_TellSeqPurged_scaffolds_final_scaff10xscaffolds.fasta >> samples.csv #path to scaffold
echo Trioza_urticae >> samples.csv #name of organism (will be used for naming only)
echo /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz,/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz >> samples.csv #Long Reads

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/treeval
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/treeval
eihic configure --samples_csv samples.csv --curation -o $OutDir

eihic run --library omni-c --curation --jobs 100 -v -np run_config


Provide sample information in tab-separated format. Please refer to
                        the sample file:
                        /ei/software/cb/eihic/0.2.1/x86_64/lib/python3.9/site-packages/eihic/etc/run_config.yaml for more information above the
                        csv format. A template is provided here
                        /ei/software/cb/eihic/0.2.1/x86_64/lib/python3.9/site-packages/eihic/etc/samples.csv. (default: None)




/ei/software/cb/eihic/0.2.1/x86_64/lib/python3.9/site-packages/eihic/etc/run_config.yaml
```





















































#### MitoHiFi
```bash
source /nbi/software/staging/RCSUPPORT-2522/stagingloader
mitohifi.py -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
```

#### HiC and PacBio Prep
```bash
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/
#urticae_286170-S3HiC_R1.fastq.gz  urticae_286170-S3HiC_R2.fastq.gz

source package 638df626-d658-40aa-80e5-14a275b7464b
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/
samtools import -@16 -r ID:urticae_286170-S3HiC -r CN:S3HiC -r PU:urticae_286170-S3HiC -r SM:urticae_286170 urticae_286170-S3HiC_R1.fastq.gz urticae_286170-S3HiC_R2.fastq.gz -o urticae_286170-S3HiC.cram
samtools index -@16 urticae_286170-S3HiC.cram
#57183981


source package /nbi/software/testing/bin/seqtk-1.2
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/
mkdir fasta
for i in *fastq.gz; do
  echo $i
  j=${i%.fastq.gz}
  echo $j
  seqtk seq -a $i > fasta/${j}.fasta
done

cd fasta
for i in *.fasta; do
  echo $i
  gzip $i
done
```












































































#### Flye
```bash
  for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@flye@g')/815m
    OutFile=T_urticae_815m
    Genomesize=815m
    PolishingIterations=1 #default
    ReadErrorRate=0.001 #default
    DataType=pacbio-hifi
    mkdir -p $OutDir
    sbatch $ProgDir/run_flye.sh $OutDir $OutFile $Genomesize $PolishingIterations $DataType $ReadErrorRate $Run1 $Run2
  done #55604458, 55609706

  for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/flye/815m/T_urticae_815m_unscaffoldeddef.fa); do
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/arthropoda_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    echo $OutFile >> logs/buscolog.txt
    sleep 30s
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/flye/815m/T_urticae_815m_unscaffoldeddef.fa); do
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/insecta_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    echo $OutFile >> logs/buscolog.txt
    sleep 30s
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/flye/815m/T_urticae_815m_unscaffoldeddef.fa); do
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    OutFile=$(basename $Genome | cut -d '.' -f1)_$(echo $Database | cut -d '/' -f7)
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    echo $OutFile >> logs/buscolog.txt
    sleep 30s
    sbatch $ProgDir/run_busco.sh $Genome $Database $OutDir $OutFile 2>&1 >> logs/buscolog.txt
done
```		
#### Canu
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi); do
    ProgDir=~/git_repos/Wrappers/NBI
    Run1=$(ls $ReadDir/urticae_hifi-reads.fastq.gz)
    Run2=$(ls $ReadDir/urticae_hifi-3rdSMRTcell.fastq.gz)
    OutDir=$(echo $ReadDir|sed 's@raw_data@assembly/genome@g'|sed 's@HiFi@canu@g')/715m
    OutFile=T_urticae_715m
    Genomesize=715m
    DataType=pacbio-hifi
    mkdir -p $OutDir
    sbatch $ProgDir/run_canu.sh $OutDir $OutFile $Genomesize $DataType $Run1 $Run2
done #57411057

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#35928   35928   8160    3075    23265   45912   86043   371196  1.242e9 T_urticae_715m.contigs.fasta

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/T_urticae_715m.contigs.fasta); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done 
#57547088

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/T_urticae_715m.contigs.fasta); do
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

#T_urticae_715m.contigs.fasta_insecta_odb10_short_summary.txt
#        C:82.6%[S:28.2%,D:54.4%],F:10.4%,M:7.0%,n:1367
#T_urticae_715m.contigs.fasta_arthropoda_odb10_short_summary.txt
#        C:83.1%[S:27.0%,D:56.1%],F:10.2%,M:6.7%,n:1013
#T_urticae_715m.contigs.fasta_hemiptera_odb10_short_summary.txt
#        C:85.0%[S:30.0%,D:55.0%],F:5.9%,M:9.1%,n:2510


Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/T_urticae_715m.contigs.fasta
T1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T502_barcoded.fastq.gz
T2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dna_qc/T_urticae/TellSeq/longranger/Turt_T504_barcoded.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.fasta@@g')_Tellseq_trimmed_2
mkdir $OutDir
sbatch $ProgDir/bwa-mem_unpaired.sh $OutDir $Outfile $Assembly $T1 $T2
#57547105, 57560718

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/T_urticae_715m.contigs.fasta
#Purge with Tellseq (Illumina) Reads:
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/bwa/T_urticae_715m.contigs_Tellseq_trimmed_2.bam
Type=short
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')_TellSeqPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57558531, 57560720, 57561619, 57583735, 57586615

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#26731   26731   5518    3075    20714   44960   91281   371196  868.9e6 T_urticae_715m.contigs_TellSeqPurged.fa

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/T_urticae_715m.contigs_TellSeqPurged.fa); do
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

#arthropoda
#        C:80.7%[S:66.2%,D:14.5%],F:11.2%,M:8.1%,n:1013
#insecta
#        C:79.8%[S:64.5%,D:15.3%],F:11.4%,M:8.8%,n:1367
#hemiptera
#        C:83.2%[S:68.9%,D:14.3%],F:6.3%,M:10.5%,n:2510


ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/T_urticae_715m.contigs.fasta
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.fasta@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2
#57624380

#Purge with HiFi reads:
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/T_urticae_715m.contigs.fasta
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/minimap2/T_urticae_715m.contigs.bam
Type=long
OutDir=$(dirname $Assembly)/purge_dups
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_dups.sh $Assembly $MappingFile $Type $OutDir $OutPrefix
#57625362

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#21439   21439   5006    5110    27685   52222   95389   371196  852.7e6 T_urticae_715m.contigssta_HiFiPurged.fa

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/T_urticae_715m.contigssta_HiFiPurged.fa); do
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

#arthropoda
#        C:81.2%[S:63.1%,D:18.1%],F:11.3%,M:7.5%,n:1013
#insecta
#        C:80.6%[S:61.8%,D:18.8%],F:11.3%,M:8.1%,n:1367
#hemiptera
#       C:83.8%[S:66.3%,D:17.5%],F:5.9%,M:10.3%,n:2510

interactive
source /nbi/software/staging/RCSUPPORT-2652/stagingloader
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/T_urticae_715m.contigssta_HiFiPurged.fa
AssemblyIndex=$(echo $Assembly).fai
MappingFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/minimap2/T_urticae_715m.contigs.bam
purge_haplotigs hist -b $MappingFile -g $Assembly -t 1
#Assess plot for low,mid,high cuttoffs
MappingIndex=$(echo $MappingFile).bai
low=3
mid=33
high=85
diploid_cuttoff=80
junk_cuttoff=80
hap_percent_cuttoff=70
rep_percent_cuttoff=250
OutDir=$(dirname $Assembly)/purge_haplotigs
OutPrefix=$(basename $Assembly | sed 's@.fa@@')_HiFiPurged
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_purge_haplotigs.sh $Assembly $AssemblyIndex $MappingFile $MappingIndex $low $mid $high $diploid_cuttoff $junk_cuttoff $hap_percent_cuttoff $rep_percent_cuttoff $OutDir $OutPrefix
#57645419

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#15804   15804   3708    5110    30860   58528   102588  371196  686.7e6 T_urticae_715m.contigssta_HiFiPurged_HiFiPurged_curated.fasta

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/purge_haplotigs/T_urticae_715m.contigssta_HiFiPurged_HiFiPurged_curated.fasta); do
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

#Arthropoda
#        C:79.9%[S:74.5%,D:5.4%],F:11.9%,M:8.2%,n:1013
#Insecta
#        C:79.4%[S:73.4%,D:6.0%],F:12.0%,M:8.6%,n:1367
#Hemiptera
#        C:82.6%[S:76.9%,D:5.7%],F:6.4%,M:11.0%,n:2510

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/purge_haplotigs/T_urticae_715m.contigssta_HiFiPurged_HiFiPurged_curated.fasta); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done 
#57650867

###############################################################################################################################

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/purge_haplotigs/T_urticae_715m.contigssta_HiFiPurged_HiFiPurged_curated.fasta
ReadDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/10x
OutDir=$(dirname $Assembly)/scaff10x
OutPrefix=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_scaff10x.sh $Assembly $ReadDir $OutDir $OutPrefix
#57674353

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#14650   14650   2554    5110    30860   58528   194157  961883  686.7e6 output_scaffolds.fasta

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/purge_haplotigs/scaff10x/T_urticae_715m.contigssta_HiFiPurged_HiFiPurged_curated_scaff10xscaffolds.fasta
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fasta@@')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $Read1 $Read2
#57681409

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/purge_haplotigs/scaff10x/T_urticae_715m.contigssta_HiFiPurged_HiFiPurged_curated_scaff10xscaffolds.fasta
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/purge_haplotigs/scaff10x/T_urticae_715m.contigssta_HiFiPurged_HiFiPurged_curated_scaff10xscaffolds_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/purge_haplotigs/scaff10x/T_urticae_715m.contigssta_HiFiPurged_HiFiPurged_curated_scaff10xscaffolds_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fasta@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#57689851

#assembly N50 (111879) too small
#n       n:500   n:N50   min     N80     N50     N20     max     sum
#14650   14650   2554    5110    30860   58528   194157  961883  686.7e6 output_scaffolds.fasta

for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/canu/715m/purge_dups/purge_haplotigs/scaff10x/yahs/T_urticae_715m.contigssta_HiFiPurged_HiFiPurged_curated_scaff10xscaffolds_scaffolds_final.fa); do
sbatch ~/git_repos/Pipelines/Trioza_merqury.sh $Assembly  
done 
#57758479
```			
#### longQC
```bash
for Reads in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/*.fastq.gz); do
Datatype=pb-hifi
OutDir=$(dirname $Reads)/longqc/$(basename $Reads | cut -d '.' -f1)_2
OutFile=$(basename $Reads | cut -d '.' -f1)
ProgDir=~/git_repos/Wrappers/NBI
echo ${OutDir}/${OutFile}
mkdir $(dirname $Reads)/longqc
sbatch $ProgDir/run_longqc.sh $Reads $OutDir $OutFile $Datatype
done #55612494,5, 55614853, 4
```

















































### T. Mathers Trioza urticae

ls /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant
Trant_hg-size_set.p_ctg.fa

Assemble with hifiasm - unzip reads first as hifiasm may not like gzipped reads. -> check assembly quality with KAT and BUSCO -> fiddle with hifiasm to improve assembly
 - HiC and Tell-seq data not utilised yet except for apicales

 default:
 n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
22201    22201    5101    5986    30211    59308    110548    75953    843665    994.1e6    Trurt_default.p_ctg.fa

more data default:
n    n:500    L50    min    N80    N50    N20    E-size    max    sum    name
20013    20013    4360    6253    34674    69264    130990    89516    877047    1.009e9    Trurt_default.bp.p_ctg.fa
