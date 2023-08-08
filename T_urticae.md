# Trioza urticae assembly
Contains commands run by T.Heaven in assembly of Trioza urticae

Unless stated otherwise commands were performed from the directory /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae

## Collect data
```bash
#HiFi reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TUF20_hifi_reads.fastq.gz raw_data/T_urticae/HiFi/urticae_hifi-reads.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TUF20_hifi_3rdSMRTcell.fastq.gz raw_data/T_urticae/HiFi/urticae_hifi-3rdSMRTcell.fastq.gz

#Tellseq reads, from 2nd tellseq run:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_I1_T502.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T502_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_R1_T502.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T502_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_R2_T502.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T502_R2.fastq.gz

ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_I1_T504.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T504_I1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_R1_T504.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T504_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_March_2022_T_urticae/220308_NB501793_0297_AHHV72BGXK/Caliber_tellseq_run2_2_psyllids_T_urticae_R2_T504.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_T504_R2.fastq.gz

#No reads here?:
ls /jic/research-groups/Saskia-Hogenhout/reads/genomic/Tellseq_2_psylid_libs_Dec_2021/211210_NB501793_0280_AHWCHWBGXH/

#Tellseq reads from 1st tellseq run
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_I1_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_I1_T505.fastq.gz
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_R1_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_R1_T505.fastq.gz
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_R2_T505.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_R2_T505.fastq.gz

ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_I1_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_I1_T507.fastq.gz
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_R1_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_R1_T507.fastq.gz
ln -s /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_tellseq_run1_2_psyllids/Full/Caliber_tellseq_run1_2_psyllids_R2_T507.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz raw_data/T_urticae/TellSeq/urticae_R2_T507.fastq.gz

#HiC reads:
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_urticae/urticae-286170_S3HiC_R1.fastq-001.gz raw_data/T_urticae/HiC/urticae_286170-S3HiC_R1.fastq.gz
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_HiC_Nov_2022/Trioza_urticae/urticae-286170_S3HiC_R2.fastq-002.gz raw_data/T_urticae/HiC/urticae_286170-S3HiC_R2.fastq.gz
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
~/git_repos/Wrappers/NBI/submit-slurm_v1.1.pl -q nbi-long -m 200000 -c 32 -t 5-00:00:00 -e -j kat_comp -i "source package 7f4fb852-f5c2-4e4b-b9a6-e648176d5543;kat comp -t 32 -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/kat/kat_comp_vs_tellseq_reads -m 31 -H 100000000 -I 100000000 '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_R1_T505.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_R1_T507.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_R2.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_R2.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_R2_T505.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_R2_T507.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T502_I1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_T504_I1.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_I1_T505.fastq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_I1_T507.fastq.gz' 'Trurt_default.bp.p_ctg.fa'"
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
F3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_R1_T505.fastq.gz
R3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_R2_T505.fastq.gz
F4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_R1_T507.fastq.gz
R4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/TellSeq/urticae_R2_T507.fastq.gz
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
```bash
#These assemblies were kept:
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm/default/Trurt_default.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/570m/32/3/10.0/0.25/T_urticae_570m_32_3_10.0_0.25.bp.p_ctg.fa

#All other assemblies were deleted - 278 in total:
for assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifi*/*m/*.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifi*/*m/*/*.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifi*/*m/*/*/*.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifi*/*m/*/*/*/*/*.fa | grep -v 'Trurt_default.bp.p_ctg.fa\|T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa\|T_urticae_570m_32_3_10.0_0.25.bp.p_ctg.fa'); do
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
done #55686243, 55713570, 55722826, 55730047
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
