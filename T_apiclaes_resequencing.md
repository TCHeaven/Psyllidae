```bash
#first batch of Psyllid resequencing (25 T. apicales individuals):
cd /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_popgen_Nov_2023/Batch1-231123
tar -xvzf X204SC23101516-Z01-F001.tar.gz

#second batch of Psyllid resequencing ():
cd /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_popgen_Nov_2023/Batch2-141223
tar -xvf X204SC23111485-Z01-F001.tar
tar -cvzf X204SC23111485-Z01-F001.tar.gz X204SC23111485-Z01-F001

ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_popgen_Nov_2023/Batch2-141223/X204SC23111485-Z01-F001/01.RawData/Dyap* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/.
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_popgen_Nov_2023/Batch2-141223/X204SC23111485-Z01-F001/01.RawData/Tap* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/.
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_popgen_Nov_2023/Batch2-141223/X204SC23111485-Z01-F001/01.RawData/TrAp* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/.
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_popgen_Nov_2023/Batch1-231123/X204SC23101516-Z01-F001/01.RawData/TrAp* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/.

ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_popgen_Nov_2023/Batch2-141223/X204SC23111485-Z01-F001/01.RawData/TrA0* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/.
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_popgen_Nov_2023/Batch2-141223/X204SC23111485-Z01-F001/01.RawData/TrAn* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/.

mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/Dyap15/Dyap15_EKDN230046453-1A_22FNL2LT3_L6_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/Dyap15/Dyap15_EKDN230046453-1A_22FNL2LT3_L6_3.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/Dyap15/Dyap15_EKDN230046453-1A_22FNL2LT3_L6_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/Dyap15/Dyap15_EKDN230046453-1A_22FNL2LT3_L6_4.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/Dyap3/Dyap3_EKDN230046444-1A_22FNL2LT3_L6_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/Dyap3/Dyap3_EKDN230046444-1A_22FNL2LT3_L6_3.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/Dyap3/Dyap3_EKDN230046444-1A_22FNL2LT3_L6_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/Dyap3/Dyap3_EKDN230046444-1A_22FNL2LT3_L6_4.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp19_3/TrAp19_3_EKDN230043615-1A_22FGM7LT3_L5_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp19_3/TrAp19_3_EKDN230043615-1A_22FGM7LT3_L5_3.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp19_3/TrAp19_3_EKDN230043615-1A_22FGM7LT3_L5_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp19_3/TrAp19_3_EKDN230043615-1A_22FGM7LT3_L5_4.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp19_7/TrAp19_7_EKDN230043618-1A_22FGM7LT3_L5_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp19_7/TrAp19_7_EKDN230043618-1A_22FGM7LT3_L5_3.fq.gz 
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp19_7/TrAp19_7_EKDN230043618-1A_22FGM7LT3_L5_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp19_7/TrAp19_7_EKDN230043618-1A_22FGM7LT3_L5_4.fq.gz 
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp13_6/TrAp13_6_EKDN230043612-1A_22FGM7LT3_L6_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp13_6/TrAp13_6_EKDN230043612-1A_22FGM7LT3_L6_3.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp13_6/TrAp13_6_EKDN230043612-1A_22FGM7LT3_L6_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp13_6/TrAp13_6_EKDN230043612-1A_22FGM7LT3_L6_4.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp13_6/TrAp13_6_EKDN230043612-1A_22FGYHLT3_L2_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp13_6/TrAp13_6_EKDN230043612-1A_22FGYHLT3_L2_5.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp13_6/TrAp13_6_EKDN230043612-1A_22FGYHLT3_L2_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp13_6/TrAp13_6_EKDN230043612-1A_22FGYHLT3_L2_6.fq.gz


mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/TrAn16_1/TrAn16_1_EKDN230046458-1A_22FNL2LT3_L6_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/TrAn16_1/TrAn16_1_EKDN230046458-1A_22FNL2LT3_L6_3.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/TrAn16_1/TrAn16_1_EKDN230046458-1A_22FNL2LT3_L6_2.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/TrAn16_1/TrAn16_1_EKDN230046458-1A_22FNL2LT3_L6_4.fq.gz
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/TrAn16_2/TrAn16_2_EKDN230046459-1A_22FNL2LT3_L6_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/TrAn16_2/TrAn16_2_EKDN230046459-1A_22FNL2LT3_L6_3.fq.gz 
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/TrAn16_2/TrAn16_2_EKDN230046459-1A_22FNL2LT3_L6_1.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/TrAn16_2/TrAn16_2_EKDN230046459-1A_22FNL2LT3_L6_4.fq.gz 
```
### QC  <a name="3"></a>
#### fastqc and qualimap <a name="4"></a>
Raw reads were assessed for quality and coverage of the clone O_v2 reference assembly:
```bash
for ReadDir in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/*/ -mindepth 1 -type d); do
    if [ ! -e ${ReadDir}/qualimap/*genome_results_gff.txt ]; then
    echo Running for:
    ls ${ReadDir}/qualimap/*genome_results_gff.txt
    Fread=$(ls ${ReadDir}/*_1.fq.gz)
    Rread=$(ls ${ReadDir}/*_2.fq.gz)
    Fread2=$(ls ${ReadDir}/*_3.fq.gz)
    Rread2=$(ls ${ReadDir}/*_4.fq.gz)
    Fread3=$(ls ${ReadDir}/*_5.fq.gz)
    Rread3=$(ls ${ReadDir}/*_6.fq.gz)
    OutDir=$(echo ${ReadDir})
    Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    Gff=
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3
    else
    echo Already Done:
    ls ${ReadDir}/qualimap/*genome_results.txt
    fi
done
#

for ReadDir in $(find /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/*/ -mindepth 1 -type d); do
    if [ ! -e ${ReadDir}/qualimap/*genome_results_gff.txt ]; then
    echo Running for:
    ls ${ReadDir}/qualimap/*genome_results_gff.txt
    Fread=$(ls ${ReadDir}/*_1.fq.gz)
    Rread=$(ls ${ReadDir}/*_2.fq.gz)
    Fread2=$(ls ${ReadDir}/*_3.fq.gz)
    Rread2=$(ls ${ReadDir}/*_4.fq.gz)
    OutDir=$(echo ${ReadDir})
    Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 
    else
    echo Already Done:
    ls ${ReadDir}/qualimap/*genome_results.txt
    fi
done
#

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/raw_data/Myzus/persicae/WGS/Archana_Feb2021/*/); do
sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
coverage=$(grep 'mean coverageData' ${ReadDir}qualimap/*_genome_results_gff.txt | rev | cut -d ' ' -f1 | rev | sed 's@X@@g')
echo $sample raw reads have average coverage of ${coverage}
done
```

### Trimming <a name="5"></a>
#### Trim galore <a name="6"></a>
Adapters and low quality regions were trimmed from raw reads via trim-galore:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/TrAp13_3); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f1 | rev)
    Fread=$(ls $ReadDir/*_1.fq.gz)
    Rread=$(ls $ReadDir/*_2.fq.gz)
    Fread2=$(ls $ReadDir/*_3.fq.gz)
    Rread2=$(ls $ReadDir/*_4.fq.gz)
    Fread3=$(ls ${ReadDir}/*_5.fq.gz)
    Rread3=$(ls ${ReadDir}/*_6.fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@dna_qc@g')/trim_galore
    OutFile=${sample}_trimmed
    Quality=20
    Length=50
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3 
done 
#59071170-97, 59168580

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/*); do
    sample=$(echo $ReadDir | rev | cut -d '/' -f1 | rev)
    Fread=$(ls $ReadDir/*_1.fq.gz)
    Rread=$(ls $ReadDir/*_2.fq.gz)
    Fread2=$(ls $ReadDir/*_3.fq.gz)
    Rread2=$(ls $ReadDir/*_4.fq.gz)
    OutDir=$(echo $ReadDir | sed 's@raw_data@dna_qc@g')/trim_galore
    OutFile=${sample}_trimmed
    Quality=20
    Length=50
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/run_trim_galore.sh $OutDir $OutFile $Quality $Length $Fread $Rread $Fread2 $Rread2 
done 
#59071018-26
```

## Alignment + QC

Psyllids
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/); do
    if [ ! -e ${ReadDir}qualimap/*genome_results_gff.txt ]; then
    echo Running for:
    ls ${ReadDir}qualimap/*genome_results_gff.txt
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo ${ReadDir})
    Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3
    else
    echo Already Done:
    ls ${ReadDir}qualimap/*genome_results.txt
    fi
done
#59168590-59168618, TrAp13_3 is missing

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/); do
    if [ ! -e ${ReadDir}qualimap/*genome_results_gff.txt ]; then
    echo Running for:
    ls ${ReadDir}qualimap/*genome_results_gff.txt
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutDir=$(echo ${ReadDir})
    Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 
    else
    echo Already Done:
    ls ${ReadDir}qualimap/*genome_results.txt
    fi
done
#59168581-59168589
```

Liberibacter:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/); do
	OutDir=$(echo ${ReadDir})Liberibacter/bwa
    if [ ! -e ${OutDir}/qualimap/*genome_results.txt ]; then
    echo Running for:
    ls ${OutDir}/qualimap/*genome_results.txt
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
	OutFile=$(echo $ReadDir | cut -d '/' -f11)
	Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
	Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.gff
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/bwa-mem.sh $OutDir $OutFile $Reference_genome $Gff $Fread $Rread 
    else
    echo Already Done:
    ls ${OutDir}/qualimap/*genome_results.txt
    fi
done
#59181595-59181622, TrAp13_3 is missing

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/Liberibacter/bwa/*.bam); do
BamFile=$file
Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
OutDir=$(dirname $file)
Gff=NA
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_qualimap.sh $BamFile $Reference_genome $OutDir $Gff
done
#59263407-34

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/Liberibacter/bwa/qualimap/*genome_results.txt); do
echo $file | cut -d '/' -f8,9,10
grep -w 'scaffold_18' $file 
grep 'scaffold_184' $file 
grep 'scaffold_1142' $file 
done
```

#### Picard and GATK <a name="8"></a>

Files were sorted by scaffold and coordinate level, duplicates were marked and removed and the files were re-indexed:
```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/bwa-mem/*_trimmed_1_sorted.bam); do
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_sort_bam.sh $file
done 
#59181934-73

source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
source package 638df626-d658-40aa-80e5-14a275b7464b
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/bwa-mem/*_trimmed_1_sorted.bam); do
    OutDir=$(dirname $file)
    OutFile=$(basename $file | sed s'@_sorted.bam@_sorted_MarkDups.bam@g')
    metrics=$(basename $file | sed s'@_trimmed_1_sorted.bam@@g')_marked_dup_metrics.txt
    java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar SortSam I=${file} O=${OutDir}/temp.bam SORT_ORDER=coordinate
    java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${OutDir}/temp.bam O=${OutDir}/${OutFile} M=${OutDir}/${metrics}
    cd ${OutDir}
    rm temp.bam
    samtools index ${OutFile}
    cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae
done
#59188706

#Input files were subsequently deleted to save space
```
Reads near detected indels realigned to remove alignment artifacts:
```bash
#The reference genome was indexed and a dictionary created:
interactive
source switch-institute ei
source package 638df626-d658-40aa-80e5-14a275b7464b
source pilon-1.22
source package /tgac/software/testing/bin/picardtools-2.1.1
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/
java -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar CreateSequenceDictionary R=T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa O=T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.dict
samtools faidx T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa 

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/
java -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar CreateSequenceDictionary R=T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa O=T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.dict
samtools faidx T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa 
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae


#GATK indel realignment:
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/bwa-mem/*MarkDups.bam); do
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa 
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_realign.sh $file $Reference 
done 
#59261801-27

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/bwa-mem/*MarkDups.bam); do
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa 
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_realign.sh $file $Reference 
done 
#59261845-53
```
### Variant calling <a name="9"></a>

Combined into a .vcf and called variants with bcftools:
```bash
source package 638df626-d658-40aa-80e5-14a275b7464b
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/bwa-mem/gatk/*realigned.bam > bamlist3.txt
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/
bcftools mpileup --threads 16 -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/bamlist3.txt --annotate AD,DP --fasta-ref /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa -O z -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/variants_1.vcf.gz

bcftools call --threads 16 --ploidy 2 -Oz -v -m -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/variants.vcf.gz  /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/variants_1.vcf.gz

ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/bwa-mem/gatk/*realigned.bam > bamlist4.txt
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/
bcftools mpileup --threads 16 -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/bamlist4.txt --annotate AD,DP --fasta-ref /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa -O z -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/variants_1.vcf.gz

bcftools call --threads 16 --ploidy 2 -Oz -v -m -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/variants.vcf.gz  /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/variants_1.vcf.gz
#59377250
```