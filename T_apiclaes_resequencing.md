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
ln -s /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_popgen_Nov_2023/Batch1-231123/X204SC23101516-Z01-F001/01.RawData/Tap* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/.

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
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/*/); do
    if ! find "${ReadDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
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
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3
    else
    echo Already Done:
    ls ${ReadDir}/qualimap/*genome_results.txt
    fi
done
#59546039-59546079

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/*/); do
    if ! find "${ReadDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
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
#59546080-59546088

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
#59071170-97, 59168580, 59546024-38

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
    if ! find "${ReadDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
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
#59168590-59168618, TrAp13_3 is missing, 59546237-59546250

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/); do
    if ! find "${ReadDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
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
    if ! find "${OutDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
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
#59181595-59181622, TrAp13_3 is missing, 59546108-21

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
echo $file | cut -d '/' -f9,10,11
grep -w 'scaffold_18' $file 
grep 'scaffold_184' $file 
grep 'scaffold_1142' $file 
done
```
All contaminants:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/); do
    OutDir=$(echo ${ReadDir})contaminants/bwa
    if ! find "${OutDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
    echo Running for:
    ls ${OutDir}/qualimap/*genome_results.txt
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutFile=$(echo $ReadDir | cut -d '/' -f11)
    Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_suspected_contaminants.fa
    Gff=NA
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/bwa-mem.sh $OutDir $OutFile $Reference_genome $Gff $Fread $Rread 
    else
    echo Already Done:
    ls ${OutDir}/qualimap/*genome_results.txt
    fi
done
#59546506-46

echo scaffold > api_conts.txt
grep 'scaffold_' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/TrAp19_8/trim_galore/contaminants/bwa/qualimap/TrAp19_8_genome_results.txt | awk '{print $1}' >> api_conts.txt
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/contaminants/bwa/qualimap/*genome_results.txt); do
echo $file | cut -d '/' -f11 > column_to_add.txt
grep 'scaffold_'  $file | awk '{print $4}' >> column_to_add.txt
paste api_conts.txt column_to_add.txt > merged_file.txt && mv merged_file.txt api_conts.txt
done
paste api_conts.txt temp_api_conts.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/contaminant_coverage.txt
awk -F'\t' '{print $45}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/contaminant_coverage.txt
sort -t$'\t' -k45n /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/contaminant_coverage.txt > temp.txt && mv temp.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/contaminant_coverage.txt


grep '>' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_suspected_contaminants.fa | sed 's@>@@g' > temp_search.txt
echo kraken > temp_ant_conts.txt
grep -wf temp_search.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/kraken2.1.3/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_kraken2nt_output.txt >> temp_ant_conts.txt
sort -t$'\t' -k4n temp_ant_conts.txt > temp.txt && mv temp.txt temp_ant_conts.txt
awk -F'\t' '{print $1, $2, $3, $4}' temp_ant_conts.txt

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/); do
    OutDir=$(echo ${ReadDir})contaminants/bwa
    if ! find "${OutDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
    echo Running for:
    ls ${OutDir}/qualimap/*genome_results.txt
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    OutFile=$(echo $ReadDir | cut -d '/' -f11)
    Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_suspected_contaminants.fa
    Gff=NA
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir -p $OutDir
    sbatch $ProgDir/bwa-mem.sh $OutDir $OutFile $Reference_genome $Gff $Fread $Rread 
    else
    echo Already Done:
    ls ${OutDir}/qualimap/*genome_results.txt
    fi
done
#59546555-63

echo scaffold > ant_conts.txt
grep 'scaffold_' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/TrAn16_5/trim_galore/contaminants/bwa/qualimap/TrAn16_5_genome_results.txt | awk '{print $1}' >> ant_conts.txt
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/contaminants/bwa/qualimap/*genome_results.txt); do
echo $file | cut -d '/' -f11 > column_to_add.txt
grep 'scaffold_'  $file | awk '{print $4}' >> column_to_add.txt
paste ant_conts.txt column_to_add.txt > merged_file.txt && mv merged_file.txt ant_conts.txt
done
paste ant_conts.txt temp_ant_conts.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/contaminant_coverage.txt


awk -F'\t' '{print $14}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/contaminant_coverage.txt
sort -t$'\t' -k14n /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/contaminant_coverage.txt > temp.txt && mv temp.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/contaminant_coverage_pallida.txt
```

#### Picard and GATK <a name="8"></a>

Files were sorted by scaffold and coordinate level, duplicates were marked and removed and the files were re-indexed:
```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/bwa-mem/*_trimmed_1_sorted.bam); do
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_sort_bam.sh $file
done 
#59181934-73, 59547731-43

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
#59188706, 59552814

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
#    if [ ! -e "$(dirname $file)/gatk" ]; then
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa 
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_realign.sh $file $Reference 
#    fi
done 
#59261801-27, 59610112-51, 59639431-70

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
#59377250,59616793,59642627
```
### Filter  <a name="10"></a>
#### Genmap filter for genome mappability <a name="11"></a>

We further refined our input files by calculating mappability of the genome with GenMap (v1.3.0) with the parameters -K 100 -E 2. This estimates k-mer uniqueness and identifies regions of the genome where Illumina reads are unable to map uniquely. We masked all regions larger than 100 bp with less than 1 i.e., max mappability. 10.1038/s41586-021-04269-6 and 10.1038/s41467-023-43383-z use k=100 for 150bp paired reads.
```bash
#Apicales
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/variants.vcf.gz | wc -l #48,883,771

VCF=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/variants.vcf.gz
OutDir=$(dirname $VCF)/genmap
OutFile=$(basename $VCF | sed 's@.vcf.gz@@g')
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
Repeatmodeller=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa.out
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_genmap_mappability_masking.sh $OutDir $OutFile $Reference $Repeatmodeller $VCF
#59650971

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable.vcf.gz | wc -l #15,554,842

#Pallida
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/variants.vcf.gz | wc -l #7,763,414

VCF=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/variants.vcf.gz
OutDir=$(dirname $VCF)/genmap
OutFile=$(basename $VCF | sed 's@.vcf.gz@@g')
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
Repeatmodeller=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa.out
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_genmap_mappability_masking.sh $OutDir $OutFile $Reference $Repeatmodeller $VCF
#59426167

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable.vcf.gz | wc -l #3,484,502

srun -p jic-short  -c 1 --mem 50G --pty bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
After filtering, kept 3941793 out of a possible 7472067 Sites

Apicales Mappable bases (genmap): 375,902,061
Apicales Masked bases (repeatmodeler): 359,790,492
Apicales Callable bases (overlap of VCF and mappable genmap regions, minus repeatmodeler masked regions): 31,003,238

Pallida Mappable bases (genmap): 361,165,258
Pallida Masked bases (repeatmodeler): 366,704,478
Pallida Callable bases (overlap of VCF and mappable genmap regions, minus repeatmodeler masked regions): 7,225,232

```python
import matplotlib.pyplot as plt

def plot_bedgraph(bedgraph_file, chromosome):
    starts, values = [], []
    with open(bedgraph_file, 'r') as f:
        for line in f:
            if line.startswith('track') or line.startswith('browser'):
                continue
            fields = line.strip().split('\t')
            if fields[0] == chromosome:
                starts.append(int(fields[1]))
                values.append(float(fields[3]))
    plt.plot(starts, values, label=f'{chromosome}')
    plt.xlabel('Position')
    plt.ylabel('Mappability')
    plt.title(f'Genmap Mappability Plot - {chromosome}')
    plt.legend()

bedgraph_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_genmap.bedgraph'
#bedgraph_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_genmap.bedgraph'
chromosomes_to_plot = ['SUPER_1', 'SUPER_2', 'SUPER_3', 'SUPER_4', 'SUPER_5', 'SUPER_6', 'SUPER_7', 'SUPER_8', 'SUPER_9', 'SUPER_10', 'SUPER_11_1', 'SUPER_12', 'SUPER_13']
plt.figure(figsize=(96, 48))

for idx, chromosome in enumerate(chromosomes_to_plot, 1):
    plt.subplot(13, 1, idx)
    plot_bedgraph(bedgraph_file, chromosome)
    plt.tight_layout()

plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/genemap_mappability_plot.png')
#plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/genemap_mappability_plot.png')
```


#### Filter samples for missingness: <a name="12"></a>
```bash
source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable.vcf.gz --missing-indv
mv out.imiss /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable.out.imiss  

vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable.vcf.gz --missing-indv
mv out.imiss /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable.out.imiss                                    
```
Plot the missingness level of samples:
```python
import matplotlib.pyplot as plt

# Read the F_MISS data from a file, skipping the header line
data = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable.out.imiss', 'r') as file:
#with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable.out.imiss', 'r') as file:
    next(file)  # Skip the header line
    for line in file:
        values = line.strip().split('\t')
        f_miss = float(values[-1])
        data.append(f_miss)

# Create the scatter plot
plt.scatter(range(len(data)), data, marker='o', color='blue')
plt.xlabel('Individual')
plt.ylabel('Missing Level (F_MISS)')
plt.title('Missing Level per Individual')
plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable.out.imiss_plot.png')
#plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable.out.imiss_plot.png')
```
Missingness for all samples is very low.

#### Filter for SNP quality: <a name="13"></a>

SNPs were filtered to keep only bi-allelic SNPs, with minimum depth of 5, minimum quality score of 30, and maximum missingness of SNPs of 10%. 
```bash
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered
bgzip -k /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf
#After filtering, kept 40 out of 40 Individuals
#After filtering, kept 11,870,982 out of a possible 21,805,206 Sites

vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered
bgzip -k /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf
#After filtering, kept 9 out of 9 Individuals
#After filtering, kept 944,135 out of a possible 3,484,111 Sites
```
#### VCF phylogeny
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/vcf2phylip.py -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf --output-folder /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/ -f -p

mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/iqtree2
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/iqtree2
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.min4.fasta
cpu=8
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $Alignment -m MF -T AUTO --threads-max $cpu
#59701773
#Best-fit model: PMB+F+R3 chosen according to BIC
#Analysis results written to:
#  IQ-TREE report:                /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.min4.fasta.iqtree
#  Tree used for ModelFinder:     /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.min4.fasta.treefile
#  Screen log file:               /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.min4.fasta.log


vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf --thin 1000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.vcf

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/vcf2phylip.py -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.vcf.recode.vcf --output-folder /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/ -f -p

mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/iqtree2
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/iqtree2
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.min4.fasta
cpu=8
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $Alignment -m MF -T AUTO --threads-max $cpu
#59702173,59727696
```
#### Network
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/splitstree
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf.gz); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/splitstree
Outfile=p_dis.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#59701785

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/mat2csv.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/splitstree/p_dis.mat
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/csv2dist.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/splitstree/p_dis_mperc.csv
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/splitstree
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/splitstree
Outfile=p_dis.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#59701786

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/mat2csv.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/splitstree/p_dis.mat
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/csv2dist.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/splitstree/p_dis_mperc.csv
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
```
#### admixture
"As a rule of thumb, we have found that 10,000 markers suffice to perform GWAS correction for continentally separated populations (for example, African, Asian, and European populations FST > .05) while more like 100,000 markers are necessary when the populations are within a continent (Europe, for instance, FST < 0.01)."
```bash
source package /tgac/software/testing/bin/admixture-1.3.0
source package /nbi/software/testing/bin/bcftools-1.8
source package /nbi/software/testing/bin/plink-1.9 
source package 01ef5a53-c149-4c9e-b07d-0b9a46176cc0

#plink cannot process more than 95 chromosome, therefore only SNPs from the 13 chromosomal scaffolds used for admixture estimation
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink

head -n 5 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.vcf
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf | grep -w '##ALT\|##INFO\|##FORMAT\|#CHROM\|##bcftools_callVersion\|##bcftools_callCommand\|##bcftools_concatVersion\|##bcftools_concatCommand\|##bcftools_filterVersion\|##bcftools_filterCommand\|##bcftools_viewVersion\|##bcftools_viewCommand\|SUPER_1\|SUPER_2\|SUPER_3\|SUPER_4\|SUPER_5\|SUPER_6\|SUPER_7\|SUPER_8\|SUPER_9\|SUPER_10\|SUPER_11_1\|SUPER_12\|SUPER_13' | sed 's@SUPER_11_1@SUPER_11@g'| sed 's@SUPER_@chr@g' >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.vcf
sed -i 's@_@-@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.vcf
bcftools sort /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.vcf -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted.vcf -Ov
counter=1
awk -v OFS='\t' '!/^#/ {$3="ID" counter; counter++} {print}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted.vcf >> temp.vcf
plink --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp.vcf --make-bed --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted
plink --bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted --indep-pairwise 50 10 0.1 --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted.pruned

#Pruned 1380987 variants from chromosome 1, leaving 152305.
#Pruned 1496916 variants from chromosome 2, leaving 157181.
#Pruned 1029480 variants from chromosome 3, leaving 100224.
#Pruned 1179591 variants from chromosome 4, leaving 128591.
#Pruned 1094794 variants from chromosome 5, leaving 110664.
#Pruned 1115119 variants from chromosome 6, leaving 110928.
#Pruned 464672 variants from chromosome 7, leaving 42312.
#Pruned 745949 variants from chromosome 8, leaving 68469.
#Pruned 635754 variants from chromosome 9, leaving 63848.
#Pruned 138209 variants from chromosome 10, leaving 14971.
#Pruned 561490 variants from chromosome 11, leaving 57352.
#Pruned 609082 variants from chromosome 12, leaving 60675.
#Pruned 300776 variants from chromosome 13, leaving 30267.
#Pruning complete.  10,752,819 of 11,850,606 variants removed.

plink --bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted --extract /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted.pruned.prune.in --make-bed --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted_pruned_set
#1,097,787 variants and 40 people pass filters and QC.
rm /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/*vcf* temp.vcf

for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40; do 
bedfile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted_pruned_set.bed
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/admixture
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_admixture_cross_validation.sh $bedfile $K $OutDir
done #59673320-358, 59679969-71, 59705950-2, 59720743-81
#Converged in 13 iterations (1590.46 sec) k=2

#59673355 = 1710 (QN/Block)

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/admixture/log*); do
CV=$(grep 'CV error' $file | sed 's@CV error (@@g'| sed 's@):@@g')
echo $CV
done

Pruned_vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted_pruned_set
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/admixture2
OutFile=apicales.genomicSNPs
Mink=2
Maxk=3
Bootstraps=200
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_admixture.sh $OutDir $OutFile $Pruned_vcf $Mink $Maxk $Bootstraps 
#59725006

###################################################################################################################################################################################

head -n 5 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.vcf
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf | grep -w '##ALT\|##INFO\|##FORMAT\|#CHROM\|##bcftools_callVersion\|##bcftools_callCommand\|##bcftools_concatVersion\|##bcftools_concatCommand\|##bcftools_filterVersion\|##bcftools_filterCommand\|##bcftools_viewVersion\|##bcftools_viewCommand\|SUPER_1\|SUPER_2\|SUPER_3\|SUPER_4\|SUPER_5\|SUPER_6\|SUPER_7\|SUPER_8\|SUPER_9\|SUPER_10\|SUPER_11\|SUPER_12\|SUPER_13' | sed 's@SUPER_@chr@g'  >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.vcf
sed -i 's@_@-@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.vcf
bgzip -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.vcf > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.vcf.gz
bcftools sort /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.vcf.gz -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted.vcf.gz -Oz
plink --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted.vcf.gz --make-bed --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted
plink --bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted --indep-pairwise 50 10 0.1 --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted.pruned

#Pruned 122236 variants from chromosome 1, leaving 2952.
#Pruned 138190 variants from chromosome 2, leaving 3473.
#Pruned 109219 variants from chromosome 3, leaving 2567.
#Pruned 89405 variants from chromosome 4, leaving 2034.
#Pruned 67343 variants from chromosome 5, leaving 1723.
#Pruned 81459 variants from chromosome 6, leaving 2055.
#Pruned 43422 variants from chromosome 7, leaving 1162.
#Pruned 69627 variants from chromosome 8, leaving 1636.
#Pruned 62067 variants from chromosome 9, leaving 1439.
#Pruned 3093 variants from chromosome 10, leaving 135.
#Pruned 67084 variants from chromosome 11, leaving 1588.
#Pruned 38037 variants from chromosome 12, leaving 987.
#Pruned 27882 variants from chromosome 13, leaving 861.
#Pruning complete.  919,064 of 941,676 variants removed.

plink --bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted --extract /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted.pruned.prune.in --make-bed --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted_pruned_set
#941,676 variants and 9 people pass filters and QC.

for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40; do 
bedfile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted_pruned_set.bed
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/admixture2
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_admixture_cross_validation.sh $bedfile $K $OutDir
done #59666987-59667025

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/admixture2/log*); do
CV=$(grep 'CV error' $file | sed 's@CV error (@@g'| sed 's@):@@g')
echo $CV
done

#Probably there is no admixture of different populations - all samples are from the same site + there are only 9 of them

Pruned_vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted_pruned_set.bed
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/admixture2
OutFile=pallida.genomicSNPs
Mink=2
Maxk=9
Bootstraps=200
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_admixture.sh $OutDir $OutFile $Pruned_vcf $Mink $Maxk $Bootstraps 
#59725003
```
#### STRUCTURE 
```bash
# Extract individual IDs and genotype data
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz --plink --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/structure/40_apicales

vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf.gz --plink --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/structure/09_pallida




# Convert PLINK files to STRUCTURE format
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/structure
source package /nbi/software/testing/bin/plink-1.9 
plink --file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/structure/193s.M_persicae --recode structure --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/filtered/structure/193s.M_persicae


srun -p jic-short  -c 1 --mem 50G --pty bash
source package /nbi/software/testing/bin/structure-2.3.4
source package /tgac/software/testing/bin/pgdspider-2.1.1.5
source jdk-1.7.0_25
java -Xmx40960m -Xms512M -jar /tgac/software/testing/bin/core/../..//pgdspider/2.1.1.5/x86_64/bin/PGDSpider2-cli.jar -inputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf -inputformat VCF -outputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.structure -outputformat STRUCTURE -spid ~/git_repos/Scripts/NBI/vcf-structure.spid
#edit population groups
awk '{print $1, $2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.structure > extracted_columns.txt
paste extracted_columns.txt <(cut -d' ' -f3- /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.structure) > temp_file && mv temp_file extracted_columns.txt
tr ' ' '\t' < extracted_columns.txt > temp_file && mv temp_file extracted_columns.txt


java -Xmx2048m -Xms512M -jar /tgac/software/testing/bin/core/../..//pgdspider/2.1.1.5/x86_64/bin/PGDSpider2-cli.jar -inputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf -inputformat VCF -outputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.structure -outputformat STRUCTURE -spid ~/git_repos/Scripts/NBI/vcf-structure.spid
#All 9 pallida are from the same location/population


```