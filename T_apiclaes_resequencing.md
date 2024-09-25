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
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/*/*.fq.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/*/*.fq.gz); do
zcat $file | echo $((`wc -l`/4))
done
```
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

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/apicales/*/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/*/); do
sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
coverage=$(grep 'mean coverageData' ${ReadDir}qualimap/*_genome_results_gff.txt | rev | cut -d ' ' -f1 | rev | sed 's@X@@g')
echo $sample raw reads have average coverage of ${coverage}
done

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/raw_data/Dyspera/pallida/*/); do
sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
coverage=$(grep 'mean coverageData' ${ReadDir}qualimap/*_genome_results.txt | rev | cut -d ' ' -f1 | rev | sed 's@X@@g')
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
#### Mitofinder
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/); do
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

source package a684a2ed-d23f-4025-aa81-b21e27e458df
seqret -sequence /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi_tant/final_mitogenome.fasta -feature -fformat gff -fopenfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi_tant/final_mitogenome.gff -osformat genbank -auto

mitofinder -j [seqid] -1 [left_reads.fastq.gz] -2 [right_reads.fastq.gz] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]   
mitofinder -j [seqid] -a [assembly.fasta] -r [genbank_reference.gb] -o [genetic_code] -p [threads] -m [memory]
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

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/); do
sample=$(echo $ReadDir | rev | cut -d '/' -f2 | rev)
coverage=$(grep 'mean coverageData' ${ReadDir}qualimap/*_genome_results_gff.txt | rev | cut -d ' ' -f1 | rev | sed 's@X@@g')
echo $sample raw reads have average coverage of ${coverage}
done

##########################################################################################################################################
#Cross species alignment
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/); do
    OutDir=$(echo ${ReadDir})/ant_map
    if ! find "${OutDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
    echo Running for:
    ls ${ReadDir}qualimap/*genome_results_gff.txt
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3
    else
    echo Already Done:
    ls ${ReadDir}qualimap/*genome_results.txt
    fi
done
#60087215-60087260

for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/); do
    OutDir=$(echo ${ReadDir})/api_map
    if ! find "${OutDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
    echo Running for:
    ls ${ReadDir}qualimap/*genome_results_gff.txt
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3/braker.gff3
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 
    else
    echo Already Done:
    ls ${ReadDir}qualimap/*genome_results.txt
    fi
done
#60087270-60087278

srun -p jic-long  -c 16 --mem 64G --pty bash
source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
source package 638df626-d658-40aa-80e5-14a275b7464b
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/ant_map/bwa-mem/*_trimmed_1_sorted.bam /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/api_map/bwa-mem/*_trimmed_1_sorted.bam); do
    OutDir=$(dirname $file)
    OutFile=$(basename $file | sed s'@_sorted.bam@_sorted_MarkDups.bam@g')
    metrics=$(basename $file | sed s'@_trimmed_1_sorted.bam@@g')_marked_dup_metrics.txt
    java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar SortSam I=${file} O=${OutDir}/temp.bam SORT_ORDER=coordinate
    java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${OutDir}/temp.bam O=${OutDir}/${OutFile} M=${OutDir}/${metrics}
    cd ${OutDir}
    rm temp.bam
    samtools index ${OutFile}
    cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/ant_map/bwa-mem/*_sorted_MarkDups.bam); do
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_realign.sh $file $Reference 
done
#60499118-157

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/api_map/bwa-mem/*_sorted_MarkDups.bam); do
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_realign.sh $file $Reference 
done
#60499159-167

source package 638df626-d658-40aa-80e5-14a275b7464b
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/bwa-mem/gatk/*realigned.bam /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/api_map/bwa-mem/gatk/*realigned.bam > bamlist3.txt
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/
bcftools mpileup --threads 16 -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/bamlist3.txt --annotate AD,DP --fasta-ref /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa -O z -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/variants_1.vcf.gz

bcftools call --threads 16 --ploidy 2 -Oz -v -m -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/variants.vcf.gz  /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/variants_1.vcf.gz

ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/bwa-mem/gatk/*realigned.bam /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/ant_map/bwa-mem/gatk/*realigned.bam > bamlist4.txt
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/
bcftools mpileup --threads 16 -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/bamlist4.txt --annotate AD,DP --fasta-ref /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa -O z -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/variants_1.vcf.gz

bcftools call --threads 16 --ploidy 2 -Oz -v -m -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/variants.vcf.gz  /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/variants_1.vcf.gz
#61923281,61923695

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/variants.vcf.gz | wc -l #49,648,480
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/variants.vcf.gz | wc -l #57,418,578

VCF=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/variants.vcf.gz
OutDir=$(dirname $VCF)/genmap
OutFile=$(basename $VCF | sed 's@.vcf.gz@@g')
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
Repeatmodeller=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa.out
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_genmap_mappability_masking.sh $OutDir $OutFile $Reference $Repeatmodeller $VCF
#63437460

VCF=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/variants.vcf.gz
OutDir=$(dirname $VCF)/genmap
OutFile=$(basename $VCF | sed 's@.vcf.gz@@g')
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
Repeatmodeller=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa.out
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_genmap_mappability_masking.sh $OutDir $OutFile $Reference $Repeatmodeller $VCF
#63437467

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/variants_callable.vcf.gz | wc -l #25,876,472
zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable.vcf.gz | wc -l #23,418,139

source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/variants_callable.vcf.gz --missing-indv
mv out.imiss /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/variants_callable.out.imiss 
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable.vcf.gz --missing-indv
mv out.imiss /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable.out.imiss 

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
import matplotlib.pyplot as plt

data = []
#with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/variants_callable.out.imiss', 'r') as file:
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable.out.imiss', 'r') as file:
    next(file)  # Skip the header line
    for line in file:
        values = line.strip().split('\t')
        f_miss = float(values[-1])
        data.append(f_miss)

plt.scatter(range(len(data)), data, marker='o', color='blue')
plt.xlabel('Individual')
plt.ylabel('Missing Level (F_MISS)')
plt.title('Missing Level per Individual')
#plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/variants_callable.out.imiss_plot.png')
plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable.out.imiss_plot.png')
exit()

vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/variants_callable.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/variants_callable_filtered
bgzip -k /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/variants_callable_filtered.recode.vcf
#After filtering, kept 49 out of 49 Individuals
#After filtering, kept 10,053,977 out of a possible 25,875,883 Sites

vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 40 --max-meanDP 40 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered
bgzip -k /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf
#After filtering, kept 49 out of 49 Individuals
#After filtering, kept 9,674,852 out of a possible 23,417,748 Sites

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/variants_callable_filtered.recode.vcf.gz); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/splitstree
mkdir $OutDir   
Outfile=p_dis.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#63454704

for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/splitstree
mkdir $OutDir   
Outfile=p_dis.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#63455422

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/mat2csv.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/splitstree/p_dis.mat
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/csv2dist.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales+/genmap/splitstree/p_dis_mperc.csv
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/mat2csv.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/splitstree/p_dis.mat
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/csv2dist.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/splitstree/p_dis_mperc.csv
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
```
Carsonella:
```bash
for ReadDir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/apicales/*/trim_galore/ /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/pallida/*/trim_galore/); do
    OutDir=$(echo ${ReadDir})/carsonella_map
    if ! find "${OutDir}/qualimap" -name "*genome_results_gff.txt" -print -quit | grep -q .; then
    echo Running for:
    ls ${ReadDir}qualimap/*genome_results_gff.txt
    Fread=$(ls ${ReadDir}*_1.fq.gz)
    Rread=$(ls ${ReadDir}*_2.fq.gz)
    Reference_genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna #this is a B.cockerelli symbiont strain
    Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/genomic.gff
    ProgDir=~/git_repos/Wrappers/NBI
    mkdir $OutDir
    sbatch $ProgDir/run_raw_read_qc.sh $OutDir $Reference_genome $Gff $Fread $Rread $Fread2 $Rread2 $Fread3 $Rread3
    else
    echo Already Done:
    ls ${ReadDir}qualimap/*genome_results.txt
    fi
done
#60088612-60088662

srun -p jic-medium  -c 4 --mem 16G --pty bash
source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
source package 638df626-d658-40aa-80e5-14a275b7464b
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/carsonella_map/bwa-mem/*_trimmed_1_sorted.bam); do
    OutDir=$(dirname $file)
    OutFile=$(basename $file | sed s'@_sorted.bam@_sorted_MarkDups.bam@g')
    metrics=$(basename $file | sed s'@_trimmed_1_sorted.bam@@g')_marked_dup_metrics.txt
    java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar SortSam I=${file} O=${OutDir}/temp.bam SORT_ORDER=coordinate
    java17 -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${OutDir}/temp.bam O=${OutDir}/${OutFile} M=${OutDir}/${metrics}
    cd ${OutDir}
    rm temp.bam
    samtools index ${OutFile}
    cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids
done

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/
java -jar /tgac/software/testing/bin/core/../..//picardtools/2.1.1/x86_64/bin/picard.jar CreateSequenceDictionary R=GCA_002009355.1_ASM200935v1_genomic.fna O=GCA_002009355.1_ASM200935v1_genomic.dict
samtools faidx T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa 
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/carsonella_map/bwa-mem/*_sorted_MarkDups.bam); do
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna 
ProgDir=/hpc-home/did23faz/git_repos/Wrappers/NBI
sbatch $ProgDir/run_realign.sh $file $Reference 
done 
#60148379-427

source package 638df626-d658-40aa-80e5-14a275b7464b
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/dna_qc/Dyspera/*/*/trim_galore/carsonella_map/bwa-mem/gatk/*realigned.bam > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/bamlist.txt
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/
bcftools mpileup --threads 16 -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/bamlist.txt --annotate AD,DP --fasta-ref /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Symbionts/Candidatus/Carsonella/ruddii/GCA_002009355.1/GCA_002009355.1_ASM200935v1_genomic.fna -O z -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants_1.vcf.gz

bcftools call --threads 16 --ploidy 2 -Oz -v -m -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants.vcf.gz /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants_1.vcf.gz
#63437768

zcat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants.vcf.gz | wc -l #19,065

source package /nbi/software/testing/bin/vcftools-0.1.15
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants.vcf.gz --missing-indv
mv out.imiss /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants.out.imiss  

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
import matplotlib.pyplot as plt

data = []
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants.out.imiss', 'r') as file:
    next(file)  # Skip the header line
    for line in file:
        values = line.strip().split('\t')
        f_miss = float(values[-1])
        data.append(f_miss)

plt.scatter(range(len(data)), data, marker='o', color='blue')
plt.xlabel('Individual')
plt.ylabel('Missing Level (F_MISS)')
plt.title('Missing Level per Individual')
plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants.out.imiss_plot.png')
exit()

vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants.vcf.gz --remove-indels --max-alleles 2 --mac 1 --max-missing 0.9 --minQ 30 --minDP 5 --min-meanDP 5 --maxDP 300 --max-meanDP 300 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants_filtered
bgzip -k /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants_filtered.recode.vcf
#After filtering, kept 49 out of 49 Individuals
#After filtering, kept 9226 out of a possible 19035 Sites

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/splitstree
for vcf in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/variants_filtered.recode.vcf.gz); do
ProgDir=~/git_repos/Wrappers/NBI
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/splitstree
Outfile=p_dis.mat
sbatch $ProgDir/run_VCF2Dis.sh $vcf $OutDir $Outfile
done
#63447010

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/mat2csv.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/splitstree/p_dis.mat
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/csv2dist.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Carsonella/ruddii/splitstree/p_dis_mperc.csv
source package 7654f72b-1692-46bb-9a56-443406d03fd9
SplitsTree
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
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.min4.fasta.* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/iqtree2/.

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/iqtree2
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.min4.fasta
cpu=12
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $Alignment -m PMB+F+R3 -B 1000 -T AUTO --threads-max $cpu
#59798091

################################################################################################################################################################################################################################
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf --thin 1000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.vcf

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/vcf2phylip.py -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.vcf.recode.vcf --output-folder /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/ -f -p

mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/iqtree2
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/iqtree2
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.min4.fasta
cpu=8
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $Alignment -m MF -T AUTO --threads-max $cpu
#59702173,59727696
#Best-fit model: PMB+F+R5 chosen according to BIC
#Analysis results written to:
#  IQ-TREE report:                /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.min4.fasta.iqtree
#  Tree used for ModelFinder:     /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.min4.fasta.treefile
#  Screen log file:               /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.min4.fasta.log
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.min4.fasta.* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/iqtree2/.

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/iqtree2
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.min4.fasta
cpu=12
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $Alignment -m PMB+F+R5 -B 1000 -T AUTO --threads-max $cpu
#59798091
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

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/plotly.sif python3
```
```python
from sklearn.manifold import MDS
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import plotly.express as px
import plotly.graph_objects as go


#Collect input
#df = pd.read_csv('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/splitstree/p_dis_mperc.csv', index_col=0)
df = pd.read_csv('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/splitstree/p_dis_mperc.csv', index_col=0)
acc =  df.index.tolist()

#MDS projection dimensionality
n_components = 3                                                                 
embedding = MDS(n_components=3,random_state=20348, dissimilarity='precomputed')  #dissimilarity parameter shows D is already known
nm_scores = embedding.fit_transform(df)

print('MDS projected coordinates matrix (nm_scores) has dimensions: ', np.shape(nm_scores))

# Collect sample info
geography = {}
head = ['ID', 'Country']
#with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/sample_info.csv') as inp:
with open('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/sample_info.csv') as inp:
    next(inp)  # skip header
    for line in inp:
        A = line.strip().split(',')
        d = dict(zip(head, A))
        geography[d['ID']] = d['Country']

colour_settings = {
    'Norway': 'Norway',
    'Finland': 'Finland',
    'Austria': 'Austria',
    'Scotland': 'Scotland'
}

# Color for geography
color_discrete_mapD = {
    'Norway': 'green',
    'Finland': 'red',
    'Austria': 'blue',
    'Scotland': 'orange'
}

# Assign colors based on geography
colour = [colour_settings[geography[x]]  for x in acc]
df3 = pd.DataFrame({
    'acc': acc, 
    'colour': colour, 
    'PC1': nm_scores[:,0], 
    'PC2': nm_scores[:,1], 
    'PC3': nm_scores[:,2]
})

# Plotting
fig = px.scatter_3d(df3, x='PC1', y='PC2', z='PC3', hover_data=['acc', 'colour'], color='colour', 
                    color_discrete_map=color_discrete_mapD, title='Country')
fig.update_traces(marker_size=4)
#fig.write_html("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/pallida/Country-python-MDS.html")
fig.write_html("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/analysis/phylogeny/Dyspera/apicales/Country-python-MDS.html")
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

# Convert PLINK files to STRUCTURE format
srun -p jic-medium  -c 16 --mem 64G --pty bash
source package /nbi/software/testing/bin/structure-2.3.4
source package /tgac/software/testing/bin/pgdspider-2.1.1.5
source jdk-1.7.0_25
java -Xmx40960m -Xms512M -jar /tgac/software/testing/bin/core/../..//pgdspider/2.1.1.5/x86_64/bin/PGDSpider2-cli.jar -inputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf -inputformat VCF -outputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.structure -outputformat STRUCTURE -spid ~/git_repos/Scripts/NBI/vcf-structure.spid
#edit population groups

#Downsample SNPs for structure analysis
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf --thin 1000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.vcf
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin1000.vcf.recode.vcf #294,085
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf --thin 2000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin2000.vcf
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin2000.vcf.recode.vcf #174810
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf --thin 4000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin4000.vcf
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin4000.vcf.recode.vcf #101299
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf --thin 10000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin10000.vcf
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin10000.vcf.recode.vcf #47376
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf --thin 50000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin50000.vcf
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin50000.vcf.recode.vcf #11640

#Convert to structure format
java -Xmx40960m -Xms512M -jar /tgac/software/testing/bin/core/../..//pgdspider/2.1.1.5/x86_64/bin/PGDSpider2-cli.jar -inputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin4000.vcf.recode.vcf -inputformat VCF -outputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin4000.structure -outputformat STRUCTURE -spid ~/git_repos/Scripts/NBI/vcf-structure.spid

java -Xmx40960m -Xms512M -jar /tgac/software/testing/bin/core/../..//pgdspider/2.1.1.5/x86_64/bin/PGDSpider2-cli.jar -inputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin50000.vcf.recode.vcf -inputformat VCF -outputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin50000.structure -outputformat STRUCTURE -spid ~/git_repos/Scripts/NBI/vcf-structure.spid
#edit population groups

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40; do 
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.thin50000.structure
Ploidy=2
Iterations=10
Burnin=100000
Reps=100000
OutDir=$(dirname $InFile)/structure
OutFile=D_apicales
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_structure.sh $InFile $Ploidy $K $Iterations $Burnin $Reps $OutDir $OutFile
done
#59748200-39

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40; do 
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted_pruned_set
Iterations=10
OutDir=$(dirname $InFile)/faststructure
OutFile=D_apicales
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_faststructure.sh $InFile $K $Iterations $OutDir $OutFile
done
#59747644-83

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40; do 
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/plink/variants_callable_filtered.recode.13.sorted
Iterations=10
OutDir=$(dirname $InFile)/faststructure-full
OutFile=D_apicales
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_faststructure.sh $InFile $K $Iterations $OutDir $OutFile
done
#59748002-42

#########################################################################################################################

# Extract individual IDs and genotype data
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf.gz --plink --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/structure/09_pallida

#Downsample SNPs for structure analysis
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf --thin 1000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin1000.vcf
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin1000.vcf.recode.vcf #170,532
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf --thin 2000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin2000.vcf
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin2000.vcf.recode.vcf #115,693
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf --thin 10000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin10000.vcf
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin10000.vcf.recode.vcf #38,473
vcftools --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf --thin 50000 --recode --recode-INFO-all --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin50000.vcf
wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin50000.vcf.recode.vcf #10,344

# Convert PLINK files to STRUCTURE format
java -Xmx2048m -Xms512M -jar /tgac/software/testing/bin/core/../..//pgdspider/2.1.1.5/x86_64/bin/PGDSpider2-cli.jar -inputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin2000.vcf.recode.vcf -inputformat VCF -outputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin2000.structure -outputformat STRUCTURE -spid ~/git_repos/Scripts/NBI/vcf-structure.spid

java -Xmx2048m -Xms512M -jar /tgac/software/testing/bin/core/../..//pgdspider/2.1.1.5/x86_64/bin/PGDSpider2-cli.jar -inputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin50000.vcf.recode.vcf -inputformat VCF -outputfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin50000.structure -outputformat STRUCTURE -spid ~/git_repos/Scripts/NBI/vcf-structure.spid
#All 9 pallida are from the same location/population

for K in 1 2 3 4 5 6 7 8 9; do 
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.thin50000.structure
Ploidy=2
Iterations=10
Burnin=100000
Reps=100000
OutDir=$(dirname $InFile)/structure
OutFile=D_pallida
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_structure.sh $InFile $Ploidy $K $Iterations $Burnin $Reps $OutDir $OutFile
done
#59748258-66

for K in 1 2 3 4 5 6 7 8 9; do 
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/plink/variants_callable_filtered.recode.13.sorted_pruned_set
Iterations=10
OutDir=$(dirname $InFile)/faststructure
OutFile=D_pallida
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_faststructure.sh $InFile $K $Iterations $OutDir $OutFile
done
#59747142-50
```
#### FST

0 to 0.05: Little to no differentiation. The populations are genetically similar.
0.05 to 0.15: Moderate differentiation.
0.15 to 0.30: High differentiation.
Above 0.30: Very high differentiation.

```bash
pwd #/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids

source package d37013e7-5691-40b6-8885-f029fe5fad54
source package c92263ec-95e5-43eb-a527-8f1496d56f1a

bcftools index /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz
bcftools index /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf.gz
bcftools index /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz

bcftools query -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz | grep 'TrAn\|TrA0' > uk.txt
bcftools query -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz | grep 'Tap' > austria.txt
bcftools query -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz | grep 'TrAp' > finland.txt
bcftools query -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz | grep 'Dyap' > norway.txt
bcftools query -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz | grep 'Dyap\|Tap' > austria+norway.txt
bcftools query -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz | grep 'Dyap\|TrAp\|Tap' > austria+norway+finland.txt

vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz --weir-fst-pop uk.txt --weir-fst-pop austria.txt --weir-fst-pop finland.txt --weir-fst-pop norway.txt --out all
#Weir and Cockerham mean Fst estimate: 0.19624
#Weir and Cockerham weighted Fst estimate: 0.40277
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz --weir-fst-pop austria+norway.txt --weir-fst-pop finland.txt --out aut+nor_v_fin
#Weir and Cockerham mean Fst estimate: 0.058454
#Weir and Cockerham weighted Fst estimate: 0.15864
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz --weir-fst-pop uk.txt --weir-fst-pop austria+norway+finland.txt --out api_v_ant
#Weir and Cockerham mean Fst estimate: 0.19313
#Weir and Cockerham weighted Fst estimate: 0.5625

vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz --weir-fst-pop austria.txt --weir-fst-pop norway.txt --out aut_v_nor
#Weir and Cockerham mean Fst estimate: 0.0099429
#Weir and Cockerham weighted Fst estimate: 0.023401

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3  ~/git_repos/Scripts/NBI/calculate_dxy.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz austria+norway.txt finland.txt
#Dxy: 0.28309154093301386
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3  ~/git_repos/Scripts/NBI/calculate_dxy.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz uk.txt austria+norway+finland.txt
#Dxy: 0.474510397938831

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3  ~/git_repos/Scripts/NBI/calculate_dxy.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz austria.txt norway.txt
#Dxy: 0.26980739608033727
```
```python
import sys
import pysam
import numpy as np

def compute_dxy(vcf_file, pop1_samples, pop2_samples):
    try:
        vcf = pysam.VariantFile(vcf_file, 'r')
    except Exception as e:
        print(f"Error opening VCF file: {e}")
        sys.exit(1)
    dxy_values = []
    for rec in vcf.fetch():
        # Extract genotypes for the two populations
        pop1_genotypes = [rec.samples[sample]['GT'] for sample in pop1_samples if sample in rec.samples]
        pop2_genotypes = [rec.samples[sample]['GT'] for sample in pop2_samples if sample in rec.samples]
        if len(pop1_genotypes) == 0 or len(pop2_genotypes) == 0:
            continue
        # Compute pairwise differences
        pairwise_differences = []
        for gt1 in pop1_genotypes:
            for gt2 in pop2_genotypes:
                if None in gt1 or None in gt2:
                    continue
                diff = np.sum(np.array(gt1) != np.array(gt2))
                pairwise_differences.append(diff)
        if pairwise_differences:
            dxy_values.append(np.mean(pairwise_differences))
    return np.mean(dxy_values) if dxy_values else float('nan')

def read_sample_ids(file_path):
    try:
        with open(file_path, 'r') as file:
            sample_ids = [line.strip() for line in file if line.strip()]
    except Exception as e:
        print(f"Error reading sample file {file_path}: {e}")
        sys.exit(1)
    return sample_ids

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <vcf_file> <pop1_samples> <pop2_samples>")
        sys.exit(1)
    vcf_file = sys.argv[1]
    group1 = sys.argv[2]
    group2 = sys.argv[3]
    pop1_samples = read_sample_ids(group1)
    pop2_samples = read_sample_ids(group2)
    dxy = compute_dxy(vcf_file, pop1_samples, pop2_samples)
    print(f"Dxy: {dxy}")

```
```bash
source package c92263ec-95e5-43eb-a527-8f1496d56f1a
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d

#Subset the groups:
bcftools view -S austria+norway.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/austria+norway.vcf
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/austria+norway.vcf
bcftools index /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/austria+norway.vcf.gz
bcftools view -S finland.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/finland.vcf
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/finland.vcf
bcftools index /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/finland.vcf.gz

#Add allele frequency data to the INFO field:
bcftools +fill-tags /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/finland.vcf.gz -- -t AF | bcftools view -Oz -o group1_af.vcf.gz
bcftools +fill-tags /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/austria+norway.vcf.gz -- -t AF | bcftools view -Oz -o group2_af.vcf.gz
bcftools index group1_af.vcf.gz
bcftools index group2_af.vcf.gz
srun -p jic-short --ntasks 1 --cpus-per-task 1 --mem 50G --pty bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/AF_difference.py group1_af.vcf.gz group2_af.vcf.gz AF_differences.tsv

awk '{if (NR > 1) print $5}' AF_differences.tsv | awk '
BEGIN {bins[0.05]=0; bins[0.1]=0; bins[0.15]=0; bins[0.2]=0; bins[0.25]=0; bins[0.3]=0; bins[0.35]=0; bins[0.4]=0; bins[0.45]=0; bins[0.5]=0}
{
  diff=$1
  if (diff <= 0.05) bins[0.05]++
  else if (diff <= 0.1) bins[0.1]++
  else if (diff <= 0.15) bins[0.15]++
  else if (diff <= 0.2) bins[0.2]++
  else if (diff <= 0.25) bins[0.25]++
  else if (diff <= 0.3) bins[0.3]++
  else if (diff <= 0.35) bins[0.35]++
  else if (diff <= 0.4) bins[0.4]++
  else if (diff <= 0.45) bins[0.45]++
  else if (diff <= 0.5) bins[0.5]++
  else if (diff <= 0.55) bins[0.55]++
  else if (diff <= 0.6) bins[0.6]++
  else if (diff <= 0.65) bins[0.65]++
  else if (diff <= 0.7) bins[0.7]++
  else if (diff <= 0.75) bins[0.75]++
  else if (diff <= 0.8) bins[0.8]++
  else if (diff <= 0.85) bins[0.85]++
  else if (diff <= 0.9) bins[0.9]++
  else if (diff <= 0.95) bins[0.95]++
  else if (diff <= 1.0) bins[1.0]++
}
END {
  for (bin in bins) print bin, bins[bin]
}' > histogram_data.txt

gnuplot << EOF
set terminal dumb size 160, 50
set title "Histogram of AF Differences"
set xlabel "Difference Bin"
set ylabel "Number of SNPs"
set boxwidth 0.03
set style fill solid
set xrange [0:1.0]
set yrange [0:8e+06]  # Adjust as necessary
set xtics 0.05
binwidth = 0.05
# Adjust the bin function to ensure it aligns correctly with the bins
bin(x, width) = width * floor(x / width + 0.5)
plot 'histogram_data.txt' using (bin(\$1, binwidth)):(\$2) with boxes title 'Number of SNPs'
EOF

#                                                                            Histogram of AF Differences
#
#                8e+06 +-------------------------------------------------------------------------------------------------------------------------------------+
#                      |      +     +      +      +      +     +      +      +     +      +      +     +      +      +      +     +      +      +     +      |
#                      |                                                                                                              Number of SNPs ******* |
#                      |                                                                                                                                     |
#                      |                                                                                                                                     |
#                7e+06 |-+                                                                                                                                 +-|
#                      |                                                                                                                                     |
#                      |                                                                                                                                     |
#                      |                                                                                                                                     |
#                      |                                                                                                                                     |
#                      |                                                                                                                                     |
#                6e+06 |-+  *****                                                                                                                          +-|
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                5e+06 |-+  *   *                                                                                                                          +-|
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                4e+06 |-+  *   *                                                                                                                          +-|
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                3e+06 |-+  *   *                                                                                                                          +-|
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                      |    *   *                                                                                                                            |
#                      |    *   * *****                                                                                                                      |
#                2e+06 |-+  *   * *   *                                                                                                                    +-|
#                      |    *   * *   *                                                                                                                      |
#                      |    *   * *   *                                                                                                                      |
#                      |    *   * *   *                                                                                                                      |
#                      |    *   * *   *                                                                                                                      |
#                      |    *   * *   *                                                                                                                      |
#                1e+06 |-+  *   * *   *  *****                                                                                                             +-|
#                      |    *   * *   *  *   *  *****                                                                                                        |
#                      |    *   * *   *  *   *  *   * ******                                                                                                 |
#                      |    *   * *   *  *   *  *   * *    * *****  *****                                                                                    |
#                      |    * + * * + *  * + *  * + * *  + * * + *  * + *  ***** *****  *****  *****   +      +      +      +     +      +      +     +      |
#                    0 +-------------------------------------------------------------------------------------------------------------------------------------+
#                      0     0.05  0.1    0.15   0.2    0.25  0.3    0.35   0.4   0.45   0.5    0.55  0.6    0.65   0.7    0.75  0.8    0.85   0.9   0.95    1
#                                                                                  Difference Bin

gnuplot << EOF
set terminal dumb size 160, 50
set title "Histogram of AF Differences"
set xlabel "Difference Bin"
set ylabel "Number of SNPs"
set logscale y 
set boxwidth 0.03
set style fill solid
set xrange [0:1.0]
set yrange [0:1e+08]  # Adjust as necessary
set xtics 0.05
binwidth = 0.05
# Adjust the bin function to ensure it aligns correctly with the bins
bin(x, width) = width * floor(x / width + 0.5)
plot 'histogram_data.txt' using (bin(\$1, binwidth)):(\$2) with boxes title 'Number of SNPs'
EOF


#                                                                            Histogram of AF Differences
#
#                 1e+08 +------------------------------------------------------------------------------------------------------------------------------------+
#                       |+     +     +      +      +     +      +      +     +      +     +      +      +     +      +      +     +      +      +     +     +|
#                       |+                                                                                                            Number of SNPs *******+|
#                       |+                                                                                                                                  +|
#                       |                                                                                                                                    |
#                 1e+07 |-+                                                                                                                                +-|
#                       |+   *****                                                                                                                          +|
#                       |+   *   *                                                                                                                          +|
#                       |+   *   * *****                                                                                                                    +|
#                       |    *   * *   *                                                                                                                     |
#                 1e+06 |-+  *   * *   *  *****  *****                                                                                                     +-|
#                       |+   *   * *   *  *   *  *   * *****                                                                                                +|
#                       |+   *   * *   *  *   *  *   * *   *  *****  *****                                                                                  +|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *****  *****                                                                     +|
#                100000 |-+  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  ****  *****                                                       +-|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *****                                                 +|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *****  *****                                    +|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *****                             +|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *****                       +|
#                 10000 |-+  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *****               +-|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *                +|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *****         +|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *****   +|
#                       |    *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *    |
#                  1000 |-+  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  +-|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  **|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |    *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  * |
#                   100 |-+  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *-|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                    10 |-+  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *-|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                     1 |-+  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *-|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |+   *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *   *  *  *  *   *  *   * *   *  *   *  *   * *   *  *   *  *   * *   *  *+|
#                       |+   * + * * + *  * + *  * + * * + *  * + *  * + * * + *  * + *  *+ *  * + *  * + * * + *  * + *  * + * * + *  * + *  * + * * + *  *+|
#                   0.1 +------------------------------------------------------------------------------------------------------------------------------------+
#                       0     0.05  0.1    0.15   0.2   0.25   0.3    0.35  0.4    0.45  0.5    0.55   0.6   0.65   0.7    0.75  0.8    0.85   0.9   0.95    1
#                                                                                  Difference Bin

awk 'NR==1 || $5 > 0.25' AF_differences.tsv | awk 'NR>1 {print $1 "\t" $2}' > positions_to_keep.txt
awk 'NR==FNR { positions[$1, $2]; next } /^#/ || ($1, $2) in positions' positions_to_keep.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered-differentiating.recode.vcf
```
```bash
awk -F'\t' '$3 > 0.3' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/aut+nor_v_fin.weir.fst > aut+nor_v_fin_meaningful.windowed.fst
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz --weir-fst-pop austria+norway.txt --weir-fst-pop finland.txt --fst-window-size 10000 --fst-window-step 5000 --out aut+nor_v_fin_10000
awk -F'\t' '$5 > 0.3' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/aut+nor_v_fin_10000.windowed.weir.fst > aut+nor_v_fin_10000_meaningful.windowed.fst
vcftools --gzvcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz --weir-fst-pop austria+norway.txt --weir-fst-pop finland.txt --fst-window-size 100000 --fst-window-step 50000 --out aut+nor_v_fin_100000
awk -F'\t' '$5 > 0.3' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/aut+nor_v_fin_100000.windowed.weir.fst > aut+nor_v_fin_100000_meaningful.windowed.fst
```
```bash
srun -p jic-medium --ntasks 1 --cpus-per-task 8 --mem 32G --pty bash
source package /nbi/software/testing/bin/plink-1.9
source package c432148d-eabc-4469-8b5f-f4755f26292b

#convert to plink
plink --double-id --allow-extra-chr --vcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf --make-bed --out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.plink

#edit phenotype data to consider the two groups
awk '
  BEGIN { OFS="\t" }
  /Dyap|Tap/ { $6 = 1 }
  /TrAp/ { $6 = 2 }
  { print }
' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.plink.fam > temp.fam && mv temp.fam /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.plink.fam

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/
gemma -bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.plink -gk 1 -o relatedness_matrix
mv output/relatedness_matrix* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/.
rm -rf output
gemma -bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.plink -k /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/relatedness_matrix.cXX.txt -lmm 4 -o gwas_output
mv output/gwas_output* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/.
rm -rf output
gemma -bfile /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.plink -k /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/relatedness_matrix.cXX.txt -lm 4 -o gwas_output_2
mv output/gwas_output_2* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/.
rm -rf output
```
When the sample size is small to moderate, the Wald test is the least reliable of the three tests. We should not trust it for such a small n as in this example (n = 10). Likelihood-ratio inference and score-test based inference are better in terms of actual error probabilities coming close to matching nominal levels. A marked divergence in the values of the three statistics indicates that the distribution of the ML estimator may be far from normality. In that case, small-sample methods are more appropriate than large-sample methods. - Agresti, A. (2007). An introduction to categorical data analysis (2nd edition). Hoboken, NJ: John Wiley & Sons.
```python
import pandas as pd
from statsmodels.stats.multitest import multipletests
import numpy as np
import matplotlib.pyplot as plt

gwas_results = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/gwas_output.assoc.txt", delim_whitespace=True)
gwas_results['fdr_p'] = multipletests(gwas_results['p_score'], method='fdr_bh')[1]
gwas_results.to_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/gwas_output_corrected.tsv", sep='\t', index=False)
gwas_results['-log10(p)'] = -np.log10(gwas_results['fdr_p'])

selected_chromosomes = ['SUPER_1', 'SUPER_2', 'SUPER_3', 'SUPER_4', 'SUPER_5',
                        'SUPER_6', 'SUPER_7', 'SUPER_8', 'SUPER_9', 'SUPER_10',
                        'SUPER_11_1', 'SUPER_12', 'SUPER_13']
gwas_results = gwas_results.sort_values(['chr', 'ps'])
fig, axes = plt.subplots(nrows=len(selected_chromosomes), ncols=1, figsize=(30, 2 * len(selected_chromosomes)), sharex=True)
colorlist = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
             '#393b79', '#FF8C00', '#637939', '#843c39', '#7b4173', '#8c6d31', '#666666']
for i, chr_id in enumerate(selected_chromosomes):
    ax = axes[i]
    subset = gwas_results[gwas_results['chr'] == chr_id]
    color = colorlist[i % len(colorlist)]  
    ax.scatter(subset['ps'], subset['-log10(p)'], color=color, s=1)  
    ax.axhline(y=1.30103, color='grey', linestyle='--', linewidth=1)  
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'Chromosome {chr_id}')

axes[-1].set_xlabel('Genomic Position')
plt.tight_layout()
plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/manhattan_plot.png', dpi=300)

############################################################################################################################

gwas_results = pd.read_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/gwas_output_2.assoc.txt", delim_whitespace=True)
gwas_results['fdr_p'] = multipletests(gwas_results['p_score'], method='fdr_bh')[1]
gwas_results.to_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/gwas_output_corrected_2.tsv", sep='\t', index=False)
gwas_results['-log10(p)'] = -np.log10(gwas_results['fdr_p'])

selected_chromosomes = ['SUPER_1', 'SUPER_2', 'SUPER_3', 'SUPER_4', 'SUPER_5',
                        'SUPER_6', 'SUPER_7', 'SUPER_8', 'SUPER_9', 'SUPER_10',
                        'SUPER_11_1', 'SUPER_12', 'SUPER_13']
gwas_results = gwas_results.sort_values(['chr', 'ps'])
fig, axes = plt.subplots(nrows=len(selected_chromosomes), ncols=1, figsize=(30, 2 * len(selected_chromosomes)), sharex=True)
colorlist = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
             '#393b79', '#FF8C00', '#637939', '#843c39', '#7b4173', '#8c6d31', '#666666']
for i, chr_id in enumerate(selected_chromosomes):
    ax = axes[i]
    subset = gwas_results[gwas_results['chr'] == chr_id]
    color = colorlist[i % len(colorlist)]  
    ax.scatter(subset['ps'], subset['-log10(p)'], color=color, s=1)  
    ax.axhline(y=1.30103, color='grey', linestyle='--', linewidth=1)  
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'Chromosome {chr_id}')

axes[-1].set_xlabel('Genomic Position')
plt.tight_layout()
plt.savefig('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/manhattan_plot_2.png', dpi=300)
```
```bash
awk 'NR == 1 || $16 < 0.05' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/gwas_output_corrected.tsv | awk 'NR>1 {print $1 "\t" $3}' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/positions_to_keep.txt
wc -l  /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/positions_to_keep.txt #114
awk 'NR==FNR { positions[$1, $2]; next } /^#/ || ($1, $2) in positions' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/positions_to_keep.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating.recode.vcf
```
```bash
srun -p jic-medium --ntasks 1 --cpus-per-task 8 --mem 32G --pty bash

source package /tgac/software/testing/bin/bedops-2.2.0
source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64
source package /tgac/software/production/bin/tabix-0.2.6

genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.gff3
vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered-differentiating.recode.vcf.gz
Species=Dyspera_apicalis
cpu=8
grep -v "#" $gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > myData.gff.gz
tabix -p gff myData.gff.gz
mkdir $(dirname $vcf)/vep
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/vep.sif vep -i $vcf --fork $cpu --species $Species --gff myData.gff.gz --fasta $genome --vcf --everything --output_file $(dirname $vcf)/vep/variant_effect_diff_output.txt


genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.gff3
vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz
Species=Dyspera_apicalis
cpu=8
grep -v "#" $gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > myData.gff.gz
tabix -p gff myData.gff.gz
mkdir $(dirname $vcf)/vep
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/vep.sif vep -i $vcf --fork $cpu --species $Species --gff myData.gff.gz --fasta $genome --vcf --everything --output_file $(dirname $vcf)/vep/variant_effect_output.txt

genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.gff3
vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida/genmap/variants_callable_filtered.recode.vcf.gz
Species=Dyspera_pallida
cpu=8
grep -v "#" $gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > myData.gff.gz
tabix -p gff myData.gff.gz
mkdir $(dirname $vcf)/vep
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/vep.sif vep -i $vcf --fork $cpu --species $Species --gff myData.gff.gz --fasta $genome --vcf --everything --output_file $(dirname $vcf)/vep/variant_effect_output.txt

genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa
gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.gff3
vcf=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/variants_callable_filtered.recode.vcf.gz
Species=Dyspera_pallida
cpu=8
grep -v "#" $gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > myData.gff.gz
tabix -p gff myData.gff.gz
mkdir $(dirname $vcf)/vep
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/vep.sif vep -i $vcf --fork $cpu --species $Species --gff myData.gff.gz --fasta $genome --vcf --everything --output_file $(dirname $vcf)/vep/variant_effect_output.txt









source package 3e7beb4d-f08b-4d6b-9b6a-f99cc91a38f9
source package 4c883633-af2d-4fac-ab67-a1574f7fe079

echo " " >> ~/git_repos/temp/snpEff.config
echo "# Dyspera pallida genome, version 1_0" >> ~/git_repos/temp/snpEff.config
echo "D_pallida_1_0.genome : Cow_parsley_psyllid" >> ~/git_repos/temp/snpEff.config
echo " " >> ~/git_repos/temp/snpEff.config
echo "# Dyspera apicalis genome, version 1_0" >> ~/git_repos/temp/snpEff.config
echo "D_apicalis_1_0.genome : Carrot_psyllid" >> ~/git_repos/temp/snpEff.config

OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/pallida+/genmap/snpEff
mkdir -p ${OutDir}/data/D_pallida_1_0
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.gff3 ${OutDir}/data/D_pallida_1_0/genes.gff
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.gtf ${OutDir}/data/D_pallida_1_0/genes.gtf
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.aa ${OutDir}/data/D_pallida_1_0/protein.fa
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.codingseq ${OutDir}/data/D_pallida_1_0/cds.fa
mkdir ${OutDir}/data/genomes
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected.fa ${OutDir}/data/genomes/D_pallida_1_0.fa
cd ${OutDir}/data/D_pallida_1_0
gzip genes.gff
gzip genes.gtf
cd ../..
cp ~/git_repos/temp/snpEff.config .
java17 -jar ~/git_repos/Scripts/NBI/snpEff.jar build -gff3 -v D_pallida_1_0
java17 -jar ~/git_repos/Scripts/NBI/snpEff.jar build -gtf22 -v D_pallida_1_0

OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/snpEff
mkdir -p ${OutDir}/data/D_apicalis_1_0
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.gff3 ${OutDir}/data/D_apicalis_1_0/genes.gff
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.gtf ${OutDir}/data/D_apicalis_1_0/genes.gtf
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.aa ${OutDir}/data/D_apicalis_1_0/protein.fa
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.codingseq ${OutDir}/data/D_apicalis_1_0/cds.fa
mkdir ${OutDir}/data/genomes
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected.fa ${OutDir}/data/genomes/D_apicalis_1_0.fa
cd ${OutDir}/data/D_apicalis_1_0
gzip genes.gff
gzip genes.gtf
cd ../..
cp ~/git_repos/temp/snpEff.config .
java17 -jar ~/git_repos/Scripts/NBI/snpEff.jar build -gff3 -v D_apicalis_1_0
java17 -jar ~/git_repos/Scripts/NBI/snpEff.jar build -gtf22 -v D_apicalis_1_0

#All:
java17 -Xmx8g -jar ~/git_repos/Scripts/NBI/snpEff.jar D_apicalis_1_0 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.vcf.gz > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.ann.vcf

awk -F'\t' 'BEGIN { OFS="\t" } /^#/ { next } { n = split($8, fields, "|"); for (i = 1; i <= n; i++) { if (fields[i] ~ /^g.*\.t1$/) { start = i - 5; end = i + 9; if (start < 1) start = 1; if (end > n) end = n; for (j = start; j <= end; j++) { if (j > start) printf "|"; printf "%s", fields[j] } printf "\n" } } }' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/variants_callable_filtered.recode.ann.vcf | sed 's@||@|.|@g' | sed 's@|@\t@g' | awk -F'\t' '$1 == "downstream_gene_variant" || $1 == "upstream_gene_variant" || $1 == "synonymous_variant" || $1 == "stop_lost&splice_region_variant" || $1 == "stop_gained&splice_region_variant" || $1 == "stop_gained" || $1 == "start_lost" || $1 == "splice_region_variant&synonymous_variant" || $1 == "splice_region_variant&stop_retained_variant" || $1 == "splice_region_variant&intron_variant" || $1 == "splice_donor_variant&intron_variant" || $1 == "splice_acceptor_variant&intron_variant" || $1 == "missense_variant&splice_region_variant" || $1 == "missense_variant" || $1 == "intron_variant" || $1 == "initiator_codon_variant"'> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/affected_genes.tsv
awk -F'\t' '{print $6}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/affected_genes.tsv | sort | uniq | wc -l #16,442 genes affected by SNPs
grep 'MODERATE\|HIGH' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/affected_genes.tsv | awk -F'\t' '{print $6}' | sort | uniq | wc -l #12,811 genes highly affected by SNPs

#GEMMA significant SNPs
bgzip /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating.recode.vcf
bcftools index /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating.recode.vcf.gz
java17 -Xmx8g -jar ~/git_repos/Scripts/NBI/snpEff.jar D_apicalis_1_0 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating.recode.vcf > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating.recode.ann.vcf
awk -F'\t' 'BEGIN { OFS="\t" } /^#/ { next } { n = split($8, fields, "|"); for (i = 1; i <= n; i++) { if (fields[i] ~ /^g.*\.t1$/) { start = i - 5; end = i + 9; if (start < 1) start = 1; if (end > n) end = n; for (j = start; j <= end; j++) { if (j > start) printf "|"; printf "%s", fields[j] } printf "\n" } } }' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating.recode.ann.vcf | sed 's@||@|.|@g' | sed 's@|@\t@g' | awk -F'\t' '$1 == "downstream_gene_variant" || $1 == "upstream_gene_variant" || $1 == "synonymous_variant" || $1 == "stop_lost&splice_region_variant" || $1 == "stop_gained&splice_region_variant" || $1 == "stop_gained" || $1 == "start_lost" || $1 == "splice_region_variant&synonymous_variant" || $1 == "splice_region_variant&stop_retained_variant" || $1 == "splice_region_variant&intron_variant" || $1 == "splice_donor_variant&intron_variant" || $1 == "splice_acceptor_variant&intron_variant" || $1 == "missense_variant&splice_region_variant" || $1 == "missense_variant" || $1 == "intron_variant" || $1 == "initiator_codon_variant"'> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating_affected_genes.tsv
awk -F'\t' '{print $6}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating_affected_genes.tsv | sort | uniq | wc -l #59 genes affected by SNPs
grep 'MODERATE\|HIGH' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating_affected_genes.tsv | awk -F'\t' '{print $6}' | sort | uniq | wc -l #8 genes highly affected by SNPs
#Missense variants:
#g10306.t1 - hypothetical protein
#g11658.t1 - calcineurin-binding
#g11670.t1 - speckle targeted PIP5k1a-regulated ploy(A) polymerase
#g14397.t1 - CLIP associating protein
#g16476.t1 - protein kintoun isoform
#g2979.t1 - zinc finger MYM-type protein 1-like
#g360.t1 - IAA leucine resistant/putative coat protein
#g6349.t1 - cell surface glycoprotein/zonadhesin isoform

#Upstream/downstream (5kb)/intron variants:
#g10307.t1 - tectonic-1 isoform
#g10346.t1 - bicaudal D-related protein homolog
#g11570.t1 - double homeobox protein A
#g11571.t1 - double homeobox protein A
#g11573.t1 - double homeobox protein 
#g11579.t1 - phospholipase A2 inhibitor/trophoblast glycoprotein
#g11640.t1 - hypothetical protein
#g11655.t1 - large proline rich protein
#g11659.t1 - hypothetical protein
#g11671.t1 - nitric oxide synthase-interacting
#g11680.t1 - importin-4-like isoform
#g11886.t1 - proton-coupled amino acid transporter-like
#g130.t1 - polymerase PA
#g13181.t1 - rho GTPase-activating protein 26 isoform
#g131.t1 - hormone sensistive lipase like/serine peptidase/rna polymerase/ribonuclease/chaplin
#g133.t1 - piggybac transposable element derived protein
#g134.t1 - hypothetical protein
#g135.t1 - hypothetical protein
#g13694.t1 - cathepsin L
#g13695.t1 - lipase
#g136.t1 - hypothetical protein
#g13957.t1 - ELYS isoform X1
#g157.t1 - hypothetical protein
#g158.t1 - hypothetical protein
#g159.t1 - glycoprotein
#g16477.t1 - odorant receptor/hypothetical protein
#g1816.t1 - transmembrane protein
#g312.t1 - dna direct rna polymerase subunit/hypothetical protein
#g324.t1 - hypothetical protein
#g361.t1 - IAA leucine resistant like/hypothetical protein
#g362.t1 - no hits
#g4074.t1 - beta-alanine transporter
#g41.t1 - putative nuclease
#g4340.t1 - hypothetical protein
#g4345.t1 - protein lifeguard 1-like isoform
#g4815.t1 - F-box/LRR-repeat protein 2/7 isoform
#g4819.t1 - vacuolar protein sorting-associated
#g4820.t1 - neuroligin-4
#g6372.t1 - peptidyl-prolyl cis-trans isomerase
#g6531.t1 - facilitated trehalose transporter
#g6532.t1 - transcription factor
#g66.t1 - hypothetical protein
#g8131.t1 - protein NDNF
#g8277.t1 - protein toll isoform
#g8278.t1 - hypothetical protein
#g8605.t1 - papilin isoform
#g8609.t1 - vanin-like protein
#g8785.t1 - unc-112-related protein
#g8871.t1 - nephrin isoform/hemicentin-1 isoform/neural cell adhesion molecule
#g8882.t1 - translation initiation factor
#g9313.t1 - colorectal mutant cancer protein


awk -F'\t' '{print $6}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating_affected_genes.tsv | sort -u > search.txt
grep -F -f search.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/interproscan/braker.aa.tsv > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating_affected_genes-annotations.tsv
awk -F'\t' '{print $14}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating_affected_genes-annotations.tsv | tr '|' '\n' | sort | uniq
#GO:0000981 DNA-binding transcription factor activity, RNA polymerase II-specific
#GO:0003676 nucleic acid binding
#GO:0003677 DNA binding
#GO:0003755 peptidyl-prolyl cis-trans isomerase activity
#GO:0003953 NAD+ nucleosidase activity
#GO:0004867 serine-type endopeptidase inhibitor activity
#GO:0005509 calcium ion binding
#GO:0005515 protein binding
#GO:0005576 extracellular region
#GO:0005737 cytoplasm
#GO:0005828 kinetochore microtubule
#GO:0005881 cytoplasmic microtubule
#GO:0006336 chromatin organization
#GO:0006355 regulation of DNA-templated transcription
#GO:0006508 proteolysis
#GO:0006606 protein import into nucleus
#GO:0006807 obsolete nitrogen compound metabolic process
#GO:0006886 intracellular protein transport
#GO:0007165 signal transduction 
#GO:0007229 integrin-mediated signaling pathway
#GO:0008234 cysteine-type peptidase activity
#GO:0016020 membrane
#GO:0016021 membrane
#GO:0016811 hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides
#GO:0019005 SCF ubiquitin ligase complex
#GO:0022857 transmembrane transporter activity
#GO:0030414 peptidase inhibitor activity
#GO:0031146 SCF-dependent proteasomal ubiquitin-dependent protein catabolic process
#GO:0031267 small GTPase binding
#GO:0043515 kinetochore binding
#GO:0046983 protein dimerization activity
#GO:0051010 microtubule plus-end binding
#GO:0055085 transmembrane transport
#GO:0061630 ubiquitin protein ligase activity
#GO:0070286 axonemal dynein complex assembly

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file search.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/braker3+helixer+tsebra_loose/braker.aa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/resequencing_psyllids/snp_calling/Dyspera/apicales/genmap/gemma/variants_callable_filtered-differentiating_affected_genes.aa

```
