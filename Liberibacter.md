# Candidatus liberibacter solanacearum assembly

Investigate scaffold 18 of the T.apicales assembly:

#### blast
```bash
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa | grep -A 1 -w '>scaffold_18' > temp_scaffold_18.fasta

Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
sbatch $ProgDir/run_blastn.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta $Database /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae temp_18
#58066323
```
Blast results suggest that scaffold 18 is C.liberibacter, however the scaffold is 1.9Mb vs the reported genome size of 1.2Mb

To investigate whether the scaffold contained missassembled apicales and liberibacter sequences merged together, hemiptera and liberibacter genes were blasted and hits plotted across the scaffold.

```bash
#Best hemiptera protein predictions were collectedf from NCBI and internally:
for file in $(ls /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/*/*/*aa.fa /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Cinara_cedri/v1/cinced3A.pep.fa /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Rhopalosiphum_maidis/v1/rmaidis_v2.gff3.prot.fa /jic/research-groups/Saskia-Hogenhout/Tom_Mathers/aphid_genomes_db/Schlechtendalia_chinensis/v1/proteins.fasta | grep -v 'Sitobion_avenae_JIC1_v2.scaffolds.braker.aa.fa\|Metopolophium_dirhodum_UK035_v1.scaffolds.braker.aa.fa\|Daktulosphaira_vitifoliae_INRAPcf7_v4.scaffolds.braker.aa.fa\|Aphis_glycines_4.v2.1.scaffolds.fa.gff.aa.fa'); do
x=$(echo $file | cut -d '/' -f8)
y=$(echo $file | cut -d '/' -f7 | sed 's@_@/@g')
cp $file ../Genomes/$y/$x/protein.faa
done

#Species IDs were added to the fasta headers:
for file in $(ls ../Genomes/*/*/*/protein.faa); do
x="[$(echo $file | cut -d '/' -f3,4 | sed 's@/@ @g')]"
echo $x
sed -i "/^>/ s/$/ $x/" $file
done

#hemiptera blast database prepared:
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_18/blastx
for file in $(ls ../Genomes/*/*/*/protein.faa | grep -v 'Frankiniella\|Thrips\|Megalurothrips'); do
echo $file
head -n 2 $file
#awk '/^>/ {split($1,a,/\[/); printf "%s|%s\n", a[1], substr($0,index($0,"[")+1,length($0)-index($0,"[")-1); next} {print}' $file > output.fasta && mv output.fasta $file
cat $file >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_18/blastx/hemip.faa
done

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_18/blastx

sed -i 's/|/_/g; s/ /_/g; s/./_/g' hemip.faa
sed -i 's/\.//g' output.fasta

source package 37f0ffda-9f66-4391-87e2-38ccd398861d
makeblastdb -in output.fasta -input_type fasta -dbtype prot -title hemip -out hemip


#blast for hemiptera genes:
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_18/blastx/hemip
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_18/blastx
Max=9999999
OutFile=hemip_18
sbatch ~/git_repos/Wrappers/NBI/run_blastx.sh $InFile $Database $OutDir $OutFile $Max
#58089671

#make liberibacter blast database:
makeblastdb -in /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/protein.faa -input_type fasta -dbtype prot -title Caliber -out Caliber

#blast for liberibacter genes:
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_18/blastx/Caliber
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_18/blastx
Max=9999999
OutFile=Caliber_18
sbatch ~/git_repos/Wrappers/NBI/run_blastx.sh $InFile $Database $OutDir $OutFile $Max
#58081242
```
Blast results reveal that liberibacter gene hits are flanked by hemiptera gene hits, however, many of the hemiptera hits have low % sequence match with the scaffold, of the 15,602 hemiptera hits only 136 have >70% match versus 998/1025 for liberibacter hits. To my eye, those hemiptera hits that have >70% sequence identity match do not look like they correspond to non-liberibacter regions in the nucmer plots.

#### BUSCO coverage

To further investigate whether scaffold 18 is tryuely a contaminant and not an apicales sequence with liberibacter regions integrated by missassembly, coverage of BUSCO genes was investiagted:
```bash
#Coverage of bacterial BUSCOs in scaffold 18:
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_18/BUSCO
Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/bacteria_odb10
OutFile=$(basename $Genome | sed 's@.fasta@@g')_$(echo $Database | cut -d '/' -f7)
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_busco_keep.sh "$Genome" "$Database" "$OutDir" "$OutFile" 
#58089630

Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_18/BUSCO/single_copy_busco_sequences.fna
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.fna@@g')
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2 #58089789

#Coverage of hemiptera BUSCOs in apicales scaffolds:
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa | grep -A 1 -w '>scaffold_1\|>scaffold_2\|>scaffold_3\|>scaffold_4\|>scaffold_5' > temp_scaffold_12345.fasta

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_12345.fasta
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_12345/BUSCO
Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
OutFile=$(basename $Genome | sed 's@.fasta@@g')_$(echo $Database | cut -d '/' -f7)
ProgDir=~/git_repos/Wrappers/NBI
mkdir -p $OutDir
sbatch $ProgDir/run_busco_keep.sh "$Genome" "$Database" "$OutDir" "$OutFile" 
#58089781

Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff_12345/BUSCO/single_copy_busco_sequences.fna
ProgDir=~/git_repos/Wrappers/NBI
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.fna@@g')
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2 #58090030
```
Coverage of BUSCOs in scaffold 18 was found to be orders of magnitude higher than coverage of BUSCOs in apicales scaffolds, corroborating blast results.

#### Nucmer

The reference C.liberibacter solanacearum assembly was aligned to scaffold_18 to investigate possilbe merging of multiple haplotypes:
```bash
#Reference vs the scaffold
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
Query=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger
OutFile=scaff18
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58078409
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18.delta

#The scaffold vs the reference
Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna
Query=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger
OutFile=scaff18-2
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58089601
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-2.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-2.delta
```
Alignments revealed that the entire C.liberibacter solanacearum reference genome is contained within scaffold 18, but ~700Mb from scaffold 18 are not present in the reference genome, ie. the increased size is not the result of merged haplotypes/duplication.

Haplotype C, these assemblies are not whole genome:
```bash
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
Query=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic.fna
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger
OutFile=scaff18-c
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58500651
source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-c.delta
mummerplot -color --layout /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-c.delta


Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
Query=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic.fna
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger
OutFile=scaff18-c2
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58500984
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-c2.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-c2.delta

#FIN114
nucmer --maxmatch --nosimplify /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta  -p temp
mummerplot -color --layout temp.delta
mummerplot -l -c temp.delta

#FIN111
nucmer --maxmatch --nosimplify /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta  -p temp
mummerplot -color --layout temp.delta
mummerplot -l -c temp.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic.fna -p temp
mummerplot -color --layout temp.delta
mummerplot -l -c temp.delta
```

#### RAxML

A phylogeny was constructed based upon bacterial BUSCO genes, to assess whether the scaffold could represent the larger genome of another Candidatus species related to C.liberibacter solanacearum.
```bash
#BUSCO
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/*.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/*.fna); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/$(basename $Genome | sed 's@.fna@@g' | sed 's@.fasta@@g')/BUSCO
    mkdir -p $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/bacteria_odb10
    OutFile=$(basename $Genome | sed 's@.fna@@g' | sed 's@.fasta@@g')_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done

#Extract complete busco IDs, keep those present in at least 3 genomes:
for file in $(ls liberibacter/*/BUSCO/bacteria_odb10/run_bacteria_odb10/full_table.tsv); do
grep -v "^#" $file | awk '$2=="Complete" {print $1}' >> liberibacter/complete_busco_ids.txt;
done

sort liberibacter/complete_busco_ids.txt |uniq -c > liberibacter/complete_busco_ids_with_counts.txt
grep -v " 2 " liberibacter/complete_busco_ids_with_counts.txt | grep -v " 1 " > liberibacter/complete_busco_ids.txt
awk '{print $2}' liberibacter/complete_busco_ids.txt > liberibacter/complete_busco_ids_3.txt

for file in $(ls liberibacter/*/BUSCO/bacteria_odb10/run_bacteria_odb10/full_table.tsv); do
ID=$(echo $file|cut -f4,5 -d '/'|sed 's@.c/BUSCO@@g')
echo $ID
rm -r gene_pred/*/*/${ID}/BUSCO/*/*/run_*/busco_sequences/single_copy_busco_sequences
done

mkdir liberibacter/busco_nt

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/*/BUSCO/bacteria_odb10/run_bacteria_odb10/busco_sequences/single_copy_busco_sequences.tar.gz); do
cd $(dirname $file)
tar -xzvf single_copy_busco_sequences.tar.gz
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae
done

#Give unique names to the complete busco genes from each assembly:
for dir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/*/BUSCO/bacteria_odb10/run_bacteria_odb10/busco_sequences/single_copy_busco_sequences); do
  sppname=$(echo $dir |cut -f 9 -d "/" | sed 's@/@_@g');
  abbrv=$(echo $sppname | sed 's@_@@g')
  echo $sppname
  echo $abbrv
  for file in ${dir}/*.fna; do
    out=$(echo $file |rev |cut -f 1 -d "/"|rev)
    cp $file liberibacter/busco_nt/${sppname}_${out}
    sed -i 's/^>/>'${abbrv}'|/g' liberibacter/busco_nt/${sppname}_${out}
  cut -f 1 -d ":" liberibacter/busco_nt/${sppname}_${out} | tr '[:lower:]' '[:upper:]' > liberibacter/busco_nt/${sppname}_${out}.1 && mv liberibacter/busco_nt/${sppname}_${out}.1 liberibacter/busco_nt/${sppname}_${out}  
  done
done

#Combine genes from each assembly into a single file per gene:
cd liberibacter/busco_nt
buscos=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/complete_busco_ids.txt
lines=$(cat $buscos)
for line in $lines; do
  for fna in $(ls *_$line.fna); do
  echo $fna
  output=$(echo $line)_nt.fasta
  echo $output
  cat $fna >> $output
  done
done
rm *.fna

#Align the gene sequences;
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/busco_nt
cd $AlignDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/sub_mafft_alignment.sh $AlignDir
#58136735

for gene in $(ls liberibacter/busco_nt/*_aligned.fasta); do
  ID=$(basename $gene |sed 's@_nt_aligned.fasta@@g')
  echo $ID
  mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/busco_nt/$ID
  cp $gene /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/busco_nt/$ID
done

#Trim the alignments:
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/busco_nt/*/*_aligned.fasta); do
  OutDir=$(dirname $Alignment)
  TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
  echo $Alignment
  echo $OutDir
  echo $TrimmedName
  singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/trimal1.4.1.sif trimal -in $Alignment -out $OutDir/$TrimmedName -keepheader -automated1
done

#Trim header names as RAxML need <60 characters in length:
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/busco_nt/*/*_aligned_trimmed.fasta); do
New=$(dirname $Alignment)/$(basename $Alignment .fasta)_edit.fasta
cat $Alignment  | cut -f1 -d '|' | sed 's@GENOMIC@@g'| sed 's@.1@@g'  > $New
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/busco_nt/*/*_edit.fasta); do
while IFS= read -r line; do
    if [[ "$line" =~ ^\>.+ ]]; then
        echo "$line" | wc -c
    fi
done < $file
done

#Run RAxML
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/busco_nt/*/*_edit.fasta); do
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/RAxML/$Prefix
ProgDir=~/git_repos/Wrappers/NBI
mkdir -p $OutDir
sbatch $ProgDir/run_RAxML_msa.sh $Alignment $Prefix $OutDir 
done

#Combine individual gene trees into a consensus tree:
source package 0351788b-6639-43fb-8e80-f28600f83cb1
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/RAxML/*.raxml.bestTree > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/tree-files.txt
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/RAxML/*.raxml.bootstraps > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/bs-files.txt
astral5 -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/tree-files.txt -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/liberibacter_phylogeny.astral.tre
astral5 -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/tree-files.txt -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/bs-files.txt -r 50 -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/liberibacter_phylogeny.bootstrapped.astral.tre
tail -n 1 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/liberibacter_phylogeny.bootstrapped.astral.tre > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/liberibacter_phylogeny.bootstrapped.consensus.astral.tre
astral5 -q /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/liberibacter_phylogeny.bootstrapped.consensus.astral.tre -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/liberibacter_phylogeny.bootstrapped.astral.tre -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/liberibacter/liberibacter_phylogeny.bootstrapped.scored.astral.tre 
```
The phylogeny confirms that scaffold 18 is closest the reference C.liberibacter solanacearum assembly.

BUSCO completeness for scaffold 18 is improved vs the reference C.liberibacter solanacearum assembly
NC_014774.1:
        C:92.7%[S:92.7%,D:0.0%],F:4.0%,M:3.3%,n:124
Scaffold_18:
        C:95.2%[S:95.2%,D:0.0%],F:1.6%,M:3.2%,n:124

### Scaffold 18 only regions

The nature of scaffold 18 regions not aligning to the reference liberibacter genome were investigated.

These unaligned regions were extracted:
```python
fasta_file = "temp_scaffold_18.fasta"
bed_file = "temp_scaffold_18.bed"

with open(fasta_file, "r") as fasta, open(bed_file, "w") as bed:
    chrom = ""
    start = 1
    for line in fasta:
        if line.startswith(">"):
            if chrom != "":
                bed.write(f"{chrom}\t{start}\t{end}\n")
            chrom = line.strip()[1:]
            start = 1
            end = 0
        else:
            end += len(line.strip())
    bed.write(f"{chrom}\t{start}\t{end}\n")
```
```bash
awk '{print "scaffold_18", $0}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18.bed | sed 's/ \+/\t/g' > scaff18.bed
bedtools subtract -a temp_scaffold_18.bed -b scaff18.bed > temp_scaffold_18_non_liberibacter.bed
bedtools getfasta -fi temp_scaffold_18.fasta -bed temp_scaffold_18_non_liberibacter.bed -fo temp_scaffold_18_non_liberibacter.fasta #774,902

awk '/^>/ {if (seqlen) print seqlen; printf "%s\t", substr($0,2); seqlen=0; next} {seqlen += length($0)} END {print seqlen}' temp_scaffold_18_non_liberibacter.fasta

```
There are 151 unaligned regions, representing 774,902 bp of the scaffold

#### Blast

The unaligned regions were classified by blast search:
```bash
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/nt_premade_02102023/nt
sbatch $ProgDir/run_blastn.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18_non_liberibacter.fasta $Database /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae temp_18 9999999
#58134218
```
```R
# Install and load necessary packages
install.packages("rentrez")
install.packages("XML")
library(rentrez)
library(XML)

accessions <- readLines("temp_out.txt")

# Function to retrieve species for a single genome accession
get_species_for_genome_accession <- function(accession) {
  summary <- assembly_summary(accession)
  species <- summary$organism_name
  return(species)
}

# Function to get species for a list of genome accessions
get_species_for_genome_accessions <- function(accessions) {
  species_list <- list()
  for (accession in accessions) {
    species <- get_species_for_genome_accession(accession)
    species_list[[accession]] <- species
  }
  return(species_list)
}

species_info <- get_species_for_genome_accessions(accessions)

for (accession in accessions) {
  print(paste("Genome Accession:", accession, "| Species:", species_info[[accession]]))
}

























# Install and load necessary packages
install.packages("rentrez")
library(rentrez)

# Function to retrieve species for a single accession (gene or genome)
get_species_for_accession <- function(accession) {
  if (grepl("^NM_", accession)) {  # Assuming gene accessions start with "NM_"
    # If it's a gene accession, query NCBI Gene database
    xml_content <- entrez_fetch(db="gene", id=accession, rettype="xml")
    doc <- xmlParse(xml_content)
    species <- xmlValue(doc[["//GBSeq_organism"]])
  } else {
    # If it's a genome accession, query NCBI Assembly database
    summary <- assembly_summary(accession)
    species <- summary$organism_name
  }
  return(species)
}

# Function to get species for a list of accessions (genes or genomes)
get_species_for_accessions <- function(accessions) {
  species_list <- list()
  for (accession in accessions) {
    species <- get_species_for_accession(accession)
    species_list[[accession]] <- species
  }
  return(species_list)
}

# Example list of gene and genome accessions
accessions <- c("NM_000014", "NM_000015", "NM_000016", "GCF_000001405.39") # Example gene and genome accessions

# Retrieve species for each accession
species_info <- get_species_for_accessions(accessions)

# Print the results
for (accession in accessions) {
  print(paste("Accession:", accession, "| Species:", species_info[[accession]]))
}

# Read the file
file_content <- readLines("nuccore_result.txt")

# Initialize an empty vector to store the extracted information
extracted_info <- c()

# Iterate through each line of the file
for (line in file_content) {
  # Check if the line starts with a number followed by a period
  if (grepl("^\\d+\\.", line)) {
    # Extract the first two words from the line
    extracted <- strsplit(line, "\\s+")[[1]][1:3]
    # Join the extracted words into a single string and append to the vector
    extracted_info <- c(extracted_info, paste(extracted, collapse = " "))
  }
}

# Print the extracted information
print(extracted_info)
```
The majority of blast hits are to Candidatus Liberibacter or Klebsiella pneumoniae, there are also some Liberibacter phage hits

None of the scaffold aligns to Klebsiella pneumoniae:
```bash
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
Query=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Klebsiella/pneumoniae/GCA_000240185.2/GCA_000240185.2_ASM24018v2_genomic.fna
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger
OutFile=scaff18-kleb
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58078409
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-kleb.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-kleb.delta
```
#### kraken

The unaligned regions were also classified via kraken:
```bash
#run kraken
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18_non_liberibacter.fasta
OutPrefix=$(basename $Assembly | sed 's@.fasta@@g')_kraken2nt
OutDir=$(dirname $Assembly)/temp_18_kraken2.1.3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/kraken/nt_14092023
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_kraken2.sh $Assembly $Database $OutDir $OutPrefix
#58128745

#regions classified as liberibacter were collected
grep 'Liberibacter' temp_18_kraken2.1.3/temp_scaffold_18_non_liberibacter_kraken2nt_output.txt | awk '{print $2}' > temp_id.txt
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file temp_id.txt --input temp_scaffold_18_non_liberibacter.fasta --output temp_scaffold_18_non_liberibacter+.fasta
grep -v '>' temp_scaffold_18_non_liberibacter+.fasta | wc -c #745,528
```
kraken corroborates the blast results classifying 85 to liberibacter, 29 to Eukaryota, 36 unknown, and 1 to viruses, with the 85 liberibacter regions representing 745,528/774,902 bp of the unaligned region.

#### Nucmer
Scaffold 18 was aligned against itself to investigate repetative regions:
```bash
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta -p scaffold_18_self
show-coords -T scaffold_18_self.delta > scaffold_18_self.coords
mummerplot -l -c scaffold_18_self.delta
mummerplot -color scaffold_18_self.delta
```
Regions of self alignments match up with regions not aligning to the reference C.liberibacter solanacearum assembly, supportive of possible missassembly by previous short read approaches.

This was compared to a self alignment of C. liberibacter asiaticus:
```bash
nucmer --maxmatch --nosimplify /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_asiaticus/GCA_002216815.1/GCA_002216815.1_ASM221681v1_genomic.fna /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_asiaticus/GCA_002216815.1/GCA_002216815.1_ASM221681v1_genomic.fna -p asi_self
show-coords -T asi_self.delta > asi_self.coords
mummerplot -l -c asi_self.delta
mummerplot -color asi_self.delta

nucmer --maxmatch --nosimplify /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_asiaticus/GCA_030585885.1/GCA_030585885.1_ASM3058588v1_genomic.fna /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_asiaticus/GCA_030585885.1/GCA_030585885.1_ASM3058588v1_genomic.fna -p asi_self2
show-coords -T asi_self2.delta > asi_self2.coords
mummerplot -l -c asi_self2.delta
mummerplot -color asi_self2.delta
```

## Reassembly

HiC contact map suggests that scaffolds 184 and 1142 are also of C. lineribacter origin, this is supported by kraken for 184.

```bash
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa | grep -A 1 -w '>scaffold_184' > temp_scaffold_184.fasta
awk 'BEGIN{RS=">"} NR>1 {sub("\n","\t",$0); gsub("\n",""); print ">"$1"\n"$2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/yahs/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_scaffolds_final.fa | grep -A 1 -w '>scaffold_1142' > temp_scaffold_1142.fasta


Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
Query=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_184.fasta
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger
OutFile=scaff18-184
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58736368
source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-184.delta
mummerplot -color --layout /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-184.delta
#ERROR: No alignment data to plot

Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
Query=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_1142.fasta
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger
OutFile=scaff18-1142
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_nucmer.sh $Reference $Query $OutDir $OutFile
#58736369
source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-1142.delta
mummerplot -color --layout /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18-1142.delta
#ERROR: No alignment data to plot

cat temp_scaffold_18.fasta temp_scaffold_184.fasta temp_scaffold_1142.fasta > Liberibacter/assembly_v1.fa
```
```bash
ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#58754389

ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
OutDir=$(dirname $Assembly)/bwa
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')_HiC
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R2.fastq.gz
Read3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R1.fastq.gz
Read4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Read1 $Read2 $Read3 $Read4
#58797861

cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter
samtools view -h minimap2/assembly_v1.fa.bam -o minimap2/assembly_v1.fa.sam
samtools fastq -@32 minimap2/assembly_v1.fa.sam > reads.fastq

ProgDir=~/git_repos/Wrappers/NBI
Reads=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/reads.fastq
OutDir=$(dirname $Reads)/hifiasm_19.5.2
OutFile=liberibacter_api
mkdir -p $OutDir
sbatch $ProgDir/run_hifiasm_default.sh $OutDir $OutFile $Reads
#58778325

#n       n:500   n:N50   min     N80     N50     N20     max     sum
#8261    8261    1530    4109    70325   151972  299874  2804342 801.5e6 liberibacter_api.bp.p_ctg.fa

source package d6092385-3a81-49d9-b044-8ffb85d0c446
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/hifiasm_19.5.2/liberibacter_api.bp.p_ctg.fa -input_type fasta -dbtype nucl -title libre  -parse_seqids -out libre
blastn -query /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic.fna -db libre -out libre_results -evalue 1e-5 -outfmt 6 -num_threads 1
awk '{print $2}' libre_results | sort | uniq > temp_id.txt 

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file temp_id.txt --input  /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/hifiasm_19.5.2/liberibacter_api.bp.p_ctg.fa --output  /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/hifiasm_19.5.2/liberibacter_api_hit.fa

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/hifiasm_19.5.2/liberibacter_api_hit.fa
Enzyme=GATC
OutDir=$(dirname $Assembly)
OutFile=$(basename $Assembly | sed 's@.fa@@')
CRead1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R1.fastq.gz
CRead2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/HiC/apicales_286172-S3HiC_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_omniHiCmap.sh $Assembly $Enzyme $OutDir $OutFile $CRead1 $CRead2
#58779970

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/hifiasm_19.5.2/liberibacter_api_hit.fa
Alignment=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/hifiasm_19.5.2/liberibacter_api_hit_mapped.PT.bam
Alignment_Index=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/hifiasm_19.5.2/liberibacter_api_hit_mapped.PT.bam.bai
Enzyme=GATC
OutDir=$(dirname $Assembly)/yahs
OutFile=$(basename $Assembly | sed 's@.fa@@')
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_yahs.sh $Assembly $Alignment $Alignment_Index $Enzyme $OutDir $OutFile
#58797844
```
#### bridging reads
```bash
cd Liberibacter/minimap2
samtools sort assembly_v1.fa.bam -o assembly_v1.fa.bam.temp && mv assembly_v1.fa.bam.temp assembly_v1.fa.bam
samtools index assembly_v1.fa.bam
cd ../..

#Extract reads overlapping with BED regions
bedtools intersect -abam Liberibacter/minimap2/assembly_v1.fa.bam -b temp_scaffold_18_non_liberibacter.bed > temp_overlapping_reads.bam

#Filter reads spanning entire BED regions
bedtools intersect -a temp_scaffold_18_non_liberibacter.bed -b temp_overlapping_reads.bam -wb -f 1.0 > spanning_reads.txt
awk '{print $0 "\t" $1 "_" $2 "_" $3}' spanning_reads.txt > output_file.txt && mv output_file.txt spanning_reads.txt
awk '{ pattern_count_all[$10]++ } $8 > 0 { pattern_count_gt_zero[$10]++ } 
     END { 
        for (pattern in pattern_count_all) 
            print pattern, pattern_count_all[pattern], (pattern_count_gt_zero[pattern] ? pattern_count_gt_zero[pattern] : 0) 
     }' spanning_reads.txt > Liberibacter/minimap2/spanning_read_depth.txt

awk '{print $1 "_" $2 "_" $3}' temp_scaffold_18_non_liberibacter.bed > temp_scaffold_18_non_liberibacter_ID.txt
awk 'NR==FNR{patterns[$10]; next} !($1 in patterns) {print $0 " 0" " 0"}' spanning_reads.txt temp_scaffold_18_non_liberibacter_ID.txt >> Liberibacter/minimap2/spanning_read_depth.txt

grep -v ">" temp_scaffold_18.fasta | tr -d '\n' | wc -c 
```
```python
import pandas as pd
import matplotlib.pyplot as plt

# Read data from tab-separated text file
df = pd.read_csv('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/minimap2/spanning_read_depth.txt', sep='\t', header=None, names=['scaffold', 'no', 'Start', 'End', 'Y', 'Y2'])

bed_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/scaff18.bed"
df_bed = pd.read_csv(bed_file, sep='\t', header=None, names=['start', 'end'])

# Define the overall range of x-values
overall_start = 1
overall_end = 1972533

# Plotting
plt.figure(figsize=(10, 6))

# Plot the ranges covered by data
for _, row in df.iterrows():
  plt.plot([row['Start'], row['End']], [row['Y'], row['Y']], color='blue')

# plot the ranges not covered by data in red at 0
for _, row in df_bed.iterrows():
  plt.plot([row['start'], row['end']], [0, 0], color='red')

plt.xlabel('Scaffold 18')
plt.ylabel('No. of reads completely spanning regions')
plt.title('Bridged regions with no alignment to Genbank C.liberibacter solanacearum')
plt.grid(True)

# Save figure to a file
plt.savefig('ranges_plot.png')

```
### other assemblies
```bash
C:\Users\did23faz\Downloads\sratoolkit.current-win64\sratoolkit.3.0.5-win64\bin>prefetch --max-size 1t --option-file C:\Users\did23faz\Documents\SRR3312978.txt --output-directory \\jic-hpc-data\Group-Scratch\Saskia-Hogenhout\tom_heaven

pwd
#/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/raw_data
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/sra-tools_3.0.9.sif fastq-dump --split-files --gzip -O SRR3312978 SRR3312978/*.sra

ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
OutDir=$(dirname $Assembly)/minimap2-2
Outfile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/raw_data/SRR3312979.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1
#58779050

ProgDir=~/git_repos/Wrappers/NBI
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
OutDir=$(dirname $Assembly)/bwa-2
Outfile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/raw_data/SRR3312978/SRR3312978_1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/raw_data/SRR3312978/SRR3312978_2.fastq.gz
mkdir $OutDir
sbatch $ProgDir/bwa-mem.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#58778991
```

#### deeplasmid
```bash
Contigs=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
OutDir=$(dirname $Contigs)/deeplasmid
OutFile=Liberibacter_v1
ProgDir=~/git_repos/Wrappers/NBI
mkdir ${OutDir}
sbatch $ProgDir/run_deeplasmid.sh $Contigs $OutDir $OutFile
#58753152
```
#### circlator
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
Read_type=illumina
OutDir=$(dirname $Assembly)/circlator
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Read1=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R1.fastq.gz
Read2=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T505_R2.fastq.gz
Read3=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R1.fastq.gz
Read4=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_apicales/TellSeq/apicales_T507_R2.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_circlator.sh $Assembly $Read_type $OutDir $Outfile $Read1 $Read2 $Read3 $Read4
#58797263

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
Read_type=pacbio-corrected
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/circlator-hifi
Outfile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_circlator.sh $Assembly $Read_type $OutDir $Outfile $Read1 $Read2
#58816003, 58817490, 58817981, 58818333, 58819283, 58819286, 58819287


```
### Phage
```bash
source package /nbi/software/testing/bin/phast-1.5

```
```bash
Liberibacter phage FP2, complete genome
38,552 bp linear DNA 
JF773396.1 GI:342674737

Liberibacter phage SC1, complete genome
40,048 bp circular DNA 
NC_019549.1 GI:423262466

Liberibacter phage SC2, complete genome
38,997 bp circular DNA 
NC_019550.1 GI:423262507

Liberibacter phage SGCA5-1 clone contig1 genomic sequence
1,465 bp linear DNA 
KX879600.1 GI:1102356448

Liberibacter phage SGCA5-1 clone contig2 genomic sequence
36,022 bp linear DNA 
KX879601.1 GI:1102356451

Liberibacter phage HHCA1-2, complete genome
38,989 bp circular DNA 
KX879602.1 GI:1103774642

Liberibacter phage P-JXGC-3, complete genome
31,449 bp circular DNA 
KY661963.1 GI:1168017351

Liberibacter phage GZQL4, complete genome
39,511 bp circular DNA 
CP124119.1 GI:2501026657

Liberibacter phage P-Myan16-2, complete sequence
36,303 bp circular DNA 
CP060690.1 GI:2076496336

Liberibacter phage P-PA19-1, complete genome
37,601 bp linear DNA 
MT899444.1 GI:1914795703

Liberibacter phage P-PA19-2, complete genome
31,505 bp linear DNA 
MT899443.1 GI:1914795665

MAG: Bacteriophage sp. isolate 2786_71427, partial genome
45,114 bp linear DNA 
OP072401.1 GI:2293579751

MAG: Bacteriophage sp. isolate 3238_77825, partial genome
45,701 bp linear DNA 
OP072565.1 GI:2293592270

MAG: Bacteriophage sp. isolate 3518_1005, partial genome
46,000 bp linear DNA 
OP072666.1 GI:2293605677

Liberibacter phage SGCA5-1 clone contig2 genomic sequence
36,022 bp linear DNA 
KX879601.1 GI:1102356451

cat Liberibacter/phages/*.fna > temp_phages.fasta
nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_phages.fasta -p phage
show-coords -T phage.delta > phage.coords
mummerplot -l -c phage.delta
mummerplot -color phage.delta
```
```bash
20. Candidatus Liberibacter asiaticus genes for hypothetical protein, phage related protein, bordetella phage Bbp38 like protein, prophage antirepressor like protein, complete cds, isolate: Iw4 (Miyako-13)
3,174 bp linear DNA 
AB623416.1 GI:352952608

21. Candidatus Liberibacter asiaticus genes for hypothetical protein, phage related protein, bordetella phage Bbp38 like protein, prophage antirepressor like protein, complete cds, isolate: V2
3,174 bp linear DNA 
AB623417.1 GI:352952613

22. Candidatus Liberibacter asiaticus genes for hypothetical protein, phage related protein, bordetella phage Bbp38 like protein, prophage antirepressor like protein, complete cds, isolate: Tw3
3,174 bp linear DNA 
AB623418.1 GI:352952618

23. Candidatus Liberibacter asiaticus genes for hypothetical protein, phage related protein, bordetella phage Bbp38 like protein, prophage antirepressor like protein, complete cds, isolate: VN50
3,015 bp linear DNA 
AB623419.1 GI:352952623

24. Candidatus Liberibacter asiaticus genes for hypothetical protein, phage related protein, bordetella phage Bbp38 like protein, prophage antirepressor like protein, complete cds, isolate: Tw1
3,175 bp linear DNA 
AB623420.1 GI:352952628

25. Candidatus Liberibacter asiaticus genes for hypothetical protein, phage related protein, bordetella phage Bbp38 like protein, prophage antirepressor like protein, complete cds, isolate: ThaiSu-1
3,172 bp linear DNA 
AB623421.1 GI:352952633

26. Candidatus Liberibacter asiaticus genes for hypothetical protein, phage related protein, bordetella phage Bbp38 like protein, prophage antirepressor like protein, complete cds, isolate: Thai-1
3,015 bp linear DNA 
AB623422.1 GI:352952638

Candidatus Liberibacter asiaticus strain FL-periwinkle isolate Type-D prophage region
11,121 bp linear DNA 
JX275492.1 GI:479012040

Candidatus Liberibacter asiaticus isolate YN-JS-11 hypothetical protein gene, partial cds; hypothetical protein gene, complete cds; and P4 family phage/plasmid primase gene, partial cds
549 bp linear DNA 
KJ777538.1 GI:662706499

67. Liberibacter phage SGCA5-1 clone contig1 genomic sequence
1,465 bp linear DNA 
KX879600.1 GI:1102356448

68. Liberibacter phage SGCA5-1 clone contig2 genomic sequence
36,022 bp linear DNA 
KX879601.1 GI:1102356451

Candidatus Liberibacter asiaticus strain TH-lime type-D prophage region genomic sequence; and hypothetical protein genes, complete cds
2,098 bp linear DNA 
MG603592.1 GI:1390637436

91. Candidatus Liberibacter asiaticus clone XC1 prophage region genomic sequence
2,651 bp linear DNA 
MN331742.1 GI:1826554177

92. Candidatus Liberibacter asiaticus clone GJ1.2 prophage region genomic sequence
2,651 bp linear DNA 
MN331743.1 GI:1826554181

93. Candidatus Liberibacter asiaticus clone MY14 prophage region genomic sequence
2,651 bp linear DNA 
MN331744.1 GI:1826554185

94. Candidatus Liberibacter asiaticus clone FC5.3 prophage region genomic sequence
2,651 bp linear DNA 
MN331745.1 GI:1826554189

95. Candidatus Liberibacter asiaticus clone FC10.1 prophage region genomic sequence
2,651 bp linear DNA 
MN331746.1 GI:1826554193

96. Candidatus Liberibacter asiaticus clone FC19.2 prophage region genomic sequence
2,650 bp linear DNA 
MN331747.1 GI:1826554197

97. Candidatus Liberibacter asiaticus clone LC10.3 prophage region genomic sequence
2,651 bp linear DNA 
MN331748.1 GI:1826554201

98. Candidatus Liberibacter asiaticus clone LC13.2 prophage region genomic sequence
2,651 bp linear DNA 
MN331749.1 GI:1826554205

99. Candidatus Liberibacter asiaticus clone LC13.1 prophage region genomic sequence
2,651 bp linear DNA 
MN331750.1 GI:1826554209

100. Candidatus Liberibacter asiaticus clone HC1.1 prophage region genomic sequence
2,651 bp linear DNA 
MN331751.1 GI:1826554213

101. Candidatus Liberibacter asiaticus clone HC1.2 prophage region genomic sequence
2,651 bp linear DNA 
MN331752.1 GI:1826554217

102. Candidatus Liberibacter asiaticus clone LC7.1 prophage region genomic sequence
2,651 bp linear DNA 
MN331753.1 GI:1826554220

103. Candidatus Liberibacter asiaticus clone QC2.1 prophage region genomic sequence
2,640 bp linear DNA 
MN331754.1 GI:1826554224

104. Candidatus Liberibacter asiaticus clone LC5.1 prophage region genomic sequence
2,640 bp linear DNA 
MN331755.1 GI:1826554228

105. Candidatus Liberibacter asiaticus clone LC5.2 prophage region genomic sequence
2,640 bp linear DNA 
MN331756.1 GI:1826554232

106. Candidatus Liberibacter asiaticus clone QC2.2 prophage region genomic sequence
2,640 bp linear DNA 
MN331757.1 GI:1826554236

107. Candidatus Liberibacter asiaticus clone MY18.1 prophage region genomic sequence
2,640 bp linear DNA 
MN331758.1 GI:1826554240

108. Candidatus Liberibacter asiaticus clone MY18.2 prophage region genomic sequence
2,639 bp linear DNA 
MN331759.1 GI:1826554244

109. Candidatus Liberibacter asiaticus clone GJ1.1 prophage region genomic sequence
2,946 bp linear DNA 
MN331760.1 GI:1826554248

110. Candidatus Liberibacter asiaticus clone GJ1.3 prophage region genomic sequence
2,946 bp linear DNA 
MN331761.1 GI:1826554252

111. Candidatus Liberibacter asiaticus clone GJ1.4 prophage region genomic sequence
2,946 bp linear DNA 
MN331762.1 GI:1826554256

112. Candidatus Liberibacter asiaticus clone NM4.1 prophage region genomic sequence
2,946 bp linear DNA 
MN331763.1 GI:1826554260

113. Candidatus Liberibacter asiaticus clone FC4.1 prophage region genomic sequence
2,946 bp linear DNA 
MN331764.1 GI:1826554263

114. Candidatus Liberibacter asiaticus clone GJ1.5 prophage region genomic sequence
2,946 bp linear DNA 
MN331765.1 GI:1826554267

Candidatus Liberibacter asiaticus genes for hypothetical protein, phage related protein, bordetella phage Bbp38 like protein, prophage antirepressor like protein, complete cds, isolate: Tw3
3,174 bp linear DNA 
AB623418.1 GI:352952618

Candidatus Liberibacter asiaticus strain FL-periwinkle isolate Type-D prophage region
11,121 bp linear DNA 
JX275492.1 GI:479012040

Candidatus Liberibacter asiaticus strain TH-lime type-D prophage region genomic sequence; and hypothetical protein genes, complete cds
2,098 bp linear DNA 
MG603592.1 GI:1390637436


```
#### Repeatmodeler
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/repeatmodeler
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmodeler.sh $Genome $OutFile $OutDir
#58762545
```
#### Repeatmasker
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
Species=bacteria
Repeat_library=$(ls $(dirname $Genome)/repeatmodeler/*/consensi.fa.classified)
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/repeatmasker
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmasker4.1.5.sh $Genome $OutFile $OutDir $Species $Repeat_library
#58763266
```
#### HiSat2
Hisat2 is the aligner used internally by braker, it is slower than star but requires less memory.
```bash
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

ERROR: the RNA comes from a different sample! Also wouldnt expect to find bacterial RNA anyway due to the Poly A sequencing method - Sam

#### DNA A box
```bash
source package 37f0ffda-9f66-4391-87e2-38ccd398861d
mkdir 18
cd 18
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/dnaa.fasta -input_type fasta -dbtype nucl -title 18 -out 18

blastn -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/18/18 -out dnaa18.vs.nt.mts1.hsp1.1e25.megablast.out -evalue 1e-25 -max_hsps 1 -max_target_seqs 99999 -num_threads 1 -task megablast -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta
```
Hit to NZ_CP019958.1:c638673-637165 Candidatus Liberibacter asiaticus, between 1,186,822-1,188,322.
```python
from Bio import SeqIO

# Function to find A box motif positions
def find_a_box_positions(genome_sequence):
    a_box_motif = "ACGTGAACGTCGT"
    positions = []
    for i in range(len(genome_sequence) - len(a_box_motif) + 1):
        if genome_sequence[i:i+len(a_box_motif)] == a_box_motif:
            positions.append(i)
    return positions

# Load the bacterial genome sequence from a FASTA file
genome_file = "bacterial_genome.fasta"
genome_record = SeqIO.read(genome_file, "fasta")
genome_sequence = str(genome_record.seq)

# Find positions of the A box motif
a_box_positions = find_a_box_positions(genome_sequence)

# Print the positions
if a_box_positions:
    print("A box motif found at positions:", a_box_positions)
else:
    print("A box motif not found in the genome sequence.")

```
NOTE: ori-finder web tools identified at different location

### Annotation
Assembly given the name JIC_LsoCDa_1.0 and locus tag WAE90 in genbank submission (JIC Candidatus Liberibacter solanacearum Dyspera apicales).
```bash
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
OutDir=$(dirname $Genome)/prokka
OutFile=JIC_LsoCDa_1.0
Locustag=WAE90
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_prokka.sh $Genome $OutDir $OutFile $Locustag
#58864242

grep '>' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.faa | wc -l #2,024
```
#### Predector
NOTE: predector is built for fungal effector prediction, however SignalP versions and tmhmm are run as part of the pipeline.

move to gruffalo
```bash
conda activate predector
for proteome in $(ls /home/theaven/scratch/uncompressed/JIC_LsoCDa_1.0.faa); do
ProgDir=/home/theaven/scratch/apps/predector
OutDir=/home/theaven/scratch/uncompressed/JIC_LsoCDa_1.0
sbatch $ProgDir/predector_singularity.sh $proteome 1.2.6
sleep 30
done
conda deactivate
#17909631
```

#### DeepTMHMM

There are five categories for each residue:
Membrane in-out
Membrane out-in
Signal peptide
Inside
Outside

Furthermore, the outputs a prediction of the protein type:
Signal peptide + globular (SP+Glob)
Signal peptide + transmembrane (SP+TM)
Transmembrane (TM)
Globular (Glob)
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/deeptmhmm
sbatch ~/git_repos/Wrappers/NBI/run_deeptmhmm.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/deeptmhmm
#58891136

#Version of the command line install is unknown, run with web client also:
awk '/^>/{n++}{print > ("subfile_"int((n-1)/500)+1".fasta")}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.faa
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/deeptmhmm/deeptmhmm*
```
#### TMHMM
```bash
source package /nbi/software/testing/bin/tmhmm-2.0c
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/tmhmm
tmhmm /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.faa > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/tmhmm/tmhmm.out
```
#### Signalp
Web clients used:

Signalp6
```bash
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp6
```
Signalp5
```bash
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp5
```
Signalp4 - split into 5 jobs
```bash
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp4.1
```
Signalp3
```bash
grep '>\|Signal peptide probability:' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp3/output.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp3/output2.txt 
```
#### Swissprot
```bash
Proteome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.faa 
OutDir=$(dirname $Proteome)/swissprot
SwissDbDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/Uniprot/swissprot_2024_March_10
SwissDbName=uniprot_sprot
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName 
#58956713
```
#### Interproscan
```bash
InFile=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.faa
OutDir=$(dirname $InFile)/interproscan
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_interproscan.sh $InFile $OutDir
#58949694
```
#### Phaster
```bash
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phaster
```
Prophage regions identified:
1 scaffold_18:3645-96338
2 scaffold_18:353952-410656
3 scaffold_18:447487-489419
4 scaffold_18:527288-565969
5 scaffold_18:564614-619688
6 scaffold_18:656708-720434
7 scaffold_18:806295-851289
8 scaffold_18:869573-951741
9 scaffold_18:937224-963457
10 scaffold_18:1053039-1093542
11 scaffold_18:1239344-1266002
12 scaffold_18:1288966-1319545
13 scaffold_18:1315885-1360196
14 scaffold_18:1357494-1384070
15 scaffold_18:1454995-1556562
16 scaffold_18:1634561-1694671
17 scaffold_18:1766081-1867628
18 scaffold_18:1867834-1883290 

##### Parse
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
```python
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.gff", "r") as infile, open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.tsv", "w") as outfile:
  # Write header to the output file
  outfile.write("Locus_Tag\tStart_Position\tStop_Position\tStrand\n")
  # Iterate over each line in the input file
  for line in infile:
    # Skip comment lines
    if line.startswith("##"):
      continue
    # Split the line into fields based on whitespace
    fields = line.strip().split("\t")
    # Ensure that the fields list has enough elements
    if len(fields) >= 3 and fields[2] == "gene":
      # Extract the relevant information
      locus_tag = "_".join(fields[8].split(";")[0].split("=")[1].split("_")[:2])
      start_position = fields[3]
      stop_position = fields[4]
      strand = fields[6]
      # Write the extracted information to the output file
      print(locus_tag)
      outfile.write(f"{locus_tag}\t{start_position}\t{stop_position}\t{strand}\n")
    else:
      print("Warning: Skipped line due to unexpected format:", line)


# Define file paths
input_file1 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.tsv"
input_file2 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp6/output.gff3"
output_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6.tsv"

# Read data from input files
data1 = {}
with open(input_file1, 'r') as file:
  header = next(file).strip().split('\t')  # Read and split header
  for line in file:
    parts = line.strip().split('\t')
    data1[parts[0]] = parts

data2 = {}
with open(input_file2, 'r') as file:
  for line in file:
    if not line.startswith('#'):
      parts = line.strip().split('\t')
      key = parts[0].split()[0]  # Extracting the first space-separated element from column 1 as the key
      data2[key] = parts[2]

# Open output file for writing
with open(output_file, 'w') as file:
  # Write header line with additional column
  file.write("\t".join(header) + "\tSignalP6\n")
  # Write data with SignalP6 values appended or '###' if no match found
  for key, values in data1.items():
    if key in data2:
      file.write("\t".join(values) + "\t" + data2[key] + "\n")
    else:
      file.write("\t".join(values) + "\t###\n")

input_file1 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6.tsv"
input_file2 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp5/output.gff3"
output_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5.tsv"

# Read data from input files
data1 = {}
with open(input_file1, 'r') as file:
  header = next(file).strip().split('\t')  # Read and split header
  for line in file:
    parts = line.strip().split('\t')
    data1[parts[0]] = parts

data2 = {}
with open(input_file2, 'r') as file:
  for line in file:
    if not line.startswith('#'):
      parts = line.strip().split('\t')
      key = parts[0].split()[0]  # Extracting the first space-separated element from column 1 as the key
      data2[key] = parts[2]

# Open output file for writing
with open(output_file, 'w') as file:
  # Write header line with additional column
  file.write("\t".join(header) + "\tSignalP5\n")
  # Write data with SignalP6 values appended or '###' if no match found
  for key, values in data1.items():
    if key in data2:
      file.write("\t".join(values) + "\t" + data2[key] + "\n")
    else:
      file.write("\t".join(values) + "\t###\n")

input_file1 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5.tsv"
input_file2 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp4.1/output.gff"
output_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4.tsv"

# Read data from input files
data1 = {}
with open(input_file1, 'r') as file:
  header = next(file).strip().split('\t')  # Read and split header
  for line in file:
    parts = line.strip().split('\t')
    data1[parts[0]] = parts

data2 = {}
with open(input_file2, 'r') as file:
  for line in file:
    if not line.startswith('#'):
      parts = line.strip().split('\t')
      key = parts[0].split()[0]  # Extracting the first space-separated element from column 1 as the key
      data2[key] = parts[2]

# Open output file for writing
with open(output_file, 'w') as file:
  # Write header line with additional column
  file.write("\t".join(header) + "\tSignalP4.1\n")
  # Write data with SignalP6 values appended or '###' if no match found
  for key, values in data1.items():
    if key in data2:
      file.write("\t".join(values) + "\t" + data2[key] + "\n")
    else:
      file.write("\t".join(values) + "\t###\n")

input_file1 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4.tsv"
input_file2 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/tmhmm/tmhmm.out"
output_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4_tmhmm.tsv"

# Read data from input files
data1 = {}
with open(input_file1, 'r') as file:
  header = next(file).strip().split('\t')  # Read and split header
  for line in file:
    parts = line.strip().split('\t')
    data1[parts[0]] = parts

data2 = {}
with open(input_file2, 'r') as file:
  for line in file:
    if line.startswith('#'):
      if 'Number of predicted TMHs:' in line:
        parts = line.strip().split()
        key = parts[1]  # Extracting the first space-separated element from column 1 as the key
        value = parts[-1]  # Extracting the last element as the value
        data2[key] = value
      else:
        continue

# Open output file for writing
with open(output_file, 'w') as file:
  # Write header line with additional column
  file.write("\t".join(header) + "\ttmhmm\n")
  # Write data with SignalP6 values appended or '###' if no match found
  for key, values in data1.items():
    if key in data2:
      file.write("\t".join(values) + "\t" + data2[key] + "\n")
    else:
      file.write("\t".join(values) + "\t###\n")


input_file1 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4_tmhmm.tsv"
input_file2 = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/deeptmhmm/TMRs.gff3"
output_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4_tmhmm_deeptmhmm.tsv"

# Read data from input files
data1 = {}
with open(input_file1, 'r') as file:
  header = next(file).strip().split('\t')  # Read and split header
  for line in file:
    parts = line.strip().split('\t')
    data1[parts[0]] = parts

data2 = {}
with open(input_file2, 'r') as file:
  for line in file:
    if line.startswith('#'):
      if 'Number of predicted TMRs:' in line:
        parts = line.strip().split()
        key = parts[1]  # Extracting the first space-separated element from column 1 as the key
        value = parts[-1]  # Extracting the last element as the value
        data2[key] = value
      else:
        continue

# Open output file for writing
with open(output_file, 'w') as file:
  # Write header line with additional column
  file.write("\t".join(header) + "\tdeeptmhmm\n")
  # Write data with SignalP6 values appended or '###' if no match found
  for key, values in data1.items():
    if key in data2:
      file.write("\t".join(values) + "\t" + data2[key] + "\n")
    else:
      file.write("\t".join(values) + "\t###\n")



# Define the regions
regions = {
    1: "scaffold_18:3645-96338",
    2: "scaffold_18:353952-410656",
    3: "scaffold_18:447487-489419",
    4: "scaffold_18:527288-565969",
    5: "scaffold_18:564614-619688",
    6: "scaffold_18:656708-720434",
    7: "scaffold_18:806295-851289",
    8: "scaffold_18:869573-951741",
    9: "scaffold_18:937224-963457",
    10: "scaffold_18:1053039-1093542",
    11: "scaffold_18:1239344-1266002",
    12: "scaffold_18:1288966-1319545",
    13: "scaffold_18:1315885-1360196",
    14: "scaffold_18:1357494-1384070",
    15: "scaffold_18:1454995-1556562",
    16: "scaffold_18:1634561-1694671",
    17: "scaffold_18:1766081-1867628",
    18: "scaffold_18:1867834-1883290"
}

# Open the input file
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4_tmhmm_deeptmhmm.tsv", "r") as f:
  lines = f.readlines()

# Add region information as an extra column
output_lines = []
for i, line in enumerate(lines):
  if i == 0:  # Header line
    output_lines.append(line.strip() + "\tRegion\n")
  else:
    fields = line.strip().split("\t")
    start = int(fields[1])
    end = int(fields[2])
    for region, region_info in regions.items():
      region_start, region_end = map(int, region_info.split(":")[1].split("-"))
      if (start >= region_start and start <= region_end) or (end >= region_start and end <= region_end):
        fields.append(str(region))
        break
    output_lines.append("\t".join(fields) + "\n")

# Write the output to a new file
with open("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4_tmhmm_deeptmhmm_phaster.tsv", "w") as f:
    f.writelines(output_lines)

import csv

# Define the input and output file paths
input_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4_tmhmm_deeptmhmm_phaster.tsv'
output_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4_tmhmm_deeptmhmm_phaster_.tsv'

# Open input and output files
with open(input_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
  reader = csv.reader(infile, delimiter='\t')
  writer = csv.writer(outfile, delimiter='\t')
  # Iterate through each row in the input file
  for row in reader:
    # Ensure the row has at least 10 columns
    while len(row) < 10:
      # Add empty columns to reach 10 columns
      row.append('')
    # Check if column 10 is empty
    if not row[9]:  # Python is zero-indexed, so column 10 corresponds to index 9
      # Add "###" to column 10
      row[9] = '###'
    # Write the modified row to the output file
    writer.writerow(row)


# Initialize an empty dictionary to store gene names and their probabilities
gene_probabilities = {}

# Open the file for reading
with open(input_file, 'r') as file:
  # Initialize variables to store gene name and signal peptide probability
  gene_name = ""
  signal_peptide_probability = ""
  # Iterate through each line in the file
  for line in file:
    # Check if the line starts with ">WAE90_"
    if line.startswith('>WAE90_'):
      # Extract the gene name and remove the '>'
      gene_name = line.strip().split()[0][1:]  # Extracting gene name from the line and removing '>'
      # Initialize the gene name in the dictionary if it's not already present
      if gene_name not in gene_probabilities:
        gene_probabilities[gene_name] = ""
    # Check if the line contains "Signal peptide probability"
    elif "Signal peptide probability" in line:
      # Extract the signal peptide probability
      signal_peptide_probability = line.strip().split(':')[-1].strip()  # Extracting probability value
      # Associate gene name with signal peptide probability in the dictionary
      gene_probabilities[gene_name] = signal_peptide_probability

# Print the dictionary containing gene names and their associated probabilities
for gene, probability in gene_probabilities.items():
  print("Gene:", gene)
  print("Signal peptide probability:", probability)

import pandas as pd
df = pd.read_csv('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4_tmhmm_deeptmhmm_phaster_.tsv', sep='\t')

# Add a new column for signal peptide probability and fill it with corresponding values
df['Signalp3'] = df['Locus_Tag'].map(gene_probabilities)

# Write the updated DataFrame back to a TSV file
df.to_csv('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4+P3_tmhmm_deeptmhmm_phaster.tsv', sep='\t', index=False)
```
```bash

```