# Candidatus liberibacter solanacearum assembly

Investigate scaffold 18 of the T.apicales assembly:
```bash
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta | grep -v '>' | wc -c #1,918,207
```

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

#### FastANI
```bash
#CLso-ZC1
fastANI -t 16 -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna -q /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/CLso-ZC1_ANI.txt
#FIN114
fastANI -t 16 -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic.fna -q /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/FIN114_ANI.txt
#FIN111
fastANI -t 16 -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic.fna -q /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/FIN111_ANI.txt

#CLso-ZC1
fastANI -t 16 --visualize -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna -q /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/CLso-ZC1_ANI.out
#FIN114
fastANI -t 16 --visualize -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic.fna -q /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/FIN114_ANI.out 
#FIN111
fastANI -t 16 --visualize -r /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic.fna -q /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/FIN111_ANI.out 

interactive
#CLso-ZC1
singularity exec --overlay /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/overlays/genoPlotR.sif:ro /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/genomation1.34.0.sif Rscript ~/git_repos/Scripts/NBI/visualize.R /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/CLso-ZC1_ANI.out.visual
#FIN114
singularity exec --overlay /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/overlays/genoPlotR.sif:ro /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/genomation1.34.0.sif Rscript ~/git_repos/Scripts/NBI/visualize.R /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/FIN114_ANI.out.visual
#FIN111
singularity exec --overlay /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/overlays/genoPlotR.sif:ro /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/genomation1.34.0.sif Rscript ~/git_repos/Scripts/NBI/visualize.R /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/FIN111_ANI.out.visual

fastANI -t 1 --visualize -r /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -q /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/CLso-ZC1_ANI2.out

singularity exec --overlay /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/overlays/genoPlotR.sif:ro /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/genomation1.34.0.sif Rscript ~/git_repos/Scripts/NBI/visualize.R /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/CLso-ZC1_ANI2.out.visual
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
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart.fasta
OutDir=$(dirname $Assembly)/minimap2
Outfile=$(basename $Assembly | sed 's@.bp.p_ctg.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
sbatch $ProgDir/run_minimap2-hifi.sh $OutDir $Outfile $Assembly $Read1 $Read2 
#58754389,6072012,6193942

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

seqkit stats /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
seqkit stats /jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
samtools flagstat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/minimap2/assembly_v1.fa.bam

sbatch $ProgDir/run_qualimap.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/minimap2/assembly_v1_corrected_fixstart.fasta.bam /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/minimap2/qualimap NA
#6193935, 6234026


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
The reference C.liberibacter solanacearum assembly was aligned to scaffold_18 to investigate possilbe merging of multiple haplotypes:
```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file temp_id.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart.fasta --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart_scaffold_18.fasta

#Reference vs the scaffold
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart_scaffold_18.fasta
Query=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/nucmer
OutFile=scaff18-1
nucmer  -p $OutDir/$OutFile $Reference $Query
show-coords -T $OutDir/$OutFile.delta > $OutDir/$OutFile.coords
awk '{print $8 "\t" $1 "\t" $2}' $OutDir/$OutFile.coords | sed '1,4d'| tr -s ' ' '\t' > $OutDir/$OutFile.bed
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/nucmer/scaff18-1.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/nucmer/scaff18-1.delta

#The scaffold vs the reference
Reference=/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna
Query=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart_scaffold_18.fasta
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/nucmer
OutFile=scaff18-2
nucmer  -p $OutDir/$OutFile $Reference $Query
show-coords -T $OutDir/$OutFile.delta > $OutDir/$OutFile.coords
awk 'NR > 4 {print $1 "\t" $2-1}' $OutDir/$OutFile.coords > $OutDir/$OutFile.bed
cat $OutDir/$OutFile.bed | tr -s ' ' '\t' > $OutDir/$OutFile.bed
mummerplot -l -c /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/nucmer/scaff18-2.delta
mummerplot -color /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/nucmer/scaff18-2.delta

/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic.fna
/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic.fna
```
```python
import os
from Bio import SeqIO

input_file = "/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic.fna"
input_file = "/jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic.fna"

#Iterate through each sequence in the multi-fasta file and write each to a separate file
input_dir = os.path.dirname(os.path.abspath(input_file))
for record in SeqIO.parse(input_file, "fasta"):
    output_file = os.path.join(input_dir, f"{record.id}.fna")
    with open(output_file, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
```
```bash
#Reference vs the scaffold
for Query in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/LVWB01000*.fna /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/LWEB01000* ); do
Reference=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart_scaffold_18.fasta
OutFile=$(echo $Query | sed 's@.fna@@g')
nucmer  -p $OutFile $Reference $Query
show-coords -T $OutFile.delta > $OutFile.coords
awk '{print $8 "\t" $1 "\t" $2}' $OutFile.coords | sed '1,4d'| tr -s ' ' '\t' > $OutFile.bed
done
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
Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)
Datatype=hifi
Correct_Datatype=pacbio-hifi
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_inspector.sh $OutFile $OutDir $Genome $Datatype $Correct_Datatype $Read1 $Read2
#59074264

circlator fixstart --verbose /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v2

source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v2.fasta -p circ
mummerplot -color -layout circ.delta

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCF_000183665.1/GCF_000183665.1_ASM18366v1_genomic.fna -p circ2
mummerplot  -l -c circ2.delta
```
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

Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected.fa
Read_type=pacbio-corrected
OutDir=$(dirname $Assembly)/circlator-hifi
Outfile=$(basename $Assembly | sed 's@.fa@@g')
Read1=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/TrAp2_hifi_reads.fastq.gz
Read2=/jic/research-groups/Saskia-Hogenhout/reads/genomic/CALIBER_PB_HIFI_July_2022/third_flow_cell/TrAp2_hifi_3rdSMRTcell.fastq.gz
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_circlator.sh $Assembly $Read_type $OutDir $Outfile $Read1 $Read2
#6054081, 6069287

circlator fixstart --verbose /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart

source package 70b0e328-5a66-4c7c-971b-b2face8a50d4
source package 09b2c824-1ef0-4879-b4d2-0a04ee1bbd6d

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart.fasta -p circl
mummerplot -color -layout circl.delta
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

cat Liberibacter/phages/*.fna | cut -d ' ' -f1 > temp_phages2.fasta
cat temp_phages2.fasta | grep '>'
>CP060690.1
>CP124119.1
>JF773396.1
>KX879601.1
>KX879602.1
>KY661963.1
>MT899443.1
>MT899444.1
>NC_019549.1
>NC_019550.1
>OP072401.1
>OP072565.1
>OP072666.1

nucmer --maxmatch --nosimplify /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta temp_phages2.fasta -p phage
mummerplot -color -layout phage.delta
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
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pgap
~/pgap.py --no-internet --no-self-update --docker singularity --container-path /hpc-home/did23faz/pgap.sif -r -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pgap -g /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa -s 'Candidatus Liberibacter solanacearum C'

./pgap.py --no-internet --no-self-update --docker singularity --container-path /home/did23faz/pgap_2023-10-03.build7061.sif -r -o mg37_results -g $HOME/.pgap/test_genomes/MG37/ASM2732v1.annotation.nucleotide.1.fasta -s 'Mycoplasmoides genitalium'

source package 1413a4f0-44e3-4b9d-b6c6-0f5c0048df88
python3 pgap.py --no-internet --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -n -o mg37_results -g ~/.pgap/test_genomes/MG37/ASM2732v1.annotation.nucleotide.1.fasta -s 'Mycoplasmoides genitalium'
python3 ~/pgap.py --debug --no-internet --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -r -o mg37_results -g ~/.pgap/test_genomes/MG37/ASM2732v1.annotation.nucleotide.1.fasta -s 'Mycoplasmoides genitalium'

python3 ~/pgap.py --no-internet --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -n -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pgap -g /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa -s 'Candidatus liberibacter'

~/pgap.py --debug --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -r -o mg37_results -g ~/.pgap/test_genomes/MG37/ASM2732v1.annotation.nucleotide.1.fasta -s 'Mycoplasmoides genitalium'

#19248434 + 19248439,19248441,19248442 (gruffalo)
#19250214, 19260710, 19260713  (gruffalo), 19261328, 19261425, 19261429, 19263765, 19264042, 19264056, 19264063, 19264068, 19264075: with 64 cpus

srun -p short  -c 64 --mem 50G --pty pgap
~/pgap.py --debug --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -r -o JIC1 -g /mnt/shared/scratch/theaven/apps/pgap/assembly_v2.fasta -s 'Candidatus liberibacter'
~/pgap.py --debug --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -r -o GCA_000496595 -g /mnt/shared/scratch/theaven/apps/pgap/GCA_000496595.1_ASM49659v1_genomic.fna -s 'Candidatus liberibacter'
~/pgap.py --debug --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -r -o GCA_003160765 -g /mnt/shared/scratch/theaven/apps/pgap/GCA_003160765.1_ASM316076v1_genomic.fna -s 'Candidatus liberibacter'
~/pgap.py --debug --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -r -o GCA_003336865 -g /mnt/shared/scratch/theaven/apps/pgap/GCA_003336865.1_ASM333686v1_genomic.fna -s 'Candidatus liberibacter'
~/pgap.py --debug --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -r -o GCA_018282155 -g /mnt/shared/scratch/theaven/apps/pgap/GCA_018282155.1_ASM1828215v1_genomic.fna -s 'Candidatus liberibacter'
~/pgap.py --debug --no-self-update --docker singularity --container-path ~/pgap_2023-10-03.build7061.sif -r -o GCA_025938115 -g /mnt/shared/scratch/theaven/apps/pgap/GCA_025938115.1_ASM2593811v1_genomic.fna -s 'Candidatus liberibacter'
```
INFO:    Using cached SIF image
PGAP version 2023-10-03.build7061 is up to date.
Output will be placed in: /mnt/shared/scratch/theaven/apps/pgap/mg37_results
WARNING: open files is less than the recommended value of 8000
PGAP failed, docker exited with rc = 1
Unable to find error in log file.

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

source package /tgac/software/testing/bin/signalp-3.0
seqkit split2 -p 2 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa -o CLsoC_JIC
#/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa.split/assembly_v2.part_001.faa
#/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa.split/assembly_v2.part_002.faa
grep '>\|Signal peptide probability:' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/signalp3/output.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp3/output2.txt 
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

#### PHASTEST
```bash
for genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/*.fasta); do
#for genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoD_ISR100/ncbi_dataset/data/GCA_002918245.2/GCA_002918245.2_ASM291824v2_genomic_fixstart.fasta); do
ID=$(echo $genome | cut -d '/' -f10)
echo $ID
OutDir=$(dirname $genome)/phastest
mkdir $OutDir
contigs=$(grep '>' $genome | wc -l)
echo $contigs
#if [ "$contigs" -eq 1 ]; then
#wget --post-file=$genome "https://phastest.ca/phastest_api" -O ${ID}_phastest
#else
#wget --post-file=$genome "https://phastest.ca/phastest_api?contigs=1" -O ${ID}_phastest
#fi
done

mkdir phastest
for genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/*.fasta); do
OutDir=$(dirname $genome)/phastest
ID=$(echo $genome | cut -d '/' -f10)
zipfile=$(ls $ID*.zip)
unzip -q $zipfile -d $OutDir
done

#from ubuntu:
for genome in $(ls phastest/*.fasta); do
ID=$(echo $genome | cut -d '/' -f2 | sed 's@.fasta@@g')
OutDir=$(dirname $genome)/phastest
contigs=$(grep '>' $genome | wc -l)
if [ "$contigs" -eq 1 ]; then
wget --post-file=$genome "https://phastest.ca/phastest_api" -O ${ID}_phastest
else
wget --post-file=$genome "https://phastest.ca/phastest_api?contigs=1" -O ${ID}_phastest
fi
done

for sub in $(ls *_phastest); do
job_id=$(cat $sub | cut -d '"' -f4)
wget "https://phastest.ca/phastest_api?acc=$job_id" -O ${ID}_phastest_results
done

for genome in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/*.fasta; do
    # Define the output directory
    OutDir=$(dirname "$genome")/phastest

    # Create the output directory if it doesn't exist
    mkdir -p "$OutDir"

    # Extract the ID from the genome path
    ID=$(echo $genome | cut -d '/' -f10)

    # Find the corresponding zip file in the current directory
    zipfile=$(ls | grep "^${ID}_.*\.zip$")

    # Check if the zip file exists and unzip it into the output directory
    if [[ -f "$zipfile" ]]; then
        unzip -o -q "$zipfile" -d "$OutDir"
        echo "Unzipped $zipfile into $OutDir"
    else
        echo "Zip file for $genome not found in current directory"
    fi
done

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/phastest/phage_regions.fna); do
OutDir=$(dirname $Genome)/prokka
ID=$(echo $Genome | cut -d '/' -f10)
Locustag=$ID
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_prokka.sh $Genome $OutDir ${ID}_phage_regions $Locustag
done
#60539805-60539872

rm /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/phastest/prokka/CLsoD_ISR100_phage_regions*
```
#### pyani
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani

for INPUT_DIR in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/phastest); do
MULTIFASTA_FILE="phage_regions.fna"
PREFIX=$(echo $INPUT_DIR | cut -d '/' -f10)
OUTPUT_DIR=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani
cd "$INPUT_DIR" || exit
awk -v prefix="$PREFIX" -v outdir="$OUTPUT_DIR" '
BEGIN { 
    RS=">"; 
    ORS=""; 
}
{
    if (NR > 1) {
        # Split the record into header and sequence
        split($0, arr, "\n", seps);
        # Remove any trailing whitespace (including tabs) from the header
        header = arr[1];
        gsub(/[ \t]+$/, "", header); 
        gsub(/[ \t]+/, "_", header); 
        clean_header = gensub(/[:\-]/, "_", "g", header);
        # Reconstruct the sequence from the remaining array elements
        sequence = "";
        for (i = 2; i <= length(arr); i++) {
            sequence = sequence arr[i] "\n";
        }
        file_name = outdir"/"prefix"_"clean_header".fna";
        # Print the sequence to the file
        print ">"header"\n"sequence > file_name;
        # Optional: Debugging output
        print "Creating file: " file_name;
    }
}' "$MULTIFASTA_FILE"
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani/${PREFIX}*); do
sed -i -e "s@>@>${PREFIX}_@g" "$file"
done
done

source package 8be892ef-2dcb-4c74-83e0-01d769d50047
average_nucleotide_identity.py -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani/ -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani/TETRA_output -m TETRA -g --gformat jpeg,jpg,pdf,png

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani2
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phages/*.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani2/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani/*.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani2/.
for file in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani2/*.fna; do
OutFile=$(basename $file | cut -d '_' -f1,2,3).fna
mv $file $(dirname $file)/$(echo $OutFile | sed 's@.fna.fna@.fna@g')
done
average_nucleotide_identity.py -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani2/ -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani2/TETRA_output -m TETRA -g --gformat jpeg,jpg,pdf,png

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3
```
```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

# Load data
data = pd.read_csv('/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani2/TETRA_output/TETRA_correlations.tab', delimiter='\t', header=0, index_col=0)

# Calculate linkage for clustering
distance_matrix = pdist(data, metric='euclidean')
linkage_matrix = linkage(distance_matrix, method='complete')

# Create a mask for values below 0.7
mask = data < 0.1

# Define a custom blue-to-red color map, with white for values below 0.7
cmap = LinearSegmentedColormap.from_list("custom_cmap", 
                                         [(0.0, "white"), 
                                          (0.7, "lightblue"), 
                                          (0.9, "purple"), 
                                          (1.0, "red")])

# Create the main figure and axes for heatmap
fig, ax = plt.subplots(figsize=(15, 15))

# Create a separate axes for the dendrogram
dendrogram_ax = fig.add_axes([0.05, 0.95, 0.9, 0.05])

# Plot the dendrogram and get the leaf ordering
dendrogram(linkage_matrix, orientation='top', labels=data.index, ax=dendrogram_ax)
dendrogram_ax.axis('off')

# Extract the order of labels from the dendrogram
dendrogram_leaves_order = dendrogram_ax.get_xticklabels()
dendrogram_leaves_order = [label.get_text() for label in dendrogram_leaves_order]

# Reorder the DataFrame based on the dendrogram leaves order
data_reordered = data.reindex(index=dendrogram_leaves_order, columns=dendrogram_leaves_order)

# Create the heatmap with the reordered data
sns.heatmap(data_reordered, mask=mask, cmap=cmap, annot=False, fmt=".2f", vmin=0.7, vmax=1, cbar=False, ax=ax)

# Customize colorbar
cbar = plt.colorbar(ax.collections[0], ax=ax, orientation='vertical', pad=0.1, fraction=0.02, aspect=40, shrink=0.8, location='left')
cbar.set_label('Color Scale', rotation=270, labelpad=20)

# Ensure a tick and label for every heatmap block on the x-axis
ax.set_xticks(np.arange(len(data_reordered)) + 0.5)  

# Directly set the x-axis labels using the index of the reordered DataFrame
ax.set_xticklabels(data_reordered.columns, rotation=90) 

ax.xaxis.set_label_position("top")

# Adjust the x-axis label size
ax.xaxis.set_tick_params(labelsize=5) 

# Ensure a tick and label for every heatmap block on the y-axis
ax.set_yticks(np.arange(len(data_reordered)) + 0.5)

# Directly set the y-axis labels using the index of the reordered DataFrame
ax.set_yticklabels(data_reordered.index, rotation=0) 

ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

# Adjust the y-axis label size
ax.yaxis.set_tick_params(labelsize=5)

# Save the figure
plt.savefig("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani2/TETRA_output/TETRA_correlations-mod2.png", dpi=300)

data_reordered.to_csv("/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/pyani2/TETRA_output/TETRA_correlations_reordered.tab", sep='\t')

```
#### Prophage hunter
```bash
source package d6092385-3a81-49d9-b044-8ffb85d0c446
source package 0dd71e29-8eb1-4512-b37c-42f7158718f4
source package 187d32ff-ff9f-4e1f-a711-3df9284d38b8
source package /tgac/software/production/bin/cufflinks-2.2.1
source package /nbi/software/testing/bin/R-3.5.2
source package /nbi/software/testing/bin/genemark-4.33_ES_ET

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prophage-hunter
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prophage-hunter
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v1.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prophage-hunter/Input.fasta
sh ~/Prophage-Hunter/prophage_hunter_RUN.sh
```
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


input_file = '/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp3/output2.txt'

gene_probabilities = {}

with open(input_file, 'r') as file:
  gene_name = ""
  signal_peptide_probability = ""
  for line in file:
    if line.startswith('>CLsoC_JIC_'):
      gene_name = line.strip().split()[0][1:] 
      if gene_name not in gene_probabilities:
        gene_probabilities[gene_name] = ""
    elif "Signal peptide probability" in line:
      signal_peptide_probability = line.strip().split(':')[-1].strip()  
      gene_probabilities[gene_name] = signal_peptide_probability

output_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp3/output.tsv"

with open(output_file, "w") as f:
    f.write("Gene:\tSignalp3:\n")
    for gene, probability in gene_probabilities.items():
        f.write(f"{gene}\t{probability}\n")
```
```bash
grep '^#' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_secreted.gff
grep -Ff <(awk -F'\t' '$11 ~ /^[0-9.]+$/ && $11 > 0.5 {print $1}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_with_SignalP6+P5+P4+P3_tmhmm_deeptmhmm_phaster.tsv) /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.gff >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_secreted.gff
sed -i 's@gnl|JIC|WAE90_1@scaffold_18@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_secreted.gff
sed -i 's@gnl|JIC|WAE90_2@scaffold_184@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_secreted.gff
sed -i 's@gnl|JIC|WAE90_3@scaffold_1142@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0_secreted.gff

sed -i 's@gnl|JIC|WAE90_1@scaffold_18@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.gff
sed -i 's@gnl|JIC|WAE90_2@scaffold_184@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.gff
sed -i 's@gnl|JIC|WAE90_3@scaffold_1142@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/JIC_LsoCDa_1.0.gff
```
#### MCScanX
```bash
#Find the contiguity of Liberibacter assemblies:
for file in $(ls Liberibacter/genomes/C*/ncbi_dataset/data/GCA_*/*.fna); do
echo $(echo $file | cut -d '/' -f3)
cat $file | grep '>' | wc -l
cat $file | grep -v '>' | wc -c
done

for file in $(ls Liberibacter/genomes/C*/ncbi_dataset/data/GCA_*/*.fna); do
Out=$(echo $file | sed 's@.fna@_fixstart@g')
circlator fixstart --verbose $file $Out
done

mkdir -p temp_download/CLsoC_JIC1
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/assembly_v2.fasta temp_download/CLsoC_JIC1/.


for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/C*/ncbi_dataset/data/GCA_*/*fixstart.fasta); do
OutDir=$(dirname $file)/prokka
OutFile=$(basename $file | sed 's@.fasta@@g')
Locustag=$(echo $file | cut -d '/' -f10)
mkdir $OutDir
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_prokka.sh $file $OutDir $OutFile $Locustag
done
#59376486, 59376512-85




for file in $(ls Liberibacter/genomes/CL*/ncbi_dataset/data/GCA_*/*.fna); do
echo $(echo $file | cut -d '/' -f3)
grep '>' $file | wc -l
cat $file | wc -c 
done

CLsoB_ZC1 #1

CLeu_ASUK1 #29
CLeu_ASNZ1 #15

CLct_Oxford #17

CLcr_BT-1 #1
CLcr_BT-0 #1

CLbr_Asol15 #25

CLas_TaiYZ2
CLas_ReuSP1
CLas_PYN
CLas_psy62 #1
CLas_PGD
CLas_JXGC #1
CLas_JRPAMB1
CLas_Ishi-1 #1
CLas_gxpsy #1
CLas_GDCZ
CLas_CoFLP
CLas_AHCA1
CLas_A4

CLam_SaoPaulo #1
CLam_PW_SP #22

CLaf_PTSAPSY #1
CLaf_Ang37 #1

nucmer --maxmatch --nosimplify Liberibacter/genomes/CLaf_Ang37/ncbi_dataset/data/GCA_017869345.1/GCA_017869345.1_ASM1786934v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/1
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLaf_PTSAPSY/ncbi_dataset/data/GCA_001021085.1/GCA_001021085.1_ASM102108v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/66
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLam_PW_SP/ncbi_dataset/data/GCA_000350385.1/GCA_000350385.1_Velvet_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/2
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLam_SaoPaulo/ncbi_dataset/data/GCA_000496595.1/cds_from_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/3
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLam_SaoPaulo/ncbi_dataset/data/GCA_000496595.1/GCA_000496595.1_ASM49659v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/4
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_9PA/ncbi_dataset/data/GCA_013778575.1/GCA_013778575.1_ASM1377857v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/5
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_A4/ncbi_dataset/data/GCA_000590865.3/GCA_000590865.3_ASM59086v3_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/6
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_AHCA17/ncbi_dataset/data/GCA_009859045.1/GCA_009859045.1_ASM985904v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/7
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_AHCA1/ncbi_dataset/data/GCA_003143875.1/GCA_003143875.1_ASM314387v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/8
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_A-SBCA19/ncbi_dataset/data/GCA_014892655.1/GCA_014892655.1_ASM1489265v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/9
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_AS-TNSK3/ncbi_dataset/data/GCA_029948395.1/GCA_029948395.1_ASM2994839v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/10
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_BCSMX/ncbi_dataset/data/GCA_025606285.1/GCA_025606285.1_ASM2560628v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/11
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_CHUC/ncbi_dataset/data/GCA_009756785.1/GCA_009756785.1_ASM975678v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/12
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_CoFLP/ncbi_dataset/data/GCA_014107775.1/GCA_014107775.1_ASM1410777v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/13
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_CRCFL16/ncbi_dataset/data/GCA_009756805.1/GCA_009756805.1_ASM975680v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/14
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_DUR1TX1/ncbi_dataset/data/GCA_009756745.1/GCA_009756745.1_ASM975674v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/15
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_DUR2TX1/ncbi_dataset/data/GCA_009756725.1/GCA_009756725.1_ASM975672v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/16
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_FL17/ncbi_dataset/data/GCA_000820625.1/GCA_000820625.1_ASM82062v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/17
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_GDCZ/ncbi_dataset/data/GCA_030585885.1/GCA_030585885.1_ASM3058588v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/18
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_GFR3TX3/ncbi_dataset/data/GCA_009756735.1/GCA_009756735.1_ASM975673v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/19
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_gxpsy/ncbi_dataset/data/GCA_000346595.1/GCA_000346595.1_ASM34659v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/21
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_HHCA16/ncbi_dataset/data/GCA_009756845.1/GCA_009756845.1_ASM975684v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/22
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_HHCA/ncbi_dataset/data/GCA_000724755.2/GCA_000724755.2_HHCA-assembly_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/23
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_Ishi-1/ncbi_dataset/data/GCA_000829355.1/GCA_000829355.1_ASM82935v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/24
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_JRPAMB1/ncbi_dataset/data/GCA_013462975.1/GCA_013462975.1_ASM1346297v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/25
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_JXGC/ncbi_dataset/data/GCA_002216815.1/GCA_002216815.1_ASM221681v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/26
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_JXGZ-1/ncbi_dataset/data/GCA_009764765.1/GCA_009764765.1_ASM976476v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/27
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_LBR19TX2/ncbi_dataset/data/GCA_009756885.1/GCA_009756885.1_ASM975688v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/28
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_LBR23TX5/ncbi_dataset/data/GCA_009756915.1/GCA_009756915.1_ASM975691v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/29
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_MAG1/ncbi_dataset/data/GCA_025938115.1/GCA_025938115.1_ASM2593811v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/30
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_Mex8/ncbi_dataset/data/GCA_009756755.1/GCA_009756755.1_ASM975675v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/31
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_MFL16/ncbi_dataset/data/GCA_009756815.1/GCA_009756815.1_ASM975681v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/32
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_PA19/ncbi_dataset/data/GCA_013309695.2/GCA_013309695.2_ASM1330969v2_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/33
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_PA20/ncbi_dataset/data/GCA_016758155.2/GCA_016758155.2_ASM1675815v2_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/34
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_PGD/ncbi_dataset/data/GCA_028473705.1/GCA_028473705.1_ASM2847370v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/35
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_psy62/ncbi_dataset/data/GCA_000023765.2/GCA_000023765.2_ASM2376v2_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/36
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_PYN/ncbi_dataset/data/GCA_028473725.1/GCA_028473725.1_ASM2847372v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/37
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_ReuSP1/ncbi_dataset/data/GCA_022220845.1/GCA_022220845.1_ASM2222084v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/38
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_SGCA16/ncbi_dataset/data/GCA_009756855.1/GCA_009756855.1_ASM975685v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/39
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_SGCA1/ncbi_dataset/data/GCA_003149415.1/GCA_003149415.1_ASM314941v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/40
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_SGCA5/ncbi_dataset/data/GCA_001430705.1/GCA_001430705.1_ASM143070v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/41
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_SGpsy/ncbi_dataset/data/GCA_003336865.1/GCA_003336865.1_ASM333686v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/42
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_Tabriz3/ncbi_dataset/data/GCA_022343665.1/GCA_022343665.1_ASM2234366v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/43
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_TaiYZ2/ncbi_dataset/data/GCA_014217975.1/GCA_014217975.1_ASM1421797v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/44
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_TX1712/ncbi_dataset/data/GCA_003160765.1/GCA_003160765.1_ASM316076v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/45
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_TX2351/ncbi_dataset/data/GCA_001969535.1/GCA_001969535.1_ASM196953v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/46
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_YCPsy/ncbi_dataset/data/GCA_001296945.1/GCA_001296945.1_ASM129694v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/47
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_YNHK-2/ncbi_dataset/data/GCA_018282155.1/GCA_018282155.1_ASM1828215v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/48
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_YNJS7C/ncbi_dataset/data/GCA_003615235.1/GCA_003615235.1_ASM361523v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/49
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_YNXP-1/ncbi_dataset/data/GCA_009764755.1/GCA_009764755.1_ASM976475v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/50
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLas_YTMX/ncbi_dataset/data/GCA_025606305.1/GCA_025606305.1_ASM2560630v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/51
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLbr_Asol15/ncbi_dataset/data/GCA_036858155.1/GCA_036858155.1_ASM3685815v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/52
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLcr_BT-0/ncbi_dataset/data/GCA_001543305.1/GCA_001543305.1_ASM154330v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/53
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLcr_BT-1/ncbi_dataset/data/GCA_000325745.1/GCA_000325745.1_ASM32574v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/54
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLct_Oxford/ncbi_dataset/data/GCA_016808295.1/GCA_016808295.1_ASM1680829v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/55
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLeu_ASNZ1/ncbi_dataset/data/GCA_003045065.1/GCA_003045065.1_ASM304506v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/56
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLeu_ASUK1/ncbi_dataset/data/GCA_019843875.1/GCA_019843875.1_ASM1984387v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/57
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLsoA_JNVH01/ncbi_dataset/data/GCA_000756225.1/GCA_000756225.1_ASM75622v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/58
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLsoA_NZ1/ncbi_dataset/data/GCA_000968085.1/GCA_000968085.1_ASM96808v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/59
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLsoA_RSTM/ncbi_dataset/data/GCA_001414235.1/GCA_001414235.1_ASM141423v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/60
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLsoB_HenneA/ncbi_dataset/data/GCA_000968075.1/GCA_000968075.1_ASM96807v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/61
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/GCA_000183665.1_ASM18366v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/62
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/63
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/64
nucmer --maxmatch --nosimplify Liberibacter/genomes/CLsoD_ISR100/ncbi_dataset/data/GCA_002918245.2/GCA_002918245.2_ASM291824v2_genomic.fna temp_phages2.fasta -p Liberibacter/phages/nucmer/65

mummerplot -color -layout Liberibacter/phages/nucmer/1.delta #nothing visible
mummerplot -color -layout Liberibacter/phages/nucmer/2.delta #nothing
mummerplot -color -layout Liberibacter/phages/nucmer/3.delta #nothing
mummerplot -color -layout Liberibacter/phages/nucmer/4.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/5.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/6.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/7.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/8.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/9.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/10.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/11.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/12.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/13.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/14.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/15.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/16.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/17.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/18.delta #1,2,3,4,5,6,7,8,9,10 #
mummerplot -color -layout Liberibacter/phages/nucmer/19.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/21.delta #1,2,3,4,5,6,7,8,9,10 ##
mummerplot -color -layout Liberibacter/phages/nucmer/22.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/23.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/24.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/25.delta #1,2,3,4,5,6,7,8,9,10 #
mummerplot -color -layout Liberibacter/phages/nucmer/26.delta #1,2,3,4,5,6,7,8,9,10 #
mummerplot -color -layout Liberibacter/phages/nucmer/27.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/28.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/29.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/30.delta #1,2,3,4,5,6,7,8,9,10 low similarity
mummerplot -color -layout Liberibacter/phages/nucmer/31.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/32.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/33.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/34.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/35.delta #1,2,3,4,5,6,7,8,9,10 #
mummerplot -color -layout Liberibacter/phages/nucmer/36.delta #1,2,3,4,5,6,7,8,9,10 #
mummerplot -color -layout Liberibacter/phages/nucmer/37.delta #1,2,3,4,5,6,7,8,9,10 #
mummerplot -color -layout Liberibacter/phages/nucmer/38.delta #1,2,3,4,5,6,7,8,9,10 #
mummerplot -color -layout Liberibacter/phages/nucmer/39.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/40.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/41.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/42.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/43.delta #nothing
mummerplot -color -layout Liberibacter/phages/nucmer/44.delta #1,2,3,4,5,6,7,8,9,10 #
mummerplot -color -layout Liberibacter/phages/nucmer/45.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/46.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/47.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/48.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/49.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/50.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/51.delta #1,2,3,4,5,6,7,8,9,10
mummerplot -color -layout Liberibacter/phages/nucmer/52.delta #nothing
mummerplot -color -layout Liberibacter/phages/nucmer/53.delta #nothing
mummerplot -color -layout Liberibacter/phages/nucmer/54.delta #nothing
mummerplot -color -layout Liberibacter/phages/nucmer/55.delta #nothing
mummerplot -color -layout Liberibacter/phages/nucmer/56.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/57.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/58.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/59.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/60.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/61.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/62.delta #nothing visible
mummerplot -color -layout Liberibacter/phages/nucmer/63.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/64.delta #small alignments
mummerplot -color -layout Liberibacter/phages/nucmer/65.delta #small hits for 6
mummerplot -color -layout Liberibacter/phages/nucmer/66.delta #nothing visible



```
```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/master

awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoD_ISR100/ncbi_dataset/data/GCA_002918245.2/prokka/GCA_002918245.2_ASM291824v2_genomic_fixstart.gff | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene;locus_tag@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/CLsoD_ISR100.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.gff | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene;locus_tag@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/CLsoC_JIC.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/prokka/GCA_001983675.1_ASM198367v1_genomic_fixstart.gff | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene;locus_tag@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/CLsoC_FIN114.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/prokka/GCA_001983655.1_ASM198365v1_genomic_fixstart.gff | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene;locus_tag@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/CLsoC_FIN111.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/prokka/GCA_000183665.1_ASM18366v1_genomic_fixstart.gff | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene;locus_tag@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/CLsoB_ZC1.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_HenneA/ncbi_dataset/data/GCA_000968075.1/prokka/GCA_000968075.1_ASM96807v1_genomic_fixstart.gff | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene;locus_tag@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/CLsoB_HenneA.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_RSTM/ncbi_dataset/data/GCA_001414235.1/prokka/GCA_001414235.1_ASM141423v1_genomic_fixstart.gff | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene;locus_tag@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/CLsoA_RSTM.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_NZ1/ncbi_dataset/data/GCA_000968085.1/prokka/GCA_000968085.1_ASM96808v1_genomic_fixstart.gff | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene;locus_tag@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/CLsoA_NZ1.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_JNVH01/ncbi_dataset/data/GCA_000756225.1/prokka/GCA_000756225.1_ASM75622v1_genomic_fixstart.gff | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene;locus_tag@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/CLsoA_JNVH01.gff

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/*.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/master/master.gff
< /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/master/master.gff tr ' ' '\t' > temp.gff && mv temp.gff /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/master/master.gff
```
```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoD_ISR100/ncbi_dataset/data/GCA_002918245.2/prokka/GCA_002918245.2_ASM291824v2_genomic_fixstart.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/CLsoD_ISR100 -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/CLsoC_JIC -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/prokka/GCA_001983675.1_ASM198367v1_genomic_fixstart.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/CLsoC_FIN114 -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/prokka/GCA_001983655.1_ASM198365v1_genomic_fixstart.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/CLsoC_FIN111 -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/prokka/GCA_000183665.1_ASM18366v1_genomic_fixstart.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/CLsoB_ZC1 -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_HenneA/ncbi_dataset/data/GCA_000968075.1/prokka/GCA_000968075.1_ASM96807v1_genomic_fixstart.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/CLsoB_HenneA -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_RSTM/ncbi_dataset/data/GCA_001414235.1/prokka/GCA_001414235.1_ASM141423v1_genomic_fixstart.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/CLsoA_RSTM -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_NZ1/ncbi_dataset/data/GCA_000968085.1/prokka/GCA_000968085.1_ASM96808v1_genomic_fixstart.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/CLsoA_NZ1 -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_JNVH01/ncbi_dataset/data/GCA_000756225.1/prokka/GCA_000756225.1_ASM75622v1_genomic_fixstart.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/CLsoA_JNVH01 -dbtype prot

for ID in CLsoD_ISR100 CLsoC_JIC CLsoC_FIN114 CLsoC_FIN111 CLsoB_ZC1 CLsoB_HenneA CLsoA_RSTM CLsoA_NZ1 CLsoA_JNVH01; do
  for ID2 in CLsoD_ISR100 CLsoC_JIC CLsoC_FIN114 CLsoC_FIN111 CLsoB_ZC1 CLsoB_HenneA CLsoA_RSTM CLsoA_NZ1 CLsoA_JNVH01; do
    DB=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/${ID}
    Query=$(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/${ID2}/*/data/*/prokka/*.faa)
    Out=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/${ID2}_v_${ID}.blast5
    blastp -db ${DB} -query ${Query} -num_threads 1 -evalue 1e-20 -num_alignments 200 -outfmt 6 -out ${Out}
  done
done
#59377072

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/intermediateData/*.blast5 > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/master/master.blast
< /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/master/master.blast tr ' ' '\t' > temp.blast && mv temp.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/master/master.blast
```
```bash
source package 038f5eb6-dc79-46b5-bb52-a86ed67aa64a
MCScanX /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/master/master
```
24525 matches imported (82310 discarded)
2469 pairwise comparisons
835 alignments generated

```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/intermediateData
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/master
for ID in CLsoB_ZC1 CLeu_ASUK1 CLeu_ASNZ1 CLct_Oxford CLcr_BT-1 CLcr_BT-0 CLbr_Asol15 CLas_psy62 CLas_JXGC CLas_Ishi-1 CLas_gxpsy CLam_SaoPaulo CLam_PW_SP CLaf_PTSAPSY CLaf_Ang37 CLsoC_JIC; do
  gff=$(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/${ID}/*/data/*/prokka/*.gff)
  awk '$3 == "gene" {split($9, parts, ";"); split(parts[1], subparts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, subparts[2]}' ${gff} | awk '{print $1, $9, $4, $5}' | sed 's@gnl|JIC|@@g' | sed 's@_gene@@g' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/intermediateData/${ID}.gff
done
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/intermediateData/*.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/master/master.gff
< /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/master/master.gff tr ' ' '\t' > temp.gff && mv temp.gff /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/master/master.gff

mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB
for ID in CLsoB_ZC1 CLeu_ASUK1 CLeu_ASNZ1 CLct_Oxford CLcr_BT-1 CLcr_BT-0 CLbr_Asol15 CLas_psy62 CLas_JXGC CLas_Ishi-1 CLas_gxpsy CLam_SaoPaulo CLam_PW_SP CLaf_PTSAPSY CLaf_Ang37 CLsoC_JIC; do
  Assembly=$(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/${ID}/*/data/*/prokka/*.faa)
  DB=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/ncbiDB/${ID}
  makeblastdb -in ${Assembly} -out ${DB} -dbtype prot
  for ID2 in CLsoB_ZC1 CLeu_ASUK1 CLeu_ASNZ1 CLct_Oxford CLcr_BT-1 CLcr_BT-0 CLbr_Asol15 CLas_psy62 CLas_JXGC CLas_Ishi-1 CLas_gxpsy CLam_SaoPaulo CLam_PW_SP CLaf_PTSAPSY CLaf_Ang37 CLsoC_JIC; do
    Query=$(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/${ID2}/*/data/*/prokka/*.faa)
    Out=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/intermediateData/${ID2}_v_${ID}.blast
    blastp -db ${DB} -query ${Query} -num_threads 32 -evalue 1e-10 -num_alignments 20 -outfmt 6 -out ${Out}
  done
done
#59379025
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/intermediateData/*.blast > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/master/master.blast
< /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/master/master.blast tr ' ' '\t' > temp.blast && mv temp.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/master/master.blast

source package 038f5eb6-dc79-46b5-bb52-a86ed67aa64a
MCScanX /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/analysis/synteny/mcscanx/interspecies/master/master
```
170176 matches imported (170385 discarded)
3208 pairwise comparisons
4152 alignments generated


To detect putative prophages, we used a combination of VirSorter2 (v2.2.2; Guo et al., 2021, p. 2), CheckV (v.0.7.0 (Nayfach et al., 2021)) and VIBRANT (v1.2.1; Kieft, Zhou & Anantharaman, 2020). First, all genomes were analyzed with VirSorter2 (settings –include-groups “dsDNAphage,ssDNA,NCLDV,laviviridae”). Resulting viral regions were retained if they scored at least 0.5 for double stranded DNA phage (dsDNA phage, n = 539). Host genome regions flanking these viral sequences were then trimmed with the CheckV ’contamination’ command. To exclude highly degraded phages, trimmed sequences were retained only if they were at least 5 kb in length (n = 462). Of these resulting sequences, a region was considered a putative prophage if it scored at least 0.9 with VirSorter2 (n = 357). Additionally, those with VirSorter scores of 0.5–0.9 were further analyzed with VIBRANT and were retained as putative prophage if VIBRANT also classified these sequences as virus (n = 74). This resulted in a final set of 431 putative prophages.
#### Virsorter1
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter1
```
#### VirSorter2
https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=1
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta); do
OutDir=$(dirname $Genome)/virsorter2
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_virsorter2.sh $Genome $OutDir
done
#59385582, 59426292, 59426508

source package 7a6ee408-8bf5-4cb5-9853-16d5ad415e8f
virsorter run --keep-original-seq -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -w /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter2-2 --include-groups dsDNAphage,ssDNA,NCLDV,laviviridae --min-length 5000 --min-score 0.5 -j 1 all

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter2-4
virsorter run --keep-original-seq -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -w /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter2-4 --include-groups dsDNAphage,ssDNA,NCLDV,laviviridae --min-length 5000 --min-score 0.5 -j 1 all --db-dir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/virsorter/db

[2024-04-04 04:04 CRITICAL] --db-dir must be provided since "template-config.yaml" has not been initialized


source package 1413a4f0-44e3-4b9d-b6c6-0f5c0048df88
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter2-3
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/virsorter2.2.4.sif virsorter run --keep-original-seq -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -w /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter2-3 --include-groups dsDNAphage,ssDNA,NCLDV,laviviridae --min-length 5000 --min-score 0.5 -j 1 all --db-dir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/virsorter/db

mkdir /home/theaven/scratch/uncompressed/hogenhoutvirsorter2-3
singularity exec /home/theaven/scratch/apps/virsorter/virsorter_2.2.4--pyhdfd78af_0 virsorter run --keep-original-seq -i /home/theaven/scratch/uncompressed/hogenhout/assembly_v2.fasta -w /home/theaven/scratch/uncompressed/hogenhoutvirsorter2-3 --include-groups dsDNAphage,ssDNA,NCLDV,laviviridae --min-length 5000 --min-score 0.5 -j 1 all

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter2-5
virsorter run --keep-original-seq -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -w /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter2-5 --include-groups dsDNAphage,ssDNA,NCLDV --min-length 5000 --min-score 0.5 -j 1 all --high-confidence-only 

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter2-6
virsorter run --keep-original-seq -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -w /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/virsorter2-6 --include-groups dsDNAphage,ssDNA,NCLDV -j 1 all 
```
```bash
conda activate viral-id-sop
```

#### VIBRANT + seperately
```bash
source package 98bd0520-6518-47be-92c1-20cf52f558e5
VIBRANT_run.py -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -folder /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/vibrant
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/vibrant/VIBRANT_integrated_prophage_coordinates_assembly_v2.tsv
awk '{print length($0)}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/vibrant/assembly_v2.phages_combined.fna
```
Prophages were identified and trimmed by removing flanking host regions using VIBRANT, geNomad, and CheckV. CheckV was run on all predicted phages and trimmed proviruses to evaluate completeness and quality. Contigs evaluated as low quality by both CheckV and VIBRANT
and had a geNomad viral score <0.9 with %1 geNomad viral hallmark gene were removed from the analysis. Contigs with eukaryotic viral taxonomies assigned by geNomad and/or CheckV were also removed. To generate a non-redundant phage reference database, all 32,401 phage scaffolds were dereplicated at 98% gANI over 85% of the phage genomes using dRep dereplicate (-sa 0.98 –ignoreGenomeQuality -l 4000 -nc 0.85 –clusterAlg single -N50W 0 -sizeW 1). Genomes with gANI R98% were classified as the same phage subspecies, and the phage genome with the highest dRep score was chosen as the representative genome from each subspecies, resulting in a total of 9,929 phage subspecies. - doi.org/10.1016/j.chom.2023.11.015
#### Genomad
```bash
srun -p jic-medium -c 8 --mem 32G --pty bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/genomad_1.8.0--pyhdfd78af_1 genomad end-to-end --restart -t 8 -v /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/genomad /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/genomad/10594875/genomad_db
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/genomad/assembly_v2_summary/assembly_v2_virus.fna
```
#### Checkv
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/vibrant/checkv
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/checkvv1.5.sif checkv end_to_end /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/vibrant/assembly_v2.phages_combined.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/vibrant/checkv -t 1 -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/checkv/checkv-db-v1.5
awk '{print length($0)}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/vibrant/checkv/proviruses.fna
awk '{print length($0)}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/vibrant/checkv/viruses.fna

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/genomad/checkv
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/checkvv1.5.sif checkv end_to_end /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/genomad/assembly_v2_summary/assembly_v2_virus.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/genomad/checkv -t 1 -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/checkv/checkv-db-v1.5

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/phastest/checkv
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/checkvv1.5.sif checkv end_to_end /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/phastest/phage_regions.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/phastest/checkv -t 1 -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/checkv/checkv-db-v1.5

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/checkv
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/checkvv1.5.sif checkv end_to_end /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/checkv -t 1 -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/checkv/checkv-db-v1.5

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/checkv/viruses.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/checkv/proviruses.fna > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/checkv/all.fna
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/checkv/checkv
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/checkvv1.5.sif checkv end_to_end /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/checkv/all.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/checkv/checkv -t 1 -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/checkv/checkv-db-v1.5
```
scaffold_18|provirus_1271814_1292265
scaffold_18|provirus_1894877_1922238

#### Prophage hunter

#### PhageBoost
```bash
source package 6ea65583-8b83-4561-8709-6ec8df0a0cb3
PhageBoost -f /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/phageboost

#processing: assembly_v2
#time after genecalls: 6.56554388999939
#time after feature calculations: 13.704559803009033
#time after predictions: 15.467231512069702
#{'phage34': 101, 'phage21': 76, 'phage0': 27, 'phage11': 26, 'phage36': 25, 'phage23': 25, 'phage31': 24, 'phage35': 23}

PhageBoost -f /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/temp_scaffold_18.fasta -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/phageboost2
#processing: temp_scaffold_18
#time after genecalls: 1.440387487411499
#time after feature calculations: 8.455031633377075
#time after predictions: 9.60255742073059
#{'phage21': 107, 'phage14': 35, 'phage30': 28, 'phage4': 26, 'phage25': 22, 'phage23': 17, 'phage19': 17, 'phage18': 13}

```
#### blast
```bash
cat phages/*.fna > phages.fna
makeblastdb -in phages.fna -dbtype nucl -title phagedb -out phagedb
sbatch $ProgDir/run_blastn.sh /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito.fa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phagedb /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter phagedb-6 999999
#61963769, 61965893

bowtie2 -x temp_scaffold_18 -f Liberibacter/phages.fna -S aligned_reads.sam -p 16
#  16 (100.00%) were unpaired; of these:
#    16 (100.00%) aligned 0 times
#    0 (0.00%) aligned exactly 1 time
#    0 (0.00%) aligned >1 times
#0.00% overall alignment rate

nucmer --prefix=alignment temp_scaffold_18.fasta Liberibacter/phages.fna
show-coords -rcl alignment.delta > alignment.coords
mummerplot -l -c alignment.delta
mummerplot -color alignment.delta
```
#### Depht
```bash
conda activate depht
```
### Orthology analysis
```bash
for assembly in $(ls -d Liberibacter/genomes/*); do
ID=$(echo $assembly | cut -d '/' -f3)
echo $ID
prokka=$(ls $assembly/*/data/*/prokka/*.faa)
grep '>' $prokka | wc -l
#proteins=$(ls $assembly/*/data/*/*.faa)
done

#pgap?:
cp Liberibacter/genomes/CLas_YNHK-2/ncbi_dataset/data/GCA_018282155.1/GCA_018282155.1_ASM1828215v1_genomic.fna lib_dwn/.
cp Liberibacter/genomes/CLas_TX1712/ncbi_dataset/data/GCA_003160765.1/GCA_003160765.1_ASM316076v1_genomic.fna lib_dwn/.
cp Liberibacter/genomes/CLas_SGpsy/ncbi_dataset/data/GCA_003336865.1/GCA_003336865.1_ASM333686v1_genomic.fna lib_dwn/.
cp Liberibacter/genomes/CLas_MAG1/ncbi_dataset/data/GCA_025938115.1/GCA_025938115.1_ASM2593811v1_genomic.fna lib_dwn/.
cp Liberibacter/genomes/CLam_SaoPaulo/ncbi_dataset/data/GCA_000496595.1/GCA_000496595.1_ASM49659v1_genomic.fna lib_dwn/.
cp Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta lib_dwn/.

ProjDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae
cd $ProjDir
IsolateAbrv=All_liberibacter
WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted
mkdir -p $WorkDir/goodProteins
mkdir -p $WorkDir/badProteins  
cd $WorkDir/formatted

Id_field=1
Taxon_code=Af01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLaf_Ang37/ncbi_dataset/data/GCA_017869345.1/prokka/GCA_017869345.1_ASM1786934v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Af02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLaf_PTSAPSY/ncbi_dataset/data/GCA_001021085.1/prokka/GCA_001021085.1_ASM102108v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Am01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLam_PW_SP/ncbi_dataset/data/GCA_000350385.1/prokka/GCA_000350385.1_Velvet_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Am02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLam_SaoPaulo/ncbi_dataset/data/GCA_000496595.1/prokka/GCA_000496595.1_ASM49659v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_9PA/ncbi_dataset/data/GCA_013778575.1/prokka/GCA_013778575.1_ASM1377857v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_A4/ncbi_dataset/data/GCA_000590865.3/prokka/GCA_000590865.3_ASM59086v3_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As03
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_AHCA17/ncbi_dataset/data/GCA_009859045.1/prokka/GCA_009859045.1_ASM985904v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As04
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_AHCA1/ncbi_dataset/data/GCA_003143875.1/prokka/GCA_003143875.1_ASM314387v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As05
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_A-SBCA19/ncbi_dataset/data/GCA_014892655.1/prokka/GCA_014892655.1_ASM1489265v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As06
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_AS-TNSK3/ncbi_dataset/data/GCA_029948395.1/prokka/GCA_029948395.1_ASM2994839v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As07
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_BCSMX/ncbi_dataset/data/GCA_025606285.1/prokka/GCA_025606285.1_ASM2560628v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As08
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_CHUC/ncbi_dataset/data/GCA_009756785.1/prokka/GCA_009756785.1_ASM975678v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As09
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_CoFLP/ncbi_dataset/data/GCA_014107775.1/prokka/GCA_014107775.1_ASM1410777v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As10
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_CRCFL16/ncbi_dataset/data/GCA_009756805.1/prokka/GCA_009756805.1_ASM975680v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As11
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_DUR1TX1/ncbi_dataset/data/GCA_009756745.1/prokka/GCA_009756745.1_ASM975674v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As12
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_DUR2TX1/ncbi_dataset/data/GCA_009756725.1/prokka/GCA_009756725.1_ASM975672v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As13
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_FL17/ncbi_dataset/data/GCA_000820625.1/prokka/GCA_000820625.1_ASM82062v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As14
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_GDCZ/ncbi_dataset/data/GCA_030585885.1/prokka/GCA_030585885.1_ASM3058588v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As15
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_GFR3TX3/ncbi_dataset/data/GCA_009756735.1/prokka/GCA_009756735.1_ASM975673v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As16
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_gxpsy/ncbi_dataset/data/GCA_000346595.1/prokka/GCA_000346595.1_ASM34659v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As17
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_HHCA16/ncbi_dataset/data/GCA_009756845.1/prokka/GCA_009756845.1_ASM975684v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As18
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_HHCA/ncbi_dataset/data/GCA_000724755.2/prokka/GCA_000724755.2_HHCA-assembly_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As19
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_Ishi-1/ncbi_dataset/data/GCA_000829355.1/prokka/GCA_000829355.1_ASM82935v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As20
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_JRPAMB1/ncbi_dataset/data/GCA_013462975.1/prokka/GCA_013462975.1_ASM1346297v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As21
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_JXGC/ncbi_dataset/data/GCA_002216815.1/prokka/GCA_002216815.1_ASM221681v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As22
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_JXGZ-1/ncbi_dataset/data/GCA_009764765.1/prokka/GCA_009764765.1_ASM976476v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As23
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_LBR19TX2/ncbi_dataset/data/GCA_009756885.1/prokka/GCA_009756885.1_ASM975688v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As24
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_LBR23TX5/ncbi_dataset/data/GCA_009756915.1/prokka/GCA_009756915.1_ASM975691v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As25
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_MAG1/ncbi_dataset/data/GCA_025938115.1/prokka/GCA_025938115.1_ASM2593811v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As26
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_Mex8/ncbi_dataset/data/GCA_009756755.1/prokka/GCA_009756755.1_ASM975675v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As27
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_MFL16/ncbi_dataset/data/GCA_009756815.1/prokka/GCA_009756815.1_ASM975681v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As28
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PA19/ncbi_dataset/data/GCA_013309695.2/prokka/GCA_013309695.2_ASM1330969v2_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As29
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PA20/ncbi_dataset/data/GCA_016758155.2/prokka/GCA_016758155.2_ASM1675815v2_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As30
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PGD/ncbi_dataset/data/GCA_028473705.1/prokka/GCA_028473705.1_ASM2847370v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As31
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_psy62/ncbi_dataset/data/GCA_000023765.2/prokka/GCA_000023765.2_ASM2376v2_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As32
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PYN/ncbi_dataset/data/GCA_028473725.1/prokka/GCA_028473725.1_ASM2847372v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As33
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_ReuSP1/ncbi_dataset/data/GCA_022220845.1/prokka/GCA_022220845.1_ASM2222084v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As34
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGCA16/ncbi_dataset/data/GCA_009756855.1/prokka/GCA_009756855.1_ASM975685v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As35
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGCA1/ncbi_dataset/data/GCA_003149415.1/prokka/GCA_003149415.1_ASM314941v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As36
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGCA5/ncbi_dataset/data/GCA_001430705.1/prokka/GCA_001430705.1_ASM143070v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As37
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGpsy/ncbi_dataset/data/GCA_003336865.1/prokka/GCA_003336865.1_ASM333686v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As38
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_Tabriz3/ncbi_dataset/data/GCA_022343665.1/prokka/GCA_022343665.1_ASM2234366v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As39
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_TaiYZ2/ncbi_dataset/data/GCA_014217975.1/prokka/GCA_014217975.1_ASM1421797v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As40
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_TX1712/ncbi_dataset/data/GCA_003160765.1/prokka/GCA_003160765.1_ASM316076v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As41
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_TX2351/ncbi_dataset/data/GCA_001969535.1/prokka/GCA_001969535.1_ASM196953v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As42
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YCPsy/ncbi_dataset/data/GCA_001296945.1/prokka/GCA_001296945.1_ASM129694v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As43
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YNHK-2/ncbi_dataset/data/GCA_018282155.1/prokka/GCA_018282155.1_ASM1828215v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As44
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YNJS7C/ncbi_dataset/data/GCA_003615235.1/prokka/GCA_003615235.1_ASM361523v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As45
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YNXP-1/ncbi_dataset/data/GCA_009764755.1/prokka/GCA_009764755.1_ASM976475v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As46
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YTMX/ncbi_dataset/data/GCA_025606305.1/prokka/GCA_025606305.1_ASM2560630v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Br01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLbr_Asol15/ncbi_dataset/data/GCA_036858155.1/prokka/GCA_036858155.1_ASM3685815v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Cr01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLcr_BT-0/ncbi_dataset/data/GCA_001543305.1/prokka/GCA_001543305.1_ASM154330v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Cr02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLcr_BT-1/ncbi_dataset/data/GCA_000325745.1/prokka/GCA_000325745.1_ASM32574v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Ct01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLct_Oxford/ncbi_dataset/data/GCA_016808295.1/prokka/GCA_016808295.1_ASM1680829v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Eu01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLeu_ASNZ1/ncbi_dataset/data/GCA_003045065.1/prokka/GCA_003045065.1_ASM304506v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Eu02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLeu_ASUK1/ncbi_dataset/data/GCA_019843875.1/prokka/GCA_019843875.1_ASM1984387v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_JNVH01/ncbi_dataset/data/GCA_000756225.1/prokka/GCA_000756225.1_ASM75622v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_NZ1/ncbi_dataset/data/GCA_000968085.1/prokka/GCA_000968085.1_ASM96808v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So03
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_RSTM/ncbi_dataset/data/GCA_001414235.1/prokka/GCA_001414235.1_ASM141423v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So04
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_HenneA/ncbi_dataset/data/GCA_000968075.1/prokka/GCA_000968075.1_ASM96807v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So05
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/prokka/GCA_000183665.1_ASM18366v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So06
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/prokka/GCA_001983655.1_ASM198365v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So07
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/prokka/GCA_001983675.1_ASM198367v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So08
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So09
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoD_ISR100/ncbi_dataset/data/GCA_002918245.2/prokka/GCA_002918245.2_ASM291824v2_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
cd $ProjDir

for WorkDir in $(ls -d analysis/orthology/orthofinder/*); do
Input_dir=$WorkDir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
~/git_repos/Scripts/NBI/orthomclFilterFasta.pl $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
done

Prefix=All_liberibacter
OutDir=orthofinder52000
sbatch ~/git_repos/Wrappers/NBI/run_orthofinder.sh $Input_dir $Prefix $OutDir
#60226609

wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups.tsv #2116
grep 'So08' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups.tsv | wc -l #1,095
grep 'So08' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups.tsv > CLsoC_JIC_OGs.tsv



####################################################################################################################################################################
ProjDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae
cd $ProjDir
IsolateAbrv=All_liberibacter_phage
WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted
mkdir -p $WorkDir/goodProteins
mkdir -p $WorkDir/badProteins  
cd $WorkDir/formatted

Id_field=1
Taxon_code=Af01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLaf_Ang37/ncbi_dataset/data/GCA_017869345.1/prokka/GCA_017869345.1_ASM1786934v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Af02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLaf_PTSAPSY/ncbi_dataset/data/GCA_001021085.1/prokka/GCA_001021085.1_ASM102108v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Am01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLam_PW_SP/ncbi_dataset/data/GCA_000350385.1/prokka/GCA_000350385.1_Velvet_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Am02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLam_SaoPaulo/ncbi_dataset/data/GCA_000496595.1/prokka/GCA_000496595.1_ASM49659v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_9PA/ncbi_dataset/data/GCA_013778575.1/prokka/GCA_013778575.1_ASM1377857v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_A4/ncbi_dataset/data/GCA_000590865.3/prokka/GCA_000590865.3_ASM59086v3_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As03
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_AHCA17/ncbi_dataset/data/GCA_009859045.1/prokka/GCA_009859045.1_ASM985904v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As04
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_AHCA1/ncbi_dataset/data/GCA_003143875.1/prokka/GCA_003143875.1_ASM314387v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As05
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_A-SBCA19/ncbi_dataset/data/GCA_014892655.1/prokka/GCA_014892655.1_ASM1489265v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As06
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_AS-TNSK3/ncbi_dataset/data/GCA_029948395.1/prokka/GCA_029948395.1_ASM2994839v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As07
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_BCSMX/ncbi_dataset/data/GCA_025606285.1/prokka/GCA_025606285.1_ASM2560628v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As08
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_CHUC/ncbi_dataset/data/GCA_009756785.1/prokka/GCA_009756785.1_ASM975678v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As09
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_CoFLP/ncbi_dataset/data/GCA_014107775.1/prokka/GCA_014107775.1_ASM1410777v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As10
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_CRCFL16/ncbi_dataset/data/GCA_009756805.1/prokka/GCA_009756805.1_ASM975680v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As11
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_DUR1TX1/ncbi_dataset/data/GCA_009756745.1/prokka/GCA_009756745.1_ASM975674v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As12
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_DUR2TX1/ncbi_dataset/data/GCA_009756725.1/prokka/GCA_009756725.1_ASM975672v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As13
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_FL17/ncbi_dataset/data/GCA_000820625.1/prokka/GCA_000820625.1_ASM82062v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As14
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_GDCZ/ncbi_dataset/data/GCA_030585885.1/prokka/GCA_030585885.1_ASM3058588v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As15
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_GFR3TX3/ncbi_dataset/data/GCA_009756735.1/prokka/GCA_009756735.1_ASM975673v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As16
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_gxpsy/ncbi_dataset/data/GCA_000346595.1/prokka/GCA_000346595.1_ASM34659v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As17
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_HHCA16/ncbi_dataset/data/GCA_009756845.1/prokka/GCA_009756845.1_ASM975684v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As18
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_HHCA/ncbi_dataset/data/GCA_000724755.2/prokka/GCA_000724755.2_HHCA-assembly_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As19
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_Ishi-1/ncbi_dataset/data/GCA_000829355.1/prokka/GCA_000829355.1_ASM82935v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As20
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_JRPAMB1/ncbi_dataset/data/GCA_013462975.1/prokka/GCA_013462975.1_ASM1346297v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As21
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_JXGC/ncbi_dataset/data/GCA_002216815.1/prokka/GCA_002216815.1_ASM221681v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As22
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_JXGZ-1/ncbi_dataset/data/GCA_009764765.1/prokka/GCA_009764765.1_ASM976476v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As23
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_LBR19TX2/ncbi_dataset/data/GCA_009756885.1/prokka/GCA_009756885.1_ASM975688v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As24
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_LBR23TX5/ncbi_dataset/data/GCA_009756915.1/prokka/GCA_009756915.1_ASM975691v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As25
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_MAG1/ncbi_dataset/data/GCA_025938115.1/prokka/GCA_025938115.1_ASM2593811v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As26
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_Mex8/ncbi_dataset/data/GCA_009756755.1/prokka/GCA_009756755.1_ASM975675v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As27
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_MFL16/ncbi_dataset/data/GCA_009756815.1/prokka/GCA_009756815.1_ASM975681v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As28
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PA19/ncbi_dataset/data/GCA_013309695.2/prokka/GCA_013309695.2_ASM1330969v2_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As29
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PA20/ncbi_dataset/data/GCA_016758155.2/prokka/GCA_016758155.2_ASM1675815v2_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As30
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PGD/ncbi_dataset/data/GCA_028473705.1/prokka/GCA_028473705.1_ASM2847370v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As31
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_psy62/ncbi_dataset/data/GCA_000023765.2/prokka/GCA_000023765.2_ASM2376v2_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As32
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PYN/ncbi_dataset/data/GCA_028473725.1/prokka/GCA_028473725.1_ASM2847372v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As33
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_ReuSP1/ncbi_dataset/data/GCA_022220845.1/prokka/GCA_022220845.1_ASM2222084v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As34
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGCA16/ncbi_dataset/data/GCA_009756855.1/prokka/GCA_009756855.1_ASM975685v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As35
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGCA1/ncbi_dataset/data/GCA_003149415.1/prokka/GCA_003149415.1_ASM314941v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As36
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGCA5/ncbi_dataset/data/GCA_001430705.1/prokka/GCA_001430705.1_ASM143070v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As37
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGpsy/ncbi_dataset/data/GCA_003336865.1/prokka/GCA_003336865.1_ASM333686v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As38
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_Tabriz3/ncbi_dataset/data/GCA_022343665.1/prokka/GCA_022343665.1_ASM2234366v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As39
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_TaiYZ2/ncbi_dataset/data/GCA_014217975.1/prokka/GCA_014217975.1_ASM1421797v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As40
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_TX1712/ncbi_dataset/data/GCA_003160765.1/prokka/GCA_003160765.1_ASM316076v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As41
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_TX2351/ncbi_dataset/data/GCA_001969535.1/prokka/GCA_001969535.1_ASM196953v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As42
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YCPsy/ncbi_dataset/data/GCA_001296945.1/prokka/GCA_001296945.1_ASM129694v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As43
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YNHK-2/ncbi_dataset/data/GCA_018282155.1/prokka/GCA_018282155.1_ASM1828215v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As44
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YNJS7C/ncbi_dataset/data/GCA_003615235.1/prokka/GCA_003615235.1_ASM361523v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As45
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YNXP-1/ncbi_dataset/data/GCA_009764755.1/prokka/GCA_009764755.1_ASM976475v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As46
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YTMX/ncbi_dataset/data/GCA_025606305.1/prokka/GCA_025606305.1_ASM2560630v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Br01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLbr_Asol15/ncbi_dataset/data/GCA_036858155.1/prokka/GCA_036858155.1_ASM3685815v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Cr01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLcr_BT-0/ncbi_dataset/data/GCA_001543305.1/prokka/GCA_001543305.1_ASM154330v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Cr02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLcr_BT-1/ncbi_dataset/data/GCA_000325745.1/prokka/GCA_000325745.1_ASM32574v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Ct01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLct_Oxford/ncbi_dataset/data/GCA_016808295.1/prokka/GCA_016808295.1_ASM1680829v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Eu01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLeu_ASNZ1/ncbi_dataset/data/GCA_003045065.1/prokka/GCA_003045065.1_ASM304506v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Eu02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLeu_ASUK1/ncbi_dataset/data/GCA_019843875.1/prokka/GCA_019843875.1_ASM1984387v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So01
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_JNVH01/ncbi_dataset/data/GCA_000756225.1/prokka/GCA_000756225.1_ASM75622v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So02
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_NZ1/ncbi_dataset/data/GCA_000968085.1/prokka/GCA_000968085.1_ASM96808v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So03
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_RSTM/ncbi_dataset/data/GCA_001414235.1/prokka/GCA_001414235.1_ASM141423v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So04
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_HenneA/ncbi_dataset/data/GCA_000968075.1/prokka/GCA_000968075.1_ASM96807v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So05
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/prokka/GCA_000183665.1_ASM18366v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So06
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/prokka/GCA_001983655.1_ASM198365v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So07
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/prokka/GCA_001983675.1_ASM198367v1_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So08
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So09
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoD_ISR100/ncbi_dataset/data/GCA_002918245.2/prokka/GCA_002918245.2_ASM291824v2_genomic_fixstart.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field

Id_field=1
Taxon_code=Af01p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLaf_Ang37/ncbi_dataset/data/GCA_017869345.1/phastest/prokka/CLaf_Ang37_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Af02p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLaf_PTSAPSY/ncbi_dataset/data/GCA_001021085.1/phastest/prokka/CLaf_PTSAPSY_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Am01p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLam_PW_SP/ncbi_dataset/data/GCA_000350385.1/phastest/prokka/CLam_PW_SP_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Am02p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLam_SaoPaulo/ncbi_dataset/data/GCA_000496595.1/phastest/prokka/CLam_SaoPaulo_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As01p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_9PA/ncbi_dataset/data/GCA_013778575.1/phastest/prokka/CLas_9PA_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As02p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_A4/ncbi_dataset/data/GCA_000590865.3/phastest/prokka/CLas_A4_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As03p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_AHCA17/ncbi_dataset/data/GCA_009859045.1/phastest/prokka/CLas_AHCA17_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As04p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_AHCA1/ncbi_dataset/data/GCA_003143875.1/phastest/prokka/CLas_AHCA1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As05p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_A-SBCA19/ncbi_dataset/data/GCA_014892655.1/phastest/prokka/CLas_A-SBCA19_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As06p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_AS-TNSK3/ncbi_dataset/data/GCA_029948395.1/phastest/prokka/CLas_AS-TNSK3_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As07p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_BCSMX/ncbi_dataset/data/GCA_025606285.1/phastest/prokka/CLas_BCSMX_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As08p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_CHUC/ncbi_dataset/data/GCA_009756785.1/phastest/prokka/CLas_CHUC_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As09p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_CoFLP/ncbi_dataset/data/GCA_014107775.1/phastest/prokka/CLas_CoFLP_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As10p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_CRCFL16/ncbi_dataset/data/GCA_009756805.1/phastest/prokka/CLas_CRCFL16_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As11p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_DUR1TX1/ncbi_dataset/data/GCA_009756745.1/phastest/prokka/CLas_DUR1TX1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As12p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_DUR2TX1/ncbi_dataset/data/GCA_009756725.1/phastest/prokka/CLas_DUR2TX1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As13p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_FL17/ncbi_dataset/data/GCA_000820625.1/phastest/prokka/CLas_FL17_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As14p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_GDCZ/ncbi_dataset/data/GCA_030585885.1/phastest/prokka/CLas_GDCZ_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As15p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_GFR3TX3/ncbi_dataset/data/GCA_009756735.1/phastest/prokka/CLas_GFR3TX3_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As16p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_gxpsy/ncbi_dataset/data/GCA_000346595.1/phastest/prokka/CLas_gxpsy_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As17p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_HHCA16/ncbi_dataset/data/GCA_009756845.1/phastest/prokka/CLas_HHCA16_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As18p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_HHCA/ncbi_dataset/data/GCA_000724755.2/phastest/prokka/CLas_HHCA_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As19p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_Ishi-1/ncbi_dataset/data/GCA_000829355.1/phastest/prokka/CLas_Ishi-1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As20p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_JRPAMB1/ncbi_dataset/data/GCA_013462975.1/phastest/prokka/CLas_JRPAMB1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As21p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_JXGC/ncbi_dataset/data/GCA_002216815.1/phastest/prokka/CLas_JXGC_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As22p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_JXGZ-1/ncbi_dataset/data/GCA_009764765.1/phastest/prokka/CLas_JXGZ-1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As23p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_LBR19TX2/ncbi_dataset/data/GCA_009756885.1/phastest/prokka/CLas_LBR19TX2_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As24p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_LBR23TX5/ncbi_dataset/data/GCA_009756915.1/phastest/prokka/CLas_LBR23TX5_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As25p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_MAG1/ncbi_dataset/data/GCA_025938115.1/phastest/prokka/CLas_MAG1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As26p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_Mex8/ncbi_dataset/data/GCA_009756755.1/phastest/prokka/CLas_Mex8_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As27p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_MFL16/ncbi_dataset/data/GCA_009756815.1/phastest/prokka/CLas_MFL16_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As28p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PA19/ncbi_dataset/data/GCA_013309695.2/phastest/prokka/CLas_PA19_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As29p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PA20/ncbi_dataset/data/GCA_016758155.2/phastest/prokka/CLas_PA20_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As30p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PGD/ncbi_dataset/data/GCA_028473705.1/phastest/prokka/CLas_PGD_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As31p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_psy62/ncbi_dataset/data/GCA_000023765.2/phastest/prokka/CLas_psy62_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As32p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_PYN/ncbi_dataset/data/GCA_028473725.1/phastest/prokka/CLas_PYN_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As33p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_ReuSP1/ncbi_dataset/data/GCA_022220845.1/phastest/prokka/CLas_ReuSP1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As34p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGCA16/ncbi_dataset/data/GCA_009756855.1/phastest/prokka/CLas_SGCA16_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As35p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGCA1/ncbi_dataset/data/GCA_003149415.1/phastest/prokka/CLas_SGCA1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As36p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGCA5/ncbi_dataset/data/GCA_001430705.1/phastest/prokka/CLas_SGCA5_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As37p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_SGpsy/ncbi_dataset/data/GCA_003336865.1/phastest/prokka/CLas_SGpsy_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As38p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_Tabriz3/ncbi_dataset/data/GCA_022343665.1/phastest/prokka/CLas_Tabriz3_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As39p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_TaiYZ2/ncbi_dataset/data/GCA_014217975.1/phastest/prokka/CLas_TaiYZ2_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As40p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_TX1712/ncbi_dataset/data/GCA_003160765.1/phastest/prokka/CLas_TX1712_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As41p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_TX2351/ncbi_dataset/data/GCA_001969535.1/phastest/prokka/CLas_TX2351_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As42p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YCPsy/ncbi_dataset/data/GCA_001296945.1/phastest/prokka/CLas_YCPsy_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As43p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YNHK-2/ncbi_dataset/data/GCA_018282155.1/phastest/prokka/CLas_YNHK-2_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As44p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YNJS7C/ncbi_dataset/data/GCA_003615235.1/phastest/prokka/CLas_YNJS7C_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As45p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YNXP-1/ncbi_dataset/data/GCA_009764755.1/phastest/prokka/CLas_YNXP-1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=As46p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_YTMX/ncbi_dataset/data/GCA_025606305.1/phastest/prokka/CLas_YTMX_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Br01p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLbr_Asol15/ncbi_dataset/data/GCA_036858155.1/phastest/prokka/CLbr_Asol15_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Cr01p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLcr_BT-0/ncbi_dataset/data/GCA_001543305.1/phastest/prokka/CLcr_BT-0_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Cr02p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLcr_BT-1/ncbi_dataset/data/GCA_000325745.1/phastest/prokka/CLcr_BT-1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Ct01p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLct_Oxford/ncbi_dataset/data/GCA_016808295.1/phastest/prokka/CLct_Oxford_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Eu01p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLeu_ASNZ1/ncbi_dataset/data/GCA_003045065.1/phastest/prokka/CLeu_ASNZ1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=Eu02p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLeu_ASUK1/ncbi_dataset/data/GCA_019843875.1/phastest/prokka/CLeu_ASUK1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So01p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_JNVH01/ncbi_dataset/data/GCA_000756225.1/phastest/prokka/CLsoA_JNVH01_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So02p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_NZ1/ncbi_dataset/data/GCA_000968085.1/phastest/prokka/CLsoA_NZ1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So03p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_RSTM/ncbi_dataset/data/GCA_001414235.1/phastest/prokka/CLsoA_RSTM_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So04p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_HenneA/ncbi_dataset/data/GCA_000968075.1/phastest/prokka/CLsoB_HenneA_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So05p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/phastest/prokka/CLsoB_ZC1_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So06p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/phastest/prokka/CLsoC_FIN111_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So07p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/phastest/prokka/CLsoC_FIN114_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So08p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/phastest/prokka/CLsoC_JIC_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
Id_field=1
Taxon_code=So09p
Fasta_file=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoD_ISR100/ncbi_dataset/data/GCA_002918245.2/phastest/prokka/CLsoD_ISR100_phage_regions.faa
~/git_repos/Scripts/NBI/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
cd $ProjDir

for WorkDir in $(ls -d analysis/orthology/orthofinder/All_liberibacter_phage); do
Input_dir=$WorkDir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
~/git_repos/Scripts/NBI/orthomclFilterFasta.pl $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
done

Prefix=liberibacter_phages
OutDir=orthofinder52000
sbatch ~/git_repos/Wrappers/NBI/run_orthofinder.sh $Input_dir $Prefix $OutDir
#60545849

wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter_phage/formatted/orthofinder52000/Results_liberibacter_phages/Orthogroups/Orthogroups.tsv #2096
grep 'So08' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter_phage/formatted/orthofinder52000/Results_liberibacter_phages/Orthogroups/Orthogroups.tsv | wc -l #1,116
grep 'So08p' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter_phage/formatted/orthofinder52000/Results_liberibacter_phages/Orthogroups/Orthogroups.tsv | wc -l #257

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/phastest/predicted_phage_regions.json); do
cd $(dirname $file)
fasta=$(ls $(echo $file | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13)/*.fasta)
gff=$(ls $(echo $file | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13)/prokka/*.gff)
faa=$(ls $(echo $file | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13)/prokka/*.faa)
ID=$(echo $file | cut -d '/' -f10)
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/extract_phage_genes.py $file $fasta $gff $faa $ID
done

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/phastest/*_phage_genes.txt > phage_genes.txt
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups.tsv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage.tsv
while IFS= read -r gene; do
  echo $gene
    sed "s@|${gene}@p|${gene}@g" /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage.tsv > temp.tsv && mv temp.tsv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage.tsv
done < "phage_genes.txt"

grep 'p|' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage.tsv

wc -l /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage.tsv #2116
grep 'So08p' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage.tsv | wc -l #225

awk -F'\t' '$2 ~ /^[0-9.]+$/ && $2 > 0.5 {print $1}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/prokka/Signalp3/output.tsv | grep -v 'Gene:' > secreted_genes.txt
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage.tsv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage_sec.tsv
while IFS= read -r gene; do
  echo $gene
    sed "s@|${gene}@s|${gene}@g" /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage_sec.tsv > temp.tsv && mv temp.tsv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage_sec.tsv
done < "secreted_genes.txt"

grep 'So08s\|So08ps' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage_sec.tsv | wc -l #104
grep 'So08ps' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter/formatted/orthofinder52000/Results_All_liberibacter/Orthogroups/Orthogroups_phage_sec.tsv | wc -l #25

grep -Ff secreted_genes.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2_secreted.gff
sed -i 's@gnl|JIC|CLsoC_JIC_2@scaffold_18@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2_secreted.gff

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.gff /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2_all.gff
sed -i 's@gnl|JIC|CLsoC_JIC_1@scaffold_1142@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2_all.gff
sed -i 's@gnl|JIC|CLsoC_JIC_2@scaffold_18@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2_all.gff
sed -i 's@gnl|JIC|CLsoC_JIC_3@scaffold_182@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2_all.gff

grep -Ff secreted_genes.txt phage_genes.txt > secreted_phage_genes.txt
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file secreted_phage_genes.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/secreted_phage.faa

proteome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/secreted_phage.faa
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/uniprot_01102023/Uniprot_01102023_reference_proteomes.dmnd
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/blast
OutFile=$(basename $proteome | sed 's@.faa@@g')
hits=10
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_diamond_blastp.sh $proteome $Database $OutDir $OutFile $hits
#61946990
```
Unique orthogroups with secreted phage region genes:
So08ps|CLsoC_JIC_00127, So08|CLsoC_JIC_00247, So08p|CLsoC_JIC_00366, So08ps|CLsoC_JIC_00498, So08ps|CLsoC_JIC_00596, So08ps|CLsoC_JIC_00658, So08p|CLsoC_JIC_01171, So08s|CLsoC_JIC_01332, So08p|CLsoC_JIC_01388, So08p|CLsoC_JIC_01530, So08p|CLsoC_JIC_01697, So08p|CLsoC_JIC_01743, So08ps|CLsoC_JIC_01855, So08ps|CLsoC_JIC_01905, So08p|CLsoC_JIC_01978

So08ps|CLsoC_JIC_00591, So08p|CLsoC_JIC_00810, So08ps|CLsoC_JIC_01115, So08ps|CLsoC_JIC_01536, So08p|CLsoC_JIC_01861, So08|CLsoC_JIC_01911

So08s|CLsoC_JIC_00325, So08s|CLsoC_JIC_01830, So08ps|CLsoC_JIC_01880

So08ps|CLsoC_JIC_01170, So08s|CLsoC_JIC_01650, So08p|CLsoC_JIC_01698

```bash
CLsoC_JIC_00127
CLsoC_JIC_00247
CLsoC_JIC_00366
CLsoC_JIC_00498
CLsoC_JIC_00596
CLsoC_JIC_00658
CLsoC_JIC_01171
CLsoC_JIC_01332
CLsoC_JIC_01388
CLsoC_JIC_01530
CLsoC_JIC_01697
CLsoC_JIC_01743
CLsoC_JIC_01855
CLsoC_JIC_01905
CLsoC_JIC_01978
grep -Ff temp.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/OG0001087.gff

CLsoC_JIC_00591
CLsoC_JIC_00810
CLsoC_JIC_01115
CLsoC_JIC_01536
CLsoC_JIC_01861
CLsoC_JIC_01911
grep -Ff temp.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/OG0001229.gff

CLsoC_JIC_00325
CLsoC_JIC_01830
CLsoC_JIC_01880
grep -Ff temp.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/OG0001407.gff

CLsoC_JIC_01170
CLsoC_JIC_01650
CLsoC_JIC_01698
grep -Ff temp.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/OG0001410.gff

sed -i 's@gnl|JIC|CLsoC_JIC_2@scaffold_18@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/OG0001087.gff
sed -i 's@gnl|JIC|CLsoC_JIC_2@scaffold_18@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/OG0001229.gff
sed -i 's@gnl|JIC|CLsoC_JIC_2@scaffold_18@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/OG0001407.gff
sed -i 's@gnl|JIC|CLsoC_JIC_2@scaffold_18@g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/OG0001410.gff

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file temp.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/unique_OG.faa

proteome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/unique_OG.faa
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/uniprot_01102023/Uniprot_01102023_reference_proteomes.dmnd
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/blast
OutFile=$(basename $proteome | sed 's@.faa@@g')
hits=10
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_diamond_blastp.sh $proteome $Database $OutDir $OutFile $hits
#61946988
```
```bash
CLsoC_JIC_00683
CLsoC_JIC_00758
CLsoC_JIC_00661
CLsoC_JIC_00739
CLsoC_JIC_01694
CLsoC_JIC_01851
CLsoC_JIC_01901
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file temp.txt --input /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa --output /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/xxx_OG.faa
```

```bash
#OA2/OI2c(modified from GCCTCGCGACTTCGCAACCCAT) - 16S
perl ~/git_repos/Scripts/NBI/in_silico_PCR.pl -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -a GCGCTTATTTTTAATAGGAGCGGCA -b GCCTCACGACTTCGCAACCCAT > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/OA2-OI2c_results.txt 2> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/OA2-OI2c_amplicons.fasta
grep 'GCGCTTATTTTTAATAGGAGCGGCA' 
#CL514F/CL514R - 50S
perl ~/git_repos/Scripts/NBI/in_silico_PCR.pl -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -a CTCTAAGATTTCGGTTGGTT -b TATATCTATCGTTGCACCAG > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CL514F-CL514R_results.txt 2> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CL514F-CL514R_amplicons.fasta
#OMB1482F/OMB2086R - Outer membrane protein
perl ~/git_repos/Scripts/NBI/in_silico_PCR.pl -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta -a GGCGTGGTTATAAGCAGAGT -b ATCTACACGCGGACCTATAC > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/OMB1482F-OMB2086_Rresults.txt 2> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/OMB1482F-OMB2086_amplicons.fasta

perl ~/git_repos/Scripts/NBI/in_silico_PCR.pl -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/GCA_000183665.1_ASM18366v1_genomic_fixstart.fasta -a GCGCTTATTTTTAATAGGAGCGGCA -b GCCTCGCGACTTCGCAACCCAT > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/OA2-OI2c_results.txt 2> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/OA2-OI2c_amplicons.fasta

perl ~/git_repos/Scripts/NBI/in_silico_PCR.pl -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic_fixstart.fasta -a GCGCTTATTTTTAATAGGAGCGGCA -b GCCTCACGACTTCGCAACCCAT > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/OA2-OI2c_results.txt 2> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/OA2-OI2c_amplicons.fasta

perl ~/git_repos/Scripts/NBI/in_silico_PCR.pl -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic_fixstart.fasta -a GCGCTTATTTTTAATAGGAGCGGCA -b GCCTCACGACTTCGCAACCCAT > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/OA2-OI2c_results.txt 2> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/OA2-OI2c_amplicons.fasta
```
ATGGGTTGCGAAGTCGCGAGGC
TGCCGCTCCTATTAAAAATAAGCGC

AACTGTTTAGAGTGATTTTTTGAAATAAGAGACTATCGTCATAATATTCTATATATTCGG
ACGATTTTTGTGAGGGCTTTCCCTTATTTGATAAAAACATTAGCAATAAATCGTTTGGGC
GGAAACACAAACATTTAATTAAAAATATTTGTGTTTCCGTTGCACAGTCGTACTATTCAG
ATATAATTTTTTCTTAATAATATTTGTAGTGTGGATTTTATTGTCTTATTCATGGAATGT
TAAATCAAAGTATTATGAATTCCATAATTTCACGTAAATAACCATAAGGTATGATAGAGG
AAAGTGCGGAATCATACTTTTTAATAGAATAACTATTTGTTATTTTAATCATATTTCAAA
TTTTTATGTATCGTTAATAAAGAGATGTTTTTTGCGAATTTTTTCCGTTATTGATGAGGG
AATTAATTCTTTTTGTTTCCATAAGAAAACGTTGTAGTAGCTCTCCTTTTAAAGTTTCTC
CCTCGGGAATTCGAATTTTCATAGGATTTACCTTGATCCCGTTGACAATTAATTCATAAT
GTAGATGTGGTCCTGTTGCCAGTCCTGTTGCTCCAATCCAACCAATGATTTGTCCTTGTT
TAACGGTTGTTTTTTCTGTGACATTTTTAGAGATTGCATCTTGATGATTATAGGAAGAAA
CAAATCCATTTGCATGACGGATAATAGTTTGTTTTCCATAACCACCTGCCCAACCAGATT
TTTCAACAATACCATCACTAACAGCTATAATAGGTGTTCCTCTTGGTGCCGCCCAGTCCA
CTCCTGTATGCATTCGCGTATAGCCTAATATTGGATGACGTCGCATACCAAATTCAGAGG
TTATACGTCCAAATGGCACGGGGGTACGGAGTAAAAAAGGTCGAGAGCTTTTTCCATTTT
CATTAAAATATTCCACTGATCCATCAAGGGGGTTTTGTAAACGATAGAAACGAGTGCGTG
TTTCACCAAAGCGAGCATGAATATATAGTAATTCGGAATCATTGGTGCTTTGGTTTTTAA
CAGGATCAACGGAGAAAAAAGTTTCAAGAAAATCAGTTGGTTTAAGATGTTCCTGTAAAT
TTAGGTTGTTTGCTAGGATTCTGATGACTAGTTGTACCAGATTTTTATTCATTCCATTAA
AAGAAGTAGCGCGCCATATTCCATCATAAATACTTGGGAATGATTCAGATGTTCGTATGT
TATCTGTATAAGAATCAATATCTATTTTAGCGGGTTCTGCGCCTAAGACATATTCATTAT
TGTCGTTAAGGGCGATAGTTAATACATGCTCTTTTTTATGGTAAATGCTGAATCTTATTA
TTGTAGATTCTGCTTCTTTTTGTAGCATGCCAATTCTGATTGTTTCGTCTTTTGTTAGCT
TATCAACGTGTATTTTGTTTTTTAAGGCTTTTACTATTTTAGAGCTTTCACTGTTAGAAT
AACCAGCTTGGATCATAGCTTTAAAAATTGTGGTATTGTGTTGTATTACAATTAAATCAT
CGGCAAATTCAGGAGATATATTTGTTGATAATTGAGGTGTTGTAATCGTTCTATTTTCTT
CGATTACTTTAACTGTTGTACTATTAATAGTGCTATTATCTTGCTCTTTATATAAAGGCG
GATCGTCGGCATAATAGAGGTTATATGATTGATTTTTAGTATAATTTAAGGGAAAATATT
GGTTTATAATTGCTTTTTTGATTGTATCATCTTTTTCTATAATATCTTTTTTAATATTAG
TAACATCTGTTGGAAAATTGAATTTTTTGGTAGTTATTTCAAGCTTAGATTCAGTGTTAT
ATAGTGTATCAATGTTGTTAATTGTATCCATCAACATTTGCGATGTTGATTCTATTTTCC
CTCCAGAAAAAATTTTTAATGGATCAAATTTAGGATAATCTCTTACATTTTGATATGGTG
TTGCGAATGTCATTCTAGCATAAGCAAATGGAATTTTTTTGATGATATCTTTATTATGAG
ATTTGATGAGAGTGGGAACATCAATAATGATTTTTTCAGGATGTTTAGATTTTGTGTTTT
TTGTAGATAATCGCGCTTTTCGCATAGAAAAAGGGTTAGGATTAAAATGAGCAGAATTTA
GATCTAATGTTTGAGGGGATAATTTTGCAGGTATCGCTACTTTTTGATGTCCATCTAAGG
CGGTCAAAAGAGCCCCTCCCATGAGAACTCCAGAAGTGACTCCAGCAAAGAATGTTGTGT
ATACCCATCGAAGAGAAACCTTGCGTCTATTTAGAAATATTTTTTCGTTATTGTCTCCTA
GAATAGGAGGGTTATCACCAAACGAAAATAGTATTTCTTTTTTATTTAATGTACTATATA
ACGTCATTTTACACGATTAACATTCAAAAGATCCATTTTTCAACAGTGCATAATTTGTAA
GGGGATATAATTTTTCTTACTAAGAGTCAAATTATAGGGTTTGGGTATGGTTAATGGTTT
GAGGGTGGTAATATTTTTTTTGATTATTTTTTATTGACTTATTTATGGCAATTGTTTTAG
TATGGCGGATGTTATTTGAGAATAGAATAGAGAAGAAAGGGAGACATGGGCGGCGATTTT
TTATATTTGTAATTAAGGATTTAATTACGAAGGTAGGAGGCGTTTTTTTCATGTTTTTCT
GGTTGTTATGACGTTGAGAGATGTTTTAGCATTAGGATTTTTTTTACATGGTGGCGCCTT
TGGCAAGAACATGAGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCAGGCTTAACAC
ATGCAAGTCGAGCGCTTATTTTTAATAGGAGCGGCAGACGGGTGAGTAACGCGTGGGAAT
CTACCTTTTTCTACGGGATAACGCACGGAAACGTGTGCTAATACCGTATACGCCCTGAGA
AGGGGAAAGATTTATTGGAGAGAGATGAGCCCGCGTTAGATTAGCTAGTTGGTGGGGTAA
ATGCCTACCAAGGCTACGATCTATAGCTGGTCTGAGAGGACGATCAGCCACACTGGGACT
GAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGGGCAA
CCCTGATCCAGCCATGCCGCGTGAGTGAAGAAGGCCTTAGGGTTGTAAAGCTCTTTCGCC
GGAGAAGATAATGACGGTATCCGGAGAAGAAGTCCCGGCTAACTTCGTGCCAGCAGCCGC
GGTAATACGAAGGGGGCGAGCGTTGTTCGGAATAACTGGGCGTAAAGGGCGCGTAGGCGG
GTAATTAAGTTAGGGGTGAAATCCCAAGGCTCAACCTTGGAACTGCCTTTAATACTGGTT
ATCTAGAGTTTAGGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGCAGATATT
CGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCCTGATACTGACGCTGAGGCGCGAAA
GCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGAGTGCTA
GCTGTTGGGTGGTTTACCATTCAGTGGCGCAGCTAACGCATTAAGCACTCCGCCTGGGGA
GTACGGTCGCAAGATTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCA
TGTGGTTTAATTCGATGCAACGCGCAGAACCTTACCAGCCCTTGACATATAGAGGACGAT
ATCAGAGATGGTATTTTCTTTTCGGAGACCTTTATACAGGTGCTGCATGGCTGTCGTCAG
CTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTGCCTCTAGTTGC
CATCAAGTTTAGATTTTATCTAGATGTTGGGTACTTTATAGGGACTGCCGGTGATAATCC
GGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGGGCTGGGCTACACACGTG
CTACAATGGTGGTTACAATGGGTTGCGAAGTCGTGAGGCGGAGCTAATCCCAAAAGGCCA
TCTCAGTTCGGATTGCACTCTGCAACTCGAGTGCATGAAGTTGGAATCGCTAGTAATCGC
GGATCAGCATGCCGCGGTGAATACGTTCTCGGGCCTTGTACACACCGCCCGTCACACCAT
GGGAGTTGGTTTTGCTTGAAGACGGTGCGCTAACCGTAAGGAGGCAGCCGGCCACGGTAG
GGTCAGCGACTGGGGTGAAGTCGTAACAAGGTAGCCGTAGGGGAACCTGTGGCTGGATCA
CCTCCTTTCTAAGGAAGATGTTTAGTATTGTTAGATTTATTTAATGATCTGAGCATTTTT
TAAAGAATTAGAGCGATAAGCTCAAGAAAAGAATTTGTTAATTTTAGCAAGTTCTAAGGG
ATCGCCGTCCATGTTTCTCTTTCTTTTTGAATAATTTTGCGATTGATGGGGTCATTTGAG
TTTATGTTAAGGGCCCATAGCTCAGGCGGTTAGAGTGCACCCCTGATAAGGGTGAGGTCG
GTAGTTCGAATCTACCTGGGCCCACCATTCAATCAGGCAAGGGGCCGTAGCTCAGCTGGG
AGAGCGCCTGCTTTGCAAGCAGGATGTCAGCGGTTCGATCCCGCTCGGCTCCACCAATTG
CGAATTTATAGTTTTTTTGTTCTAGGGGATTTTTTTTTAGAGCAATAGTTTTTTGAAAAT
TGAATAGAAGGTAGATTTTTTTGTATTTTTCATATTGGCATTGTATGCGATATGGGAGGT
ACCGACGTTGTATAACCGCACGTTGAAGATTTATCTCAGGAAATTGGTCTATTGAAAGAG
CATAATTTATTTATGTTTTTTTAATTAAGAAACGTTTGTAATGAACTTTATGACGTATTG
ACAATGAGAGTGATCAAGCGCGATAAGGGCATTTGGTGGATGCCTTGGCATGCACAGGCG
ATGAAGGACGTAATACGCTGCGATAAGCTACGGGGAGCTGCAAATGAGCATTGATCCGTA
GATTTCCGAATGGGGCAACCCACCTTAGGTGTCTAGGAAAGTATACTATTAAGGTTTAAT
TTTCTAGGTACTTGAAGGTATCTTTACCTGAATAAAATAGGGTAAAAGAAGCGAACGCAG
GGAACTGAAACATCTAAGTACCTGTAGGAAAGGACATCAATTGAGACTCCGTTAGTAGTG
GCGAGCGAACGCGGATCAGGCCAGTGGTAGGGAAGATTTAAGTAGAATTATCTGGGAAGG
TAAGCCATAGAAGGTGATAGCCCCGTACACGTAATAATTTTTTCTATCCTTGAGTAGGGC
GGGACACGTGAAATCCTGTTTGAAGATGGGGCGACCACGCTCCAAGCCTAAGTACTCGTG
CATGACCGATAGTGAACTAGTACCGTGAGGGAAAGGCGAAAAGAACCCCTACTAGGGGAG
TGAAATAGACCCTGAAACCGAATGCCTACAAACAGTCGGAGGCTGTAAAGCTGACGGCGT
ACCTTTTGTATAATGGGTCAACGACTTAGTGTGGCAAGCGAGCTTAAGCCGATAGGTGTA
GGCGCAGCGAAAGCGAGTCTGAATAGGGCGTTTAGTTTGTTGCATTAGACCCGAAACCGA
GTGATCTAGCCATGAGCAGGTTGAAGGTTGGGTAACACCAATTGGAGGACCGAACCCGTA
TCTGTTGCAATAGATTGGGATGACTTGTGGCTAGGGGTGAAAGGCCAATCAAACTCGGAG
ATAGCTGGTTCTCCGCGAAATCTATTTAGGTAGAGCGTTAACTGAATACCCTCGGGGGTA
GAGCACTGGATAGGCTATGGGGGCTTACCGCCTTACTGATCCTAACCAAACTCCGAATAC
CGAGGAGTAATAGTTGGCAGACACACAGTGGGTGCTAACGTCCATTGTGGAGAGGGAAAC
AACCCTGACCTCCAGCTAAGGTCCCGAAGTCATGGCTAAGTGGGAAAGGAAGTGAAAATC
CCATAACAACCAGGATGTTGGCTTAGAAGCAGCCATCATTTAAAGAAAGCGTAGCAGCTC
ACTGGTCTAAAATTAAGGATTTTTGCGCCGAAAATGTAACGGGGCTCAAGCCATGCACCG
AAGCTGAGGATTTGTTTATTTTTTTTAAGATAAGCGAGTGGTAGCGGAGCGTTCCGTAAG
CTGATGAAGGAGGATCTGTGAGGACTTCTGGAGGTATCGGAAGTGAGAATGTTGACATGA
GTAACGATAAAGAGGGTGAGAAACCCTCTCGCCGAAAGACCAAGGGTTCCTGCTTAAAGT
TAATCTGAGCAGGGTTAGCCGGCCCCTAAGGTGAGGCGGAAACGCGTAGCTGATGGGAAC
CACATTAATATTTGTGGGCCTGGTGTAAGTGACGGATTAAGTATATTGTACATTTTTATT
GGATTAGATGTGCTTTGGATTAATTCCAGGAAATAGCTTCACCGTATAGACCGTACCCGA
AACCGACACAGGTGGTCAGGTAGAGTATACTAAGGCGCTTGAGAGAACTGCGTTGAAGGA
ACTCGGCAAATTGCACGCGTAACTTCGGGATAAGCGTGACCTTTTTTTGGGCAACCAGGA
GGAGGTGTCACAGATTAGGGGGTAGCGACTGTTTACCAAAAACACAGGGCTCTGCGAAGT
CGTAAGACGAAGTATAGGGCCTGACGCCTGCCCGGTGCTGGAAGGTTAATAGGAGGGGTG
AGAGCTCTGAATTGAAGCCCCAGTAAACGGCGGCCGTAACTATAACGGTCCTAAGGTAGC
GAAATTCCTTGTCGGGTAAGTTCCGACCTGCACGAATGGCGTAACGACTTCCCCACTGTC
TCCAACGCAGACTCAGTGAAATTGAATTCCCCGTGAAGATGCGGGGTTCCTGCGGTTAGA
CGGAAAGACCCCGTGCACCTTTACTATAGCTTTACATTGGCGTTTGCGTTGATATGTGTA
GGATAGGTGGTAGGCATTGAAGCAGAGGCGTTAGTTTTTGTGGAGCCATCCTTGAAATAC
CACCCTTATCCATGTGGACGTCTAACTGCGCTCCGTTATCCGGGGCCGGGACATTGTATG
GTGGGTAGTTTGACTGGGGCGGTCGCCTCCGAAAGAGTAACGGAGGCGCGCGATGGTAGG
CTCAGAGCGGTCAGAAATCGCTTGTTGAGTGCAATGGCATAAGCCTGCCTGACTGTGAGA
CTGACAAGTCGAGCAGAGACGAAAGTCGGTCATAGTGATCCGGTGGTTCCGAGTGGAAGG
GCCATCGCTCAACGGATAAAAGGTACGCCGGGGATAACAGGCTGATGACCCCCAAGAGCT
CATATCGACGGGGTTGTTTGGCACCTCGATGTCGGCTCATCGCATCCTGGGGCTGTAGAA
GGTCCCAAGGGTTTGGCTGTTCGCCAATTAAAGCGGTACGTGAGCTGGGTTCAGAACGTC
GTGAGACAGTTCGGTCCCTATCTGCCGTGGGTGTAGGAATATTGACAGGATCTTTCCCTA
GTACGAGAGGACCGGGATGGACGTATCTCTGGTGGACCTGTTGTTATGCCAATAGCATAG
CAGGGTAGCTAAATACGGAATGGATAACCGCTGAAAGCATCTAAGTGGGAAACCAACCTG
AAAACGAGTATTCCCTATCAGAGCTGTGGGAGACTACCACGTTGATAGGCTGGATGTGGA
AGCTAGGTAACTAGTGAAGCTGACCAGTACTAATAGCTCGATTGGCTTGATTGCTCTTAT
TGTCCATAGTCATATAATGTTTTAGACCATTTTCCATTCGTTTTTTATAGACTTGGTGGC
TTTTGCGGGGTTTCTGCACTCGTTCCCATTTCGAACACGGCCGTTAAATGCCCTAGCGCC
TATGGTACTTCATCTTAAGATGCGGGAGAGTCGGTCGCTGCCAGGTCTATAAAAATCGGG
TTCTTATTTTATTGTTATTTGTATATTGATCAAAATTATTGTCGCGGGGTGGAGCAGTCT
GGTAGCTCGTCAGGCTCATAACCTGAAGGTCGTGGGTTCAAATCCTACCCCCGCAACCAT
TGTATTATGGAATCTTATTACAATGAACCAAAACATAGAGCAACAAATTTTGGATGCGTT
GAAAATCCTGTACATACCAGGAGAAACAATTAATATCGTGGATATGAAAAGACTATCCAA
TATCTGCATAGTCCAAAATACAGCCTACCTATCCATAACAGTCCCTCATAATCTAGAAAA
ACAATTGCAATCCTTGCGCTTGAATGCACAACAAATTGTCCAAAATATCCCGCAAATTAA
AAATGCTGTCGTTACTCTAACTGAAAATAAAAGTAAACCTATCCTAGATCCAATAATTGA
AAATAAATTAAAAATTAATGCACTCATCGCAATAGCTTCCGGAAAAGGTGGCGTTGGAAA
ATCAACTACAGCCGTTAATCTCGCTTGTGCCTTGAAAAATAAAAATAAAAATGTCGCCAT
ACTCGATGCAGATATTTATGGTCCTTCCATACCAAAACTGCTACAATTATCTGGAAAAGC
AGAAATATTAGAAAAAAAATTCCTAAAGCCTATGGAAAATTATGGAATAAAAATAATGTC
CATGGCTTCTCTTGTAGATGATAATGTGGCCATGATATGGCGTGGTCCAATGGTACAATC
CGCTATAATGCATATGTTTCAGAATGTTTCTTGGGGTCAATTGGATTTTTTATTGATAGA
TATGCCTCCTGGAACAGGAGATGCTCATTTAACTGTTGCGCAAAAAATTCCTCTTTCTGG
TGTTGTGATCGTTTCTACTCCGCAAGATTTGGCTCTAATTGATGTGAAAAGAGCGATTAA
TATGTATCAAAAGATGAAGGTGCCTATCATTGGTATAATCGAAAATATGAGTTATTTCGT
AACTTCTGATACCGGGAAAAGATATGATTTATTTGGAAATGGTGGAGTACGTGCCGAAGC
AGAAAAGATGGGTATACCTTTTCTAGAATCAATTCCTTTTGATATGGATGTTAGAATTTT
ATCAGATCTAGGAATCCCAATTATAATAGATAATCCAAATTCGGTCGTTTCAAAAATGTA
TCAAAAAATTTCGGATCGTATCCAAGAATATCTTTTATCAAAAGTATAAATATTATGAAA
AAATAAATACATATATTCACTATTATCCGTTTTGATTATTGTAAAAATTTTTTTTTTAGT
TCAAGAGGGGTCTCACTAGAGAAGTGGTTAATTTGTAGAGGCAAATGTTTGTATCATGAA
TACTTGTAGAGAAGATCTTCAAAAGAATAATAGAAATCCTTTTTATCTAATGTTTGAATC
CATAGGTGTCGTTTATGGAGATATTGGTACTAGTGTTTTGTATGCTTTTAAAGAAGCATT
AAAAACAATGAATCATACTTCTCTCGTAGGAAGGGTAGAGATTATAGGTTTGGTTTCCTT
GATGATTTGGATCCTTACTATAATTGTAACAATAAAATATGTGTTATTATTATTAAGAGC
TGATAATGGTGGTGAAGGGGGAATATTATCTCTTCTTGCTTTGTTATTGAAAAAAATTCC
CCAACATTCAACTATATTAATAGTATTAGGTTTGATTGGTTTTGCGTTATTTATTGGAGA
TACTATGGTAACTCCAGCACTTTCCGTTCTTTCTGCAGTAGAAGGTATAAGATATATAAC
GCCGGAATTGGATGGTTTTATTATTCTTATTGCGCTAGGAATATTAATTTTATTGTTTAT

```bash
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/orthology/orthofinder/All_liberibacter_phage/formatted/orthofinder52000/Results_liberibacter_phages/Orthogroup_Sequences/*
```
#### Orthovenn3
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_JNVH01/ncbi_dataset/data/GCA_000756225.1/prokka/GCA_000756225.1_ASM75622v1_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/CLsoA_JNVH01.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_NZ1/ncbi_dataset/data/GCA_000968085.1/prokka/GCA_000968085.1_ASM96808v1_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/CLsoA_NZ1.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_RSTM/ncbi_dataset/data/GCA_001414235.1/prokka/GCA_001414235.1_ASM141423v1_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/CLsoA_RSTM.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_HenneA/ncbi_dataset/data/GCA_000968075.1/prokka/GCA_000968075.1_ASM96807v1_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/CLsoB_HenneA.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/prokka/GCA_000183665.1_ASM18366v1_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/CLsoB_ZC1.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/prokka/GCA_001983655.1_ASM198365v1_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/CLsoC_FIN111.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/prokka/GCA_001983675.1_ASM198367v1_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/CLsoC_FIN114.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/CLsoC_JIC.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoD_ISR100/ncbi_dataset/data/GCA_002918245.2/prokka/GCA_002918245.2_ASM291824v2_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/CLsoD_ISR100.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLas_A4/ncbi_dataset/data/GCA_000590865.3/prokka/GCA_000590865.3_ASM59086v3_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/Clas_A4.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLcr_BT-0/ncbi_dataset/data/GCA_001543305.1/prokka/GCA_001543305.1_ASM154330v1_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/Clcr_BT-0.faa

cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLam_SaoPaulo/ncbi_dataset/data/GCA_000496595.1/prokka/GCA_000496595.1_ASM49659v1_genomic_fixstart.faa /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/orthovenn3/Clam_SaoPaulo.faa
```
#### Circos
```bash
Assembly=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart.fasta
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos
mkdir $OutDir

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/fasta-to-karyotype.py $Assembly > $OutDir/karyotype.txt

samtools depth -a /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/minimap2/assembly_v1_corrected_fixstart_sorted.bam > ${OutDir}/coverage.txt
samtools faidx $Assembly
cut -f1,2 ${Assembly}.fai > ${OutDir}/genome.txt
bedtools makewindows -g ${Assembly}.fai -w 2000 -s 2000 > $OutDir/genome.windows
bedtools multicov -bams /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/minimap2/assembly_v1_corrected_fixstart_sorted.bam -bed $OutDir/genome.windows > $OutDir/genome.cov.histogram

awk '{print $0 "\t0.5"}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/nucmer/scaff18-1.bed | awk '{$1=$1; print}' OFS="\t" > $OutDir/scaff18-1.bed

cp /hpc-home/did23faz/git_repos/temp/bands.conf $OutDir/.
cp /hpc-home/did23faz/git_repos/temp/ideogram.conf $OutDir/.
cp /hpc-home/did23faz/git_repos/temp/ideogram.label.conf $OutDir/.
cp /hpc-home/did23faz/git_repos/temp/ideogram.position.conf $OutDir/.
cp /hpc-home/did23faz/git_repos/temp/ticks.conf $OutDir/.

source package 22619b3e-43a5-4546-ab14-4561f701f247
cd $OutDir
circos -conf /hpc-home/did23faz/git_repos/temp/liberibacter_circos.conf
```
```bash
# Define a parameter that will be used in <rules> block to
# globally toggle rules on/off. This can be convenient if you
# have many <rules> blocks and want to turn them all off using
# one parameter.
use_rules = yes

<<include etc/colors_fonts_patterns.conf>>

<<include ideogram.conf>>
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
</image>

karyotype = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/karyotype.txt

<ideogram>

show_label = no
label_font = default
label_radius = 1r + 75p
label_size = 50
label_parallel = yes

<spacing>
default = 0.005r
</spacing>

radius    = 0.9r
thickness = 0p
fill      = no

</ideogram>

<plots>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/genome.cov.histogram
extend_bin = yes
thickness = 1
r0 = 0.70r
r1 = 0.88r
orientation = out
min = 0
max = 1000
fill = yes
fill_color = blue
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/scaff18-1.bed
extend_bin = no
thickness = 1
r0 = 0.89r
r1 = 1.04r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = red
fill_transparency = 0.5
</plot>

</plots>

chromosomes_units   = 100000
chromosomes_display_default = yes


################################################################
# The remaining content is standard and required. It is imported 
# from default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files, 
# look in etc/ in the Circos distribution.

<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>
```
```bash
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos

for file in $(ls /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983655.1/*.bed /jic/research-groups/Saskia-Hogenhout/TCHeaven/Genomes/Candidatus/Liberibacter_solanacearum/GCA_001983675.1/*.bed); do 
awk '{print $0 "\t0.5"}' $file | awk '{$1=$1; print}' OFS="\t" > $OutDir/$(basename $file)
done

source package 22619b3e-43a5-4546-ab14-4561f701f247
cd $OutDir
circos -conf /hpc-home/did23faz/git_repos/temp/liberibacter_circos2.conf
circos -conf /hpc-home/did23faz/git_repos/temp/liberibacter_circos3.conf
```
```bash
# Define a parameter that will be used in <rules> block to
# globally toggle rules on/off. This can be convenient if you
# have many <rules> blocks and want to turn them all off using
# one parameter.
use_rules = yes

<<include etc/colors_fonts_patterns.conf>>

<<include ideogram.conf>>
<<include ticks.conf>>

<image>
file = circos2
<<include etc/image.conf>>
</image>

karyotype = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/karyotype.txt

<ideogram>

show_label = no
label_font = default
label_radius = 1r + 75p
label_size = 50
label_parallel = yes

<spacing>
default = 0.005r
</spacing>

radius    = 0.9r
thickness = 0p
fill      = no

</ideogram>

<plots>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/genome.cov.histogram
extend_bin = yes
thickness = 1
r0 = 0.70r
r1 = 0.88r
orientation = out
min = 0
max = 1000
fill = yes
fill_color = blue
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LWEB01000001.1.bed
extend_bin = no
thickness = 1
r0 = 0.89r
r1 = 1.04r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = red
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LWEB01000002.1.bed
extend_bin = no
thickness = 1
r0 = 0.89r
r1 = 1.04r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = green
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LWEB01000003.1.bed
extend_bin = no
thickness = 1
r0 = 0.89r
r1 = 1.04r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = yellow
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LWEB01000004.1.bed
extend_bin = no
thickness = 1
r0 = 0.89r
r1 = 1.04r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = purple
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LWEB01000005.1.bed
extend_bin = no
thickness = 1
r0 = 0.89r
r1 = 1.04r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = black
fill_transparency = 0.5
</plot>

</plots>

chromosomes_units   = 100000
chromosomes_display_default = yes


################################################################
# The remaining content is standard and required. It is imported 
# from default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files, 
# look in etc/ in the Circos distribution.

<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>
```
```bash
# Define a parameter that will be used in <rules> block to
# globally toggle rules on/off. This can be convenient if you
# have many <rules> blocks and want to turn them all off using
# one parameter.
use_rules = yes

<<include etc/colors_fonts_patterns.conf>>

<<include ideogram.conf>>
<<include ticks.conf>>

<image>
file = circos3
<<include etc/image.conf>>
</image>

karyotype = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/karyotype.txt

<ideogram>

show_label = no
label_font = default
label_radius = 1r + 75p
label_size = 50
label_parallel = yes

<spacing>
default = 0.005r
</spacing>

radius    = 0.9r
thickness = 0p
fill      = no

</ideogram>

<plots>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/genome.cov.histogram
extend_bin = yes
thickness = 1
r0 = 0.70r
r1 = 0.88r
orientation = out
min = 0
max = 1000
fill = yes
fill_color = blue
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000001.1.bed
extend_bin = no
thickness = 1
r0 = 0.89r
r1 = 1.04r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = red
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000002.1.bed
extend_bin = no
thickness = 1
r0 = 0.95r
r1 = 1.04r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = green
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000003.1.bed
extend_bin = no
thickness = 1
r0 = 0.90r
r1 = 1.0r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = yellow
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000004.1.bed
extend_bin = no
thickness = 1
r0 = 0.85r
r1 = 0.95r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = purple
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000005.1.bed
extend_bin = no
thickness = 1
r0 = 0.8r
r1 = 0.9r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = black
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000006.1.bed
extend_bin = no
thickness = 1
r0 = 0.75r
r1 = 0.85r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = cyan
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000007.1.bed
extend_bin = no
thickness = 1
r0 = 0.7r
r1 = 0.8r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = magneta
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000008.1.bed
extend_bin = no
thickness = 1
r0 = 0.65r
r1 = 0.75r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = pink
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000009.1.bed
extend_bin = no
thickness = 1
r0 = 0.6r
r1 = 0.7r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = white
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000010.1.bed
extend_bin = no
thickness = 1
r0 = 0.55r
r1 = 0.65r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = gray
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000011.1.bed
extend_bin = no
thickness = 1
r0 = 0.5r
r1 = 0.6r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = brown
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000012.1.bed
extend_bin = no
thickness = 1
r0 = 0.45r
r1 = 0.55r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = blue_5
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000013.1.bed
extend_bin = no
thickness = 1
r0 = 0.4r
r1 = 0.5r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = red_5
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000014.1.bed
extend_bin = no
thickness = 1
r0 = 0.35r
r1 = 0.45r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = green_5
fill_transparency = 0.5
</plot>

<plot>
type = histogram
file = /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/circos/LVWB01000015.1.bed
extend_bin = no
thickness = 1
r0 = 0.3r
r1 = 0.4r
orientation = out
min = 0
max = 1.4
fill = yes
fill_color = yellow_5
fill_transparency = 0.5
</plot>
</plots>

chromosomes_units   = 100000
chromosomes_display_default = yes


################################################################
# The remaining content is standard and required. It is imported 
# from default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files, 
# look in etc/ in the Circos distribution.

<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>
```
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/prokka/*fna); do
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $Genome)/BUSCO
    mkdir -p $OutDir 
    Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/bacteria_odb10
    OutFile=$(echo $Genome | cut -d '/' -f10)_$(echo $Database | cut -d '/' -f7)
    sbatch $ProgDir/run_busco_keep.sh $Genome $Database $OutDir $OutFile 
done
#6253820-885

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/prokka/BUSCO/bacteria_odb10/run_bacteria_odb10/full_table.tsv); do
grep -v "^#" $file | awk '$2=="Complete" {print $1}' >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/bacteria_complete_busco_ids.txt;
done

sort /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/bacteria_complete_busco_ids.txt |uniq -c > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/bacteria_complete_busco_ids_with_counts.txt
grep -v " 2 " /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/bacteria_complete_busco_ids_with_counts.txt | grep -v " 1 " > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/bacteria_complete_busco_ids_3.txt
awk '{print $2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/bacteria_complete_busco_ids_3.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/bacteria_complete_busco_ids.txt


mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/prokka/BUSCO/bacteria_odb10/run_bacteria_odb10/busco_sequences/single_copy_busco_sequences.tar.gz); do
cd $(dirname $file)
tar -xzvf single_copy_busco_sequences.tar.gz
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae
done

#Give unique names to the complete busco genes from each assembly:
for dir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/*/*/data/*/prokka/BUSCO/bacteria_odb10/run_bacteria_odb10/busco_sequences/single_copy_busco_sequences); do
  sppname=$(echo $dir |cut -f10 -d "/" | sed 's@/@_@g')
  abbrv=$(echo $dir |cut -f10 -d "/" | sed 's@/@_@g')
  echo $sppname
  echo $abbrv
  for file in ${dir}/*.fna; do
    out=$(echo $file |rev |cut -f 1 -d "/"|rev)
    cp $file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/${sppname}_${out}
    sed -i 's/^>/>'${abbrv}'|/g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/${sppname}_${out}
  cut -f 1 -d ":" /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/${sppname}_${out} | tr '[:lower:]' '[:upper:]' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/${sppname}_${out}.1 && mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/${sppname}_${out}.1 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/${sppname}_${out}  
  done
done

#Combine genes from each assembly into a single file per gene:
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt
buscos=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/bacteria_complete_busco_ids.txt
lines=$(cat $buscos)
for line in $lines; do
  for fna in $(ls *_$line.fna); do
  output=$(echo $line)_nt.fasta
  cat $fna >> $output
  done
done
rm *.fna

#Align the gene sequences;
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt
for file in $(ls ${AlignDir}/*_nt.fasta); do
OutFile=$(basename $file | sed 's@_nt.fasta@_nt_aligned.fasta@g')
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt
ProgDir=~/git_repos/Wrappers/NBI
echo "$file" >> mafft_log.txt
sbatch $ProgDir/sub_mafft_alignment.sh $file $OutDir $OutFile 2>&1 >> mafft_log.txt
done

for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/*_aligned.fasta); do
  ID=$(basename $gene |sed 's@_nt_aligned.fasta@@g')
  echo $ID
  mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/$ID
  cp $gene /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/$ID
done

#Trim the alignments:
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/*/*_aligned.fasta); do
  OutDir=$(dirname $Alignment)
  TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
  echo $Alignment
  echo $OutDir
  echo $TrimmedName
  sed -i '/^>/! s/[yrmwn]/N/g' "$Alignment" #trimal does not like all IUPAC symbols, just atgcn
  singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/trimal1.4.1.sif trimal -in $Alignment -out $OutDir/$TrimmedName -keepheader -automated1
done

#Trim header names
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/*/*_aligned_trimmed.fasta); do
New=$(dirname $Alignment)/$(basename $Alignment .fasta)_edit.fasta
cat $Alignment  | cut -f1 -d '|'  > $New
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/*/*_edit.fasta); do
while IFS= read -r line; do
    if [[ "$line" =~ ^\>.+ ]]; then
        echo "$line" | wc -c
    fi
done < $file
done

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir
for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/*/*_edit.fasta); do
ln -s $gene /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir/.
done

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/iqtree2
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/iqtree2
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir
cpu=4
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $AlignDir -m MFA -T AUTO --threads-max $cpu
#6256372, 6274514

#Akaike Information Criterion:           GTR+F+R5
#Corrected Akaike Information Criterion: GTR+F+R5
#Bayesian Information Criterion:         GTR+F+I+G4
#Best-fit model: GTR+F+I+G4 chosen according to BIC

AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir
cpu=4
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $AlignDir -m GTR+F+I+G4 -B 1000 -T AUTO --threads-max $cpu -redo
#6289787

#Analysis results written to:
#  IQ-TREE report:                /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir.iqtree
#  Maximum-likelihood tree:       /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir.treefile
#  Likelihood distances:          /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir.mldist

#Ultrafast bootstrap approximation results written to:
#  Split support values:          /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir.splits.nex
#  Consensus tree:                /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir.contree
#  Screen log file:               /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/phylogeny/bac_busco_nt/AlignDir.log

```
```bash
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_JNVH01/ncbi_dataset/data/GCA_000756225.1/GCA_000756225.1_ASM75622v1_genomic_fixstart.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_NZ1/ncbi_dataset/data/GCA_000968085.1/GCA_000968085.1_ASM96808v1_genomic_fixstart.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoA_RSTM/ncbi_dataset/data/GCA_001414235.1/GCA_001414235.1_ASM141423v1_genomic_fixstart.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_HenneA/ncbi_dataset/data/GCA_000968075.1/GCA_000968075.1_ASM96807v1_genomic_fixstart.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoB_ZC1/ncbi_dataset/data/GCA_000183665.1/GCA_000183665.1_ASM18366v1_genomic_fixstart.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN111/ncbi_dataset/data/GCA_001983655.1/GCA_001983655.1_ASM198365v1_genomic_fixstart.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_FIN114/ncbi_dataset/data/GCA_001983675.1/GCA_001983675.1_ASM198367v1_genomic_fixstart.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/assembly_v2.fasta /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoD_ISR100/ncbi_dataset/data/GCA_002918245.2/GCA_002918245.2_ASM291824v2_genomic_fixstart.fasta > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso.fasta

makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso.fasta -input_type fasta -dbtype nucl -title CLso -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso

blastn -query 16s+23s.fasta -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso -out CLso_results_16s+23s.out -evalue 1e-21 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -num_threads 1 -max_target_seqs 99
blastn -query 50s.fasta -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso -out CLso_results_50s.out -evalue 1e-21 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -num_threads 1 -max_target_seqs 99

awk '{print ">" $2 "\n" $NF}' CLso_results_16s+23s.out > CLso_results_16s+23s.fasta
awk '{print ">" $2 "\n" $NF}' CLso_results_50s.out > CLso_results_50s.fasta
```
CKC_RS06180 WP_013461446.1
CKC_RS02790 WP_013461982.1
CKC_RS02920 WP_013462002.1 
CKC_RS04235 WP_013462290.1
CKC_RS04955 WP_013462425.1
CKC_RS05515 WP_013462539.1
  
CKC_04420 ADR52633.1

DJ66_RS03635  WP_013461553.1
DJ66_RS00115  WP_034441411.1
DJ66_RS00995  WP_034441527.1
DJ66_RS01365  WP_034441657.1
DJ66_RS01465  WP_034442013.1
DJ66_RS03455  WP_034442419.1
DJ66_RS00825  WP_034442854.1
DJ66_RS04250  WP_034442983.1
DJ66_RS05455  WP_034443047.1
DJ66_RS00075  WP_045960327.1
DJ66_RS00565  WP_045960390.1
DJ66_RS00625  WP_045960409.1
DJ66_RS00645  WP_045960415.1
DJ66_RS00660  WP_045960420.1
DJ66_RS01460  WP_045960579.1
DJ66_RS01570  WP_045960600.1
DJ66_RS01575  WP_045960602.1
DJ66_RS02195  WP_045960718.1
DJ66_RS04460  WP_045960882.1
DJ66_RS05290  WP_045960949.1
DJ66_RS05295  WP_045960950.1
DJ66_RS05330  WP_045960955.1
DJ66_RS05400  WP_045960971.1
DJ66_RS05405  WP_045960972.1
DJ66_RS05900  WP_045960993.1
DJ66_RS04795  WP_052691167.1
DJ66_RS01090  WP_161802672.1
DJ66_RS04950  WP_244463173.1
DJ66_RS06360  WP_244463216.1
DJ66_RS05360  WP_244463237.1
DJ66_RS00260  
DJ66_RS05280  
DJ66_RS05415 
```bash
makeblastdb -in path_genes_Clso.fasta -input_type fasta -dbtype prot -title CLso_eff -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso_eff

blastx -query path_genes.fasta -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso_eff -out saberi_check.out -evalue 1e-5 -outfmt 6 -num_threads 1 -max_target_seqs 99

ls path_genes_Clso.faa
ls path_genes_Clso.fasta

makeblastdb -in path_genes_Clso.faa -input_type fasta -dbtype prot -title CLso_effectors -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso_effectors
blastp -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/genomes/CLsoC_JIC/1/data/JIC1/prokka/assembly_v2.faa -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso_effectors -out blastp_effectors_Clso.out -evalue 1e-5 -outfmt 6 -num_threads 1 -max_target_seqs 99 

makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/inspector/assembly_v1_corrected_fixstart.fasta -input_type fasta -dbtype nucl -title CLso_JIC -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso_JIC
blastn -query path_genes_Clso.fasta -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/Liberibacter/CLso_JIC -out blastn_effectors_Clso.out -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -num_threads 1 -max_target_seqs 99

```

