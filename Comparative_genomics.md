
### Phylogeny
```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny

pwd
#/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes

for file in $(ls */*/*/*.fna.gz); do
dir=$(dirname $file)
cd $dir
gunzip *.fna.gz
cd ../../..
done

#BUSCO - keeping output files
for Genome in $(ls */*/*/*.fna | grep -v 'cds_from_genomic'); do
OutDir=$(dirname $Genome)/BUSCO
Database=/jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5/hemiptera_odb10
OutFile=$(basename $Genome | sed 's@.fna@@g')_$(echo $Database | cut -d '/' -f7)
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir 
if [ ! -e "$OutDir/hemiptera_odb10/run_hemiptera_odb10/busco_sequences/single_copy_busco_sequences.tar.gz" ]; then
sbatch $ProgDir/run_busco_keep.sh "$Genome" "$Database" "$OutDir" "$OutFile" 
echo "Not done for $Genome - Failed"
else
echo "done for $Genome"
fi
done
#58077709-58072790, 58089658, 58089689, 58089850, 58736464 - 58736698 (medium), 58746250

#Extract complete busco IDs, keep those present in at least 3 genomes:
for file in $(ls */*/*/BUSCO/hemiptera_odb10/run_hemiptera_odb10/full_table.tsv); do
grep -v "^#" $file | awk '$2=="Complete" {print $1}' >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/complete_busco_ids.txt;
done

sort /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/complete_busco_ids.txt |uniq -c > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/complete_busco_ids_with_counts.txt
grep -v " 2 " /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/complete_busco_ids_with_counts.txt | grep -v " 1 " > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/complete_busco_ids.txt
awk '{print $2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/complete_busco_ids.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/complete_busco_ids_3.txt

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt

for file in $(ls */*/*/BUSCO/hemiptera_odb10/run_hemiptera_odb10/busco_sequences/single_copy_busco_sequences.tar.gz); do
cd $(dirname $file)
tar -xzvf single_copy_busco_sequences.tar.gz
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes
done

#Give unique names to the complete busco genes from each assembly:
for dir in $(ls -d */*/*/BUSCO/hemiptera_odb10/run_hemiptera_odb10/busco_sequences/single_copy_busco_sequences); do
  sppname=$(echo $dir |cut -f1,2,3 -d "/" | sed 's@/@_@g');
  abbrv=$(echo $dir | cut -d '/' -f1 | cut -c 1-3)_$(echo $dir | cut -d '/' -f2 | cut -c 1-3)_$(echo $dir | cut -d '/' -f3)
  echo $sppname
  echo $abbrv
  for file in ${dir}/*.fna; do
    out=$(echo $file |rev |cut -f 1 -d "/"|rev)
    cp $file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/${sppname}_${out}
    sed -i 's/^>/>'${abbrv}'|/g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/${sppname}_${out}
  cut -f 1 -d ":" /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/${sppname}_${out} | tr '[:lower:]' '[:upper:]' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/${sppname}_${out}.1 && mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/${sppname}_${out}.1 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/${sppname}_${out}  
  done
done

#Combine genes from each assembly into a single file per gene:
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt
buscos=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/complete_busco_ids.txt
lines=$(cat $buscos)
for line in $lines; do
  for fna in $(ls *_$line.fna); do
  output=$(echo $line)_nt.fasta
  cat $fna >> $output
  done
done
rm *.fna

#Align the gene sequences;
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt
for file in $(ls ${AlignDir}/*_nt.fasta); do
OutFile=$(basename $file | sed 's@_nt.fasta@_nt_aligned.fasta@g')
Jobs=$(squeue -u did23faz| grep 'mafft'  | wc -l)
while [ $Jobs -gt 100 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'mafft'| wc -l)
done
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt
ProgDir=~/git_repos/Wrappers/NBI
echo "$file" >> mafft_log.txt
sbatch $ProgDir/sub_mafft_alignment.sh $file $OutDir $OutFile 2>&1 >> mafft_log.txt
done
#58797311

for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/*_aligned.fasta); do
  ID=$(basename $gene |sed 's@_nt_aligned.fasta@@g')
  echo $ID
  mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/$ID
  cp $gene /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/$ID
done

#Trim the alignments:
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/*/*_aligned.fasta); do
  OutDir=$(dirname $Alignment)
  TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
  echo $Alignment
  echo $OutDir
  echo $TrimmedName
  singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/trimal1.4.1.sif trimal -in $Alignment -out $OutDir/$TrimmedName -keepheader -automated1
done

#Diuraphis_noxia biotype_2_v1 has numerous ambiguous nucleotides...

#Trim header names as RAxML need <60 characters in length:
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/*/*_aligned_trimmed.fasta); do
New=$(dirname $Alignment)/$(basename $Alignment .fasta)_edit.fasta
cat $Alignment  | cut -f1 -d '|'  > $New
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/*/*_edit.fasta); do
while IFS= read -r line; do
    if [[ "$line" =~ ^\>.+ ]]; then
        echo "$line" | wc -c
    fi
done < $file
done

Jobs=$(squeue -u did23faz| grep 'RAxML'  | wc -l)
while [ $Jobs -gt 100 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'RAxML'| wc -l)
done
#Run RAxML
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/*/*_edit.fasta); do
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/RAxML/$Prefix
ProgDir=~/git_repos/Wrappers/NBI
mkdir -p $OutDir
Jobs=$(squeue -u did23faz| grep 'RAxML'  | wc -l)
while [ $Jobs -gt 190 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'RAxML'| wc -l)
done
sbatch $ProgDir/run_RAxML_msa.sh $Alignment $Prefix $OutDir 
done

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/RAxML
count=0
for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/*/*_edit.fasta); do
ID=$(echo $gene | cut -d '/' -f11)
echo $ID 
x=$(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/RAxML/${ID}.log)
            if [[ -f ${x} ]]; then
               ((count++))
            fi
#mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/RAxML/${ID}* 
done
echo "$count"
```

```python
with open("hemip.faa", "r") as f:
    lines = f.readlines()

with open("output.fasta", "w") as f:
    for line in lines:
        if line.startswith(">"):
            line = line.replace("|", "_").replace(" ", "_").replace(".", "_")
        f.write(line)
```
```bash
singularity exec ~/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --model-filepath ../databases/helixer/invertebrate_v0.3_a_0600/invertebrate_v0.3_a_0600.h5 --subsequence-length 213840 --overlap-offset 106920 --overlap-core-length 160380 --fasta-path Arabidopsis_lyrata.v.1.0.dna.chromosome.8.fa  \
  --species Arabidopsis_lyrata --gff-output-path Arabidopsis_lyrata_chromosome8_helixer.gff3

```
## Synteny
#### Helixer
```bash
mkdir -p analysis/synteny/helixer

for fasta in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Diaphorina/citri/GCA_030643865.1/GCA_030643865.1_ASM3064386v1_genomic.fna); do
model_filepath=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/helixer/invertebrate_v0.3_m_0100/invertebrate_v0.3_m_0100.h5
lineage=invertebrate
species=$(echo $fasta | cut -d '/' -f8 | cut -c 1)_$(echo $fasta | cut -d '/' -f9)
outfile=$species
outdir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/${species}
ProgDir=~/git_repos/Wrappers/NBI
mkdir $outdir
sbatch $ProgDir/run_helixer.sh $fasta $model_filepath $lineage $species $outfile $outdir
done
#58856194-6, 58857547, 58885166, 58902576, 58967085

grep 'gene' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.gff | wc -l #21,036
grep 'gene' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.gff | wc -l #18,217
grep 'gene' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.gff | wc -l #34,557
```
#### AGAT
```bash
source package 4c883633-af2d-4fac-ab67-a1574f7fe079

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Diaphorina/citri/GCA_000475195.1/GCA_000475195.1_Diaci_psyllid_genome_assembly_version_1.1_genomic.fna
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.gff
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Bactericera/cockerelli/GCA_024516035.1/GCA_024516035.1_ASM2451603v1_genomic.fna
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.gff
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Pachypsylla/venusta/GCA_012654025.1/GCA_012654025.1_Pven_dovetail_genomic.fna
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.gff
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.gff
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.gff
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_HiFiPurged_curated_break_scaffolds_final_nomito_filtered_corrected_softmasked.fa
Gff=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.gff
OutFile=$(dirname $Gff)/$(basename $Gff | sed 's@.gff@.faa@g')
agat_sp_extract_sequences.pl -g $Gff -f $Genome -t cds --output $OutFile --clean_final_stop --protein
```
#### Genespace
```bash
source package 03380c15-2730-4b19-b17a-5a435e152681

mkdir -p analysis/synteny/genespace/raw_genomes/D_citri
mkdir -p analysis/synteny/genespace/raw_genomes/T_apicales
mkdir -p analysis/synteny/genespace/raw_genomes/T_anthrisci
mkdir -p analysis/synteny/genespace/raw_genomes/T_urticae # highly fragmented
mkdir -p analysis/synteny/genespace/raw_genomes/B_cockerelli
mkdir -p analysis/synteny/genespace/raw_genomes/P_venusta

ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.gff analysis/synteny/genespace/raw_genomes/B_cockerelli/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.faa analysis/synteny/genespace/raw_genomes/B_cockerelli/.

ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.gff analysis/synteny/genespace/raw_genomes/P_venusta/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.faa analysis/synteny/genespace/raw_genomes/P_venusta/.

ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.gff analysis/synteny/genespace/raw_genomes/D_citri/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.faa analysis/synteny/genespace/raw_genomes/D_citri/.

ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.gff analysis/synteny/genespace/raw_genomes/T_anthrisci/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.faa analysis/synteny/genespace/raw_genomes/T_anthrisci/.

ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.gff analysis/synteny/genespace/raw_genomes/T_apicales/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa analysis/synteny/genespace/raw_genomes/T_apicales/.

ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.gff analysis/synteny/genespace/raw_genomes/T_urticae/.
ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.faa analysis/synteny/genespace/raw_genomes/T_urticae/.

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/genespace/raw_genomes/*/*.faa); do
Out=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/genespace/peptide/$(basename $file | sed 's@.faa@.fa@g')
cat $file | cut -d ' ' -f1 > $Out
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/genespace/raw_genomes/*/*.gff); do
Out=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/genespace/bed/$(basename $file | sed 's@.gff@.bed@g')
cat "$file" | awk '$3 == "gene"' | cut -f1,4,5,9 | awk -F'\t' -v OFS='\t' '{ $4 = $4 ".1"; gsub("ID=", "", $4); print }' > $Out
done

source package 03380c15-2730-4b19-b17a-5a435e152681
R
```
```R
library(GENESPACE)
gpar <- init_genespace(
  wd = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/genespace", 
  path2mcscanx = "/software/03380c15-2730-4b19-b17a-5a435e152681/bin/MCScanX/")
gpar <- run_genespace(gsParam = gpar)
```
>D_citri_KI478532.1_000002.1
D_citri_KI472552.1_000002.1
>D_citri_KI474893.1_000002.1