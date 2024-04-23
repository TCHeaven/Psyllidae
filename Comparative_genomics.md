### BUSCO plot
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Trioza/urticae/v1/BUSCO/hemiptera_odb10/short_summary.specific.hemiptera_odb10.hemiptera_odb10.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco/short_summary.specific.hemiptera_odb10.Tri_urt.txt
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Dyspera/apicales/v1/BUSCO/hemiptera_odb10/short_summary.specific.hemiptera_odb10.hemiptera_odb10.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco/short_summary.specific.hemiptera_odb10.Dys_api.txt
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Dyspera/pallida/v1/BUSCO/hemiptera_odb10/short_summary.specific.hemiptera_odb10.hemiptera_odb10.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco/short_summary.specific.hemiptera_odb10.Dys_pal.txt
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Diaphorina/citri/GCA_000475195.1/BUSCO/hemiptera_odb10/short_summary.specific.hemiptera_odb10.hemiptera_odb10.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco/short_summary.specific.hemiptera_odb10.Dia_cit.txt
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Bactericera/cockerelli/GCA_024516035.1/BUSCO/hemiptera_odb10/short_summary.specific.hemiptera_odb10.hemiptera_odb10.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco/short_summary.specific.hemiptera_odb10.Bac_coc.txt
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Pachypsylla/venusta/GCA_012654025.1/BUSCO/hemiptera_odb10/short_summary.specific.hemiptera_odb10.hemiptera_odb10.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco/short_summary.specific.hemiptera_odb10.Pac_ven.txt
source package 97a7f391-0f5c-4bdb-b951-14d0d5ee4576
generate_plot.py -wd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco
```
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
sbatch $ProgDir/run_RAxML_msa.sh $Alignment $OutDir $Prefix
done
#ERROR: $Prefix and $OutDir were the wrong way around resulting in evaluations being deleted, results are GTR+G

#Tested propberly for 9997at7524:
#9997at7524.raxml.log:AIC score: 88062.537364 / AICc score: 88476.378774 / BIC score: 89480.735172
#GTR+G+FO.raxml.log:AIC score: 88062.550680 / AICc score: 88476.392090 / BIC score: 89480.748487
#GTR+G.raxml.log:AIC score: 88062.550680 / AICc score: 88476.392090 / BIC score: 89480.748487
#GTR+R4+FO.raxml.log:AIC score: 88049.098703 / AICc score: 88481.312512 / BIC score: 89490.469678
#GTR.raxml.log:AIC score: 97980.036779 / AICc score: 98390.278538 / BIC score: 99393.599953
#JC+G.raxml.log:AIC score: 90954.927524 / AICc score: 91340.650467 / BIC score: 92336.048264
#JC.raxml.log:AIC score: 100426.479828 / AICc score: 100808.795163 / BIC score: 101802.965936
#Lower is better, score for GTR+G is ~ as good as any other method, therefore not worth retesting all the trees.

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/RAxML

count=0
for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/*/*_edit.fasta); do
ID=$(echo $gene | cut -d '/' -f11)
echo $ID 
x=$(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/RAxML/${ID}.log)
            if [[ -f ${x} ]]; then
               ((count++))
            fi
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/RAxML/${ID}* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/RAxML/.
done
echo "$count"


#Combine individual gene trees into a consensus tree:
source package 0351788b-6639-43fb-8e80-f28600f83cb1
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/RAxML/*.raxml.bestTree > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/tree-files.txt
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/RAxML/*.raxml.bootstraps > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/bs-files.txt
astral5 -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/tree-files.txt -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.astral.tre
astral5 -t 3 -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/tree-files.txt -b /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/bs-files.txt -r 100 -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped2.astral.tre
tail -n 1 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped2.astral.tre > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped.consensus2.astral.tre
head -n 100 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped2.astral.tre > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstraps.astral.tre


astral5 -t 3 -q /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped.consensus2.astral.tre -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstraps.astral.tre -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped.scored2.astral.tre 
#59183148,59185037,59260253, 59315880 
#no .csv was produced

astral5 -t 3 -i /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/tree-files.txt -q /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped.consensus2.astral.tre -o /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped.localposteriorprobabilityscored2.astral.tre 2> astral-scored.log

/home/theaven/scratch/apps/ASTER/bin/astral4 -i /home/theaven/scratch/apps/ASTER/tree-files.txt -C -c  /home/theaven/scratch/apps/ASTER/hemiptera_phylogeny.bootstrapped.consensus2.astral.tre -o /home/theaven/scratch/apps/ASTER/hemiptera_phylogeny.bootstrapped.consensus2.astral.castle.tre --root API_MEL_GCF_003254395.2 --genelength 1245


mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir
for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/*_aligned.fasta); do
ln -s $gene /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir/.
done

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir2
for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir/*_aligned.fasta); do
name=$(basename $gene)
cat $gene | cut -d '|' -f1 > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir2/${name}
done

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3
for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/*/*_edit.fasta); do
ln -s $gene /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3/.
done

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3/
cpu=4
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $AlignDir -m MFP -mtree -T AUTO --threads-max $cpu -b 1000 --verbose
#59196096, 59251934, 59260066

#Not sure if is running correctly, find model seperately:
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2-test
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2-test
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3/
cpu=32
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $AlignDir -m MF -T 32 --threads-max $cpu
#59316762, more cpus: 59322062 - runs out of memory after 4 days only tested 12 models

#Randomly select 100 BUSCO genes and use these for model selection, also select the GTR model as almost certainly best base model
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3/ | shuf -n 100 | xargs -I {} ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3/{} /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir-100/{}
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2-test
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2-test
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir-100/
cpu=40
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $AlignDir -m MF -T AUTO -mset GTR --threads-max $cpu
#59400442
#Best-fit model: GTR+F+I+R10 chosen according to BIC

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3/
cpu=12
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $AlignDir -m GTR+F+I+R10 -B 1000 -T AUTO --threads-max $cpu
#59425290
```
Input data: 161 sequences with 3121748 nucleotide sites
Number of constant sites: 472042 (= 15.1211% of all sites)
Number of invariant (constant or ambiguous constant) sites: 472042 (= 15.1211% of all sites)
Number of parsimony informative sites: 2416927
Number of distinct site patterns: 2748603

Analysis results written to:
  IQ-TREE report:                /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3.iqtree
  Maximum-likelihood tree:       /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3.treefile
  Likelihood distances:          /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3.mldist

Ultrafast bootstrap approximation results written to:
  Split support values:          /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3.splits.nex
  Consensus tree:                /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3.contree
  Screen log file:               /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3.log

```bash
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2/.

ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/RAxML/*.raxml.bestTree > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/tree-files2.txt
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3/*_edit.fasta > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/align-files.txt
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -t /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3.contree --gcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/tree-files.txt --prefix concord -T 1
#Tree with concordance factors written to concord.cf.tree
#Annotated tree (best viewed in FigTree) written to concord.cf.tree.nex
#Tree with branch IDs written to concord.cf.branch
#Concordance factors per branch printed to concord.cf.stat

singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -te /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3.contree -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3/ --scfl 100 --prefix concord2 -T 8
#59526546, 59645092?

#Partition information was printed to concord2.best_model.nex
#Tree with site concordance factors written to concord2.cf.tree
#Annotated tree (best viewed in FigTree) written to concord2.cf.tree.nex
#Tree with branch IDs written to concord2.cf.branch
#Concordance factors per branch printed to concord2.cf.stat

#Analysis results written to:
#  IQ-TREE report:                concord2.iqtree
#  Maximum-likelihood tree:       concord2.treefile
#  Screen log file:               concord2.log

mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/concord2* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree2/.
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/concord* /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/.


singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -t /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped.consensus2.astral.tre --gcf /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/tree-files.txt --prefix concord3 -T 8
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -te /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped.consensus2.astral.tre -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3/ --scfl 100 --prefix concord4 -T 8
#59645173
```
For every branch of a reference tree, gCF is defined as the percentage of “decisive” gene trees containing that branch. 

sCF is defined as the percentage of decisive alignment sites supporting a branch in the reference tree.

gCF and sCF complement classical measures of branch support (e.g. bootstrap) in phylogenetics by providing a full description of underlying disagreement among loci and sites.

```bash
#Run RAxML
conda activate raxml
for Alignment in $(ls /mnt/shared/scratch/theaven/uncompressed/hogenhout/phylogeny/9997at7524_nt_aligned_trimmed_edit.fasta); do
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=/mnt/shared/scratch/theaven/uncompressed/hogenhout/phylogeny/RAxML/$Prefix
ProgDir=/home/theaven/scratch/apps/phylogeny
mkdir -p $OutDir
Jobs=$(squeue -u theaven| grep 'RAxML'  | wc -l)
while [ $Jobs -gt 50 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -u theaven| grep 'RAxML'| wc -l)
done
sbatch $ProgDir/run_RAxML_msa.sh $Alignment $OutDir $Prefix
done
#19226105
conda deactivate
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
#### R8s
Making the species tree ultrametric CAFE requires a tree that is ultramatric. There are many ways to obtain ultrametric trees (also known as timetrees, these are phylogenetic trees scaled to time, where all paths from root to tips have the same length). Here, we use a fast program called r8s. You will need to know the number of sites in the alignment used to estimate the species tree (the one you want to make ultrametric), and then you can specify one or more calibration points (ideally, the age or age window of a documented fossil) to scale branch lengths into time units. We provide you with a script that prepares the control file for running r8s on the species tree above (the number of sites is 35157236, and the calibration point for cats and humans is 94). In your shell, type:

```bash
source package /nbi/software/testing/bin/r8s-1.80

echo "#NEXUS" > r8s_ctl_file3.txt
echo "begin trees; " >> r8s_ctl_file3.txt
echo "tree hemiptera_tree = [&R] $(cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/AlignDir3.contree)" >> r8s_ctl_file3.txt
echo "End;" >> r8s_ctl_file3.txt
echo "begin rates;" >> r8s_ctl_file3.txt
echo "blformat nsites=3121748 lengths=persite ultrametric=no;" >> r8s_ctl_file3.txt
echo "collapse;" >> r8s_ctl_file3.txt
#echo "reroot DRO_MEL_GCF_000001215.4;" >> r8s_ctl_file.txt
echo "reroot BOM_MOR_GCF_014905235.1;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_1 THR_PAL_GCA_012932325.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_2 BOM_MOR_GCF_014905235.1 ANO_STE_GCF_013141755.1;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_3 BOM_MOR_GCF_014905235.1 TRI_CAS_GCF_000002335.3;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_4 BOM_MOR_GCF_014905235.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_5 BAC_COC_GCA_024516035.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_6 CIM_LEC_GCA_000648675.3 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_7 DAK_VIT_INRAPCF7_V5 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_8 PSE_VIB_GCA_033439095.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_9 CIM_LEC_GCA_000648675.3 PHI_SPU_GCA_018207615.1;" >> r8s_ctl_file3.txt
echo "mrca callibration_node_10 NIL_LUG_GCA_014356525.1 PHI_SPU_GCA_018207615.1;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_1 min_age=206.1 max_age=404.6;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_2 min_age=151.9 max_age=344.7;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_3 min_age=195 max_age=361.6;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_4 min_age=330.4 max_age=376;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_5 min_age=211.9 max_age=351;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_6 min_age=112.5 max_age=391.7;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_7 min_age=87.1 max_age=162;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_8 min_age=168 max_age=286.4;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_9 min_age=234.9 max_age=366.2;" >> r8s_ctl_file3.txt
echo "constrain taxon=callibration_node_10 min_age=169.6 max_age=344.5;" >> r8s_ctl_file3.txt
echo "divtime method=pl algorithm=tn cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;" >> r8s_ctl_file3.txt
echo "describe plot=chronogram;" >> r8s_ctl_file3.txt
echo "describe plot=tree_description;" >> r8s_ctl_file3.txt
echo "end;" >> r8s_ctl_file3.txt


r8s -b -f r8s_ctl_file.txt > temp_r8s3.txt
tail -n 1 temp_r8s3.txt | cut -c 16- > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree_r8s2.txt

r8s -b -f r8s_ctl_file3.txt > temp_r8s4.txt
tail -n 1 temp_r8s4.txt | cut -c 16- > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree_r8s3.txt

echo "#NEXUS" > r8s_ctl_file22.txt
echo "begin trees; " >> r8s_ctl_file22.txt
echo "tree hemiptera_tree = [&R] $(cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/astralreroot.txt)" >> r8s_ctl_file22.txt
echo "End;" >> r8s_ctl_file22.txt
echo "begin rates;" >> r8s_ctl_file22.txt
echo "blformat nsites=3121748 lengths=persite ultrametric=no;" >> r8s_ctl_file22.txt
echo "collapse;" >> r8s_ctl_file22.txt
#echo "reroot DRO_MEL_GCF_000001215.4;" >> r8s_ctl_file2.txt
#echo "reroot API_MEL_GCF_003254395.2;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_1 THR_PAL_GCA_012932325.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_2 BOM_MOR_GCF_014905235.1 ANO_STE_GCF_013141755.1;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_3 BOM_MOR_GCF_014905235.1 TRI_CAS_GCF_000002335.3;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_4 BOM_MOR_GCF_014905235.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_5 BAC_COC_GCA_024516035.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_6 CIM_LEC_GCA_000648675.3 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_7 DAK_VIT_INRAPCF7_V5 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_8 PSE_VIB_GCA_033439095.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_9 CIM_LEC_GCA_000648675.3 PHI_SPU_GCA_018207615.1;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_10 NIL_LUG_GCA_014356525.1 PHI_SPU_GCA_018207615.1;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_1 min_age=206.1 max_age=404.6;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_2 min_age=151.9 max_age=344.7;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_3 min_age=195 max_age=361.6;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_4 min_age=330.4 max_age=376;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_5 min_age=211.9 max_age=351;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_6 min_age=112.5 max_age=391.7;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_7 min_age=87.1 max_age=162;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_8 min_age=168 max_age=286.4;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_9 min_age=234.9 max_age=366.2;" >> r8s_ctl_file22.txt
echo "constrain taxon=callibration_node_10 min_age=169.6 max_age=344.5;" >> r8s_ctl_file22.txt
echo "divtime method=pl algorithm=tn cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;" >> r8s_ctl_file22.txt
echo "describe plot=chronogram;" >> r8s_ctl_file22.txt
echo "describe plot=tree_description;" >> r8s_ctl_file22.txt
echo "end;" >> r8s_ctl_file22.txt


r8s -b -f r8s_ctl_file22.txt > temp_r8s2.txt
grep 'tree hemiptera_tree' temp_r8s2.txt | cut -c 16- > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/astral_r8s.txt
tail -n 17 temp_r8s2.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/astral_r8s2.txt
#tail -n 1 temp_r8s2.txt | cut -c 16- > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/astral_r8s.txt

r8s -b -f r8s_ctl_file22.txt > temp_r8s22.txt
tail -n 1 temp_r8s22.txt | cut -c 16- > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/astral_r8s2.txt

echo "#NEXUS" > r8s_ctl_file222.txt
echo "begin trees; " >> r8s_ctl_file222.txt
echo "tree hemiptera_tree = [&R] $(cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/iqreroot.txt)" >> r8s_ctl_file222.txt
echo "End;" >> r8s_ctl_file222.txt
echo "begin rates;" >> r8s_ctl_file222.txt
echo "blformat nsites=3121748 lengths=persite ultrametric=no;" >> r8s_ctl_file222.txt
echo "collapse;" >> r8s_ctl_file222.txt
#echo "reroot DRO_MEL_GCF_000001215.4;" >> r8s_ctl_file2.txt
#echo "reroot API_MEL_GCF_003254395.2;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_1 THR_PAL_GCA_012932325.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_2 BOM_MOR_GCF_014905235.1 ANO_STE_GCF_013141755.1;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_3 BOM_MOR_GCF_014905235.1 TRI_CAS_GCF_000002335.3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_4 BOM_MOR_GCF_014905235.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_5 BAC_COC_GCA_024516035.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_6 CIM_LEC_GCA_000648675.3 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_7 DAK_VIT_INRAPCF7_V5 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_8 PSE_VIB_GCA_033439095.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_9 CIM_LEC_GCA_000648675.3 PHI_SPU_GCA_018207615.1;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_10 NIL_LUG_GCA_014356525.1 PHI_SPU_GCA_018207615.1;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_1 min_age=206.1 max_age=404.6;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_2 min_age=151.9 max_age=344.7;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_3 min_age=195 max_age=361.6;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_4 min_age=330.4 max_age=376;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_5 min_age=211.9 max_age=351;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_6 min_age=112.5 max_age=391.7;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_7 min_age=87.1 max_age=162;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_8 min_age=168 max_age=286.4;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_9 min_age=234.9 max_age=366.2;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_10 min_age=169.6 max_age=344.5;" >> r8s_ctl_file222.txt
echo "divtime method=pl algorithm=tn cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;" >> r8s_ctl_file222.txt
echo "describe plot=chronogram;" >> r8s_ctl_file222.txt
echo "describe plot=tree_description;" >> r8s_ctl_file222.txt
echo "end;" >> r8s_ctl_file222.txt

r8s -b -f r8s_ctl_file222.txt > temp_r8s222.txt
tail -n 1 temp_r8s222.txt | cut -c 16- > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/busco_nt/iqtree_r8s222.txt

echo "#NEXUS" > r8s_ctl_file222.txt
echo "begin trees; " >> r8s_ctl_file222.txt
echo "tree hemiptera_tree = [&R] $(cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/hemiptera_phylogeny.bootstrapped.consensus2.astral.castle.tre)" >> r8s_ctl_file222.txt
echo "End;" >> r8s_ctl_file222.txt
echo "begin rates;" >> r8s_ctl_file222.txt
echo "blformat nsites=3121748 lengths=persite ultrametric=no;" >> r8s_ctl_file222.txt
echo "collapse;" >> r8s_ctl_file222.txt
#echo "reroot DRO_MEL_GCF_000001215.4;" >> r8s_ctl_file2.txt
#echo "reroot API_MEL_GCF_003254395.2;" >> r8s_ctl_file22.txt
echo "mrca callibration_node_1 THR_PAL_GCA_012932325.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_2 BOM_MOR_GCF_014905235.1 ANO_STE_GCF_013141755.1;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_3 BOM_MOR_GCF_014905235.1 TRI_CAS_GCF_000002335.3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_4 BOM_MOR_GCF_014905235.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_5 BAC_COC_GCA_024516035.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_6 CIM_LEC_GCA_000648675.3 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_7 DAK_VIT_INRAPCF7_V5 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_8 PSE_VIB_GCA_033439095.1 ACY_PIS_LSR1_AL4_V3;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_9 CIM_LEC_GCA_000648675.3 PHI_SPU_GCA_018207615.1;" >> r8s_ctl_file222.txt
echo "mrca callibration_node_10 NIL_LUG_GCA_014356525.1 PHI_SPU_GCA_018207615.1;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_1 min_age=206.1 max_age=404.6;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_2 min_age=151.9 max_age=344.7;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_3 min_age=195 max_age=361.6;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_4 min_age=330.4 max_age=376;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_5 min_age=211.9 max_age=351;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_6 min_age=112.5 max_age=391.7;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_7 min_age=87.1 max_age=162;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_8 min_age=168 max_age=286.4;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_9 min_age=234.9 max_age=366.2;" >> r8s_ctl_file222.txt
echo "constrain taxon=callibration_node_10 min_age=169.6 max_age=344.5;" >> r8s_ctl_file222.txt
echo "divtime method=pl algorithm=tn cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;" >> r8s_ctl_file222.txt
echo "describe plot=chronogram;" >> r8s_ctl_file222.txt
echo "describe plot=tree_description;" >> r8s_ctl_file222.txt
echo "end;" >> r8s_ctl_file222.txt

r8s -b -f r8s_ctl_file222.txt > temp_r8s-castle.txt

tail -n 1 temp_r8s-castle.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/phylogeny/astral_r8s-castle.txt
```
Running r8s on rerooted trees from ITOL (as the reroot command in r8s doesnt seem to have any effect), this removes the bootstrap values for some reason. R8s non-rerooted trees keep the bootstraps and rerooted trees prior to r8s have bootstrap, I can't be bothered trying to troubleshoot this so am just manually adding the bootstrap values back into the raw newick files.
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
grep 'gene' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.gff | wc -l #27,628
```
#### AGAT
```bash
source package 4c883633-af2d-4fac-ab67-a1574f7fe079

Genome=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Diaphorina/citri/GCA_030643865.1/GCA_030643865.1_ASM3064386v1_genomic.fna
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
```bash
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.gff | awk '{print $1, $9, $4, $5}' | awk 'NF{print "Dc" $0}' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/Dc.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.gff | awk '{print $1, $9, $4, $5}' | awk 'NF{print "Bc" $0}' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/Bc.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.gff | awk '{print $1, $9, $4, $5}' | awk 'NF{print "Pv" $0}' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/Pv.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.gff | awk '{print $1, $9, $4, $5}' | awk 'NF{print "Tp" $0}' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/Tp.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.gff | awk '{print $1, $9, $4, $5}' | awk 'NF{print "Tn" $0}' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/Tn.gff
awk '$3 == "gene" {split($9, parts, "="); printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, parts[2]}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.gff | awk '{print $1, $9, $4, $5}' | awk 'NF{print "Tu" $0}' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/Tu.gff

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/*.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.gff
#.gff variant names do not have .1 variant information, blast file gene names do
```


```bash
mkdir -p /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/D_citri -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/B_cockerelli -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/P_venusta -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_apicales -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_anthrisci -dbtype prot
makeblastdb -in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.faa -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_urticae -dbtype prot

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/D_citri -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/B_cockerelli_v_D_citri.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/D_citri -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/P_venusta_v_D_citri.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/D_citri -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_D_citri.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/D_citri -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_D_citri.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/D_citri -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_urticae_v_D_citri.blast
#59190238

blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/B_cockerelli -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/P_venusta_v_B_cockerelli.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/B_cockerelli -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_B_cockerelli.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/B_cockerelli -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_B_cockerelli.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/B_cockerelli -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_urticae_v_B_cockerelli.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/B_cockerelli -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/D_citri_v_B_cockerelli.blast
#59190267

blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/P_venusta -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/B_cockerelli_v_P_venusta.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/P_venusta -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_P_venusta.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/P_venusta -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_P_venusta.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/P_venusta -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_urticae_v_P_venusta.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/P_venusta -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/D_citri_v_P_venusta.blast
#59190270

blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_apicales -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/B_cockerelli_v_T_apicales.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_apicales -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/P_venusta_v_T_apicales.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_apicales -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_T_apicales.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_apicales -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_urticae_v_T_apicales.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_apicales -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/D_citri_v_T_apicales.blast
#59190272

blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_anthrisci -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/B_cockerelli_v_T_anthrisci.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_anthrisci -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/P_venusta_v_T_anthrisci.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_anthrisci -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_T_anthrisci.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_anthrisci -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_urticae_v_T_anthrisci.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_anthrisci -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/D_citri_v_T_anthrisci.blast
#59190274

blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_urticae -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/B_cockerelli_v_T_urticae.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_urticae -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/P_venusta_v_T_urticae.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_urticae -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_T_urticae.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_urticae -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_T_urticae.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_urticae -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/D_citri_v_T_urticae.blast
#59190275

blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_urticae -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/helixer/T_urticae.faa -num_threads 32 -evalue 1e-10 -num_alignments 6 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_urticae_v_T_urticae.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_anthrisci -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_anthrisci.faa -num_threads 32 -evalue 1e-10 -num_alignments 6 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_T_anthrisci.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/T_apicales -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/helixer/T_apicales.faa -num_threads 32 -evalue 1e-10 -num_alignments 6 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_T_apicales.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/P_venusta -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/P_venusta/P_venusta.faa -num_threads 32 -evalue 1e-10 -num_alignments 6 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/P_venusta_v_P_venusta.blast
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/B_cockerelli -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/B_cockerelli/B_cockerelli.faa -num_threads 32 -evalue 1e-10 -num_alignments 6 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/B_cockerelli_v_B_cockerelli.blast 
blastp -db /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/ncbiDB/D_citri -query /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.faa -num_threads 32 -evalue 1e-10 -num_alignments 6 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/helixer/D_citri/D_citri.faa -num_threads 32 -evalue 1e-10 -num_alignments 5 -outfmt 6 -out /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/D_citri_v_D_citri.blast
#59196288

cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/*.blast > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.blast
awk 'BEGIN {FS=OFS="\t"} {sub(/\.1$/, "", $1); sub(/\.1$/, "", $2); print}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.blast > temp.blast && mv temp.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.blast
< /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.blast tr ' ' '\t' > temp.blast && mv temp.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.blast
< /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.gff tr ' ' '\t' > temp.gff && mv temp.gff /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.gff

source package 038f5eb6-dc79-46b5-bb52-a86ed67aa64a
MCScanX /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master

#1180980 matches imported (1343327 discarded)
#188870 pairwise comparisons
#

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_T_anthrisci.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_T_apicales.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_T_anthrisci.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_T_apicales.blast > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant/master.blast
cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/Tn.gff /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/Tp.gff > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant/master.gff
awk 'BEGIN {FS=OFS="\t"} {sub(/\.1$/, "", $1); sub(/\.1$/, "", $2); print}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant/master.blast > temp.blast && mv temp.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant/master.blast
< /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant/master.blast tr ' ' '\t' > temp.blast && mv temp.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant/master.blast
< /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant/master.gff tr ' ' '\t' > temp.gff && mv temp.gff /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant/master.gff
MCScanX /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/tapitant/master
```
```python
from Bio import SeqIO

# Path to your multi-FASTA file
fasta_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Diaphorina/citri/GCA_030643865.1/GCA_030643865.1_ASM3064386v1_genomic.fna"
fasta_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Bactericera/cockerelli/GCA_024516035.1/GCA_024516035.1_ASM2451603v1_genomic.fna"
fasta_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes/Pachypsylla/venusta/GCA_012654025.1/GCA_012654025.1_Pven_dovetail_genomic.fna"
fasta_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_apicales/hifiasm_19.5/880m/29/3/3.0/0.75/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_apicales_880m_29_3_3.0_0.75_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa"
fasta_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_anthrisci/hifiasm_19.5/820m/48/1/10.0/0.25/break10x/purge_dups/sanger/MitoHifi/filtered/inspector/repeatmasker/softmask/T_anthrisci_820m_48_1_10.0_0.25_break_TellSeqPurged_curated_nomito_filtered_corrected_softmasked.fa"
fasta_file = "/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/filtered/purge_dups/purge_haplotigs/break10x/yahs/filtered/inspector/repeatmasker/softmask/T_urticae_715m_12_2_3.0_0.5_filtered_HiFiPurged_HiFiPurged_curated_break_scaffolds_final_nomito_filtered_corrected_softmasked.fa"

# Iterate over each sequence in the multi-FASTA file
for record in SeqIO.parse(fasta_file, "fasta"):
    # Print the ID and length of each sequence
    print(f"Contig ID: {record.id}, Length: {len(record.seq)}")
```
```bash

awk -F'\t' '{print $1}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.gff | grep '^Tp' | sort -u


cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.blast /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.blast
cp /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master.gff /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff

sed -i 's/TnSUPER_1/Tn1/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_2/Tn2/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_3/Tn3/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_4/Tn4/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_5/Tn5/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_6/Tn6/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_7/Tn7/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_8/Tn8/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_9/Tn9/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_10/Tn10/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_11/Tn11/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_12/Tn12/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TnSUPER_13/Tn13/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff

sed -i 's/TpSUPER_1/Tp1/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_2/Tp2/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_3/Tp3/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_4/Tp4/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_5/Tp5/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_6/Tp6/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_7/Tp7/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_8/Tp8/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_9/Tp9/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_10/Tp10/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's@TpSUPER_11_1@Tp11@' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_12/Tp12/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/TpSUPER_13/Tp13/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/Tp11_1/Tp11/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff

sed -i 's/BcCM044955.1/Bc1/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044956.1/Bc2/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044957.1/Bc3/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044958.1/Bc4/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044959.1/Bc5/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044960.1/Bc6/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044961.1/Bc7/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044962.1/Bc8/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044963.1/Bc9/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044964.1/Bc10/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044965.1/Bc11/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044966.1/Bc12/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/BcCM044967.1/Bc13/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff

sed -i 's/DcCP125321.1/Dc9/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125322.1/Dc2/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125323.1/Dc12/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125324.1/Dc7/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125325.1/Dc5/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125326.1/Dc13/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125327.1/Dc11/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125328.1/Dc8/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125329.1/Dc6/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125330.1/Dc10/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125331.1/Dc4/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125332.1/Dc3/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/DcCP125333.1/Dc1/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff

sed -i 's/PvCM022874.1/Pv1/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022875.1/Pv2/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022876.1/Pv3/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022877.1/Pv4/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022878.1/Pv5/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022879.1/Pv6/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022880.1/Pv7/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022881.1/Pv8/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022882.1/Pv9/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022883.1/Pv10/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022884.1/Pv11/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvCM022885.1/Pv12/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff
sed -i 's/PvJAAEEI010039819.1/Pv13/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited.gff

source package 038f5eb6-dc79-46b5-bb52-a86ed67aa64a
MCScanX /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/master/master-edited
```

```bash
awk '$1 ~ /^T_apicales_SUPER_4/ && $2 ~ /^T_anthrisci_SUPER_9/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_T_anthrisci.blast
awk '$1 ~ /^T_apicales_SUPER_9/ && $2 ~ /^T_anthrisci_SUPER_4/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_apicales_v_T_anthrisci.blast

awk '$1 ~ /^T_anthrisci_SUPER_4/ && $2 ~ /^T_apicales_SUPER_9/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_T_apicales.blast
awk '$1 ~ /^T_anthrisci_SUPER_9/ && $2 ~ /^T_apicales_SUPER_4/' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/analysis/synteny/mcscanx/intermediateData/T_anthrisci_v_T_apicales.blast
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

source package 038f5eb6-dc79-46b5-bb52-a86ed67aa64a
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

#### MEGAN
We modified a previous pipeline (Nikoh et al., 2010) to detect genes of bacterial origin expressed in our RNA-seq data. First, the transcriptome assemblies from both RNA-seq samples were searched by BLASTX v2.2.25+ (-evalue 1 × 10−3 -outfmt 7) against the nr database (posted January 2012) and the blast results were visualized in Megan v4 (Huson et al., 2011) as metatranscriptomic data. By plotting number of sequences assigned to distinct taxonomic units, we identified several HGT candidates and contaminant bacteria.
```bash

```