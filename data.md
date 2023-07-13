
```bash
screen -S size
srun -J jellyfish -p long --mem-per-cpu 8G --cpus-per-task 8 --pty bash
conda activate jellyfish
cd /projects/nano_diagnostics
jellyfish count -t 8 -C -m 19 -s 10G -o ${OutFile}_19mer_out2 --min-qual-char=? <(zcat $4) <(zcat $5) <(zcat $6) <(zcat $7) <(zcat $8) <(zcat $9) <(zcat $10) <(zcat $11) <(zcat $12) <(zcat $13) <(zcat $14) <(zcat $15) <(zcat $16) <(zcat $17)  

jellyfish count -t 8 -C -m 21 -s 10G -o ${OutFile}_21mer_out --min-qual-char=? <(zcat $4) <(zcat $5) <(zcat $6) <(zcat $7) <(zcat $8) <(zcat $9) <(zcat $10) <(zcat $11) <(zcat $12) <(zcat $13) <(zcat $14) <(zcat $15) <(zcat $16) <(zcat $17) 

jellyfish count -t 8 -C -m 25 -s 10G -o ${OutFile}_25mer_out --min-qual-char=? <(zcat $4) <(zcat $5) <(zcat $6) <(zcat $7) <(zcat $8) <(zcat $9) <(zcat $10) <(zcat $11) <(zcat $12) <(zcat $13) <(zcat $14) <(zcat $15) <(zcat $16) <(zcat $17) 

jellyfish count -t 8 -C -m 31 -s 10G -o ${OutFile}_31mer_out --min-qual-char=? <(zcat $4) <(zcat $5) <(zcat $6) <(zcat $7) <(zcat $8) <(zcat $9) <(zcat $10) <(zcat $11) <(zcat $12) <(zcat $13) <(zcat $14) <(zcat $15) <(zcat $16) <(zcat $17) 

jellyfish histo -h 10000000 -o ${OutFile}_19mer_out.histo ${OutFile}_19mer_out
jellyfish histo -h 10000000 -o ${OutFile}_21mer_out.histo ${OutFile}_21mer_out
jellyfish histo -h 10000000 -o ${OutFile}_25mer_out.histo ${OutFile}_21mer_out
jellyfish histo -h 10000000 -o ${OutFile}_31mer_out.histo ${OutFile}_21mer_out
```

```bash
source package 7f4fb852-f5c2-4e4b-b9a6-e648176d5543
source package 3fe68588-ae0b-4935-b029-7a2dfbf1c4f3

kat count -m 21 --hashing-method jellyfish reads.fastq -o counts.jf
kat plot spectra counts.jf -o kmer_distribution.png
kat comp size kmer_distribution.png

```
Assemblies with default hifiasm settings:
```bash
ls /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trant/Trant_default.p_ctg.fa
ls /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trapi/Trapi_default.p_ctg.fa
ls /hpc-home/tmathers/JIC_TM_scratch_DIR/Caliber_pb_HiFi_assembly/with_third_flow_cell/Trurt/Trurt_default.bp.p_ctg.fa
```
assembly/genome/T_anthrisci/hifiasm/default

curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux.tar.xz --output canu/canu-2.2.Linux.tar.xz 