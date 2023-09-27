image: gitpod/workspace-full-vnc
# Update 7th September to reflect code base changes
ports:
- name: JBrowseWeb
  description: The JBrowse Webserver port
  port: 3000
  onOpen: open-browser
  visibility: public

- name: HiGlass
  description: The HiGlass port
  port: 8989
  onOpen: open-browser
  visibility: public

tasks:
  - name: Install Singularity 3.11.4
  # https://docs.sylabs.io/guides/3.0/user-guide/installation.html
    init: |
      cd /workspace/treeval-curation/

      sudo apt-get update && sudo apt-get install -y \
      build-essential \
      libssl-dev \
      uuid-dev \
      libgpgme11-dev \
      squashfs-tools \
      libseccomp-dev \
      pkg-config

      mkdir -p $GOPATH/src/github.com/sylabs && \
      cd $GOPATH/src/github.com/sylabs && \
      wget https://github.com/sylabs/singularity/releases/download/v3.11.4/singularity-ce-3.11.4.tar.gz && \
      tar -xzf singularity-ce-3.11.4.tar.gz && \
      cd ./singularity-ce-3.11.4 && \
      ./mconfig

      ./mconfig && \
      make -C ./builddir && \
      sudo make -C ./builddir install
  
  - name: Install Nextflow
  # https://www.nextflow.io/docs/latest/getstarted.html
    init: |
      cd /workspace/treeval-curation/
      
      wget -qO- https://get.nextflow.io | bash

      chmod +x nextflow

      nextflow self-update

  - name: Install JBrowse2
  # https://jbrowse.org/jb2/download/#jbrowse-cli-tools
    command: |
      cd /workspace/treeval-curation/
      
      npm install -g @jbrowse/cli

      jbrowse create jbrowse2

      cd jbrowse2/

      npx serve . -l 3000

  - name: Install TreeVal Pipeline
  # https://github.com/sanger-tol/treeval
    init: |
      cd /workspace/treeval-curation/
        
      git clone -b pre-tag https://github.com/sanger-tol/treeval.git

  - name: Install Curtation Pretext
  # https://github.com/sanger-tol/curationpretext
    init: |
      cd /workspace/treeval-curation/

      git clone -b dev https://github.com/sanger-tol/curationpretext.git

  - name: Install HiGlass
  # https://docs.higlass.io/tutorial.html
    init: |
      cd /workspace/treeval-curation/
      
      pip install higlass-manage

      higlass-manage start
  
  - name: Alias Nextflow
    init: |
      cd /workspace/treeval-curation/
      
      echo "alias nextflow_cmd='/workspace/treeval-curation/nextflow'" >> ~/.bashrc

      source ~/.bashrc

  - name: Download busco for nematode
    init: |
      cd /workspace/treeval-curation/

      curl https://dp24.cog.sanger.ac.uk/Busco.tar.gz | tar xzf -

  - name: Download Nematode Test data and make synteny
    init: |
      cd /workspace/treeval-curation/

      curl https://dp24.cog.sanger.ac.uk/Nematode.tar.gz | tar xzf -

      mkdir -p /workspace/treeval-curation/synteny/nematode/

      cp /workspace/treeval-curation/Oscheius_DF5033/genomic_data/Oscheius_DF5033.fa /workspace/treeval-curation/synteny/nematode/SuperNematode.fa

  - name: Download Lepidoptera data
    init: |
      cd /workspace/treeval-curation/

      curl https://dp24.cog.sanger.ac.uk/ilTorViri5.tar.gz | tar xzf -

  - name: Download Genomic Alignment data
    init: |
      cd /workspace/treeval-curation/

      curl https://dp24.cog.sanger.ac.uk/AlignmentData.tar.gz | tar xzf -  

  - name: Open Tutorial Page
    init: |
      gp preview https://bga23.org/treeval-curation/Tutorial/

github:
  prebuilds:
    # enable for the master/default branch (defaults to true)
    master: true
    # add a "Review in Gitpod" button as a comment to pull requests (defaults to true)
    addComment: true
    # add a "Review in Gitpod" button to pull requests (defaults to false)
    addBadge: true
    # add a label once the prebuild is ready to pull requests (defaults to false)
    addLabel: prebuilt-in-gitpod

vscode:
  extensions:                               # based on nf-core.nf-core-extensionpack
    - codezombiech.gitignore                # Language support for .gitignore files
    - esbenp.prettier-vscode                # Markdown/CommonMark linting and style checking for Visual Studio Code
    - EditorConfig.EditorConfig             # override user/workspace settings with settings found in .editorconfig files
    - mechatroner.rainbow-csv               # Highlight columns in csv files in different colors
    - nextflow.nextflow                     # Nextflow syntax highlighting
    - oderwat.indent-rainbow                # Highlight indentation level
    - streetsidesoftware.code-spell-checker # Spelling checker for source code

```bash
cd /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/raw_fasta
for i in *.fasta.gz; do
gunzip $i;
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/scripts/GA_data_prep.py ${i/.gz} ncbi 10;
done
mv GallusGallus/ ../../gene_alignment_data/bird/
cd ../
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 scripts/GA_csv_gen.py /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_data
```
```bash
mv GCA_012654025.1_Pven_dovetail_genomic.fna P_venustra.fasta gene_alignment_prep/raw_data/
mv GCA_024516035.1_ASM2451603v1_genomic.fna B_cockerel.fasta

cp /hpc-home/did23faz/treeval/treeval-resources/synteny/insects/P_venustra.fasta /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/raw_fasta/PachypsyllaVenustra-P_venustra.cdna.fasta

cp /hpc-home/did23faz/treeval/treeval-resources/synteny/insects/B_cockerel.fasta /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/raw_fasta/BactericeraCockerelli-B_cockerel.cdna.fasta

mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/protein.faa /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/raw_fasta/DiaphorinaCitri-D_citri000475195.pep.fasta
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/cds_from_genomic.fna /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/raw_fasta/DiaphorinaCitri-D_citri000475195.cds.fasta
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/rna.fna /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/raw_fasta/DiaphorinaCitri-D_citri000475195.rna.fasta
mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/GCF_000475195.1_Diaci_psyllid_genome_assembly_version_1.1_genomic.fna /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/raw_fasta/DiaphorinaCitri-D_citri000475195.cdna.fasta

cd /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/raw_fasta
for i in *.fasta; do
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_prep/scripts/GA_data_prep.py ${i} ncbi 100;
done
#NOTE: files must be named in x-x.x.fasta format

mv DiaphorinaCitri/ ../../gene_alignment_data/insects/
mv PachypsyllaVenustra/ ../../gene_alignment_data/insects/
mv BactericeraCockerelli/ ../../gene_alignment_data/insects/

cd ../
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 scripts/GA_csv_gen.py /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_data

nano /hpc-home/did23faz/treeval/treeval-resources/treeval_yaml/T_urticae.yaml
```
Create yaml file:
```bash
assembly:
  level: contig
  sample_id: Turt7151223005
  latin_name: Trioza_urticae
  classT: insects
  asmVersion: 1
  gevalType: DTOL
reference_file: /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/T_urticae_715m_12_2_3.0_0.5.bp.p_ctg.fa
assem_reads:
  pacbio: /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiFi/fasta
  hic: /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/raw_data/T_urticae/HiC/cram
  supplementary: path
alignment:
  data_dir: /hpc-home/did23faz/treeval/treeval-resources/gene_alignment_data
  common_name: "" # For future implementation (adding bee, wasp, ant etc)
  geneset: "DiaphorinaCitri-D_citri000475195,BactericeraCockerelli-B_cockerel,PachypsyllaVenustra-P_venustra"
  #Path should end up looking like "{data_dir}{classT}/{common_name}/csv_data/{geneset}-data.csv"
self_comp:
  motif_len: 0
  mummer_chunk: 10
intron:
  size: "50k"
telomere:
  teloseq: TTAGGG
synteny:
  synteny_genome_path: /hpc-home/did23faz/treeval/treeval-resources/synteny
busco:
  lineages_path: /jic/research-groups/Saskia-Hogenhout/BUSCO_sets/v5
  lineage: hemiptera_odb10
```
```bash
nextflow run main.nf -profile singularity --input /hpc-home/did23faz/treeval/treeval-resources/treeval_yaml/T_urticae.yaml -entry FULL --outdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Psyllidae/assembly/genome/T_urticae/hifiasm_19.5/715m/12/2/3.0/0.5/treeval
```