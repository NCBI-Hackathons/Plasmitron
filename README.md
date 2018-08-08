# Plasmitron
![My image](https://github.com/NCBI-Hackathons/Plasmitron/blob/master/plasmids.png)

Identification of Plasmids from Pacbio Long Read Bacterial Sequences

**Contributors**

 * Laraib Malik
 * Rick Copin
 * John Didion

 ## Table of Contents
* [Intro](https://github.com/NCBI-Hackathons/Plasmitron.git#intro)
* [Quickstart](https://github.com/NCBI-Hackathons/Plasmitron.git#quickstart)
* [Help](https://github.com/NCBI-Hackathons/Plasmitron.git#help)
* [Inputs](https://github.com/NCBI-Hackathons/Plasmitron.git#inputs)
* [Outputs](https://github.com/NCBI-Hackathons/Plasmitron.git#outputs)
* [WorkFlow](https://github.com/NCBI-Hackathons/Plasmitron.git#workFlow)

## Intro
### Goal
Retreive and visualize plasmid sequences from  dataset generated using Long Read sequencing technology.
### Challenges: 
Plasmid hunting with PacBio can be tricky especilly for small size plasmids. One reason is that size selection done as part of library preparation tends to remove fragments < 8kb. The other reason is that they are harder to detect in the assemblies because there will be a subset of reads that has multiple tandem copies of the same plasmid. That is somewhat mitigated by using e.g. circlator, but it still remains an issue.

### Why study plasmid diversity?
Plasmid discovery and characterization are important goals for clinical genomic analysis because almost all bacterial strains harbor at least one plasmid with potentially syndrome- and tissue-specific functions. 
 
### Plasmid diversity in NCBI.
As of August 2018,  18702 plasmid sequences have been identified and deposited in the US National Center for Biotechnology and Information (NCBI) Refseq database.  With the ease and speed of whole genome sequencing, new plasmid are discovered every day, highlighting the impressive breadth of plasmid diversity.

### The Advantage to study Plasmid diversity using PacBio. 
The extent and importance of plasmid contribution to pathogenesis is largely under-appreciated. This is mainly due to the complications inherent to their identification, analysis and characterization. Plasmid are rich in repetitive sequences. As such and in theory, long read sequencing technology is the most appropriate technique to study plasmid genomics.
  
## Quickstart
   - Open terminal or connect to server
    
    git clone https://github.com/NCBI-Hackathons/Plasmitron.git
  
## Help
### Dependencies
1. wget
3. Blast
4. conda
5. minimap2
6. fastq-dump
7. samtools
8. Python
9. Canu
10. Circulator
    
**wget**
   - For mac use Homebrew
   
    brew install wget

**BLAST**
   - Download system compatible BLAST version.
  For MAC:
  
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+.dmg
    export PATH=$PATH:$<path_to>/ncbi-blast-2.2.29+/bin
    
  Download BLAST databases: Follow instrucitons at : https://www.ncbi.nlm.nih.gov/books/NBK52640/
  
    $export BLASTDB=$<path_to>/blastdb
 
**conda**
       
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh

**minimap2**
   
    conda install -c bioconda minimap2  

**fastq-dump**

    wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -vxzf sratoolkit.tar.gz

**samtools**

    wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2
    tar -xjvf samtools.tar.bz2
    cd samtools-1.3.1
    make
    make install

**Canu**

    If you are on linux or macOS:
    * Download a release for your system from https://github.com/marbl/canu/releases
    * Unpack: tar -xJf canu-1.7.1.Darwin-amd64.tar.xz
    * Add the Canu bin folder to your path: export PATH=/path/to/canu/Darwin-amd64/bin:$PATH
    
    Otherwise, you must build from source:
    * Download the source release from https://github.com/marbl/canu/releases
    * Unpack: tar -xzf v1.7.1.tar.gz
    * cd /path/to/canu/src
    * make -j <threads> (where threads is an int >= 1)

**Circlator**

    The installation instructions for Circlator are here: http://sanger-pathogens.github.io/circlator/
    But by far the easiest way to use it is via Docker
    * Install Docker: https://www.docker.com/
    * docker pull sangerpathogens/circlator
    * docker run --rm -it -v /path/to/my/data:/data sangerpathogens/circlator \
    circlator all /data/assembly.fasta /data/reads.fasta /data/output_directory \
    (where 'assembly.fasta' is the name of the Canu assembly and 'reads.fasta' is \
    the fasta file with raw reads)

## Inputs

1. Reference sequence for the bacterial genome. Download the fasta file for the sequenced species using NCBI. 
2. Database for reference sequences for plasmids from the bacteria of interest.
3. SRA number for the sequencing dataset or the relevant fasta/fastq file with the sequenced reads.

## Outputs

 1. A list of potential plasmid sequences present in the dataset, along with the coverage percentage.

## Workflow
 
 1. Obtain SRA dataset and extract fasta files using fastq-dump
        
        fastq-dump -fasta 0 SRR<number>
        
   Or
    
        wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR744/SRR7445584/SRR7445584.sra
        fastq-dump -fasta 0 SRR<number>.sra
       
 2. Align against the bacterial chromosome using minimap2 where sequence.fasta is the reference:
    
        minimap2-2.11_x64-linux/minimap2 -ax map-pb reference_sequence.fasta SRR<number>.fasta > chromAlign.sam

 3. Obtain fasta sequence for reads that are unaligned against the bacterial genome using provided scripts under the "script" folder. 
 
        python getUnmappedSAM.py ../minimap_output/chromAlign.sam ../sra/SRR<number>.fasta ../minimap_output/<name_output>.fasta
 
 4. Align the filtered reads against the reference file made by concatinating all the plasmid reference sequences. Do this alignment using a similar minimap2 command as before. You can use the option "--secondary=no" to avoid multi-mapping.
 
 5. Check for plasmids that have a significant number of reads aligning against them. This will tell you which plasmids are of interest for further analysis. The cut-off for our dataset was '20'. Note that the detected plasmids may have high sequence similarity between them.
 
        samtools view -S -F 4 alignPlasmid.sam | awk '{print $3}' | sort | uniq -c | awk '{if ($1 > 20) print}'
        
 6. Now to check for the number of bases covered by the reads within each plasmid, run the following steps to obtain a coverage estimate file using samtools:
  
        samtools view -Sb alignPlasmid.sam > alignPlasmid.bam
        samtools sort alignPlasmid.bam alignPlasmid.sorted
        samtools depth alignPlasmid.sorted.bam > alignPlasmid.depth
        
 7. Given the depth file and the alignment samfile, run the script provided in folder "analyse" to obtain the final table:
 
        python analyseDepth.py alignPlasmid.depth alignPlasmid.sam out_table.txt
 
## WorkFlow

 1.  Create customized blast databases.

### Step 1.
# ------------------
    
      # Retrieve plasmid fasta.
      # creates plasmid fasta sequences in plasmids/fasta

      # from NCBI website, go to Nucleotide database:
      # plasmid[title] AND staphylococcus[title]
      # Filters: Species = bacteria
      # Filters: Molecule types = genomic RNA/DNA
      # Filters: Genetic compartments = Plasmid
      # download accession table
      # save as "accession_plasmids.txt" in /plasmids/assembly_ID/


### Step 2. 
# ------------------
      # Create customized blast databases.
      # use makeblastdb

      
### Step 3. 
# ------------------


### Step 4. 
# ------------------


### Step 5. 
# ------------------
   

# ==============


## Case Study
### Staphylococcus aureus genomes and mobile genetic elements (MGE).
The proliferation of S. aureus genomic sequences in public databases reflects strong interest in understanding S. aureus genome diversity and evolution. With an average of 2,800 coding sequences, it is estimated that 44% of S. aureus genes are NOT shared by all S. aureus strains. These genes constitute the  ‘accessory genome’, which is variable between strains and mostly made of mobile genetic elements (MGE) enriched in hypothetical proteins. 

Staphylococcal MGE encompass any intra- or extra-chromosomal DNA segment that can be independently mobilized within or between S. aureus cells. It includes plasmids, transposons, integrons, genomic islands, S. aureus pathogenicity islands (SaPIs), integrative conjugative elements, staphylococcal chromosome cassettes, and phages. Together, phages and plasmids are the main source of MGE diversity among S. aureus strains. 
