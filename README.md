# Plasmitron
![My image](https://github.com/NCBI-Hackathons/Plasmitron/blob/master/plasmids.png)

Identification of Plasmids from Pacbio Long Read Bacterial Sequences

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
  
## Quickstart
   - Open terminal or connect to server
    
    git clone https://github.com/NCBI-Hackathons/Plasmitron.git
  
## Help
### Dependencies
1. wget
2. magicBlast
3. Blast
4. R
5. minimap2
    
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
 
**magicBLAST**
   - Download system compatible magicBLAST version.
  For MAC:
    
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/LATEST/ncbi-magicblast-1.3.0.dmg
    export PATH=$PATH:$<Path_to>/ncbi-magicblast-1.3.0/bin/

**minimap2**
   - Download precompiled binaries from the github respository using:
    
    curl -L https://github.com/lh3/minimap2/releases/download/v2.12/minimap2-2.12_x64-linux.tar.bz2 | tar -jxvf -
    ./minimap2-2.12_x64-linux/minimap2

**R**
   - Install R locally
   
    cd Path_to/Got_plasmid/
    wget http://cran.rstudio.com/src/base/R-3/R-3.4.1.tar.gz 
    tar -xzvf R-3.4.1.tar.gz 
    cd R-3.4.1 
     
    ./configure --prefix=$HOME/R 
    make && make install 
    export PATH=$PATH:$HOME/R/bin
    
**R_modules**

    R
    source("https://bioconductor.org/biocLite.R")
    biocLite("GenomicRanges")
    biocLite("IRanges")
    biocLite("Biostrings")
    install.packages("data.table")
    install.packages("reutils")
    install.packages("devtools")
    install.packages("biofiles")
    devtools::install_github("gschofl/biofiles")

## Inputs

1. Reference sequence for the bacterial genome. Download the fasta file for the sequenced species using NCBI. 
2. Database for reference sequences for plasmids from the bacteria of interest.
3. SRA number for the sequencing dataset or the relevant fasta/fastq file with the sequenced reads.

## Outputs
**file.txt**

## Workflow
 
 1. Obtain SRA (SRR7445584) dataset and extract fasta files using fasterq-dump
 2. Align against the bacterial chromosome using minimap2 where sequence.fasta is the reference:
    
        minimap2-2.11_x64-linux/minimap2 -ax map-pb sequence.fasta SRR7445584.fasta > chromAlign.sam

 3. Obtain fasta sequence for reads that are unaligned against the bacterial genome using provided scripts under the "filter" folder. 
 
 4. Align the filtered reads against the reference file made by concatinating all the plasmid reference sequences. Do this alignment using the minimap2 command as before. 
 
 5. Check for plasmids that have a significant number of reads aligning against them. The cut-off for our dataset was '20':
 
        samtools view -S -F 4 <samfile> | awk '{print $3}' | sort | uniq -c | awk '{if ($1 > 20) print}'
 
## WorkFlow

 1.  Create customized blast databases.
 2.

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
      # Use magicBlast on one SRRA versus individual plasmid databases.
      # creates magicBlast output files in plasmids/magic_output/


### Step 4. 
# ------------------
      # Generate a table with the % of plasmid sequences covered by contigs.



### Step 5. 
# ------------------
      # Generate plasmid visualization





# ==============


## Case Study
### Staphylococcus aureus genomes and mobile genetic elements (MGE).
The proliferation of S. aureus genomic sequences in public databases reflects strong interest in understanding S. aureus genome diversity and evolution. With an average of 2,800 coding sequences, it is estimated that 44% of S. aureus genes are NOT shared by all S. aureus strains. These genes constitute the  ‘accessory genome’, which is variable between strains and mostly made of mobile genetic elements (MGE) enriched in hypothetical proteins. 

Staphylococcal MGE encompass any intra- or extra-chromosomal DNA segment that can be independently mobilized within or between S. aureus cells. It includes plasmids, transposons, integrons, genomic islands, S. aureus pathogenicity islands (SaPIs), integrative conjugative elements, staphylococcal chromosome cassettes, and phages. Together, phages and plasmids are the main source of MGE diversity among S. aureus strains. 

### Why study Staphylococcus aureus plasmid diversity?
MGE discovery and characterization are important goals for clinical genomic analysis because almost all S. aureus strains harbor at least one plasmid with potentially syndrome- and tissue-specific functions. 
 
### Plasmid diversity in NCBI.
 As of August 2017,  327 unique plasmid sequences have been identified and deposited in the US National Center for Biotechnology and Information (NCBI) database.  With the ease and speed of whole genome sequencing, new MGE are discovered every day, highlighting the impressive breadth of S. aureus plasmid diversity. There are ~40,000 primary (unannotated) S. aureus datasets in the Sequence Read Archive (SRA). 

### The challenge to study MGE diversity. 
The extent and importance of plasmid contribution to S. aureus pathogenesis is largely under-appreciated. This is mainly due to the complications inherent to their identification, analysis and characterization. 
 
Plasmid are rich in repetitive sequences. High level of sequence identity facilitates genetic recombination and contributes to the emergence of mosaic MGE. As such, MGE are ever-changing and can be hard to identify. Moreover, nomenclature in public databases is constantly evolving and inconsistency in annotation among MGE is common and complicates functional inter- and intra-species analyses. 
