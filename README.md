MultiMotifMaker
=========

Overview
--------
MultiMotifMaker is a multi-thread tool for identifying DNA methylation motifs from 
Pacbio reads. It is an efficient, multi threads implementation of MotifMaker 
([https://github.com/PacificBiosciences/MotifMaker](https://github.com/PacificBiosciences/MotifMaker)) to take full advantage of 
multi-processors to achieve significant acceleration. 

DNA methylation is the most common form of DNA modification in the genomes of 
prokaryotes and eukaryotes and plays a vital role in many critical biological processes.
The methylation of DNA bases is catalyzed by DNA methyltransferases, which bind to a special DNA sequence motif.
Recently, the third generation sequencing technologies, such as Pacbio SMRT, 
provide a new way to identify base methylation in the genome. For each methylated site, 
the SMRT pipeline outputs a sequence of 41 bases centered by the methylated base. 
The number of methylated site ranges from ten thousands for E. coli to ten millions for human.
Identifying methylation motif from the output of SMRT pipeline differs from previous de novo motif finding algorithms.
See the publication for more background on modification detection: ([http://nar.oxfordjournals.org/content/early/2011/12/07/nar.gkr1146.full](http://nar.oxfordjournals.org/content/early/2011/12/07/nar.gkr1146.full))

Algorithm
---------
Existing motif finding algorithms such as MEME, Gibbs motif Sampler and MEpigram. 
None of those methods considers the base centralized sequences as input. 
In addition, none of those methods includes the base modification signals in their algorithms. 

In order to find methylation motifs from SMRT output methylation sequences, 
the PacBio developed a tool, MotifMaker.  
Given a list of modification detections and a genome sequence, MotifMaker searches 
all possible motifs using a motif score. The search is gradually performed from short 
to longer motifs using a branch-and-bound method. However, MotifMaker generally executes 
in single-threaded and the search process is very time consuming(MotifMaker).

Here, we give a rough overview of the algorithm used by MultiMotifMaker.
The branch-and-bound search step, which is designed to search motifs from 
modification sequences through a series of iterations, is the most time-consuming procedure of 
overall workflow of MotifMaker. Since every expansion node of the solution space tree 
can calculate its subtree independently, multiple expansion nodes can be computed at the same time. 
Therefore, according to the branching rule in the branch-and-bound method, we may first calculate the candidate living nodes 
of the first *k* layers of the invisible tree with *n* son nodes for every expansion node. 
Then, branch-and-bound search method will be applied to these candidate living nodes respectively. 
Consequently, it is possible to submit these computing tasks to a thread pool to achieve parallel computation. 
Every thread will search their local solution space trees to obtain local optimal solutions independently. 
Additionally, we set two global variables to represent the current optimal solution and the paths searched. 
After every thread ends, two global variables should be updated. The following task threads 
will adopt these two values as parameters to prune for their search processes. Finally, 
when all executions finish, we got the final optimal motif. We take full advantage of multi-core
concurrent computation, which has obviously accelerated.


Usage
----

The jar supplied in artifacts/MultiMotifMaker.jar bundles all
dependencies and should be runnable on most systems. The sample datasets can be found in ./resources.
(The original datasets are as follows, *Geobacter metallireducens* data: [https://github.com/PacificBiosciences/MotifMaker/tree/master/src/test/resources](https://github.com/PacificBiosciences/MotifMaker/tree/master/src/test/resources) , *E.coli* data: [https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly](https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly), *Arabidopsis* data: NCBI(PRJNA314706),[https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA314706&go=go](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA314706&go=go).)

For command-line motif finding, run the 'find' sub-command, and pass
the reference fasta and the modifications.gff(.gz) file emitted by the
PacBio modification detection workflow (SMRT Analysis: [https://www.pacb.com/products-and-services/analytical-software/smrt-analysis/](https://www.pacb.com/products-and-services/analytical-software/smrt-analysis/)).   
The reprocess command annotates the gff with motif information for better genome browsing.

```
$ java -jar artifacts/MultiMotifMaker.jar

Usage: MultiMotifMaker [options] [command] [command options]
  Options:
    -h, --help
                 Default: false
  Commands:
    find      Run motif finding
      Usage: find [options]
        Options:
        * -f, --fasta      Reference fasta file
        * -g, --gff        modifications.gff or .gff.gz file
          -m, --minScore   Minimum Qmod score to use in motif finding
                           Default: 30.0
        * -o, --output     Output motifs csv file
          -x, --xml        Output motifs xml file
          -l, --layer      Search depth used to Parallelize motif finder
                           Default: 1
          -t, --thread     The concurrency to parallelize motif finder,this parameter should be the number of CPU
                           Default: 16

    reprocess      Reprocess gff file with motif information
      Usage: reprocess [options]
        Options:
          -c, --csv           Raw modifications.csv file
        * -f, --fasta         Reference fasta file
        * -g, --gff           original modifications.gff or .gff.gz file
              --minFraction   Only use motifs above this methylated fraction
                              Default: 0.75
        * -m, --motifs        motifs csv
        * -o, --output        Reprocessed modifications.gff file
```
  For example: java -jar ./artifacts/MultiMotifMaker.jar find --fasta ./resources/lambda/lambda.fasta --gff ./resources/lambda/lambda_modifications.gff
  --minScore 30.0 --output ./resources/lambda/lambda_motifs.csv --xml ./resources/lambda/lambda_motifs.xml --layer 1 --thread 8
  or java -jar ./artifacts/MultiMotifMaker.jar reprocess --fasta ./resources/lambda/lambda.fasta --gff ./resources/lambda/lambda_modifications.gff --motifs ./resources/lambda/lambda_motifs.csv --output ./resources/lambda/lambda_modifications_new.gff

Output file descriptions:
-------------------------

Using the ``find`` command:

- Output csv file: This file follows the same format as the standard
  "Fields included in motif_summary.csv" described in the Methylome
  Analysis White Paper
  ([https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Methylome-Analysis-Technical-Note](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Methylome-Analysis-Technical-Note)).

- Output xml file: This is an output used by SMRT Portal and is not
  necessary using the command line. Simply do not include the -x
  command option. The information contained in this file is used to
  fill in the standard motif report table in SMRT Portal and is
  redundant with the CSV output file.

Using the ``reprocess`` command:
(Reprocessing will update a modifications.gff file with information based on new Modification QV thresholds)

- Output gff file: The format of the output file is the same as the
  input file, and is described in the Methylome Analysis White Paper
  ([https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Methylome-Analysis-White-Paper](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Methylome-Analysis-White-Paper))
  under "Fields included in the modifications.gff file".

DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
