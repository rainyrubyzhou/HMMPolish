# HMMPolish: A coding region polishing tool for TGS-sequenced RNA viruses
=======================================================================

**HMMPolish** uses **Sequence Alignment Graph** 
 to represent multiple sequences and utilizes **Profile Hidden Markov Models** (pHMMs) to polish coding regions of target RNA virus. 
 Users are welcome to read our published paper at Briefings in Bioinformatics [here](https://doi.org/10.1093/bib/bbad264) and try our tool as guided below.

HMMPolish requires the following as input:
+ Reads file for graph construction.
+ Seed sequence file for graph construction (Can be either draft assembly from assembly tools including Canu/Shasta/Flye etc, or the raw sequence).
+ pHMMs of targeted protein family/domain.
 
 
## Installation
#### Dependencies
- Biopython
- HMMER
- pandas
- BLASR

#### Installation
```console
git clone --recursive https://github.com/rainyrubyzhou/HMMPolish HMMPolish
cd HMMPolish/src
python HMMPolish.py -h
```
Successful installation will end with usage information using the above commands.

## Usage of HMMPolish: 
>**Command Usage:**
```console
python HMMPolish.py --read READFILE --seed SEEDFILE 
                    --hmm HMMFILE -o OUTFILE   --wei WEI      
```

>**Mandatory args:**
```console
--read  <str, e.g. "raw_HIV.fasta">
Reads file for graph construction (in fasta format). 

--seed | <str, e.g. "canu_contig.fasta">
Backbone sequence file for graph construction (in fasta format). 

--hmm | <str, e.g. "canu_contig.fasta">
```
>**Optional args:**
```console
--wei WEI | <float, e.g. 0.8> 
Weight of Viterbi score in the recursive function (default: 0.9).

-v, --verbose | increase output verbosity
Will output the path score and the HMMER results of polished sequences under this mode.

-h | Print the usage information. 
```


>**Example usage:** 
```console
python HMMPolish.py --read test/test_noro.fa --seed test/canu_ass.fa --hmm test/7_noro_profile.hmm -o polished.fa
```
+ `../data/test_noro.fa ` example reads file 

>**Example Outputs:** 
+ `polished.fa` containing polished reads of coding regions covered by each pHMM.

## Data Availability

Users can download the profiles of some common RNA viruses following the referred accession numbers we provided [here](src/common_profiles.md).


## Contact
Other than raising issues on Github, you can also contact YU Runzhou (runzhouyu2-c@my.cityu.edu.hk) for help in installation/usage or any other related query.


