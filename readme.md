# HMMPolish: A coding region polishing tool for TGS-sequenced RNA viruses
=======================================================================

**HMMPolish** uses **Sequence Alignment Graph** 
 to represent multiple sequences and utilizes **Profile Hidden Markov Models** (pHMMs) to polish coding regions of target RNA virus. 

HMMPolish requires the following as input:
+ Reads file for graph construction.
+ Seed sequence file for graph construction (Can be either draft assembly from aseembly tools including Canu/Shasta/Flye etc, or the raw sequence).
+ pHMMs of targeted protein family/domain.
 
 
## Installation
#### Dependencies
- Conda
- python
- Required python package: **Biopython**>=1.70, networkx >= 2.5.1, pandas >= 1.1.3

#### Installation
```console
git clone --recursive https://github.com/rainyrubyzhou/HMMPolish HMMPolish
cd HMMPolish/src
python HMMPolish.py -h
```
Successfull installation will end with usage information using above commands.

