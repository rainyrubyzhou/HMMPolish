# HMMPolish: A coding region polishing tool for TGS-sequenced RNA viruses
=======================================================================

**HMMPolish** uses **Sequence Alignment Graph** 
 to represent multiple sequences and utilizes **Profile Hidden Markov Models** (pHMMs) to polish coding regions of target RNA virus. 

HMMPolish requires the following as input:
+ Reads file for graph construction.
+ Backbone sequence file for graph construction (Can be either draft assembly from aseembly tools including Canu/Shasta/Flye etc, or the raw sequence).
+ pHMMs of targeted protein family/domain.
 
 
