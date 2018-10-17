## Bacterial WGS training : Exercise 2

|**Title**| Sequence quality and assambly.|
|---------|-------------------------------------------|
|**Training dataset:**|  
|**Questions:**| <ul><li>How do I know if my data was correctly sequenced?</li><li>How can I improve my data quality?</li><li>Which steps lead to an assambly?</li><li>How do I map?</li><li>How do I fix errors in the mapping?</li><li>How do I assamble the reads?</li></ul>|
|**Objectives**:|<ul><li>Check quality of sequenced data.</li><li>Trimm low quality segments and adapters.</li><li>Introduction to the process of assembly.</li><li>Mapp raw reads against reference.</li><li>Remove duplicated reads.</li><li>Assamble mapped reads.</li></ul>|  
|**Time estimation**:| 1 h 30 min |
|**Key points**:|<ul><li>Analysis of sequence quality.</li><li>Assembly steps.</li></ul>|

## How do I know if my data was correctly sequenced?

Despite the improvement of sequencing methods, there is no error-free technique. The Phred quality score (Ewing et al., 2005) has been used since the late 90s as a measure of the quality of each sequenced nucleotide. Phred quality scores not only allow us to determine the accuracy of sequencing and of each individual position in an assembled consensus sequence, but it is also used to compare the efficiency of the sequencing methods.

Phred quality scores Q are defined (Ewing et al., 2005) as a property which is logarithmically related to the base-calling error probabilities P . The Phred quality score is the negative ratio of the error probability to the reference level of P = 1 expressed in Decibel (dB): 

<p align="center">𝑄 = −10 log<sub>10</sub> P</p>

A correct measuring of the sequencing quality is essential for identifying problems in the sequencing and removal of low-quality sequences or sub sequences. Conversion of typical Phred scores used for quality thresholds into accuracy can be ead in the following table:

|**Phred score**| Error probability | Accuracy|
|----------------|--------------------|---------|
|10|1/10|90%|
|20|1/100|99%|
|30|1/1000|99.9%|
|40|1/10000|99.99%|
|50|1/100000|99.999%|
|60|1/1000000|99.99999%|
|70|1/10000000|99.999999%|

There are multiple software to read and generate statistics to help with the interpretation of the quality of a sequence. One of the most commonly used methods for this task is FastQC (Andrews, 2010), a java program that run on any system and has both command line and graphic interface.

