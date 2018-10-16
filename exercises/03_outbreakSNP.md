# Bacterial WGS training : Exercise 3

|**Title**| SNP-based bacterial outbreak investigation.|
|---------|-------------------------------------------|
|**Training dataset:**|                                |
|**Questions:**| <ul><li>Do I have the needed depth of coverage?</li><li>Have I chosen the correct reference?</li><li>How do I create a SNP matrix? How many SNPs do I have?</li><li>How can I visualize my phylogenetic tree? Which problems can I encounter?</li><li>Which strains belong to the outbreak?</li></ul>|
|**Objectives**:|<ul><li>Trimming and quality control of raw reads.</li><li>Mapping against genome reference and duplicate filter.</li><li>Variant Calling.</li><li>SNP matrix creation.</li><li>Maximum Likelihood phylogeny.</li><li>Visualization of results.</li></ul>|  
|**Time estimation**:| 1 h|
|**Key points**:|<ul><li>Importance of reference selecion in SNP-based tipification.</li><li>Variant calling and SNP reconstruction is a key step in the process.</li><li>Interpretation of results is case, species and epidemiology dependant.</li></ul>|
  

## Introduction
Although scientific community efforts have been focused on assembly-based methods and the optimization of reconstructing complete genomes, variant calling is a essential procedure that allows per base comparison between different genomes ([Olson et al 2015](https://www.frontiersin.org/articles/10.3389/fgene.2015.00235/full)).

SNP-based strain typing using WGS can be performed via reference-based mapping of either reads or assembled contigs. There are many available microbial SNP pipelines, such as [Snippy](https://github.com/tseemann/snippy), NASP, SNVphyl, CFSAN SNP Pipeline, or Lyve-SET
