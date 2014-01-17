#!/bin/bash


### REFERENCES ### 

EXONS=


### HARD CODED PATHS ###

### GTF PATH ###
#GTF=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/gencode18/gencode.v18.annotation.gtf 
#CHRS=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/GRCh37_ensemble19/chrs

### GENOME IDXS ###

#CHRS_100_CLIP=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/GRCh37_ensemble19/indexes/GRCh37.p11.genome_25_F1.index
#CHRS_100_F2=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/GRCh37_ensemble19/indexes/GRCh37.p11.genome_50_F2.index
#CHRS_100_F3=/export/uec-gs1/knowles/analysis/tade/references_and_annotation/human/GRCh37_ensemble19/indexes/GRCh37.p11.genome_50_F3.index



GTF=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/mini.gtf 
CHRS=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/chrs 
CHR_LIST=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gtfar-100-GENOME.txt

EXONS_100_REF=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gtfar-100-EXONS.fa 
INTRONS_100_REF=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gtfar-100-INTRONS.fa 
GENES_100_REF=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gtfar-100-GENES.fa 

#### TRANSCRIPTOME REFERENCES ####

#EXONS_100=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gc18_100_exonSeqs.MINI.fa
#GENES_100=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gc18_100_geneSeqs.MINI.fa
#INTRONS_100=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gc18_100_intronSeqs.MINI.fa
#EXONS_50=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gc18_50_exonSeqs.MINI.fa
#GENES_50=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gc18_50_geneSeqs.MINI.fa
#INTRONS_50=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gc18_50_intronSeqs.MINI.fa
#EXONS_75=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gc18_75_exonSeqs.MINI.fa
#GENES_75=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gc18_75_geneSeqs.MINI.fa
#INTRONS_75=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/refs/gc18_75_intronSeqs.MINI.fa


#### TRANSCRIPTOME INDEXES ####_
#EXONS_100_F1=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_exonSeqsMINI.index.F1
#EXONS_100_F2=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_exonSeqsMINI.F2.index
#EXONS_100_F3=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_exonSeqsMINI.index.F3
#EXONS_100_F4=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_exonSeqsMINI.index.F4
#GENES_100_F0=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_geneSeqsMINI.21gapIndex.F0
#GENES_100_F1=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_geneSeqsMINI.21gapIndex.F1
#GENES_100_IDX=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_geneSeqsMINI.25gap.F1.index
#INTRONS_100_F1=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_intronSeqsMINI.index.F1
#INTRONS_100_F2=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_intronSeqsMINI.F2.index
#INTRONS_100_F3=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_intronSeqsMINI.index.F3
#INTRONS_100_F4=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_100_intronSeqsMINI.index.F4
#EXONS_50_F1=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_exonSeqsMINI.index.F1
#EXONS_50_F2=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_exonSeqsMINI.index.F2
#EXONS_50_F3=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_exonSeqsMINI.index.F3
#EXONS_50_F4=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_exonSeqsMINI.index.F4
#GENES_50_F0=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_geneSeqsMINI.21gapIndex.F0
#GENES_50_F1=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_geneSeqsMINI.21gapIndex.F1
#GENES_50_F2=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_geneSeqsMINI.25gapIndex.F2
#INTRONS_50_F1=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_intronSeqsMINI.index.F1
#INTRONS_50_F2=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_intronSeqsMINI.index.F2
#INTRONS_50_F3=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_intronSeqsMINI.index.F3
#INTRONS_50_F4=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_50_intronSeqsMINI.index.F4
#EXONS_75_F1=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_exonSeqsMINI.index.F1
#EXONS_75_F2=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_exonSeqsMINI.index.F2
#EXONS_75_F3=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_exonSeqsMINI.index.F3
#EXONS_75_F4=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_exonSeqsMINI.index.F4
#GENES_75_F0=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_geneSeqsMINI.21gapIndex.F0
#GENES_75_F1=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_geneSeqsMINI.21gapIndex.F1
#GENES_75_F2=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_geneSeqsMINI.25gapIndex.F2
#INTRONS_75_F1=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_intronSeqsMINI.index.F1
#INTRONS_75_F2=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_intronSeqsMINI.index.F2
#INTRONS_75_F3=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_intronSeqsMINI.index.F3
#INTRONS_75_F4=/export/uec-gs1/knowles/analysis/tade/gtfar_tests/2014_TESTS/mini_refs/indexes/gc18_75_intronSeqsMINI.index.F4
