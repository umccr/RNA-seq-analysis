**STAR command parameters used in bcbio-WTS workflow** 

Argument | Value | Description
---------|-------|-------
--genomeDir | /bcbio/genomes/Hsapiens/GRCh37/star/ | Genome file
--readFilesIn | fastq's | Input data
--runThreadN | 21 | Threads
--outFileNamePrefix | CR180160_VPT-WH031 | Preference
--outReadsUnmapped | Fastx | Output of unmapped reads in separate fasta/fastq files
--outFilterMultimapNmax | 10 | Max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
--outStd | BAM_Unsorted | Which output will be directed to stdout (standard out) - alignments in BAM format, unsorted. Requires –outSAMtype BAM Unsorted
--limitOutSJcollapsed | 2000000 | Max number of collapsed junctions
--outSAMtype | BAM Unsorted | Output Unsorted
--outSAMmapqUnique | 60 | Integer0to255 The number of loci a read maps to
--outSAMunmapped | Within | Output of unmapped reads in the SAM format - output unmapped reads within the main SAM file (i.e. Aligned.out.sam)
--outSAMattributes | NH HI NM MD AS | The SAM attributes
--sjdbGTFfile | ref-transcripts.gtf | Path to annotation file
--sjdbOverhang | 150 | Specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina 2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the ideal value is max(ReadLength)-1. In most cases, a generic value of 100 will work as well as the ideal value.
--outSAMattrRGline | ID:CR180160_VPT-WH031 PL:illumina PU:CR180160_VPT-WH031 SM:CR180160_VPT-WH031 | The read group line
--chimSegmentMin | 12 |  The minimum mapped length of the two segments that is allowed.
--chimJunctionOverhangMin | 12 | Minimum overhang for a chimeric junction - By default, it would require 20b of read sequence on each side of a chimeric junction
--chimScoreDropMax | 30 |  Max drop (difference) of chimeric score (the sum of scores of all chimeric segements) from the read length. The higher value (>=37) will allow chimeras with poorer alignment scores to be output 
--chimSegmentReadGapMax | 5 | 
--chimScoreSeparation | 5 | Minimum difference (separation) between the best chimeric score and the next one
--chimOutType | withinBam | Include chimeric alignments together with normal alignments
--outSAMstrandField | intronMotif | For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute, which STAR will generate with --outSAMstrandField intronMotif option. As required, the XS strand attribute will be generated for all alignments that contain splice junctions. The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed.
--quantMode | TranscriptomeSAM | Outputs alignments translated into transcript coordinates in the Aligned.toTranscriptome.out.bam file (in addition to alignments in genomic coordinates in Aligned.*.sam/bam files). 

**Notes:**

- According to the discussion [here](https://github.com/alexdobin/STAR/issues/333) setting `-alignIntronMax to 50000 --chimScoreDropMax to 37 --alignSplicedMateMapLminOverLmate 0` helps rescue few chimeric alignments.
	- The default for `-alignIntronMax` is `1000000` - The smaller value would prohibit long gaps between mates.


- `--chimSegmentMin` can be further reduced to 5 may be. Needs testing.

- Reading different resources online, seems we are using the most _optimal_ params to capture fusion genes including `--chimSegmentMin`, `--chimJunctionOverhangMin` and `chimSegmentReadGapMax`

- [Best practices for star](https://pdfs.semanticscholar.org/43f4/f2087eb12a5ea1dc578387253c89f07c9df3.pdf) specifes few additional parameters that can be explored 

```
STAR --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax 3
--alignSJDBoverhangMin 10 --alignMatesGapMax 200000 (default 0) --alignIntronMax 200000
(default 0) --alignSJstitchMismatchNmax 5 -1 5 5 
--twopassMode Basic
```

- Other combination of [options](https://arriba.readthedocs.io/en/v0.12.0/execution/) for chimeric reads to consider are:

```
--chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 
--chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 
-1 5 5 --chimSegmentReadGapMax 3 --chimMainSegmentMultNmax 10
```

