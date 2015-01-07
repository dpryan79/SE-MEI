SE-MEI
======

This repository contains a number of scripts and programs useful for finding mobile element insertion (MEI) sites from single-end datasets. In short, these programs align soft-clipped regions of alignments against repeat sequences and processes the results to find their locations.

Requirements
------------

 * The pysam python module
 * bowtie1 or another short read aligner.
 * A local aligner, like STAR.

Compilation
-----------

The C files in this repository can be compiled with the supplied Makefile by simply typeing `make` in the source directory.

General pipeline
----------------

1. Map with `STAR` or another local aligner against the genome.
2. Use `extractSoftclipped` on the resulting BAM file to create fastq files of the soft-clipped regions.
3. Download and filter the output from RepeatMasker from UCSC or another source. Filtering could be done as follows (piping to `sed` is only needed if you aligned to an Ensembl reference):

        head -n 2 1/chr1.fa.out
        for chrom in ls */
        do
        awk 'BEGIN{OFS="\t"}{if(NR>3) print $5,$6,$7,$11}' $chrom/*.out | sed 's/chr//g' >> rmask.bed
        done

4. Extract the repeat sequences with `bed2fasta.py` using the just BED file from step 3. Note that this is a VERY slow program and should probably just be rewritten. One could alternatively just use `bedtools getfasta -name`, which is probably faster.
5. Align those results with bowtie1 or another short-read aligner to the repeat sequences.
6. Run `compactRepeats` on the resulting BAM file to create a quality controlled list in BED format of putative insertion sites.
7. Use `bedtools cluster` to group insertion sites into clusters. Ensure that the input is sorted (bedtools can be used for this).
8. Use `awk` to filter each sample's insertion sites by repeat type.
9. Use `filter.py` to produce a BED file of only sites with at least two supporting alignments. A repeat type must be specified.

A more detailed example
-----------------------

Below are the exact commands used to process a single example sample. The fastq file isn't provided as this is meant to simply be an illustrative example.

1. Align with `STAR` (note that multiple samples would be processed in this manner, so the index isn't unloaded until completion):

    STAR --genomeDir indexes --genomeLoad LoadAndKeep --outSAMstrandField intronMotif --outFilterMultimapNmax 2 --outSAMattributes Standard --readFilesCommand zcat --readFilesIn foo.fq.gz --outFileNamePrefix foo --outStd SAM | samtools view -Sbo foo.bam -

2. Extract the soft-clipped sequences in fastq format:

    

Comparisons between groups or samples
-------------------------------------

A common goal is to find insertion sites unique to a given group of samples. This can be done as follows:

    compareGroups Sample1.bed,Sample2.bed Sample3.bed,Sample4.bed All.bed > unique.bed

Note that `compareGroups` requires that the BED files are sorted in lexicographic order. Comparisons between samples can be done in a similar way, by specifying only a single file in one or both groups.
