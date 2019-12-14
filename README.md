SE-MEI
======

This repository contains a number of scripts and programs useful for finding mobile element insertion (MEI) sites from single-end datasets. In short, these programs align soft-clipped regions of alignments against repeat sequences and processes the results to find their locations.

Requirements
------------

 * The pysam python module
 * bowtie1 or another short read aligner.
 * A local aligner, like STAR.
 * htslib 1.1, which is a submodule

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
8. Use `filter.py` to produce a BED file of only sites with at least two supporting alignments. A repeat type must be specified, since the output only contains a single repeat type.

A more detailed example
-----------------------

Below are the exact commands used to process a single example sample. The fastq file isn't provided as this is meant to simply be an illustrative example.

1. Align with `STAR` (note that multiple samples would be processed in this manner, so the index, which we assume has already been created, isn't unloaded until completion):

    ```
    STAR --genomeDir indexes --genomeLoad LoadAndKeep --outSAMstrandField intronMotif --outFilterMultimapNmax 2 --outSAMattributes Standard --readFilesCommand zcat --readFilesIn foo.fq.gz --outFileNamePrefix foo --outStd SAM | samtools view -Sbo foo.bam -
    ```

2. Extract the soft-clipped sequences in fastq format:

    ```
    extractSoftclipped foo.bam > foo.fq.gz
    ```

3. Download, filter and merge RepeatMasker sequences from UCSC (this example uses a mouse dataset and extract only LINE/L1 and SINE/Alu regions):

    ```
    wget --quiet http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromOut.tar.gz
    tar xf chromOut.tar.gz
    for d in `ls -d */`
    do
        awk '{if(NR>3 && ($11=="LINE/L1" || $11=="SINE/Alu")) printf("%s\t%s\t%s\t%s\n",$5,$6,$7,$11)}' $d/*.out | sed -e 's/chr//g' >> rmask.bed #Convert UCSC to Ensembl chromosome names
    done
    ```

4. Extract the genomic sequence of the repeat regions:

    ```
    bed2fasta.py rmask.bed reference.fa > rmask.fa
    ```

   Alternatively:

    ```
    bedtools getfasta -name -fi reference.fa -bed rmask.bed -fo /dev/stdout | fold -w 60 > rmask.fa #This is faster than `bed2fasta.py`
    ```

5. Index and align to the repeat regions (with bowtie1 in this example):

    ```
    bowtie-build rmask.fa rmask
    bowtie -a --best --strata -p 6 rmask <(zcat foo.fq.gz) -S | samtools view -Sbo foo.bam - #You can ignore the "duplicated sequence" warnings
    ```

6. Perform a bit of QC and sort the sites:

    ```
    compactRepeats foo.bam | bedtools sort > foo.bed #This could be merged with step 5 by specifying "-" as the input file.
    ```

7. Group reads into clusters (with multiple samples, concatenate the output from step 6 of each sample and then use `bedtools sort` before clustering):

    ```
    bedtools cluster -d 50 foo.bed > foo.clustered.bed
    ```

8. Produce unique insertions of a single type

    ```
    filter.py foo.clustered.bed "LINE/L1" foo.LINE.bed
    ```

Comparisons between groups or samples
-------------------------------------

A common goal is to find insertion sites unique to a given group of samples. This can be done as follows:

    compareGroups Sample1.bed,Sample2.bed Sample3.bed,Sample4.bed All.bed > unique.bed

Note that `compareGroups` requires that the BED files are sorted in lexicographic order. Comparisons between samples can be done in a similar way, by specifying only a single file in one or both groups.
