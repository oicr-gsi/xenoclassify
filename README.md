# xenoclassify
Xenoclassify is a tool that provides a means by which to classify short read sequencing data generated from xenograft samples. It requires alignment to the reference genomes of the graft and host species using bwa mem. Once aligned, reads (or read pairs) are assessed to identify the likely source of the cells from which the DNA was extracted. The output is a bam file with reads tagged to indicate the source species.

