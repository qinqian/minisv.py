## <a name ="start"></a>Getting started

```sh

mamba create -n msvpy python==3.12 poetry seqtk samtools
conda activate msvpy
git clone https://github.com/qinqian/minisv.py
cd minisv.py && make

# extract candidate SVs from hg38- and pangenome-based long read alignments
minisv extract -n COLO829T tests/demo/COLO829T.hs38l.paf.gz | gzip - > COLO829T_hs38l.rsv.gz
minisv extract -x 5 -n COLO829T tests/demo/COLO829T.chm13g.paf.gz | gzip - > COLO829T_chm13g.rsv.gz
minisv extract -n COLO829BL tests/demo/COLO829BL.hs38l.paf.gz | gzip - > COLO829BL.rsv.gz

# joint filtering analysis
minisv isec COLO829T_hs38l.rsv.gz COLO829T_chm13g.rsv.gz | gzip > COLO829T.rsv.gz

# single-sample germline calling
zcat COLO829T_hs38l.rsv.gz | sort -k1,1 -k2,2n -S4G | minisv merge - > COLO829T_germline.msv

# single-sample mosaic calling
zcat COLO829T.rsv.gz | sort -k1,1 -k2,2n -S4G | minisv merge - > COLO829T_mosaic.msv

# tumor-normal paired somatic calling
zcat COLO829{T,BL}.rsv.gz | sort -k1,1 -k2,2n -S4G | minisv merge - | grep COLO829T | grep -v COLO829BL > COLO829T_somatic.msv

# convert to vcf
minisv genvcf COLO829T_germline.msv > COLO829T_germline.vcf
minisv genvcf COLO829T_mosaic.msv > COLO829T_mosaic.vcf
minisv genvcf COLO829T_somatic.msv > COLO829T_somatic.vcf

```


## <a name="intro"></a>Introduction

Minisv is a lightweight mosaic/somatic structural variation (SV) caller for
long genomic reads. Different from other SV callers, it prefers to combine read
alignments against multiple reference genomes. Minisv retains an SV on a read
only if the SV is observed on the read alignments against all references. This
simple strategy reduces alignment errors and filters out germline SVs in the
sample (when used with the assembly of input reads) or in the population (when
used with pangenomes).

Given PacBio HiFi reads at high coverage, minisv
achieves higher specificity for mosaic SV calling, has comparable accuracy to
tumor-normal paired SV callers, and is the only tool accurate enough for
calling large chromosomal alterations (CAs) from a single tumor sample without
matched normal.

`minisv.py` is a equivalent implementation as [minisv.js][minisvjs], and add extra functions.

## Table of Contents

<!--<img aligh="right" width="278" src="doc/example1.png">-->

- [Getting Started](#start)
- [Introduction](#intro)
- [Data preprocessing](#process)
- [Design](#design)
- [Calling SVs](#call-sv)
  - [Germline SVs](#call-germline)
  - [Somatic SVs in tumor-normal pairs](#call-pair)
  - [Large somatic SVs in tumor-only samples](#call-tonly)
  - [Mosaic SVs](#call-mosaic)
- [Generic filtering for other caller](#filter)
- [Trio filtering for sniffles2](#trio_filtering)
- [Ensembling filtered caller results](#ensemble)
- [Comparing SVs](#compare)
- [Filtering interface](#filtering)
- [Limitations](#limit)


## <a name="process"></a>Data preprocessing

Minisv seamlessly works with PAF or GAF (graph PAF) formats. It
**requires** the `ds:Z` tag outputted by [minigraph][mg]-0.21+ or
[minimap2][mm2]-2.28+. For minigraph, use `-cxlr` for long reads. For minimap2,
use `-cxmap-hifi -s50 --ds` for HiFi reads and `-cxlr:hq` for ONT reads. The minimap2 option for Nanopore
reads varies with the read error rate.

Getting the reference genomes, 

```sh
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
wget -c https://zenodo.org/records/6983934/files/chm13-90c.r518.gfa.gz?download=1 -O chm13-90c.r518.gfa.gz
```

### Optional step: de novo assembly of paired normal sample

Generate denovo assembly from the paired normal sample using [hifiasm][hifiasm].

```sh
hifiasm -t32 -o normal_asm normal.fastq.gz
# for ONT
hifiasm --ont -t32 -o normal_asm normal.fastq.gz

gfatools gfa2fa normal_asm.bp.hap1.p_ctg.gfa | bgzip -c - > normal_asm.bp.hap1.p_ctg.fa.gz
gfatools gfa2fa normal_asm.bp.hap2.p_ctg.gfa | bgzip -c - > normal_asm.bp.hap2.p_ctg.fa.gz

```


### Hifi read alignment

```sh
# hg38
minimap2 -cx map-hifi -s50 --ds -t 24 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz tumor.fastq.gz | gzip - > tumor_hg38.paf.gz
minimap2 -cx map-hifi -s50 --ds -t 24 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz normal.fastq.gz | gzip - > normal_hg38.paf.gz
# T2T
minimap2 -cx map-hifi -s50 --ds -t 24 chm13v2.0.fa.gz tumor.fastq.gz | gzip - > tumor_chm13.paf.gz
# denovo assembly
minimap2 -cx map-hifi -s50 --ds -t 24 <(zcat normal_asm.bp.hap1.p_ctg.fa.gz normal_asm.bp.hap2.p_ctg.fa.gz) -I100g --secondary=no tumor.fastq.gz | gzip - > tumor_self.paf.gz

# pangenome graph genome
minigraph -cxlr -t 24 --ds chm13-90c.r518.gfa.gz tumor.fastq.gz | gzip - > tumor_chm13g.paf.gz
```

### ONT read alignment

```sh
# hg38
minimap2 -cxlr:hq --ds -t 24 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz tumor.fastq.gz | gzip - > tumor_hg38.paf.gz
minimap2 -cxlr:hq -s50 --ds -t 24 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz normal.fastq.gz | gzip - > normal_hg38.paf.gz
# T2T
minimap2 -cxlr:hq -s50 --ds -t 24 chm13v2.0.fa.gz tumor.fastq.gz | gzip - > tumor_chm13.paf.gz
# denovo assembly
minimap2 -cxlr:hq -s50 --ds -t 24 <(zcat normal_asm.bp.hap1.p_ctg.fa.gz normal_asm.bp.hap2.p_ctg.fa.gz) -I100g --secondary=no tumor.fastq.gz | gzip - > tumor_self.paf.gz

# pangenome graph genome
minigraph -cxlr -t 24 --ds chm13-90c.r518.gfa.gz tumor.fastq.gz | gzip - > tumor_chm13g.paf.gz
```

## <a name="design"></a>Minisv standalone mode design

Minisv can call 1) germline SVs, 2) somatic SVs in tumor-normal pairs, 3) large
somatic SVs in a single tumor, 4) mosaic SVs and 5) de novo SVs in a trio; the
exact use also depends on the input data types. Minisv achieves this variety of
calling modes in three steps:

1. Extract SVs from read alignment with `extract`. This command processes one
   read at a time, extracts long INDELs or breakends and outputs in a BED-like
   minisv format. If you have tumor-normal pairs, remember to use `-n` to
   specify samples such that you can distinguish the samples in step 3.

2. Intersect candidate SVs extracted from multiple alignment files with `isec`.
   You can skip this step if you use one reference genome only, but you would
   lose the key advantage of minisv. The first two steps are shared by all
   calling modes, though the input alignments may be different. 

3. Call SVs supported by multiple reads using `merge`. This command counts the
   number of supporting reads for each sample. You can filter out SVs supported
   by normal reads to call somatic SVs - this is how minisv works with multiple
   samples jointly.


## <a name="call-sv"></a>Calling SVs

### <a name="call-germline"></a>Germline SVs

```sh
minisv extract -b data/hs38.cen-mask.bed tumor_hg38.paf.gz > sv.hg38l.rsv
cat sv.hg38l.rsv | sort -k1,1 -k2,2n -S4g | minisv merge - > sv.hg38l.msv
minisv genvcf sv.hg38l.msv > sv.hg38l.vcf
```
For calling germline SVs, you only need one linear reference genome. This is the
simplest use case. However, minisv does not infer genotypes. It is not the best
tool for germline SV calling.

### <a name="call-pair"></a>Somatic SVs in tumor-normal pairs

```sh
minisv extract -n TUMOR -b data/hs38.cen-mask.bed tumor_hg38.paf.gz | gzip - > tumor_hs38l.rsv.gz
minisv extract -n TUMOR tumor_chm13.paf.gz | gzip - > tumor_chm13l.rsv.gz
minisv extract -x 5 -n TUMOR tumor_chm13g.paf.gz | gzip - > tumor_chm13g.rsv.gz
minisv extract -q 0 -x 0 -n TUMOR tumor_self.paf.gz | gzip - > tumor_self.rsv.gz

minisv extract -n NORMAL normal_hg38.paf.gz | gzip - > normal.rsv.gz

# joint filtering analysis
minisv isec tumor_hs38l.rsv.gz tumor_chm13l.rsv.gz tumor_chm13g.rsv.gz tumor_self.rsv.gz | gzip > tumor_tgs_filter.rsv.gz

zcat tumor_tgs_filter.rsv.gz normal.rsv.gz | sort -k1,1 -k2,2 -S4g \
  | minisv merge - | grep TUMOR | grep -v NORMAL > sv-paired.hg38l+tgs.msv
```

The last command selects SVs only present in TUMOR but not in NORMAL.

Graph alignment greatly reduces alignment errors when the normal assembly is
not available. When you have the normal assembly, intersecting with graph
alignment can be optional. While graph alignment still improves specificity, it
may affect sensitivity a little - a classical sensitivity vs specificity
problem.

Calling *de novo* SVs will be similar.

### <a name="call-tonly"></a>Large somatic SVs in tumor-only samples

```sh
minisv extract -n TUMOR -b data/hs38.cen-mask.bed tumor_hg38.paf.gz | gzip - > tumor_hs38l.rsv.gz
minisv extract -n TUMOR tumor_chm13.paf.gz | gzip - > tumor_chm13l.rsv.gz
minisv extract -x 5 -n TUMOR tumor_chm13g.paf.gz | gzip - > tumor_chm13g.rsv.gz

# joint filtering analysis
minisv isec tumor_hs38l.rsv.gz tumor_chm13l.rsv.gz tumor_chm13g.rsv.gz | gzip > tumor_tg_filter.rsv.gz

zcat tumor_tg_filter.rsv.gz | sort -k1,1 -k2,2 -S4g | minisv merge - > sv.hg38l+tg.msv
minisv.js view -IC sv.hg38l+tg.msv
```

In the lack of the normal assembly in this case, pangenome graph alignment is
critical for reducing alignment errors and for filtering out common germline
SVs. There may be several hundred SVs called with this procedure. Most of the
called SVs below 10 kb are rare germline events, not somatic. Nevertheless,
most large chromosomal alterations, such as translocations, foldback
inversions (tagged by `foldback` in the output) and events longer 100 kb, are
likely somatic becase such large alterations rarely occur to germline.

### <a name="call-mosaic"></a>Mosaic SVs

Mosaic calling can be applied to a single sample from either normal or tumor sample, take tumor sample as an example below.

```sh
minisv extract -n TUMOR -b data/hs38.cen-mask.bed tumor_hg38.paf.gz | gzip - > tumor_hs38l.rsv.gz
minisv extract -n TUMOR tumor_chm13.paf.gz | gzip - > tumor_chm13l.rsv.gz
minisv extract -x 5 -n TUMOR tumor_chm13g.paf.gz | gzip - > tumor_chm13g.rsv.gz
minisv extract -q 0 -x 0 -n TUMOR tumor_self.paf.gz | gzip - > tumor_self.rsv.gz

minisv isec tumor_hs38l.rsv.gz tumor_chm13l.rsv.gz tumor_chm13g.rsv.gz tumor_self.rsv.gz | gzip > tumor_tgs_filter.rsv.gz
zcat tumor_tgs_filter.rsv.gz | sort -k1,1 -k2,2 -S4g | minisv.js merge -c2 -s0 - > sv.hg38l+tgs.msv
```

Having the phased sample assembly is critical to the calling of small mosaic
SVs. If you do not have the assembly, please perform graph alignment. However,
because there are more small rare SVs than small mosaic SVs, only large mosaic
chromosomal alteration calls are reliable. If you have the assembly, graph
alignment may still help specificity but may hurt sensitivity around [VNTRs][vntr-wiki].

## <a name="filter"></a>Generic filtering for other callers

Downloading example dataset from the [zenodo][msv-zenodo], use `COLO829_hifi1.tar.gz` as an example.

```sh
tar xvfz COLO829_hifi1.tar.gz
cd COLO829_hifi1

# de novo assembly-based filtering of severus
minisv filterasm -c 5 -a -b minisv.py/data/hs38.cen-mask.bed severus/read_ids.csv minisv/COLO829T.self.Q0.gsv.gz severus_filter.stat severus/severus_somatic.vcf > severus_asm.vcf

# de novo assembly-based filtering of sniffles2
minisv filterasm -c 5 -a -b minisv.py/data/hs38.reg.bed sniffles2/grch38_snf_somatic.vcf minisv/COLO829T.self.Q0.gsv.gz snf_filter.stat sniffles2/grch38_snf_somatic.vcf > snf_asm.vcf

# de novo assembly-b5sed filtering of savana
minisv filterasm -c 5 -a -b minisv.py/data/hs38.reg.bed savana/grch38_T_tag.sv_breakpoints_read_support.tsv minisv/COLO829T.self.Q0.gsv.gz savana_filter.stat savana/grch38_T_tag.classified.somatic.vcf > sanava_asm.vcf

# de novo assembly-b5sed filtering of nanomonsv
minisv filterasm -c 5 -a -b minisv.py/data/hs38.reg.bed nanomonsv/grch38_parse.nanomonsv.supporting_read.txt minisv/COLO829T.self.Q0.gsv.gz nanomonsv_filter.stat nanomonsv/grch38_tnpair.vcf > nanomonsv_asm.vcf

# de novo assembly-b5sed filtering of minisv ltg
minisv filterasm -c 5 -a -b minisv.py/data/hs38.reg.bed minisv/COLO829_hifi1T.hg38l+tg.pair-c2s1.msv minisv/COLO829T.self.Q0.gsv.gz msv_filter.stat minisv/COLO829_hifi1T.hg38l+tg.pair-c2s1.msv > msv_filter.msv

```


## <a name="ensemble"></a>Ensembling filtered caller results

The ensemble depends on the `union` function from minisv.js. TODO: add union in minsv.py. This function takes the filtered SV call sets from each caller after denovo assembly-based filtering, then ensemble the SV sets into a in silico truth set. `ensembleunion` collapsed the truth set by selecting one SV per group in the `union` output.

```sh
minisv.js union -p -c 5 -b minisv.py/data/hs38.reg.bed severus_asm.vcf msv_filter.msv sanava_asm.vcf nanomonsv_asm.vcf > ensemble.vcf
minisv ensembleunion ensemble.vcf > ensemble_collapsed.vcf
```


## <a name="compare"></a>Comparing SVs

The `eval` command of minisv compares two or multiple SV callsets. To compare
two callsets:
```sh
minisv eval -b data/hs38.reg.bed -l 100 call1.vcf call2.msv
```
where `-l` specifies the minimum SV length and `-b` specifies confident regions.
The command line outputs TP, FN and FP. Minisv considers two SVs, *S1* and
*S2*, to be the same if both ends of *S1* are within 500 bp from ends of *S2*
and the INDEL types of *S1* and *S2* are the same. Minisv compares all types of
SVs that can be associated with two ends. You can also specify the minimum read
support (`-c`) and the minimum SV length (`-l`) on the command line.

If three or more callsets are given on the command line, minisv will generate an
output like:
```txt
SN  980     0.6091  0.8580  0.6135  0.9104  0.8754  C1
SN  0.0673  110     0.0762  0.1319  0.1119  0.0665  C2
SN  0.8408  0.6727  958     0.6411  0.9216  0.8469  C3
SN  0.1980  0.4000  0.2119  326     0.2313  0.2154  C4
SN  0.7071  0.5818  0.7359  0.6074  536     0.7012  C5
SN  0.8459  0.5636  0.8361  0.6442  0.8731  947     C6
```
where the diagonal gives the count of SVs and the number at row *R* and column
*C* equals to the fraction of calls in *R* found in *C*. In this example,
84.59% of calls in C6 were found in C1 and 87.54% of calls in C1 found in C6.
Generally, higher fraction on a row is correlated with higher sensitivity;
higher fraction on a column is correlated with higher specificity for the
caller on the column.

If you have 3+ callsets, you may also use option `-M` for evaluation in the
consensus mode. In this mode, an FP is a call in callset *C* that is not found
in another *other* callsets. Conversely, an FN is a call that is supported by
two or more *other* callsets but is not called in *C*. Here is sample output:
```txt
RN  845  98   0.1160  C1
RP  980  49   0.0500  C1
RN  994  918  0.9235  C2
RP  110  18   0.1636  C2
...
```
Here the false negative rate (RN) for *C1* is 11.6% and the false positive rate
(RP) is 5.0% based on the definition above.

Another way of SV evaluation is to count the true positive and false negative using `union` result of ensembled SV call set, 

```sh
# k8 javascript
minisv.js union -c 5 -b minisv.py/data/hs38.reg.bed severus_asm.vcf msv_filter.msv sanava_asm.vcf nanomonsv_asm.vcf > ensemble.vcf
# python 
minisv union -c 5 -b minisv.py/data/hs38.reg.bed severus_asm.vcf msv_filter.msv sanava_asm.vcf nanomonsv_asm.vcf > ensemble.vcf
```


This will output the supported SV number for each combination of SV callers, the first column is the binary vector for the appearance of
each caller as in the same order of `union` input, the second caller show the number of SVs, e.g., 1000 row means Severus uniquely reports 
16 SVs.

```txt
1000    16
0100    4
1100    3
0010    10
1010    1
0110    1
1110    2
0001    4
1001    1
0101    2
1101    3
0011    0
1011    2
0111    3
1111    38
```


Minisv can also integrate and ensemble sv from the vcfs and read id files directly:

```sh
minisv advunion -p -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed \
        -i1 output/minisv_puretumor_somatic_asm/snf_HCC1954_hifi1_somatic_generation2.vcf \
        -i1 ../1a.alignment_sv_tools/output/savana12/HCC1954_hifi1/grch38/grch38_T_tag.sv_breakpoints_read_support.tsv \
        -i1 ../1a.alignment_sv_tools/output/severus/HCC1954_hifi1/grch38_cutoff2_read_ids/read_ids.csv \
        -i1 ../1a.alignment_sv_tools/output/nanomonsv/HCC1954_hifi1/T/grch38_parse.nanomonsv.supporting_read.txt \
        -i1 output/msv_somatic/HCC1954_hifi1T.hg38l+tg.pair-c2s1.msv \
        -i2 output/minisv_puretumor_somatic_asm/snf_HCC1954_hifi1_somatic_generation2.vcf \
        -i2 ../1a.alignment_sv_tools/output/savana12/HCC1954_hifi1/grch38/grch38_T_tag.classified.somatic.vcf \
        -i2 ../1a.alignment_sv_tools/output/severus/HCC1954_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf \
        -i2 ../1a.alignment_sv_tools/output/nanomonsv/HCC1954_hifi1/grch38_tnpair.vcf \
        -i2 output/msv_somatic/HCC1954_hifi1T.hg38l+tg.pair-c2s1.msv \
        HCC1954T.self.Q0.gsv.gz > HCC1954_withasmreadids.msv

# keep the SV call with maximum read counts
minisv advunion -u -p -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed \
        -i1 output/minisv_puretumor_somatic_asm/snf_HCC1954_hifi1_somatic_generation2.vcf \
        -i1 ../1a.alignment_sv_tools/output/savana12/HCC1954_hifi1/grch38/grch38_T_tag.sv_breakpoints_read_support.tsv \
        -i1 ../1a.alignment_sv_tools/output/severus/HCC1954_hifi1/grch38_cutoff2_read_ids/read_ids.csv \
        -i1 ../1a.alignment_sv_tools/output/nanomonsv/HCC1954_hifi1/T/grch38_parse.nanomonsv.supporting_read.txt \
        -i1 output/msv_somatic/HCC1954_hifi1T.hg38l+tg.pair-c2s1.msv \
        -i2 output/minisv_puretumor_somatic_asm/snf_HCC1954_hifi1_somatic_generation2.vcf \
        -i2 ../1a.alignment_sv_tools/output/savana12/HCC1954_hifi1/grch38/grch38_T_tag.classified.somatic.vcf \
        -i2 ../1a.alignment_sv_tools/output/severus/HCC1954_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf \
        -i2 ../1a.alignment_sv_tools/output/nanomonsv/HCC1954_hifi1/grch38_tnpair.vcf \
        -i2 output/msv_somatic/HCC1954_hifi1T.hg38l+tg.pair-c2s1.msv \
        HCC1954T.self.Q0.gsv.gz > HCC1954_withasmreadids_collapsed.msv
```


Minisv seamlessly parses the VCF format and the minisv format. It has been
tested with Severus, Sniffles2, cuteSV, SAVANA, SVision-Pro, nanomonsv, SvABA
and GRIPSS.


## <a name="filtering"></a>Filtering interface alone

We first extracted the somatic SVs from the Sniffles2,

```sh
minisv snfpair -n 2 -t 1 grch38_multi.vcf.gz  > snf2.vcf
```

Then, given the GRCh38-based vcf and their corresponding readnames files from the Severus, SAVANA, nanomonsv, Sniffles2, the filtering pipeline sequentially extracts somatic reads from the bam, identify SV signal from realignment to de novo assembly[hifiasm][hifiasm] and pangenome[mg-zenodo][mg-zenodo], and generate consensus SV calls from the Severus, SAVANA, nanomonsv after filtering by multiple references. 

```sh
minisv sv-cross-ref-filter --maskb minisv/data/hg38.cen-mask.bed  --mm2 minimap2_path --mg minigraph_path -b minisv/data/hs38.reg.bed --vcf severus_somatic.vcf savana_grch38_T_tag.classified.somatic.vcf nanomonsv_grch38_tnpair.vcf snf2.vcf --readid_tsv severus_read_ids.csv savana_grch38_T_tag.sv_breakpoints_read_support.tsv nanomonsv_grch38_parse.nanomonsv.supporting_read.txt snf2.vcf grch38.cram grch38.fa asm.bp.hap1.fa.gz asm.bp.hap2.fa.gz CHM13-464.gfa.gz output_folder
```

A typical output folder contains the following output files, the `l+g+s_union_dedup.msv` is the final filtered consensus SV callset, the rest of the intermediate files includes the filtered SV callset and filtering statistics for each caller at different cutoffs, extracted reads(`som_reads.fq.gz`) and intermediate SV signals (`gsv.gz`), realignment results (`paf.gz`, `gaf.gz`):

```sh
denovo_aligned.gaf.gz  l+g_union_dedup.msv    minisv_timings.tsv               nanomonsv_l+g+s_3_filtered.vcf   nanomonsv_l+s_4_filtered.vcf   savana_l+g_5_filtered.stat    savana_l+s_3_filtered.stat   severus_l+g_4_filtered.stat    severus_l+g+s_5_filtered.stat  sniffles2_l+g_3_filtered.stat    sniffles2_l+g+s_4_filtered.stat  sniffles2_l+s_5_filtered.stat
denovo_aligned.paf.gz  l+g_union.msv          nanomonsv_l+g_3_filtered.stat    nanomonsv_l+g+s_4_filtered.stat  nanomonsv_l+s_5_filtered.stat  savana_l+g_5_filtered.vcf     savana_l+s_3_filtered.vcf    severus_l+g_4_filtered.vcf     severus_l+g+s_5_filtered.vcf   sniffles2_l+g_3_filtered.vcf     sniffles2_l+g+s_4_filtered.vcf   sniffles2_l+s_5_filtered.vcf
denovo.gsv.gz          l+g_union_stat.msv     nanomonsv_l+g_3_filtered.vcf     nanomonsv_l+g+s_4_filtered.vcf   nanomonsv_l+s_5_filtered.vcf   savana_l+g+s_3_filtered.stat  savana_l+s_4_filtered.stat   severus_l+g_5_filtered.stat    severus_l+s_3_filtered.stat    sniffles2_l+g_4_filtered.stat    sniffles2_l+g+s_5_filtered.stat  som_reads.fq.gz
graph.gsv.gz           l_only_union.msv       nanomonsv_l+g_4_filtered.stat    nanomonsv_l+g+s_5_filtered.stat  readid.names                   savana_l+g+s_3_filtered.vcf   savana_l+s_4_filtered.vcf    severus_l+g_5_filtered.vcf     severus_l+s_3_filtered.vcf     sniffles2_l+g_4_filtered.vcf     sniffles2_l+g+s_5_filtered.vcf
gs.gsv.gz              l_only_union_stat.msv  nanomonsv_l+g_4_filtered.vcf     nanomonsv_l+g+s_5_filtered.vcf   savana_l+g_3_filtered.stat     savana_l+g+s_4_filtered.stat  savana_l+s_5_filtered.stat   severus_l+g+s_3_filtered.stat  severus_l+s_4_filtered.stat    sniffles2_l+g_5_filtered.stat    sniffles2_l+s_3_filtered.stat
l+g+s_union_dedup.msv  l+s_union_dedup.msv    nanomonsv_l+g_5_filtered.stat    nanomonsv_l+s_3_filtered.stat    savana_l+g_3_filtered.vcf      savana_l+g+s_4_filtered.vcf   savana_l+s_5_filtered.vcf    severus_l+g+s_3_filtered.vcf   severus_l+s_4_filtered.vcf     sniffles2_l+g_5_filtered.vcf     sniffles2_l+s_3_filtered.vcf
l+g+s_union.msv        l+s_union.msv          nanomonsv_l+g_5_filtered.vcf     nanomonsv_l+s_3_filtered.vcf     savana_l+g_4_filtered.stat     savana_l+g+s_5_filtered.stat  severus_l+g_3_filtered.stat  severus_l+g+s_4_filtered.stat  severus_l+s_5_filtered.stat    sniffles2_l+g+s_3_filtered.stat  sniffles2_l+s_4_filtered.stat
l+g+s_union_stat.msv   l+s_union_stat.msv     nanomonsv_l+g+s_3_filtered.stat  nanomonsv_l+s_4_filtered.stat    savana_l+g_4_filtered.vcf      savana_l+g+s_5_filtered.vcf   severus_l+g_3_filtered.vcf   severus_l+g+s_4_filtered.vcf   severus_l+s_5_filtered.vcf     sniffles2_l+g+s_3_filtered.vcf   sniffles2_l+s_4_filtered.vcf
```

## <a name="trio_filtering"></a>Trio WGS filtering interface alone

Suppose we have a trio of family long read WGS, and we generated the diplod denovo assembly for the parents, we could use this interface to reduce false positives based on the dad and mom's denovo assembly

```sh
minisv sv-trio-filter --mm2 minimap2 --mg minigraph --vcf test.vcf --readid_tsv test.vcf test.bam standard_ref.fna dad_hap1.fasta dad_hap2.fasta mom_hap1.fasta mom_hap2.fasta trio_output_filtered
```


## <a name="limit"></a>Limitations

1. Minisv simply counts supporting reads to call an SV. More complex
   algorithms, such as phasing, VNTR reposition and machine learning, may
   improve its accuracy. On the other hand, matching the accuracy of other
   callers with such a simple algorithm implies the potential of our
   multi-reference approach.

2. Minisv does not output genotypes. This makes it less useful for germline SV
   calling. On the HG002 benchmark data, minisv is close to but not as good as
   the best germline SV caller.

3. For the best accuracy, minisv needs the alignment of all reads against
   multiple reference genomes or pangenome graphs. This greatly increases the
   running time. A better strategy is to only align reads with SVs to reduce
   alignment time. It is possible to build a pipeline on top of minisv but this
   has not been implemented.

[minisvjs]: https://github.com/lh3/minisv
[msv-zenodo]: https://zenodo.org/uploads/14715665
[mg-zenodo]: https://zenodo.org/records/6286521
[vntr-wiki]: https://en.wikipedia.org/wiki/Variable_number_tandem_repeat
[mg]: https://github.com/lh3/minigraph
[mm2]: https://github.com/lh3/minimap2
[hifiasm]: https://github.com/chhylp123/hifiasm
