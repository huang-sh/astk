### Introduction

**ASTK** is a command line software for comprehensive alternative splicing analysis(AS) analyses including AS event analysis, AS splicing sites sequence feature extraction, AS gene function analysis, potential regulatory mechanism analysis of AS.

---

### Installation

```bash
## create a new conda environment for astk and install python and R
$ conda create -n astk -c conda-forge r-base=4.1 python=3.8 -y
## activate conda environment
$ conda activate astk

## install the development version from github
$ pip install git+https://github.com/huang-sh/astk.git@dev


```

After installed **astk**, you should install **astk**' dependent R packages with:

```bash
$ astk install -r 
...
```

### ASTK docker image

It is recommended to use ASTK docker image because all dependencies have been installed in ASTK docker images.

```
$ docker pull huangshing/astk:latest
```

### Command

**ASTK** works with a command/subcommand structure:

```bash
astk subcommand options
```

**ASTK** provides multiple groups sub-commands for comprehensive AS analysis:

**AS differential splicing analysis**

* **meta**: generate metadata of AS differential splicing analysis contrast groups; it's helpful when you have multiple condition for anlysis
* **generateEvent**: generate AS events using genome GTF annotation
* **generatePsi**: calculates AS events PSI values
* **diffSplice**: run AS differential splicing analysis
* **dsflow**: wrapper of **generateEvent**, **generatePsi** and **diffSplice**

**PSI/dPSI analysis and plot**

* **sigFilter**: select significant AS event
* **psiFilter**: filter AS event with PSI value
* **pca**: PSI value PCA ploting
* **heatmap**: PSI value heatmap ploting
* **volcano**: dPSI value volcano ploting
* **upset**: AS event upset ploting

**alternative exon/intron lenght analysis**

* **lenCluster**: AS event clustering based on alternative exon/intron length
* **lenDist**: alternative exon/intron length distribution plotting
* **lenPick**: selecting specfic exon/intron length of AS event

**gene function enrichment analysis**

* **enrich**: AS gene function enrichment
* **enrichCompare**: AS gene function comparsion
* **gsea**: gene set enrichment analysis
* **nease**: AS events analysis using NEASE

**motif analysis**

* **motifEnrich**：motif enrichment
* **motifFind**：motif discovery
* **motifPlot**: motif plot
* **motifMap**：motit RNA map
* **getmeme**: extract motif from a meme motif file
* **seqlogo**: draw seqLogo figure

**chromatin analysis**

* **signalProfile**：profile chromatin signal of splicing sites

**Eukaryotic Linear Motif**

* **elms**： search Eukaryotic Linear Motifs within amino acid sequence coding by alternative exon DNA sequence.

**AS sites coordinate extract**

* **getcoor**, extract AS site coordiante and generate BED file and fasta file, it also can set AS site upstream or downstream width.

**useful utilities**

* **install**: install other dependent software
* **list**: list 20 orgnism annotation OrgDb

### ASKT docker iamge Usage

We could create a shortcut for the docker command with alias command. It is convenient for us to run the docker version of ASTK multiple times.

```bash
$ alias astkdocker="docker run --rm -v /home/test/project:/project -e MY_USER=$(id -u) huangshing/astk"
```

Please replace `/home/test/project`  with your path. This directory should contain some reference files and all files you need to analyze.

```
$ ll -h | cut -d " " -f 5-

  446M Aug 12 21:08 ATAC.e16.5.fb.bigwig
    27 Aug  8 16:40 data
  1.4G Aug  8 22:03 gencode.v38.annotation.gtf
  847M Aug  8 17:35 gencode.vM25.annotation.gtf
  2.6G Aug  8 17:35 GRCm38.primary_assembly.genome.fa

```

And then we just run astk like:

```
$ astkdocker astk meta -h
Usage: astk meta [OPTIONS]

  generate metadata template file

Options:
  -p1, --control PATH           file path for condtion 1  [required]
  -p2, --treatment PATH         file path for condtion 2  [required]
  -gn, --groupName TEXT         group name
  -repN, --replicate INTEGER    replicate, number
  -o, --output PATH             metadata output path  [required]
  -repN1, --replicate1 INTEGER  replicate1, number
  -repN2, --replicate2 INTEGER  replicate2, number
  --condition TEXT              condition name
  -fn, --filename               file name
  --split TEXT                  name split symbol and index
  -h, --help                    Show this message and exit.


```

### Usage

#### Prepare

transcript TPM quantitative file is required for **astk** input. In the example analysis , we use RNA-seq data of mouse forebrain, midbrain and hindbrain RNA-Seq data from embryonic 11.5-16.5 day and p0. And **salmon** was used for transcript TPM.

you can download data from [encode forbrain](https://www.encodeproject.org/search/?type=Experiment&control_type!=*&assay_term_name=polyA%20plus%20RNA-seq&replicates.library.biosample.donor.organism.scientific_name=Mus%20musculus&biosample_ontology.term_name=forebrain&status=released)

**exammple**
According to the [tutorial](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)

```
$ salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
$ salmon quant -i salmon_index --gcBias --useVBOpt --seqBias  -r ENCFF329ACL.fastq.gz --validateMappings -o fb_e11.5_rep1 -p 30
$ ls fb_e11.5_rep1
aux_info  cmd_info.json  lib_format_counts.json  libParams  logs  quant.sf
```

> NOTE: --gencode is required.

And we get following output including quantification data of fb embryonic development stages

```
$ ls data/quant/fb*/*quant.sf

data/quant/fb_e10.5_rep1/quant.sf  data/quant/fb_e12.5_rep1/quant.sf  data/quant/fb_e14.5_rep1/quant.sf  data/quant/fb_e16.5_rep1/quant.sf
data/quant/fb_e10.5_rep2/quant.sf  data/quant/fb_e12.5_rep2/quant.sf  data/quant/fb_e14.5_rep2/quant.sf  data/quant/fb_e16.5_rep2/quant.sf
data/quant/fb_e11.5_rep1/quant.sf  data/quant/fb_e13.5_rep1/quant.sf  data/quant/fb_e15.5_rep1/quant.sf  data/quant/fb_p0_rep1/quant.sf
data/quant/fb_e11.5_rep2/quant.sf  data/quant/fb_e13.5_rep2/quant.sf  data/quant/fb_e15.5_rep2/quant.sf  data/quant/fb_p0_rep2/quant.sf

```

#### meta

**meta** is used to generate contrast group metadata table for AS differential splicing analysis.

For example, we can generate multiple developmental stage contrast groups with e11.5 stage as control:

```bash
$ mkdir metadata -p
$ astk meta -o metadata/fb_e11_based -repN 2 \
    -c1 data/quant/fb_e11.5_rep*/quant.sf \
    -c2 data/quant/fb_e1[2-6].5_rep*/quant.sf  data/quant/fb_p0_rep*/quant.sf \
    -gn fb_e11_12 fb_e11_13 fb_e11_14 fb_e11_15 fb_e11_16 fb_e11_p0

```

**meta** arguments

* -o:  output file path
* -repN:  number of replicate samples
* -c1: condition 1(ctrl) sample transcript quantification files path
* -c2: condition 2(case) sample transcript quantification files path
* -gn: group names

the output of **meta** is a CSV file and JSON file. CSV file is convenient for viewing in excel, and JSON file will be used in other sub-commands.

![fb_e11_meta.csv](demo/img/fb_e11_based_metadata.png)

#### dsflow

**dsflow** is wrapper of **generateEvent**, **generatePsi** and **diffSplice**. It's used to simplify the differential splicing analysis workflow. And above 7 contrast groups differential splicing analysis could be performed with one command line:

```bash
$ mkdir result
$ astk dsflow -od result/fb_e11_based -md metadata/fb_e11_based.json \
    -gtf gencode.vM25.annotation.gtf -t ALL &

$ ls result/fb_e11_based
dpsi  psi  ref  sig01  tpm
```

**dsflow** arguments:

* -od: output directory
* -md: meta data, the meta output json file
* -gtf: genome annotation GTF file
* -t: alternative splicing type, ALL is for all supported types.

the output of **dsflow**  contain four directories:

* ref is the directory including AS event reference annotation files
* tpm is the directory including sample TPM files
* psi is the directory including AS event PSI file
* dpsi  is the directory including differential splicing result
* sig01 is the directory including result that filter using pval < 0.05 and |PSI| > 0.1

#### psiFilter

**psiFilter** is used to filter AS event with PSI value.

```bash
$ astk pf -i result/fb_e11_based/psi/fb_e11_p0_SE_c2.psi \
    -psi 0.8 -o result/fb_e11_based/psi/fb_e11_p0_SE_c2_08.psi
$ astk pf -i result/fb_e11_based/psi/fb_e11_p0_SE_c2.psi \
    -psi -0.2 -o result/fb_e11_based/psi/fb_e11_p0_SE_c2_02.psi
```

**psiFilter** arguments

* FILE: input psi file
* -psi: when option value > 0, it denotes that select AS events that PSI > option value; however, when option value < 0, it denotes that select AS events that PSI < abs(option value)
* -o: output file

#### pca

**pca** sub-commands is using for PCA analysis of PSI.

```bash
$ astk pca -i result/fb_e11_based/psi/fb_e11_12_SE_c1.psi \
    result/fb_e11_based/psi/fb_e11_1[2-6]_SE_c2.psi \
    result/fb_e11_based/psi/fb_e11_p0_SE_c2.psi \
    -o img/fb_pca.pdf -w 6 --height 4 \
    -gn fb_e11.5 fb_e12.5 fb_e13.5 fb_e14.5 fb_e15.5 fb_e16.5 fb_p0

```

**pca** arguments

* -i : input PSI files
* -o: output figure path
* --width: figure width
* --height: figure height
* -gn: group names for different file data labe

![fb_pca.png](demo/img/fb_pca.png)

#### heatmap

**heatmap** is used for ploting heatmap of PSI. **hm** is short alias of **heatmap**.

```bash
$ astk hm -i result/fb_e11_based/sig01/psi/fb_e11_12_SE_c1.sig.psi \
    result/fb_e11_based/sig01/psi/fb_e11_1*_SE_c2.sig.psi \
    -o img/fb_hm.png -fmt png

```

**heatmap** arguments

* -i : PSI files
* -o: output figure
* -o : output path

![fb_hm.png](demo/img/fb_hm.png)

#### barplot

**barplot** is used to draw barplot figure that show AS events counts distribution withing different condition.

```bash
astk barplot -i result/fb_e11_based/sig01/fb_e11_p0_*.sig.dpsi \
    -o img/fb_e11_p0_bar.png -dg -xl A3 A5 AF AL MX RI SE

```

**barplot** arguments:

* -i : input dpsi files, support multiple files
* -o: output figure
* -dg: AS events can be divided into two groups based on dPSI values  (group +: dPSI > 0, group -: dPSI < 0)
* -xl: x labels

![fb_e11_p0_bar.png](demo/img/fb_e11_p0_bar.png)

#### volcano

**volcano** is used for dPSI volcano ploting.  **vol** is short alias of **volcano**.

```bash
$ astk volcano -i result/fb_e11_based/dpsi/fb_e11_p0_SE.dpsi \
    -o img/fb_e11_p0_SE.vol.png 

```

**vol** arguments:

* -i : input dpsi files, support multiple files
* -o: output directory

![fb_hm.png](demo/img/fb_e11_p0_SE.vol.png)

#### upset

**upset** is used for dPSI upset ploting.  **vol** is short alias of **upset**.

```bash
$ astk upset -i result/fb_e11_based/sig01/fb_e11_12_SE.sig.dpsi \
    result/fb_e11_based/sig01/fb_e11_14_SE.sig.dpsi \
    result/fb_e11_based/sig01/fb_e11_16_SE.sig.dpsi \
    -o img/fb_upset.png -xl e11_12 e11_14 e11_16 

```

**upset** arguments:

* -i : input dpsi files, support multiple files
* -o: output file
* -xl: x labels

![fb_upset.png](demo/img/fb_upset.png)

#### lenCluster

**lenCluster** provide a function for cluster AS events based on alternative exon/intron  length. **lc** is short alias of **lenCluster**.

```bash
$ astk lc -i result/fb_e11_based/sig01/*dpsi -lr 1 51 251 1001 \
    -od result/fb_e11_based/lenc

```

**lc** arguments:

* -i: input dpsi files
* -lr: length range
* -od: output directory
* -fmt: figure format

#### enrich

**enrich** is used for genes GO term enrichment enrichment. GO term enrichment map networks and enrichment clustering are provided.

```bash
$ mkdir img/enrich -p
$ astk enrich -i result/fb_e11_based/sig01/fb_e11_13_SE.sig.dpsi \
    -ont BP -qval 0.1 -orgdb mm  -fmt png \
    -od img/enrich/fb_e11_13_SE 

```

**enrich** arguments:

* -i: dpsi file
* -ont:  ontology
* -od : output directory
* -qval : q-value
* -org : organism, for example, ‘hs’ for human, ‘mm’ for mouse

GO terms enrichment result and enrichment clustering  have figure and text formats.

#### enrichCompare

**enrichCompare** is used for  gene functional characteristics comparison of different AS genes, short alias: **ecmp**.

Comparison between The dpsi > 0.1 and dpsi < 0.1 in the 7 group (fb_16.5 vs fb_p0)

```bash
$ mkdir img/ecmp
$ astk ecmp -i result/fb_e11_based/lenc/*/fb_e11_12_SE.sig.dpsi \
     -ont BP -org mm  -fmt png \
     -od img/enrich/fb_e11_12_SE_lc
```

**enrichCompare** arguments:

* -i : dpsi files
* -od : output directory
* -ont : ontology
* -qval : q-value
* -org : organism, for example, ‘hs’ for human, ‘mm’ for mouse

![GO.cmp.BP.png](demo/img/enrich/fb_e11_12_SE_lc/GO.cmp.BP.qval0.1_pval0.1.png)

#### motifEnrich

**motifEnrich** is used for performing motif enrichment within splicing sites flanking sequence using RBP motif database. me is short alias.

```bash
$ astk me -te result/fb_e11_based/psi/fb_e11_p0_SE_c2_08.psi \
    -ce result/fb_e11_based/psi/fb_e11_p0_SE_c2_02.psi \
    -od img/motif/fb_e11_p0_SE_me -org mm \
    -fi GRCm38.primary_assembly.genome.fa

```

Arguments:

* -te: input treatment event file
* -ce: input control event file
* -od: output directory
* -org: organism
* -fi: genome fasta, need index

#### motifFind

**motifFind** is used for performing motif discovery and the compared to known RBP motif. mf is short alias.

```bash
astk mf -te result/fb_e11_based/psi/fb_e11_p0_SE_c2_08.psi \
    -od img/motif/fb_e11_p0_SE_mf -org mm \
    -fi GRCm38.primary_assembly.genome.fa
```

Arguments:

* -te: input treatment event file
* -od: output directory
* -org: organism
* -fi: genome fasta, need index

#### getmeme

**getmeme** is used for querying ASTK built-in motif data.

```bash
astk getmeme M316_0.6 M083_0.6 -db CISBP-RNA \
    -org mm -o img/motif/query.meme
```

Arguments:

* MOTIFID...: input motif IDs
* -db: motif database
* -org: organism
* -o: output file path
* 

#### motifPlot

**motifPlot** is used for drawing motif figure using motif meme data

```bash
astk mp -mi M083_0.6 -db CISBP-RNA -org mm \
 -o img/motif/M083_0.6_plot.png -w 10

```

Arguments:

* MOTIFID...: input motif IDs
* -db: motif database
* -org: organism
* -o: output file path

![motif_plot.png](demo/img/motif/M083_0.6_plot.png)

#### mmap

**mmap** is used for generating motif map to show motif distribution.

```bash
astk mmap -e result/fb_e11_based/psi/fb_e11_p0_SE_c2_08.psi \
    -n a1 a2 a3 a4 -c 150 150 150 150 \
    -mm img/motif/query.meme -od img/motif/motif_map \
    -fi GRCm38.primary_assembly.genome.fa
```

Arguments:

* -fa: input fasta files
* -n: fasta files names
* -c: center positions
* -mm: motif meme file
* -od: output directory
* -fi: genome fasta, need index

#### signalProfile

**signalProfile** is used to profile chromatin signal of splicing sites flank.

```bash
astk pf -i result/fb_e11_based/psi/fb_e11_16_AF_c2.psi \
    -psi 0.8 -o result/fb_e11_based/psi/fb_16_AF_08.psi

astk signalProfile -o img/fb_16_AF_low_ATAC.png \
    -e result/fb_e11_based/psi/fb_16_AF_08.psi \
    -bw ATAC.e16.5.fb.bigwig \
    -ssl A1 A2 A3 A4 A5 -fmt png

```

Arguments:

* -o: output file
* -e: AS event file that including  AS event ID
* -bw: bigwig file
* -ssl: splicing site labels
* -fmt: figure format

![fb_16_AF_low_ATAC.png](demo/img/fb_16_AF_low_ATAC.png)

#### elms

**elms** is using for searching Eukaryotic Linear Motifs within amino acid sequence coding by alternative exon sequence.

```bash
astk elms -i result/fb_e11_based/sig01/fb_e11_p0_SE.sig.dpsi -g mm10 -o img/elm.csv
```

Arguments:

* -i: dpsi file path
* -g: genome assembly
* -o: output

### FAQs

**ERROR: compilation failed for package ‘magick’**

```bash
$ astk install -r
ERROR: compilation failed for package ‘magick’
* removing ‘/home/user/software/anaconda/envs/astk/lib/R/library/magick’
```

You could install the software manually,  and the re-run `astk install -r`. If you use conda, you can do it like:

```bash
$ conda install -c conda-forge r-magick -y
...
$ astk install -r
```

**namespace ‘rlang’ 0.4.12 is already loaded, but >= 1.0.0 is required**

```bash
Error in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]) : 
  namespace ‘rlang’ 0.4.12 is already loaded, but >= 1.0.0 is required
Calls: source ... asNamespace -> loadNamespace -> namespaceImport -> loadNamespace
Execution halted
```

You could install the rlang to higher version manually. For example:

```
# R console
> install.packages("rlang")
```
