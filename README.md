# Introduction

NECAT is an error correction and *de-novo* assembly tool for Nanopore long noisy reads.

# Installation

We have sucessfully tested `NECAT` on

* Ubuntu 16.04 (GCC 5.4.0, Perl v5.22.1)
* CentOS 7.3.1611 (GCC 4.8.5, Perl v5.26.2)

If you meet problems in running `NECAT` like
```shell
Syntax error at NECAT/Linuax-amd64/bin/Plgd/Project.pm line 46, near "${cfg{"
```
Please update your `perl` to a newer version (such as v5.26).

There are two ways to install `NECAT`.

### Install from executable binaries

```shell
$ wget https://github.com/xiaochuanle/NECAT/releases/download/SourceCodes20200119/necat_20200119_Linux-amd64.tar.gz
$ tar xzvf necat_20190307_linux_amd64.tar.gz
$ export PATH=$PATH:$(pwd)/NECAT/Linux-amd64/bin
```


### Build from source codes

```shell
$ wget https://github.com/xiaochuanle/NECAT/releases/download/v0.01/necat_20190307_linux_amd64.tar.gz
$ tar xzvf necat_20190307_linux_amd64.tar.gz
$ export PATH=$PATH:$(pwd)/NECAT/Linux-amd64/bin
```

After installation, all the executable files can be found in `NECAT/Linux-amd64/bin`.  The third line above is used for adding `NECAT/Linux-amd64/bin` to the system `PATH`.


# Quick Start

### Step 1: Create a config file

Create a config file template using the following command:

```shell
$ necat.pl config ecoli_config.txt
```

The template looks like

``` shell
PROJECT=
ONT_READ_LIST=
GENOME_SIZE=
THREADS=4
MIN_READ_LENGTH=3000
PREP_OUTPUT_COVERAGE=40
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=30
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true
```
Filling and modifying the relative information, we have

``` shell
PROJECT=ecoli
ONT_READ_LIST=read_list.txt
GENOME_SIZE=4600000
THREADS=20
MIN_READ_LENGTH=3000
  ......
```

`read_list.txt` in the second line above contains the ***full paths*** of all read files. It looks like

``` shell
$ cat read_list.txt
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161027_Spenn_001_001_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161101_Spenn_002_002_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161103_Spenn_003_003_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161108_Spenn_004_004_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161108_Spenn_004_005_all.fastq
```

Please note that files in `read_list.txt` need not be the same format. Each file can independently be either `FASTA` or `FASTQ`, and can further be compressed in [GNU Zip (gzip) format](https://www.gnu.org/software/gzip/manual/gzip.html).

### Step 2: Correct raw reads
Correct the raw noisy reads using the following command:
``` Shell
$ necat.pl correct ecoli_config.txt
```
The pipeline only corrects longest 40X (`PREP_OUTPUT_COVERAGE`) raw reads. The corrected reads are in the files `./ecoli/1-consensus/cns_iter${NUM_ITER}/cns.fasta`.   
The longest 30X (`CNS_OUTPUT_COVERAGE`) corrected reads are extracted for assembly, which are in the file `./ecoli/1-consensus/cns_final.fasta`

### Step 3: Assemble contigs

After correcting the raw reads, we assemble the contigs using the following command. If the correcting-step is not done, the command  automatically runs the correcting-step first.

```Shell
$ necat.pl assemble ecoli_config.txt
```
The assembled contigs are in the file `./ecoli/4-fsa/contigs.fasta`.

### Step 4: Bridge contigs

After assembling the contigs, we run the bridging-step using the following command. The command  checks and runs the preceding steps first.

```Shell
$ necat.pl bridge ecoli_config.txt
```
The bridged contigs are in the file  `./ecoli/6-bridge_contigs/bridged_contigs.fasta`.

If `POLISH_CONTIGS` is set, the pipeline uses the corrected reads to polish the bridged contigs. The polished contigs are in the file `./ecoli/6-bridge_contigs/polished_contigs.fasta`

# Running with multiple computation nodes

On PBS and SGE systems, users may plan to run `NECAT` with multiple computation nodes. This is done by setting the config file (Step 1 of Quick Start) like
```shell
USE_GRID=true
GRID_NODE=4
```
In the above example, `4` computation nodes will be used and each computation node will run with `THREADS` CPU threads.

# Contact

* Chuan-Le Xiao, xiaochuanle@126.com
