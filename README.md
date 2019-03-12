# Introduction

NECAT is an error correction and assembly tool for Nanopore long noisy reads.

# Installation

We currently only provide executable binaries, which can be downloaded and installed in the following way.
```shell
$ wget https://github.com/xiaochuanle/NECAT/releases/download/v0.01/necat_20190307_linux_amd64.tar.gz
$ tar xzvf necat_20190307_linux_amd64.tar.gz
$ export PATH=$PATH:$(pwd)/NECAT/Linux-amd64/bin
```
After decompression, all the executable files can be found in `NECAT/Linux-amd64/bin`.  The third line above is used for adding `NECAT/Linux-amd64/bin` to the system `PATH`.

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
OVLP_FAST_OPTIONS="-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000"
OVLP_SENSITIVE_OPTIONS="-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000"
CNS_FAST_OPTIONS="-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0"
CNS_SENSITIVE_OPTIONS="-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0"
TRIM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400"
ASM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400"
NUM_ITER=2
CNS_OUTPUT_COVERAGE=45
CLEANUP=0
USE_GRID=false
GRID_NODE=0
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1 --coverage=40"
FSA_ASSEMBLE_OPTIONS=""
FSA_CTG_BRIDGE_OPTIONS="--dump --read2ctg_min_identity=80 --read2ctg_min_coverage=4 --read2ctg_max_overhang=500 --read_min_length=5000 --ctg_min_length=1000 --read2ctg_min_aligned_length=5000 --select_branch=best"
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

`read_list.txt` in the second line above contains the full paths of all read files. It looks like

``` shell
$ cat read_list.txt
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161027_Spenn_001_001_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161101_Spenn_002_002_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161103_Spenn_003_003_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161108_Spenn_004_004_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161108_Spenn_004_005_all.fastq
```

### Step 2: Correct raw reads
Correct the raw noisy reads using the following command:
``` Shell
$ necat.pl correct ecoli_config.txt
```
The corrected reads is `./ecoli/1-consensus/cns_iter${NUM_ITER}/cns.fasta`.   
The extracted longest 45x corrected reads is `./ecoli/1-consensus/cns_final.fasta`

### Step 3: Assemble contigs

After correcting the raw reads, we assemble the contigs using the following command. If the correcting-step is not done, the command  automatically runs the correcting-step first.

```Shell
$ necat.pl assemble ecoli_config.txt
```
The assembled contigs is `./ecoli/4-fsa/contigs.fasta`.

### Step 4: Bridge contigs

After assembling the contigs, we run the bridging-step using the following command. The command  checks and runs the preceding steps first.

```Shell
$ necat.pl bridge ecoli_config.txt
```
The bridged contigs is `./ecoli/6-bridge_contigs/bridged_contigs.fasta`.