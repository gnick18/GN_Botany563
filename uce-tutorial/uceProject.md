# Installation

First I installed the program by using conda and the command conda install phyluce.

- This gave an error because my default python environent is 3+
  -  program requires 2.7 to work

Created a new conda environment using the following command

```shell
conda create -n "phyluce" python=2.7 ipython
```

I then installed conda with `conda install phyluce` and it worked perfectly.

# Tutorial 1 Walkthrough

I followed their code line for line to install and unzip the data

```shell
# create a project directory
mkdir uce-tutorial

# change to that directory
cd uce-tutorial

# download the data into a file names fastq.zip
wget -O fastq.zip https://ndownloader.figshare.com/articles/1284521/versions/1

# make a directory to hold the data
mkdir raw-fastq

# move the zip file into that directory
mv fastq.zip raw-fastq

# move into the directory we just created
cd raw-fastq

# unzip the fastq data
unzip fastq.zip

# delete the zip file
rm fastq.zip

# you should see 6 files in this directory now
ls -l
```

1. Counting the read data

`for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done`

This produced the following.

```shell
Alligator_mississippiensis_GGAGCTATGG_L001_R1_001.fastq.gz
1750000
Anolis_carolinensis_GGCGAAGGTT_L001_R1_001.fastq.gz
1874362
Gallus_gallus_TTCTCCTTCA_L001_R1_001.fastq.gz
376559
Mus_musculus_CTACAACGGC_L001_R1_001.fastq.gz
1298196

```

2. Cleaning the read data

The data is just raw untrimmed fastq data. I need to remove the adapter contamination and low quality bases. 

Following the tutorial I inslaled illumiprocessor with `conda install illumiprocessor` and ran the following command. I also created a conf file to pass into the program.

```shell
# this is the section where you list the adapters you used.  the asterisk
# will be replaced with the appropriate index for the sample.
[adapters]
i7:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

# this is the list of indexes we used
[tag sequences]
BFIDT-166:GGAGCTATGG
BFIDT-016:GGCGAAGGTT
BFIDT-045:TTCTCCTTCA
BFIDT-011:CTACAACGGC

# this is how each index maps to each set of reads
[tag map]
Alligator_mississippiensis_GGAGCTATGG:BFIDT-166
Anolis_carolinensis_GGCGAAGGTT:BFIDT-016
Gallus_gallus_TTCTCCTTCA:BFIDT-045
Mus_musculus_CTACAACGGC:BFIDT-011

# we want to rename our read files something a bit more nice - so we will
# rename Alligator_mississippiensis_GGAGCTATGG to alligator_mississippiensis
[names]
Alligator_mississippiensis_GGAGCTATGG:alligator_mississippiensis
Anolis_carolinensis_GGCGAAGGTT:anolis_carolinensis
Gallus_gallus_TTCTCCTTCA:gallus_gallus
Mus_musculus_CTACAACGGC:mus_musculus
```



```shell
illumiprocessor     --input raw-fastq/     --output clean-fastq     --config illumiprocessor.conf     --cores 2
```

This created a log file, a directory with the clean data. 

3. Quality control

I ran the same command as described in the tutorial and got the following output.

```shell
# run this script against all directories of reads

for i in *;
do
    phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/ --csv;
done

#OUTPUT
All files in dir with alligator_mississippiensis-READ2.fastq.gz,3279362,294890805,89.9232243955,0.0100399681299,40,100,100.0
All files in dir with anolis_carolinensis-READ-singleton.fastq.gz,3456457,314839345,91.0873026917,0.00799863728967,40,100,100.0
All files in dir with gallus_gallus-READ2.fastq.gz,749026,159690692,213.197795537,0.0588973605565,40,251,250.0
All files in dir with mus_musculus-READ-singleton.fastq.gz,2332785,211828511,90.8049867433,0.0102813002698,40,100,100.0
```

4. Assemble the data 

We will be using trinity to assemble the data.

First I need to make a configuration file for the assembly.

* After trying to do trinity as the tutorial laid out I noticed that this is no longer workable in the Mac OS. I am attempting to use an alternative assembly program.

After fighting with it I noticed that the assembly config file will work with other assembly programs in phyluce. 

My config file:

```shell
[samples]
alligator_mississippiensis:/Users/gnickles/Documents/GN_Botany563/uce-tutorial/clean-fastq/alligator_mississippiensis/split-adapter-quality-trimmed/
anolis_carolinensis:/Users/gnickles/Documents/GN_Botany563/uce-tutorial/clean-fastq/anolis_carolinensis/split-adapter-quality-trimmed/
gallus_gallus:/Users/gnickles/Documents/GN_Botany563/uce-tutorial/clean-fastq/gallus_gallus/split-adapter-quality-trimmed/
mus_musculus:/Users/gnickles/Documents/GN_Botany563/uce-tutorial/clean-fastq/mus_musculus/split-adapter-quality-trimmed/
```

The command I used to run abyss:

```shell
phyluce_assembly_assemblo_abyss     --conf assembly.conf     --output trinity-assemblies     --clean
```

- This is going to take a bit to run as it is on my personal computer


