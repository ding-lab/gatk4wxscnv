# gatk4wxscnv
Pipeline for WXS CNV using GATK4, based on https://gatkforums.broadinstitute.org/gatk/discussion/9143/how-to-call-somatic-copy-number-variants-using-gatk4-cnv

## Requirements
-Python3.6
-R and packages(`naturalsort`, `getopt`, `optparse`, `DNAcopy`)
-Java1.8
-GATK4
-Bedtools2

## Simple run
```bash
python3 gatk4wxscnv.py --config config.yml
```

`config.yml` file contains all the absolute paths to the softwares and input files. 
```
$ less config.yml
JAVAPATH: '/usr/bin/java'
GATKPATH: '/home/software/gatk-4.beta.5/gatk-package-4.beta.5-local.jar'
referencePath: 'PATH_TO_REFERENCE.fa'
exomeBedPath: 'PATH_TO_EXOME_CAPTURE.bed'
outputDIR: 'PATH_TO_OUTPUTDIR'
normalBamPaths: 'PATH_TO_WXS_NORMAL_BAM_PATH.txt'
cancerBamPaths: 'PATH_TO_WXS_CANCER_BAM_PATH.txt'
```
You can modify the paths in `config.yml` file. 
`normalBamPaths` and `cancerBamPaths` will each contain absolute paths to the BAM files, where each line correspond to a single BAM file. For example,

```
$ less normalBamPaths.txt
/usr/normalbams/normal1.bam
/usr/normalbams/normal2.bam
...
```


## Docker version
```bash
docker build -t gatk4wxscnv .
docker run -it gatk4wxscnv
```

