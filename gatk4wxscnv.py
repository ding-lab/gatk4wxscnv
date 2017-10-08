import os
import sys
import subprocess
import shlex
import configparser
import argparse
import re
import pysam
import yaml

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', default='config.yml')
    args = vars(parser.parse_args())

    config = args['config']
    return config

def filePathParser(inputFile):
    '''input file is a list of either cancer or normal BAM Paths
    reads this file and creates a list of BAM paths to be processed'''
    fileList = []
    with open(inputFile, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                filePath = line.strip('\n')
                fileList.append(filePath)


    return fileList

def sampleNameBam(bamFile):
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name

def execute(cmdString):
    '''uses subprocess to execute a command, waits until command is done'''
    execution = subprocess.Popen(shlex.split(cmdString))
    execution.wait()
    return 0

def formatBedChromosome(bedPath):
    '''this either removes 'chr' string from bed file chromosome or adds'''
    newBedFile = re.sub(string=bedPath, pattern=r'.bed$', repl='.new.bed')
    with open(newBedFile, 'w') as g:
        linecount=0
        with open(bedPath, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    chromosome, start, end, *args = line.strip('\n').split('\t')
                    chromosomeString = re.sub(string=chromosome, pattern=r'chr', repl='')
                    g.write('\t'.join([chromosomeString, start, end]) + '\t' + str(linecount) + '\n')
                    linecount+=1
    return newBedFile

def mergeBed(exomeBedPath):
    '''convertBedToTargetFile does not work if BED file has overlapping regions
    thus using BEDTOOL's mergeBed command create a mergedBed file'''
    mergedExomeBedPath = re.sub(string=exomeBedPath, pattern=r'.bed$', repl='.merged.bed')
    if os.path.isfile(mergedExomeBedPath):
        pass
    else:
        with open(mergedExomeBedPath, 'w') as f:
            mergeCmd = f'mergeBed -i {exomeBedPath}'
            print(mergeCmd)
            subprocess.call(shlex.split(mergeCmd), stdout=f)

        f.close()
    return mergedExomeBedPath

def convertBedToTargetFile(exomeMergedBedPath, JAVAPATH, GATKPATH):
    '''converts exome capture BED file (0-based) into GATK target format (1-based)'''
    exomeGatkTarget = exomeMergedBedPath + '.gatkTarget.tsv'
    if os.path.isfile(exomeGatkTarget):
        pass
    else:
        cmd = f'{JAVAPATH} -jar -Xmx4g {GATKPATH} ConvertBedToTargetFile --input {exomeMergedBedPath} --output {exomeGatkTarget}'
        execute(cmd)

    return exomeGatkTarget


def configure(configPath):
    '''configures the path to softwares and input files from a yaml file'''
    with open(configPath, 'r') as f:
        pathMaps = yaml.safe_load(f)

    try:
        JAVAPATH = pathMaps['JAVAPATH']
        GATKPATH = pathMaps['GATKPATH']
        referencePath= pathMaps['referencePath']
        exomeBedPath = pathMaps['exomeBedPath']
        outputDIR = pathMaps['outputDIR']
        normalBamPaths = pathMaps['normalBamPaths']
        cancerBamPaths =pathMaps['cancerBamPaths']
    except KeyError:
        print('one of the required input path not specified. Exiting...')

    return JAVAPATH, GATKPATH, referencePath, exomeBedPath, outputDIR, normalBamPaths, cancerBamPaths


def calculateTargetCoverage(bamPath, exomeGatkTargetPath, outputDIR, JAVAPATH, GATKPATH):
    '''this step may require parallelization if access to LSF is possible'''
    sampleName = sampleNameBam(bamPath)
    outputCovPath = outputDIR + sampleName + '.pcov'
    if os.path.isfile(outputCovPath):
        print(f'{outputCovPath} exists...skipping...')
    else:
        cmd = f'{JAVAPATH} -jar -Xmx4g {GATKPATH} CalculateTargetCoverage -I {bamPath} -T {exomeGatkTargetPath} -transform PCOV -groupBy SAMPLE -targetInfo FULL -O {outputCovPath}'
        execute(cmd)
    return outputCovPath

def panelOfNormalList(normalOutputDir):
    '''writes the absolute path of output from calculateTargetCoverage of normal samples '''
    normal_wxs_pcovPaths = normalOutputDir + 'normal_wxs_pcovPaths.txt'
    print(normal_wxs_pcovPaths)
    if os.path.isfile(normal_wxs_pcovPaths):
        pass
    else:
        with open(normal_wxs_pcovPaths, 'w') as f:
            for file in os.listdir(normalOutputDir):
                if re.search(string=file, pattern=r'.pcov$') != None:
                    f.write(normalOutputDir + file + '\n')

    return normal_wxs_pcovPaths

def combineReadCounts(normal_wxs_pcovPaths, outputDIR, JAVAPATH, GATKPATH):
    '''combines results of calculateTargetCoverage of normal samples to create PoN'''
    combined_normals = outputDIR + 'combined-normals.txt'
    if os.path.isfile(combined_normals):
        pass
    else:
        cmd = f'{JAVAPATH} -jar -Xmx4g {GATKPATH} CombineReadCounts -inputList {normal_wxs_pcovPaths} -O {combined_normals}'
        execute(cmd)
    return combined_normals

def createPanelOfNormals(combinedNormal, outputPoN, JAVAPATH, GATKPATH, minimumTargetFactorPercentileThreshold=5):
    '''creates PoN
    this step requires increasing memory with increasing number of normal bams
    example: with ~800 bam -> -R "rusage[mem=24576]" -M 30000000 and -Xmx24g'''
    if os.path.isfile(outputPoN):
        pass
    else:
        cmd = f'{JAVAPATH} -jar -Xmx24g {GATKPATH} CreatePanelOfNormals -I {combinedNormal} -O {outputPoN} --disableSpark --minimumTargetFactorPercentileThreshold {minimumTargetFactorPercentileThreshold} --outputFailedSamples failedSamples.txt --noQC false --maximumColumnZerosPercentage 10'
        execute(cmd)
    return outputPoN

def normalizeSomaticReadCounts(cancerCoverage, pon, ptn, tn, JAVAPATH, GATKPATH):
    '''create normalized somatic readcounts'''
    if os.path.isfile(ptn) and os.path.isfile(tn):
        pass
    else:
        cmd = f'{JAVAPATH} -jar -Xmx4g {GATKPATH} NormalizeSomaticReadCounts -I {cancerCoverage} -PON {pon} -PTN {ptn} -TN {tn}'
        print(cmd)
        execute(cmd)
    return ptn, tn

def performSegmentation(tn, seg, JAVAPATH, GATKPATH):
    '''perform segmantation'''
    if os.path.isfile(seg):
        pass
    else:
        cmd = f'{JAVAPATH} -jar -Xmx4g {GATKPATH} PerformSegmentation -TN {tn} -O {seg} -LOG'
        print(cmd)
        execute(cmd)
    return seg

def createSequenceDictionary(referencePath, JAVAPATH, GATKPATH):
    '''creates sequence dictionary of reference file if it does not exist'''
    seqDictionary = referencePath + '.dict'
    if os.path.isfile(seqDictionary):
        print(f'{seqDictionary} exists already, pass...')
    else:
        cmd = f'{JAVAPATH} -jar -Xmx4g {GATKPATH} CreateSequenceDictionary -R {referencePath} -O {seqDictionary}'
        execute(cmd)
    return seqDictionary

def plotSegmentedCopyRatio(ptn, tn, seg, sequenceDictionary, outputDIR, samplePrefix, JAVAPATH, GATKPATH):
    '''creates optional plot'''

    cmd = f'{JAVAPATH} -jar -Xmx4g {GATKPATH} PlotSegmentedCopyRatio -PTN {ptn} -TN {tn} -S {seg} -SD {sequenceDictionary} -O {outputDIR} -pre {samplePrefix} -LOG '
    print(cmd)
    execute(cmd)

    return 0

def callSegments(tn, seg, cnvcall, JAVAPATH, GATKPATH):
    '''creates final text output of segmented outputfile'''
    if os.path.isfile(cnvcall):
        pass
    else:
        cmd = f'{JAVAPATH} -jar -Xmx4g {GATKPATH} CallSegments -TN {tn} -S {seg} -O {cnvcall}'
        print(cmd)
        execute(cmd)
    return cnvcall


def somatic_CNV(cancerCoverage, pon, sequenceDictionary, outputDIR, JAVAPATH, GATKPATH):
    '''high level wrapper to call somatic CNV once panel of normal has been created'''
    samplePrefix = re.sub(string=os.path.basename(cancerCoverage), pattern=r'.pcov', repl='')
    ptn = outputDIR + samplePrefix + '.ptn'
    tn = outputDIR + samplePrefix + '.tn'
    seg = outputDIR + samplePrefix + '.seg'
    cnvcall = outputDIR + samplePrefix + '.cnv'
    if os.path.isfile(cnvcall):
        pass
    else:
        normalizeSomaticReadCounts(cancerCoverage, pon, ptn, tn, JAVAPATH, GATKPATH)
        performSegmentation(tn, seg, JAVAPATH, GATKPATH)
        plotSegmentedCopyRatio(ptn, tn, seg, sequenceDictionary, outputDIR, samplePrefix, JAVAPATH, GATKPATH)
        callSegments(tn, seg, cnvcall, JAVAPATH, GATKPATH)

    return cnvcall

def main():
    configPath = argument_parser()

    # configure file paths

    JAVAPATH, GATKPATH, referencePath, exomeBedPath, outputDIR, normalBamPaths, cancerBamPaths = configure(configPath)
    # configure output paths
    tmpCancerDir = outputDIR + 'cancer/'
    tmpNormalDir = outputDIR + 'normal/'

    os.makedirs(tmpCancerDir, exist_ok=True)
    os.makedirs(tmpNormalDir, exist_ok=True)

    # read in all the bam files to be processed into a list
    wxs_normal_absPaths = filePathParser(normalBamPaths)
    wxs_cancer_absPaths = filePathParser(cancerBamPaths)

    # process BED file into proper GATK target tsv file
    mergedExomeBedPath = mergeBed(exomeBedPath)
    newMergedExomeBedPath = formatBedChromosome(mergedExomeBedPath)
    exomeGatkTargetPath = convertBedToTargetFile(newMergedExomeBedPath,JAVAPATH, GATKPATH)
    print(exomeGatkTargetPath)

    # create GATK dictionary file of reference fasta, if not already present
    sequenceDictionary = createSequenceDictionary(referencePath, JAVAPATH, GATKPATH)

    normal_pcov_list = []
    cancer_pcov_list = []

    # calculate target coverage for normal
    for normalBam in wxs_normal_absPaths:
        normal_pcov_list.append(calculateTargetCoverage(normalBam, exomeGatkTargetPath, tmpNormalDir, JAVAPATH, GATKPATH))

    # calculate target coverage for cancer
    for cancerBam in wxs_cancer_absPaths:
        cancer_pcov_list.append(calculateTargetCoverage(cancerBam, exomeGatkTargetPath, tmpCancerDir, JAVAPATH, GATKPATH))


    # create panel of normal
    normal_wxs_pcovPaths = panelOfNormalList(tmpNormalDir)
    combinedNormal = combineReadCounts(normal_wxs_pcovPaths, tmpNormalDir, JAVAPATH, GATKPATH)
    pon = createPanelOfNormals(combinedNormal, outputDIR + 'normals.pon', JAVAPATH, GATKPATH, 5)

    # Now calculate somatic CNV for all cancer coverages against panel of normal
    for cancerCoverage in cancer_pcov_list:
        somatic_CNV(cancerCoverage, pon, sequenceDictionary, outputDIR, JAVAPATH, GATKPATH)

    return 0


if __name__ == '__main__':
    # argument parsing
    main()
