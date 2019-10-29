#!/usr/bin/python
# Copyright @2001-2019 Python Software Foundation

# Author: Jasper Boom

# Prequisites:
# - sudo apt-get install python
# - sudo apt-get install python-pip
# - sudo pip install pandas

# Imports
import os
import sys
import argparse
import re
import pandas as pd

# The getNewCsv function.
# This function creates a tabular output file which contains 4 columns. These
# columns are a umi number, a unique umi, a count of the number of reads
# belonging to that unique umi and a representing umi read.
# The function creates a dataframe with the desired column name. It then adds
# a row for every unique umi based on the above mentioned columns. The
# "UMI Number" column is used as index and the tabular file is created.
def getNewCsv(dicUniqueUmi, dicUniqueUmiCount, flTabular):
    dfOutput = pd.DataFrame(columns=["UMI Number", "UMI", "Read Count",
                                     "UMI Read"])
    intCount = 0
    intNumber = 1
    for strCode in dicUniqueUmi:
        dfOutput.loc[intCount] = [intNumber, strCode,
                                  dicUniqueUmiCount[strCode],
                                  dicUniqueUmi[strCode]]
        intCount += 1
        intNumber += 1
    dfOutput = dfOutput.set_index("UMI Number")
    dfOutput.to_csv(flTabular, sep="\t", encoding="utf-8")

# The getFastaFiles function.
# This function creates separate fasta files for every unique umi. The function
# creates a unique name for every umi file and combines that with the desired
# output path. A file is opened or created based on this combination and a
# read header and the read itself are appended to it. 
def getFastaFiles(strHeader, strRead, strCode, flZip):
    strFileIdentifier = "UMI_" + strCode + ".fasta"
    strFileName = flZip + strFileIdentifier
    with open(strFileName, "a") as flOutput:
        flOutput.write(strHeader)
        flOutput.write(strRead)

# The getTargetZero function.
# This function checks if both the forward and reverse primer can be found, if
# that succeeds, the forward or reverse (when working with single UMIs) or
# forward and reverse (when working with double UMIs) nucleotides are isolated
# based on the length of the UMI. This isolation is done from the first
# nucleotide at the 5'-end of a read and the 3'-end of a read.
def getTargetZero(strRead, intUmiLength, strSearch, strForward, strReverse):
    tplCheckForward = re.search(strForward, strRead)
    if tplCheckForward != None:
        tplCheckReverse = re.search(strReverse, strRead)
        if tplCheckReverse != None:
            if strSearch == "umi5":
                return strRead[0:int(intUmiLength)]
            elif strSearch == "umidouble":
                return (strRead[0:int(intUmiLength)],
                        strRead[-int(intUmiLength):])
            elif strSearch == "umi3":
                return strRead[-int(intUmiLength):]
            else:
                pass
        else:
            pass
    else:
        pass

# The getTargetFront function.
# This function searches for a regex string in the provided read. It will
# isolate either a forward, reverse or double UMIs. The isolation is based on
# this read structure UMI-PRIMERF-PRODUCT-PRIMERR-UMI.
# When looking for the forward UMI, the first position of PRIMERF is used, when
# looking for the reverse UMI, the last position of PRIMERR is used, when
# looking for double UMIs both positions are used.
# The mentioned positions + or - the UMI length result in a UMI code.
def getTargetFront(strRead, intUmiLength, strSearch, strForward, strReverse):
    if strSearch == "umi5" or strSearch == "umidouble":
        intPositionForward = re.search(strForward, strRead).end()
        intPositionUmiForward = intPositionForward + int(intUmiLength)
        strUmiCodeForward = strRead[intPositionForward:intPositionUmiForward]
        if strSearch == "umi5":
            tplCheckReverse = re.search(strReverse, strRead)
            if tplCheckReverse != None:
                return strUmiCodeForward
            else:
                pass
        elif strSearch == "umidouble":
            intPositionReverse = re.search(strReverse, strRead).start()
            intPositionUmiReverse = intPositionReverse - int(intUmiLength)
            strUmiCodeReverse = strRead[intPositionUmiReverse:intPositionReverse]
            return strUmiCodeForward, strUmiCodeReverse
        else:
            pass
    elif strSearch == "umi3":
        tplCheckForward = re.search(strForward, strRead)
        if tplCheckForward != None:
            intPositionReverse = re.search(strReverse, strRead).start()
            intPositionUmiReverse = intPositionReverse - int(intUmiLength)
            strUmiCodeReverse = strRead[intPositionUmiReverse:intPositionReverse]
            return strUmiCodeReverse
        else:
            pass
    else:
        pass

# The getTargetBehind function.
# This function searches for a regex string in the provided read. It will
# isolate either a forward, reverse or double UMIs. The isolation is based on
# this read structure UMI-PRIMERF-PRODUCT-PRIMERR-UMI.
# When looking for the forward UMI, the first position of PRIMERF is used, when
# looking for the reverse UMI, the last position of PRIMERR is used, when
# looking for double UMIs both positions are used.
# The mentioned positions + or - the UMI length result in a UMI code.
def getTargetBehind(strRead, intUmiLength, strSearch, strForward, strReverse):
    if strSearch == "umi5" or strSearch == "umidouble":
        intPositionForward = re.search(strForward, strRead).start()
        intPositionUmiForward = intPositionForward - int(intUmiLength)
        strUmiCodeForward = strRead[intPositionUmiForward:intPositionForward]
        if strSearch == "umi5":
            tplCheckReverse = re.search(strReverse, strRead)
            if tplCheckReverse != None:
                return strUmiCodeForward
            else:
                pass
        elif strSearch == "umidouble":
            intPositionReverse = re.search(strReverse, strRead).end()
            intPositionUmiReverse = intPositionReverse + int(intUmiLength)
            strUmiCodeReverse = strRead[intPositionReverse:intPositionUmiReverse]
            return strUmiCodeForward, strUmiCodeReverse
        else:
            pass
    elif strSearch == "umi3":
        tplCheckForward = re.search(strForward, strRead)
        if tplCheckForward != None:
            intPositionReverse = re.search(strReverse, strRead).end()
            intPositionUmiReverse = intPositionReverse + int(intUmiLength)
            strUmiCodeReverse = strRead[intPositionReverse:intPositionUmiReverse]
            return strUmiCodeReverse
        else:
            pass
    else:
        pass

# The getReverseComplement function.
# This function creates a complementary string using a sequence as input. The
# function loops through a list version of the sequence and checks and changes
# every character. It then returns a joined string.
def getReverseComplement(strLine):
    dicComplementCodes = {"A": "T", "T": "A", "G": "C", "C": "G", "M": "K",
                          "R": "Y", "W": "W", "S": "S", "Y": "R", "K": "M",
                          "V": "B", "H": "D", "D": "H", "B": "V", "N": "N"}
    lstLine = list(strLine)
    for intPosition in range(len(lstLine)):
        lstLine[intPosition] = dicComplementCodes[lstLine[intPosition]]
    return "".join(lstLine)

# The getRegex function.
# This function creates a regex string using a sequence as input. This regex
# string is based on the IUPAC ambiguity codes. The function loops through
# a list version of the sequence and checks per character if it is a
# ambiguous character, and changes it into regex code if so. It then returns
# a joined string. 
def getRegex(strLine):
    dicAmbiguityCodes = {"M": "[AC]", "R": "[AG]", "W": "[AT]", "S": "[CG]",
                         "Y": "[CT]", "K": "[GT]", "V": "[ACG]", "H": "[ACT]",
                         "D": "[AGT]", "B": "[CGT]", "N": "[GATC]"}
    lstLine = list(strLine)
    for intPosition in range(len(lstLine)):
        if lstLine[intPosition] != "A" and lstLine[intPosition] != "T" and\
           lstLine[intPosition] != "G" and lstLine[intPosition] != "C":
            lstLine[intPosition] = dicAmbiguityCodes[lstLine[intPosition]]
        else:
            pass
    return "".join(lstLine)

# The getUmiCode function.
# This function controls the umi searching process. It first uses the functions
# getRegex and getReverseComplement to create regex strings of both the
# standard forward and reverse primers/scaffolds and regex strings of the
# reverse complement of the herefore mentioned primers/scaffolds. These regex
# strings are then directed to the desired functions, this depends on the
# search method choice [primer/scaffold/zero].
def getUmiCode(strRead, strProcess, intUmiLength, strSearch, strForward,
               strReverse):
    strRead = strRead.strip("\n")
    strRegexForward = getRegex(strForward)
    #strRegexComplementForward = getRegex(getReverseComplement(\
    #                                     strForward[::-1]))
    #strRegexReverse = getRegex(strReverse)
    strRegexComplementReverse = getRegex(getReverseComplement(\
                                         strReverse[::-1]))
    if strProcess == "primer":
        try:
            return getTargetBehind(strRead, intUmiLength, strSearch,
                                   strRegexForward, strRegexComplementReverse)
        except AttributeError:
            pass
            #try:
            #    return getTargetBehind(strRead, intUmiLength, strSearch,
            #                           strRegexReverse,
            #                           strRegexComplementForward)
            #except AttributeError:
            #    pass
    elif strProcess == "scaffold":
        try:
            return getTargetFront(strRead, intUmiLength, strSearch,
                                  strRegexForward, strRegexComplementReverse)
        except AttributeError:
            pass
            #try:
            #    return getTargetFront(strRead, intUmiLength, strSearch,
            #                          strRegexReverse,
            #                          strRegexComplementForward)
            #except AttributeError:
            #    pass
    elif strProcess == "zero":
        try:
            return getTargetZero(strRead, intUmiLength, strSearch,
                                 strRegexForward, strRegexComplementReverse)
        except AttributeError:
            pass
            #try:
            #    pass
            #except AttributeError:
            #    pass
    else:
        pass

# The getUmi function.
# This function opens the input file and loops through it. It isolates the
# read headers and reads. For every read the getUmiCode function is called
# which outputs a UMI. This UMI is then used to separate and count their
# occurrences, which is kept track of in dicUniqueUmi and dicUniqueUmiCount
# respectively. Every UMI, read header and read is used to create separate
# fasta files with the getFastaFiles function. After every read has been
# processed the getNewCsv function is called.
def getUmiCollection(flInput, flTabular, flZip, flBlast, strProcess,
                     intUmiLength, strSearch, strForward, strReverse,
                     strFormat, strOperand):
    dicUniqueUmi = {}
    dicUniqueUmiCount = {}
    with open(flInput) as oisInput:
        for strLine in oisInput:
            if strLine[0] == strOperand and bool(re.match("[A-Za-z0-9]",
               strLine[1])) == True:
                strHeader = strLine
                strRead = next(oisInput)
                try:
                    strUmiCode = getUmiCode(strRead.upper(), strProcess,
                                            intUmiLength, strSearch, strForward,
                                            strReverse)
                except UnboundLocalError:
                    pass
                try:
                    if strUmiCode != None:
                        if strSearch == "umi5" or strSearch == "umi3":
                            strCode = strUmiCode
                        elif strSearch == "umidouble":
                            strCode = strUmiCode[0] + strUmiCode[1]
                        else:
                            pass
                        if strCode not in dicUniqueUmi:
                            dicUniqueUmi[strCode] = strRead.strip("\n")
                            dicUniqueUmiCount[strCode] = 1
                            with open(flBlast, "a") as flOutput:
                                flOutput.write(">" + strCode + "\n")
                                flOutput.write(strRead.strip("\n") + "\n")
                        elif strCode in dicUniqueUmi:
                            dicUniqueUmiCount[strCode] += 1
                        else:
                            pass
                    else:
                        pass
                except UnboundLocalError:
                    pass
                try:
                    getFastaFiles(strHeader, strRead, strCode, flZip)
                except UnboundLocalError:
                    pass
    getNewCsv(dicUniqueUmi, dicUniqueUmiCount, flTabular)

# The setFormat function.
# This function specifies the first character of the read headers based on the
# input file format. It then calls the getUmi function.
def setFormat(flInput, flTabular, flZip, flBlast, strProcess, strFormat,
              intUmiLength, strSearch, strForward, strReverse):
    if strFormat == "fasta":
        strOperand = ">"
    elif strFormat == "fastq":
        strOperand = "@"
    else:
        pass
    getUmiCollection(flInput, flTabular, flZip, flBlast, strProcess, 
                     intUmiLength, strSearch, strForward, strReverse,
                     strFormat, strOperand)

# The argvs function.
def parseArgvs():
    parser = argparse.ArgumentParser(description="Use a python script to\
                                                  accumulate all UMIs and\
                                                  output a tabular, ZIP and\
                                                  BLAST file.")
    parser.add_argument("-v", action="version", version="%(prog)s [0.1.0]")
    parser.add_argument("-i", action="store", dest="fisInput",
                        help="The location of the input file(s)")
    parser.add_argument("-o", action="store", dest="fosOutput",
                        help="The location of the tabular output file(s)")
    parser.add_argument("-z", action="store", dest="fosOutputZip",
                        help="The location of the zip output file(s)")
    parser.add_argument("-q", action="store", dest="fosBlastFile",
                        help="The location of the BLAST output file(s)")
    parser.add_argument("-p", action="store", dest="disProcess",
                        help="The UMI search approach [primer/scaffold(adapter)/zero]")
    parser.add_argument("-f", action="store", dest="disFormat",
                        help="The format of the input file(s) [fasta/fastq]")
    parser.add_argument("-l", action="store", dest="disUmiLength",
                        help="The length of the UMI sequences")
    parser.add_argument("-s", action="store", dest="disSearch",
                        help="Search UMIs at 5'-end [umi5], 3'-end [umi3] or at\
                              5'-end and 3'-end [umidouble]")
    parser.add_argument("-a", action="store", dest="disForward",
                        help="The 5'-end search nucleotides")
    parser.add_argument("-b", action="store", dest="disReverse",
                        help="The 3'-end search nucleotides")
    argvs = parser.parse_args()
    return argvs

# The main function.
def main():
    argvs = parseArgvs()
    setFormat(argvs.fisInput, argvs.fosOutput, argvs.fosOutputZip,
              argvs.fosBlastFile, argvs.disProcess, argvs.disFormat,
              argvs.disUmiLength, argvs.disSearch, argvs.disForward,
              argvs.disReverse)

if __name__ == "__main__":
    main()
