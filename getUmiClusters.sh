#!/bin/bash
# Copyright @2007-2019 Free Software Foundation

# Author: Jasper Boom

# Prequisites:
# - sudo apt-get install python
# - sudo apt-get install python-pip
# - sudo pip install pandas

# - sudo apt-get install libargtable2-dev
# - sudo ./configure && make && make install

# Galaxy prequisites:
# - sudo ln -s /path/to/folder/galaxy-tool-umi-clustering/getUmiClusters.py
#              /usr/local/bin/getUmiClusters.py

# The getFormatFlow function.
# This function creates a temporary storage directory in the output directory.
# It then calls the getUmiClusters.py script with the correct input values.
# The output tabular file is copied to the expected location and removed.
# The output BLAST file is copies to the expected lacation and removed.
# A ZIP file is created of the generated fastA file and copied to the expected
# location. After that the temporary storage directory is removed.

SCRIPTDIR=$(dirname "$(readlink -f "$0")")

getFormatFlow() {
    strDirectory=${fosOutput::-4}
    mkdir -p "${strDirectory}_temp"
    python $SCRIPTDIR"/getUmiClusters.py" -i ${fisInput} -o ${strDirectory}_temp/flTempCsv.csv \
                      -z ${strDirectory}_temp/ \
                      -q ${strDirectory}_temp/flTempBlast.fasta \
                      -f ${disFormat} -p ${disProcess} -l ${disUmiLength} \
                      -s ${disSearch} -a ${disForward} -b ${disReverse}
    cat ${strDirectory}_temp/flTempCsv.csv > ${fosOutputTabular}
    rm ${strDirectory}_temp/flTempCsv.csv
    cat ${strDirectory}_temp/flTempBlast.fasta > ${fosBlastFile}
    rm ${strDirectory}_temp/flTempBlast.fasta
    zip -jqr ${strDirectory}_temp/flTempZip.zip ${strDirectory}_temp/*
    cat ${strDirectory}_temp/flTempZip.zip > ${fosOutputZip}
    rm -rf ${strDirectory}_temp
}

# The main function.
main() {
    getFormatFlow
}

# The getopts function.
while getopts ":i:o:z:q:p:f:l:s:a:b:vh" opt; do
    case ${opt} in
        i)
            fisInput=${OPTARG}
            ;;
        o)
            fosOutputTabular=${OPTARG}
            ;;
        z)
            fosOutputZip=${OPTARG}
            ;;
        q)  
            fosBlastFile=${OPTARG}
            ;;
        p)
            disProcess=${OPTARG}
            ;;
        f)
            disFormat=${OPTARG}
            ;;
        l)
            disUmiLength=${OPTARG}
            ;;
        s)
            disSearch=${OPTARG}
            ;;
        a)
            disForward=${OPTARG}
            ;;
        b)
            disReverse=${OPTARG}
            ;;
        v)
            echo ""
            echo "getUmiClusters.sh [0.1.0]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: getUmiClusters.sh [-h] [-v] [-i INPUT] [-o TABULAR]"
            echo "                         [-z ZIP] [-q BLAST] [-p PROCESS]"
            echo "                         [-f FORMAT] [-l LENGTH] [-s SEARCH]"
            echo "                         [-a FORWARD] [-b REVERSE]"
            echo ""
            echo "Optional arguments:"
            echo " -h                    Show this help page and exit"
            echo " -v                    Show the software's version number"
            echo "                       and exit"
            echo " -i                    The location of the input file(s)"
            echo " -o                    The location of the tabular output file(s)"
            echo " -z                    The location of the zip output file(s)"
            echo " -q                    The location of the BLAST output file(s)"
            echo " -p                    The UMI search approach"
            echo "                       [primer/scaffold(adapter)/zero]"
            echo " -f                    The format of the input file(s)"
            echo "                       [fasta/fastq]"
            echo " -l                    The length of the UMI sequences"
            echo " -s                    Search UMIs at 5'-end [umi5],"
            echo "                       3'-end [umi3] or at 5'-end and" 
            echo "                       3'-end [umidouble]"
            echo " -a                    The 5'-end search nucleotides"
            echo " -b                    The 3'-end search nucleotides"
            echo ""
            echo "Use a python script to accumulate all UMIs and output a"
            echo "tabular file, a BLAST file and a zip file."
            echo "The tabular file will contain all unique UMI nucleotides, a"
            echo "count of the number of reads that are associated with that"
            echo "umi and a read that represents the sequences associated with"
            echo "that umi."
            echo "The BLAST file can be used to identify all UMI clusters."
            echo "The zip file will contain fastA files for every unique UMI"
            echo "and contain reads associated with that UMI."
            echo ""

            exit
            ;;
        \?)
            echo ""
            echo "You've entered an invalid option: -${OPTARG}."
            echo "Please use the -h option for correct formatting information."
            echo ""

            exit
            ;;
        :)
            echo ""
            echo "You've entered an invalid option: -${OPTARG}."
            echo "Please use the -h option for correct formatting information."
            echo ""

            exit
            ;;
    esac
done

main
