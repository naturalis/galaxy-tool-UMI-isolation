# galaxy-tool-UMI-isolation
Use a python script to accumulate all UMIs and output a tabular file, a BLAST file and a zip file.\
The tabular file will contain all unique UMI nucleotides, a count of the number of reads that are\
associated with that UMI and a read that represents the sequences associated with that umi.\
The BLAST file can be used to identify all UMI clusters.  The zip file will contain fastA files for every\
unique UMI and contain reads associated with that UMI.

# About this repository
This repo is kept as a backup to guarantee continuity for Naturalis.\
An author maintained verion of this tool can be found here:\
[https://github.com/JasperBoom/galaxy-tools-umi-isolation](https://github.com/JasperBoom/galaxy-tools-umi-isolation)\
Active development on this tool is done in the Caltha package:\
[https://github.com/JasperBoom/caltha](https://github.com/JasperBoom/caltha)


# Getting started
### Installing
Install the tool for use in Galaxy  
(user: **galaxy**)  
```
cd /home/galaxy/Tools
```
```
git clone https://github.com/naturalis/galaxy-tool-UMI-isolation
chmod -R 755 galaxy-tool-UMI-isolation
```
Add the following line to /home/galaxy/galaxy/config/tool_conf.xml
```
<tool file="/home/galaxy/Tools/galaxy-tool-UMI-isolation/getUmiIsolation.xml"/>
```

## Source(s)
* __Giardine B, Riemer C, Hardison RC, Burhans R, Elnitski L, Shah P__,  
  Galaxy: A platform for interactive large-scale genome analysis.  
  Genome Research. 2005; 15(10) 1451-1455 __doi: 10.1101/gr.4086505__  
  [GALAXY](https://www.galaxyproject.org/)

## Author(s)
* [Jasper Boom](https://github.com/JasperBoom)
* [Rutger Vos](https://github.com/rvosa)

## Citation
* __Boom J__ & __Vos RA__, galaxy-tool-UMI-isolation. 2019.  
  Github repository: https://github.com/naturalis/galaxy-tool-UMI-isolation
