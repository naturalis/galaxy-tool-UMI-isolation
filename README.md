# galaxy-tool-UMI-isolation
Use a python script to accumulate all UMIs and output a tabular file, a BLAST file and a zip file. 
The tabular file will contain all unique UMI nucleotides, a count of the number of reads that are 
associated with that UMI and a read that represents the sequences associated with that umi. The BLAST 
file can be used to identify all UMI clusters.  The zip file will contain fastA files for every 
unique UMI and contain reads associated with that UMI.

# Getting started
### Installing
Download and install the tool according to the following steps.
```
cd /home/galaxy/Tools
git clone https://github.com/naturalis/galaxy-tool-UMI-isolation
chmod -R 755 galaxy-tool-UMI-isolation
```
Continue with the tool installation
```
sudo mkdir -m 755 /home/galaxy/tools/directoryname
sudo cp /home/Tools/galaxy-tool-UMI-isolation/getUmiClusters.sh /home/galaxy/tools/directoryname/getUmiClusters.sh
sudo cp /home/Tools/galaxy-tool-UMI-isolation/getUmiClusters.xml /home/galaxy/tools/directoryname/getUmiClusters.xml
```
Edit the following file in order to make galaxy display the tool.
```
/home/galaxy/config/tool_conf.xml
```
```
<tool file="airdentification/getUmiClusters.xml"/>
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
