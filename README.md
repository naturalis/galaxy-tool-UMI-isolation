# galaxy-tool-UMI-isolation
Use a python script to cluster all UMIs and output a tabular file, a BLAST file and a zip file.  
Use a python script to accumulate all UMIs and output a tabular file, a BLAST file and a zip file. The tabular file will contain all unique UMI nucleotides, a count of the number of reads that are associated with that umi and a read that represents the sequences associated with that umi.  
The BLAST file can be used to identify all UMI clusters.  
The zip file will contain fastA files for every unique UMI and contain reads associated with that UMI.

# Getting started

### Prerequisites
Download and install the following software according to the following steps.
```
sudo apt-get install python-pip
sudo pip install pandas
```

### Installing
Download and install the tool according to the following steps.
```
sudo mkdir -m 755 /home/Tools
cd /home/Tools
sudo git clone https://github.com/JasperBoom/galaxy-tool-UMI-isolation
sudo chmod -R 755 galaxy-tool-UMI-isolation
```
The following file in the galaxy-tool-metadata folder should be made avaible from any location.
```
sudo ln -s /home/Tools/galaxy-tool-UMI-isolation/getUmiClusters.py /usr/local/bin/getUmiClusters.py
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

## Citation
* __Boom J__, galaxy-tool-UMI-isolation. 2018.  
  Github repositry: https://github.com/JasperBoom/galaxy-tool-UMI-isolation
