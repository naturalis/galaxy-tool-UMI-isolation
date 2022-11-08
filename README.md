# galaxy-tool-UMI-isolation  
Accumulate all UMIs and output a tabular file, a BLAST file (fasta) and a zip file. 

**The tabular file** will contain all unique UMI nucleotides, a count of the number of reads that are  
associated with that UMI and a read that represents the sequences associated with that umi.  
**The BLAST file** can be used to identify all UMI clusters.  
**The zip file** will contain fasta files for each UMI and the reads associated to it.  

## Installation
### Manual
Clone this repo in your Galaxy ***Tools*** directory:  
`git clone https://github.com/naturalis/galaxy-tool-UMI-isolation`  

Make the python script executable:  
`chmod 755 galaxy-tool-UMI-isolation/getUmiIsolation.sh`  
`chmod 755 galaxy-tool-UMI-isolation/getUmiIsolation.py` 

Append the file ***tool_conf.xml***:    
`<tool file="/path/to/Tools/galaxy-tool-UMI-isolation/getUmiIsolation.sh" />`  

### Ansible
Depending on your setup the [ansible.builtin.git](https://docs.ansible.com/ansible/latest/collections/ansible/builtin/git_module.html) module could be used.  
[Install the tool](https://docs.ansible.com/ansible/latest/collections/ansible/builtin/git_module.html#examples) 
by including the following in your dedicated ***.yml** file:  

`  - repo: https://github.com/naturalis/galaxy-tool-UMI-isolation`  
&ensp;&ensp;`file: getUmiIsolation.xml`  
&ensp;&ensp;`version: master`  

## About this repository
This repo is kept as a backup to guarantee continuity for Naturalis.\
An author maintained verion of this tool can be found here:\
[https://github.com/JasperBoom/galaxy-tools-umi-isolation](https://github.com/JasperBoom/galaxy-tools-umi-isolation)\
Active development on this tool is done in the Caltha package:\
[https://github.com/JasperBoom/caltha](https://github.com/JasperBoom/caltha)

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
