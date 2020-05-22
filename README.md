# galaxy-tool-UMI-isolation
This repository contains a Galaxy specific xml file that can be used to install
Caltha on the Naturalis Galaxy instance.  
A description and further text about Caltha can be found on the
main [repository](https://github.com/JasperBoom/caltha).  

## Old description (before Caltha)
Use a python script to accumulate all UMIs and output a tabular file, a BLAST
file and a zip file. The tabular file will contain all unique UMI nucleotides,
a count of the number of reads that are associated with that UMI and a read
that represents the sequences associated with that umi. The BLAST file can be
used to identify all UMI clusters. The zip file will contain fastA files for
every unique UMI and contain reads associated with
that UMI.

About this repository:

This repo is kept as a backup to guarantee continuity for Naturalis.
Active development on this tool is done in the Caltha package:
[https://github.com/JasperBoom/caltha](https://github.com/JasperBoom/caltha)

Getting started:

Installing:

```
Install the tool for use in Galaxy  
(user: **galaxy**)  

cd /home/galaxy/Tools

git clone https://github.com/naturalis/galaxy-tool-UMI-isolation
chmod -R 755 galaxy-tool-UMI-isolation

Add the following line to /home/galaxy/galaxy/config/tool_conf.xml

<tool file="/home/galaxy/Tools/galaxy-tool-UMI-isolation/getUmiIsolation.xml"/>
```

## Source(s)
* __Giardine B, Riemer C, Hardison RC, Burhans R, Elnitski L, Shah P__,  
  Galaxy: A platform for interactive large-scale genome analysis.  
  Genome Research. 2005; 15(10) 1451-1455. __doi: 10.1101/gr.4086505__  
  [Galaxy](https://www.galaxyproject.org/)
* __Python Software Foundation__,  
  Python 3.8+. 2019.  
  [Python](https://www.python.org/)
* __Rognes T, Flouri T, Nichols B, Quince C, Mahe F__,  
  VSEARCH: A versatile open source tool for metagenomics.  
  PeerJ. 2016. __doi: 10.7717/peerj.2584__  
  [Vsearch](https://github.com/torognes/vsearch)
* __Augspurger T, Ayd W, Bartak C, Battiston P, Cloud P, Garcia M__,  
  Python Data Analysis Library.  
  [Pandas](https://pandas.pydata.org/)
* __Langa L, Willing C, Meyer C, Zijlstra J, Naylor M, Dollenstein Z__,  
  The uncompromising Python code formatter.  
  [Black](https://black.readthedocs.io/en/stable/)
* __Ziadé T, Cordasco I__,  
  Your tool for style guide enforcement.  
  [Flake8](http://flake8.pycqa.org/en/latest/index.html)
* __Sottile A, Struys K, Kuehl C, Finkle M__,  
  A framework for managing and maintaining multi-language pre-commit hooks.  
  [Pre-commit](https://pre-commit.com/)
* __Python Software Foundation__,  
  The Python Package index.  
  [PyPI](https://pypi.org/)

## Author(s)
* [Jasper Boom](https://github.com/JasperBoom)

```
Copyright (C) 2018 Jasper Boom

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
```
