# MindTheGap 

| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# What is MindTheGap ?

MindTheGap  performs detection and assembly of DNA insertion variants in NGS read datasets with respect to a reference genome. It takes as input a set of reads and a reference genome. It outputs two sets of FASTA sequences: one is the set of breakpoints of detected insertion sites, the other is the set of assembled insertions for each breakpoint. For each breakpoint, MindTheGap either returns  1) single insertion sequence;  2)  a set of candidate insertion sequences; 3) nothing at all (when the insertion is too complex to be assembled).

# Getting the latest source code

## Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

## Instructions

    # get a local copy of MindTheGap source code
    git clone --recursive https://github.com/GATB/MindTheGap.git
    
    # compile the code an run a simple test on your computer
    cd MindTheGap
    sh INSTALL

# User manual	 

Type `MindTheGap ` without any arguments for usage instructions.

Short manual: http://mindthegap.genouest.org/files/user_manual.pdf

#Contact

To contact a developer, request help, etc: https://gatb.inria.fr/contact/