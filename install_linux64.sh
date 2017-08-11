#!/bin/bash

## This script will install the tools required for the pipeline
## It will fetched each tool from the web and place it into the tools/ subdirectory.
## Paths to all installed tools can be found in the file tools.groovy at the
## end of execution of this script. These paths can be changed if a different
## version of software is required. Note that R must be installed manually
##
## Last Modified: 22nd October by Nadia Davidson

mkdir -p tools/bin 
cd tools 

#a list of which programs need to be installed
commands="bpipe hisat2 stringtie gffread Trinity blat cluster python lace gtf2flatgtf"
#samtools bowtie2 dedupe reformat"

#installation method
function bpipe_install {
   wget -O bpipe-0.9.9.2.tar.gz https://github.com/ssadedin/bpipe/releases/download/0.9.9.2/bpipe-0.9.9.2.tar.gz
   tar -zxvf bpipe-0.9.9.2.tar.gz ; rm bpipe-0.9.9.2.tar.gz
   ln -s $PWD/bpipe-0.9.9.2/bin/* $PWD/bin/
}

function Trinity_install {
    wget https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.0.5.tar.gz
    tar -zxvf v2.0.5 ; rm v2.0.5
    make -C trinityrnaseq-2.0.5
    echo "$PWD/trinityrnaseq-2.0.5/Trinity \$@" > $PWD/bin/Trinity
    chmod +x $PWD/bin/Trinity
}

function hisat2_install {
    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.5-Linux_x86_64.zip
    unzip hisat2-2.0.5-Linux_x86_64.zip ; rm hisat2-2.0.5-Linux_x86_64.zip
    ln -s $PWD/hisat2-2.0.5/* $PWD/bin/
}

function stringtie_install {
   wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz
   tar xvfz stringtie-1.3.3b.Linux_x86_64.tar.gz
   rm stringtie-1.3.3b.Linux_x86_64.tar.gz
   ln -s $PWD/stringtie-1.3.3b.Linux_x86_64/stringtie $PWD/bin/
}

function gffread_install {
   wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.9.8c.Linux_x86_64.tar.gz
   tar xvfz gffread-0.9.8c.Linux_x86_64.tar.gz 
   rm gffread-0.9.8c.Linux_x86_64.tar.gz
   ln -s $PWD/gffread-0.9.8c.Linux_x86_64/gffread $PWD/bin/ 
}

function blat_install {
   wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
   mv blat $PWD/bin
   chmod +x $PWD/bin/blat
}

function fasta_formatter_install {
    wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    tar -jxvf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    rm fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
}

function cluster_install {
    g++ -o bin/cluster ../cluster.c
}

function gtf2flatgtf_install {
    g++ -o bin/gtf2flatgtf ../gtf2flatgtf.c
}

function python_install {
   wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash ./Miniconda3-latest-Linux-x86_64.sh -b -p $PWD/miniconda
   rm ./Miniconda3-latest-Linux-x86_64.sh
   ln -s $PWD/miniconda/bin/* $PWD/bin/
   bin/conda install -y pandas networkx numpy matplotlib
}

function lace_install {
    wget https://github.com/Oshlack/Lace/releases/download/v0.99/Lace-0.99.tar.gz -O Lace-0.99.tar.gz
    tar -xvf Lace-0.99.tar.gz ; rm Lace-0.99.tar.gz
    cd Lace-0.99
    ../bin/conda env create environment.yml
    cd ../
    echo "source $PWD/bin/activate lace ; $PWD/bin/python $PWD/Lace-0.99/Lace.py \$@" > bin/lace
    chmod +x bin/lace
    bin/conda install pandas
}

#function dedupe_install {
#    wget --no-check-certificate https://sourceforge.net/projects/bbmap/files/BBMap_36.59.tar.gz
#    tar -zxvf BBMap_36.59.tar.gz
#    rm BBMap_36.59.tar.gz
#    for script in `ls $PWD/bbmap/*.sh` ; do
#	s=`basename $script`
#	s_pre=`echo $s | sed 's/.sh//g'`
#	echo "$PWD/bbmap/$s \$@" > $PWD/bin/$s_pre
#	chmod +x $PWD/bin/$s_pre
#    done
#}


echo "// Path to tools used by the pipeline" > ../tools.groovy

for c in $commands ; do 
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then 
	echo "$c not found, fetching it"
	${c}_install
	c_path=`which $PWD/bin/$c 2>/dev/null`
    fi
    echo "$c=\"$c_path\"" >> ../tools.groovy
done

#finally check that R is install
R_path=`which R 2>/dev/null`
if [ -z $R_path ] ; then
    echo "R not found!"
    echo "Please go to http://www.r-project.org/ and follow the installation instructions."
    echo "Note that the IRanges R package must be installed."
fi
echo "R=\"$R_path\"" >> ../tools.groovy

#loop through commands to check they are all installed
echo "Checking that all required tools were installed:"
Final_message="All commands installed successfully!"
for c in $commands ; do
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then
	echo -n "WARNING: $c could not be found!!!! " 
	echo "You will need to download and install $c manually, then add its path to tools.groovy"
	Final_message="WARNING: One or more command did not install successfully. See warning messages above. \
                       You will need to correct this before running the pipeline."
    else 
        echo "$c looks like it has been installed"
    fi
done
echo "**********************************************************"
echo $Final_message



