/***********************************************************
 ** Stages to preform genome-guided assembly with Stringties
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 28th August 2017
 *********************************************************/

//Output directory
genome_guided_assembly_dir="genome_guided_assembly"

//User specified genome assembly
genome_guided_assembly_file=""

build_genome_index = {
	output.dir=genome_guided_assembly_dir
	produce("genome.1.ht2"){
	  exec "${hisat2}-build $genome $output.prefix.prefix"
        }
}

gtf_to_splice_sites = {
       output.dir=genome_guided_assembly_dir
       produce("splicesites.txt"){
          exec "cat $annotation | ${hisat2}_extract_splice_sites.py - > $output"
       }
}

map_reads_to_genome = {
	output.dir=genome_guided_assembly_dir
	produce("genome_mapped.bam","mapped2genome.sum"){
	   exec """
	   $hisat2 $hisat2_options 
	   	     --known-splicesite-infile $input.txt 
	   	     --dta
	             --summary-file $output2
	             -x $input.ht2.prefix.prefix 
		     -1 $reads_R1 -2 $reads_R2 |
		     $samtools view -Su - | $samtools sort - -o $output1
	   """
	  }
}

genome_assembly = {
	output.dir=genome_guided_assembly_dir
	produce("genome_assembly.gtf"){
	   if(genome_guided_assembly_file!=""){
	      exec "cp $genome_guided_assembly_file $output"
	   } else {
	      exec " $stringtie $input.bam -o $output $stringtie_options"
	   }
	}
}

build_genome_guided_assembly = segment { 
			      	      	build_genome_index +
			                gtf_to_splice_sites +
					map_reads_to_genome +
					genome_assembly 
					}

