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
          exec "cat $annotation | $python ${hisat2}_extract_splice_sites.py - > $output"
       }
}

map_reads_to_genome = {
	def input_reads_option=""
	if(reads_R2=="") 
	     input_reads_option = "-U "+input
	else 
	     input_reads_option = "-1 "+input1+" -2 "+input2
	doc "Aligning reads to genome using HISAT2"
	output.dir=genome_guided_assembly_dir
	produce(branch.name+".bam",branch.name+".summary"){
	   exec """
	   $hisat2 $hisat2_options 
	   	     --known-splicesite-infile $input.txt 
	   	     --dta
	             --summary-file $output2
	             -x $input.ht2.prefix.prefix 
		     $input_reads_option |
		     $samtools view -Su - | $samtools sort - -o $output1
	   """
	  }
}

genome_assembly = {
	output.dir=genome_guided_assembly_dir
	produce(branch.name+".gtf"){
	   if(genome_guided_assembly_file!=""){
	      exec "cp $genome_guided_assembly_file $output"
	   } else {
	      exec "$stringtie $input.bam -o $output $stringtie_options"
	   }
	}
}


build_genome_guided_assembly = segment { 
			      	      	build_genome_index +
			                gtf_to_splice_sites +
					fastqInputFormat *
					[ map_reads_to_genome +
					genome_assembly ] 
					}

