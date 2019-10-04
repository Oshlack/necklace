/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 
 *********************************************************/
ST_mapped_dir="mapped_reads"

build_ST_index = {
        doc "Building HISAT2 index for the superTranscriptome"
        output.dir=ST_mapped_dir
        produce("ST.1.ht2"){
          exec "${hisat2}-build $input.fasta $output.prefix.prefix"
        }
}

map_reads_to_ST = {
	def final_map_input_option=""
        if(reads_R2=="")
             final_map_input_option = "-U "+input
        else
             final_map_input_option = "-1 "+input1+" -2 "+input2
	def thrds=nthreads/reads_R1.split(",").size()
	def map_threads=Math.max(thrds.intValue(),1)
	println "Using "+map_threads+" threads"
	doc "Aligning reads back to the superTranscriptome using HISAT2"
        output.dir=ST_mapped_dir
        produce(branch.name+".bam",branch.name+".splice.sites",branch.name+".summary"){
           exec """
           $hisat2 $hisat2_options --dta
	      --threads $map_threads
	      --summary-file $output3
	      --pen-noncansplice 0
	      --novel-splicesite-outfile $output2
              -x $input.ht2.prefix.prefix
              $final_map_input_option |
            $samtools view -Su - | $samtools sort - -o $output1 ;
	      $samtools index $output1
           """
           }
}

map_reads = segment { 
  build_ST_index + fastqInputFormat * [ map_reads_to_ST ]
}
