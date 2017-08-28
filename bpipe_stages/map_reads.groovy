/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 
 *********************************************************/
ST_mapped_dir="mapped_reads"

build_ST_index = {
        output.dir=ST_mapped_dir
        from("SuperDuper.fasta"){ produce("ST.1.ht2"){
          exec "${hisat2}-build $input $output.prefix.prefix"
        }}
}

map_reads_to_ST = {
        output.dir=ST_mapped_dir
        produce(branch.name+".bam",branch.name+".splice.sites",branch.name+".summary"){
           exec """
           ${hisat2} --dta
	      --summary-file $output3
	      --pen-noncansplice 0
	      --novel-splicesite-outfile $output2
              -x $input.ht2.prefix.prefix
              -1 $input1 -2 $input2 |
            $samtools view -Su - | $samtools sort - -o $output1 ;
	      $samtools index $output1
           """
           }
}

set_input = {
   def files=reads_R1.split(",")+reads_R2.split(",")
   forward files
}

fastqInputFormat="%_R*.gz"
map_reads = segment { build_ST_index +
	    	      set_input + fastqInputFormat * [ map_reads_to_ST ]
}





