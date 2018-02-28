/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 28/2/2018
 *********************************************************/

de_novo_assembly_file=""
trinity_options="--max_memory 50G --normalize_reads"

de_novo_assemble = {
    output.dir="de_novo_assembly"
    produce("de_novo_assembly.fasta"){
	if(de_novo_assembly_file!=""){
	doc "Using user's de novo assembly"
	    exec "cp $de_novo_assembly_file $output"
    	} else {
	doc "De novo assembling the reads with Trinity"
       	exec """
          ${Trinity} --seqType fq $trinity_options
       	    --left $reads_R1  --right $reads_R2 --CPU $threads
	    --full_cleanup ;
	    mv trinity_out_dir.Trinity.fasta $output
	    ""","trinity"
	}
    }
}

