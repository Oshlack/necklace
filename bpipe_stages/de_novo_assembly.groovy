/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 28/2/2018
 *********************************************************/

de_novo_assemble = {
    def assembly_input_reads_option="--left "+reads_R1+" --right "+reads_R2
    //single end case
    if(reads_R2==""){assembly_input_reads_option="--single "+reads_R1}
    output.dir="de_novo_assembly"
    produce("de_novo_assembly.fasta"){
	if(de_novo_assembly_file!=""){
	doc "Using user's de novo assembly"
	    exec "cp $de_novo_assembly_file $output"
    	} else {
	doc "De novo assembling the reads with Trinity"
       	exec """
          ${Trinity} --seqType fq $trinity_options
       	    $assembly_input_reads_option --CPU $threads --full_cleanup ;
	    mv trinity_out_dir.Trinity.fasta $output
	    ""","trinity"
	}
    }
}

