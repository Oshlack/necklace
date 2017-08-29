/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 
 *********************************************************/

max_memory="50G"
de_novo_assemble = {
    doc "De novo assembling the reads with Trinity"
    output.dir="de_novo_assembly"
    produce("de_novo_assembly.fasta"){
       exec """
          ${Trinity} --seqType fq --max_memory $max_memory 
       	    --left $reads_R1  --right $reads_R2 --CPU $threads
	    --normalize_reads --full_cleanup ;
	    mv trinity_out_dir.Trinity.fasta $output
	    ""","trinity"
    }
}

