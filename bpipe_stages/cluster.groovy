/***********************************************************
 ** Stages to cluster de novo assembled contigs prior to running lace
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 18th August 2017
 *********************************************************/


//config
blat_minIdentity="98"
blat_minScore="200"
rel_blat_minScore="200"

//output directory
cluster_dir="cluster_files"

blat_relST_genomeST = {
    output.dir=cluster_dir
    from("genome_superT.relative.fasta","genome_superT.fasta") produce("relST_genomeST.psl"){
       exec "$blat $input1 $input2 -t=dnax -q=dnax -minScore=$rel_blat_minScore $output.psl"
    }
}

blat_relST_denovo = {
    output.dir=cluster_dir
    from("genome_superT.relative.fasta","de_novo_assembly.fasta") produce("relST_denovo.psl"){
       exec "$blat $input1 $input2 -t=dnax -q=dnax -minScore=$rel_blat_minScore $output.psl"
    }
}

//Combine these three blat commands?
blat_genomeST_denovo = {
    output.dir=cluster_dir
    from("genome_superT.fasta","de_novo_assembly.fasta") produce("genomeST_denovo.psl"){
       exec "$blat $input1 $input2 -minScore=$blat_minScore -minIdentity=$blat_minIdentity $output.psl"
    }
}

get_genome_ST_names = {
    output.dir=cluster_dir
    from("genome_superT.fasta") produce ("genome_superT.names"){
    exec "grep \"^>\" $input | sed 's/>//g' | cut -f1 -d \" \" > $output "
    }
}

make_cluster_files_for_lace = {
   output.dir=cluster_dir
   from("genomeST_denovo.psl","relST_denovo.psl","relST_genomeST.psl",
        "genome_superT.names","genome_superT.fasta","de_novo_assembly.fasta") 
	produce ("clusters.txt","sequences.fasta"){ 
      exec """
      	   $cluster $inputs.psl $input.names > $output1 ;
	   cat $inputs.fasta > $output2
	   """
    }
}

cluster_files = segment { [blat_relST_genomeST, 
	            blat_relST_denovo,
		    blat_genomeST_denovo,
		    get_genome_ST_names ] + make_cluster_files_for_lace  }



