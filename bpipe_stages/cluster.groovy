/***********************************************************
 ** Stages to cluster de novo assembled contigs prior to running lace
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 18th August 2017
 *********************************************************/

//output directory
cluster_dir="cluster_files"

get_genome_ST_names = {
    output.dir=cluster_dir
    from("genome_superT.fasta") produce ("genome_superT.names"){
    exec "grep \"^>\" $input | sed 's/>//g' | cut -f1 -d \" \" > $output "
    }
}


blat_denovo_ref = {
    output.dir=cluster_dir
    from("de_novo_assembly.fasta") produce("denovo_ref.psl"){
       exec "$blat $input $proteins_related_species $blat_related_options $output.psl"
    }
}

cut_chimeras = {
    output.dir=cluster_dir
    from("de_novo_assembly.fasta","denovo_ref.psl") produce("denovo.clusters","denovo.fasta"){
       exec "$chimera_breaker -c $output.clusters -f $output.fasta $input.psl $input.fasta"
    }
}

blat_genomeST_denovo = {
    output.dir=cluster_dir
    from("genome_superT.fasta","denovo.fasta") produce("denovo_genomeST.psl"){
       exec """
       cat $input2 $input1 > ${output.dir}/temp.fasta ; 
       $blat $blat_options ${output.dir}/temp.fasta $input1 $output.psl ;
       """ 
    }
}

ref_cluster = {
    output.dir=cluster_dir
    from("denovo_genomeST.psl") produce("ref.clusters"){
       if(!scaffold_genes){ //remove all denovo contigs that align to more than one gene
          exec """ 
	  sed 1,5d $input | cut -f 10,14 | sort -k2 | uniq -u -f1 | awk '{printf("%s\\t%s\\n", \$2, \$1)}' > $output
          """
	} else { //otherwise allow the reference genes to be clustered together.
          exec """	
          $chimera_breaker -n -m 0.4 -c $output.clusters $input ${output.dir}/temp.fasta ;
          rm ${output.dir}/temp.fasta
       	  """
	}
    }
}

make_cluster_files_for_lace = {
   output.dir=cluster_dir
   from("ref.clusters",
	"denovo.clusters",
	"denovo_genomeST.psl",
        "genome_superT.names",
	"genome_superT.fasta",
	"denovo.fasta") 
	produce ("clusters.txt","sequences.fasta"){ 
      exec """
      	   cut -f 14 $input.psl | sort -u > ${output.dir}/temp.cl ;
	   $remove_clusters_match ${output.dir}/temp.cl $input2 > ${output.dir}/denovo_final.clusters;
	   rm ${output.dir}/temp.cl ;
	   $final_cluster $input.names $input1 ${output.dir}/denovo_final.clusters > $output1;
	   cat $inputs.fasta > $output2
	   """
    }
}



// OLD STUFF

/**blat_genomeST_ref = {
    output.dir=cluster_dir
    from("genome_superT.fasta") produce("genomeST_ref.psl"){
       exec "$blat $input $proteins_related_species $blat_related_options $output.psl"
    }
}

//sed 1,5d all_ref2.psl | cut -f 10,14 | sort -k2 | uniq -u -f1 | awk '{printf("%s\t%s\n", $2, $1)}'

 cluster_files = segment { [blat_genomeST_ref, 
	            blat_denovo_ref,
		    blat_genomeST_denovo,
		    get_genome_ST_names ] + make_cluster_files_for_lace  } 
**/

cluster_files = segment { 
	      blat_denovo_ref + cut_chimeras + 
	      blat_genomeST_denovo + ref_cluster + 
	      get_genome_ST_names + 
	      make_cluster_files_for_lace }

