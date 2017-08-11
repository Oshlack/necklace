/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 
 *********************************************************/

VERSION="0.80"

codeBase = file(bpipe.Config.config.script).parentFile.absolutePath
load codeBase+"/tools.groovy"
load args[0]

/********** Configuration ********/

//output_dir=""

/******************* Here are the pipeline stages **********************/

run_check = {
    doc "check that the data files exist"
    if (output_dir) {
        output.dir=output_dir
    }
    produce("checks") {
        exec """
            echo "Running XXXX? version $VERSION" ;
            echo "Checking for the data files..." ;
	    for i in $genome $annotation $reads_R1 $reads_R2 $prot_related_species ; 
                 do ls $i 2>/dev/null || { echo "CAN'T FIND ${i}..." ;
		 echo "PLEASE FIX PATH... STOPPING NOW" ; exit 1  ; } ; 
	    done ;
            echo "All looking good" ;
            echo "running  XXXXX? version $VERSION.. checks passed" > $output
        ""","checks"
    }
}

// genome-guided assembly 
build_genome_index = {
	output.dir=genome_guided_assembly_dir
	produce("genome.1.ht2"){
	  exec "${hisat2}-build $genome $output.prefix.prefix"
        }
}

build_ST_index = {
        output.dir=ST_mapped
        from("SuperDuper.fasta"){ produce("ST.1.ht2"){
          exec "${hisat2}-build $input $output.prefix.prefix"
        }}
}

gtf_to_splice_sites = {
       output.dir=genome_guided_assembly_dir
       produce("splicesites.txt"){
          exec "cat $annotation | ${hisat2}_extract_splice_sites.py - > $output"
       }
}

map_reads_to_genome = {
	output.dir=genome_guided_assembly_dir
	produce("genome_mapped.bam"){
	   exec """
	   ${hisat2} --known-splicesite-infile $input.txt --dta 
	             -x $input.ht2.prefix.prefix 
		     -1 $reads_R1 -2 $reads_R2 |
		     samtools view -Su - | samtools sort - $output.prefix
	   """
	   }
}

map_reads_to_ST = {
        output.dir=ST_mapped
        produce("mapped_reads.bam"){
           exec """
           time ${hisat2} --dta
                     -x $input.ht2.prefix.prefix
                     -1 $reads_R1 -2 $reads_R2 |
                     samtools view -Su - | samtools sort - $output.prefix
           """
           }
}

genome_assembly = {
	output.dir=genome_guided_assembly_dir
	produce("genome_assembly.gtf"){
	   exec " ${stringtie} $input.bam -o $output"
	}
}

merge_genome_annotations = {
      output.dir=genome_superTranscriptome_dir
      produce("genome_merged.gft"){
	exec """
	     cat $annotation > $output.dir/ref_annotations_combined.gtf ;
	     ${stringtie} --merge -G $output.dir/ref_annotations_combined.gtf -o $output $input $annotation
	     """
      }
}

flatten_gtf = {
     output.dir=genome_superTranscriptome_dir
     produce("genome_merged.flattened.gft"){
	exec "$gtf2flatgtf $input $output" 
     }
}

extract_exons_from_fasta = {
     output.dir=genome_superTranscriptome_dir
     produce("genome_superT.fasta"){
        exec "${gffread} $input -g $genome -w $output"
     }
}


max_memory="50G"
threads=6
de_novo_assemble = {
    output.dir="de_novo_assembly"
    produce("de_novo_assembly.fasta"){
       exec """
          ${Trinity} --seqType fq --max_memory $max_memory 
       	    --left $reads_R1  --right $reads_R2 --CPU $threads
	    --normalize_reads --full_cleanup ;
	    mv trinity_out_dir.Trinity.fasta $output
	    """
    }
}

make_relatives_superTranscriptome = {
    output.dir="relatives_superTranscriptome"
    produce("genome_superT.relative.fasta"){
        exec """
	   cat $annotation_related_species > $output.dir/related_species_annotations_combined.gtf ;
	   $stringtie --merge -G $output.dir/related_species_annotations_combined.gtf 
	   	      -o $output.dir/annotation_related_species.merged.gtf 
		      $annotation_related_species ; 
           $gtf2flatgtf $output.dir/annotation_related_species.merged.gtf 
	   		$output.dir/annotation_related_species.flattened.gtf ;
	   $gffread $output.dir/annotation_related_species.flattened.gtf -g $genome_related_species -w $output
	"""
     }
}

blat_relST_genomeST = {
    output.dir=cluster_dir
    from("genome_superT.relative.fasta","genome_superT.fasta") produce("relST_genomeST.psl"){
       exec "$blat $input1 $input2 -t=dnax -q=dnax -minScore=200 $output.psl"
    }
}

blat_relST_denovo = {
    output.dir=cluster_dir
    from("genome_superT.relative.fasta","de_novo_assembly.fasta") produce("relST_denovo.psl"){
       exec "$blat $input1 $input2 -t=dnax -q=dnax -minScore=200 $output.psl"
    }
}

//Combine these three blat commands?
blat_genomeST_denovo = {
    output.dir=cluster_dir
    from("genome_superT.fasta","de_novo_assembly.fasta") produce("genomeST_denovo.psl"){
       exec "$blat $input1 $input2 -t=dnax -q=dnax -minScore=200 $output.psl"
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

run_lace = {
   from("clusters.txt","sequences.fasta") 
   produce ("SuperDuper.fasta"){
   exec """
      $lace --cores $threads --tidy
             --outputDir lace_files
             $input2 $input1 ;
      mv lace_files/SuperDuper.fasta ../
	 """
   }
}

//instead of doing this align the superTranscripts to the genome and then use liftOver??
make_fasta_from_gtf = { 
    from("genome_assembly.gtf") produce ("genome.gtf.fasta"){
    exec "cat $annotation $input | $gffread - -g $genome -w $output"
    }
}

blat_de_novo_to_sT = {
    from("SuperDuper.fasta","de_novo_assembly.fasta"){
	produce("de_novo_to_sT.psl"){
	   exec "$blat $input1 $input2 -minScore=200 -minIdentity=98 $output"
	}
    }
}

blat_sT_to_genome = {
    from("SuperDuper.fasta") produce("sT_to_genome.psl"){
	   exec "$blat $genome $input -minScore=100 -minIdentity=98 $output"
    }
}

//blat SuperDuper.fasta ../analysis/cuff_GM2.fasta -minScore=200 -minIdentity=98 cuff_GM2.psl &
//blat /mnt/storage/shared/genomes/galGal5/fasta/galGal5.fa SuperDuper.fasta -minScore=100 -minIdentity=98 SuperDuper.genome.gal5.psl

/*** Pipeline segments ****/
genome_guided_assembly_dir="genome_guided_assembly"
build_genome_guided_assembly = segment { 
			      	      	build_genome_index +
			                gtf_to_splice_sites +
					map_reads_to_genome +
					genome_assembly }

genome_superTranscriptome_dir="genome_superTranscriptome"
make_genome_superTranscriptome = segment { 
			       	 	merge_genome_annotations +
			       	 	flatten_gtf +
					extract_exons_from_fasta
					 }

cluster_dir="cluster_files"
cluster_files = segment { [blat_relST_genomeST, 
	            blat_relST_denovo,
		    blat_genomeST_denovo,
		    get_genome_ST_names ] + make_cluster_files_for_lace  }


ST_mapped="mapped_reads"
map_reads = segment { build_ST_index +
	    	      map_reads_to_ST
}

//annotate_sT = segment { [blat_de_novo_to_sT,blat_sT_to_genome]}

run { [ build_genome_guided_assembly + make_genome_superTranscriptome, 
    de_novo_assemble, make_relatives_superTranscriptome ] + cluster_files +
    run_lace + map_reads }

//annotate_sT }

