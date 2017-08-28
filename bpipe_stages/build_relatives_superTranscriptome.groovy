/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 
 *********************************************************/


relatives_superTranscriptome_dir="relatives_superTranscriptome"


build_relatives_superTranscriptome = {
    output.dir=relatives_superTranscriptome_dir
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


