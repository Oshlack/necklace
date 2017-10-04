/***********************************************************
 ** Stages to construct a superTranscriptome using genome
 ** annotation and genome-guided assembly
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 28th August 2017
 *********************************************************/

genome_superTranscriptome_dir="genome_superTranscriptome"

merge_genome_annotations = {
      output.dir=genome_superTranscriptome_dir
      produce("ref_annotations_combined.gtf","genome_merged.gft"){
	exec """
	     cat $annotation > $output1 ;
	     ${stringtie} --merge -G $output1 -o $output2 $input $annotation
	     """
      }
}

flatten_gtf = {
     output.dir=genome_superTranscriptome_dir
     from("genome_merged.gft") { produce("genome_merged.flattened.gft"){
	exec "$gtf2flatgtf $input $output" 
     }}
}

extract_exons_from_fasta = {
     output.dir=genome_superTranscriptome_dir
     produce("genome_superT.fasta"){
        exec "${gffread} $input -g $genome -w $output"
     }
}

build_genome_superTranscriptome = segment {
			       	 	merge_genome_annotations +
			       	 	flatten_gtf +
					extract_exons_from_fasta
					}
