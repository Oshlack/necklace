/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 
 *********************************************************/

counts_dir="counts"

get_splice_blocks = {
   output.dir=counts_dir
   produce("blocks.gtf"){
   exec """
    	$samtools faidx $input.fasta ;
	cut -f1,2 ${input.fasta}.fai > $output.dir/gene.sizes ;
	$make_blocks $output.dir/gene.sizes $inputs.sites > $output
	"""
    }
}

count_reads = {
    output.dir=counts_dir
    from("*.bam","blocks.gtf") produce ("gene.counts","block.counts"){
           exec """
		$featureCounts -T $threads --primary -p -t exon -g gene_id 
			       -a $input.gtf -o $output1 $inputs.bam ;
		$featureCounts -T $threads --primary -p -t exon -g gene_id 
			       --fraction -O -f -a $input.gtf -o $output2 $inputs.bam ;
		"""
    }
}

get_counts = segment { get_splice_blocks + count_reads }

