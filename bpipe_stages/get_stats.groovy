/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 
 *********************************************************/

stats_dir="stats"

get_basic_info = {
   output.dir=stats_dir
   produce("basic_info.txt"){
      from("clusters.txt","SuperDuper.fasta",
	   "genome_merged.flattened.gft",
	   "ref_annotations_combined.gtf",
	   "*.summary"
	   ){

      exec """
         ann_text="Annotation" ; 
	 gg_text=`echo "\$ann_text + Genome guided assembly"` ;
	 dn_text=`echo "\$gg_text + De novo assembly"` ;

	 genes_dn_only=`cut -f2 $input1 | grep "^REL-" | sort -u  | wc -l` ;
	 genes_gg_only=`cut -f2 $input1 | grep "^MSTRG" | sort -u  | wc -l` ;
   	 genes_total=`cut -f2 $input1 | sort -u  | wc -l` ;
	 genes_gg=\$((\$genes_total-\$genes_dn_only)) ;
	 genes_ann=\$((\$genes_gg-\$genes_gg_only)) ;
	 echo "" >> $output ;
	 echo "Genes found:" >> $output ;
         echo "-------------" >> $output ;
	 echo "\$genes_ann\t\$ann_text" >> $output ;
	 echo "\$genes_gg\t\$gg_text" >> $output ;
	 echo "\$genes_total\t\$dn_text" >> $output ;
	 echo "" >> $output ;

	 bp_all=`grep -v "^>" $input2 | wc -c` ;
	 bp_all=\$(($bp_all-\$genes_total)) ;
	 bp_gg=`awk 'END { print s } { s += \$5 - \$4 + 1 }' $input3` ;
	 bp_ann=`$stringtie --merge $input4 | $gtf2flatgtf /dev/stdin /dev/stdout | awk 'END { print s } { s += \$5 - \$4 + 1 }' - ` ; 
	 echo "SuperTranscriptome Size (bp):" >> $output ;
	 echo "-----------------------------" >> $output ;
	 echo "\$bp_ann\t\$ann_text" >> $output ;
	 echo "\$bp_gg\t\$gg_text" >> $output ;
	 echo \$bp_all"\t\$dn_text" >> $output ;


	 reads_all=`cut -f1 -d " " $input5 | head -n1`
	 echo "Read pairs:" >> $output ;
	 echo "-----------" >> $output ;
	 echo "\$reads_all\tTotal >> $output ;
	 echo $inputs.summary >> $output ;
	 echo "\$bp_gg\t(\%)Map to genome" >> $output ;
	 echo \$bp_all"\t(\%)Map to superTranscriptome" >> $output ;
	 

      """
      }
   }
}

get_stats = segment { get_basic_info }

