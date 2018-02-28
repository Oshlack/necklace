/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 28/2/2018 
 *********************************************************/

do_stats=true
stats_dir="stats"

get_gene_info = {
   output.dir=stats_dir
   produce("gene_info.txt"){
      from("clusters.txt"){
      exec """
         ann_text="Annotation" ; 
	 gg_text=`echo "\$ann_text + Genome guided assembly"` ;
	 dn_text=`echo "\$gg_text + De novo assembly"` ;
	 genes_dn_only=`cut -f2 $input1 | grep "^REL-" | sort -u  | wc -l` ;
	 genes_gg_only=`cut -f2 $input1 | grep "^MSTRG" | sort -u  | wc -l` ;
   	 genes_total=`cut -f2 $input1 | sort -u  | wc -l` ;
	 genes_gg=\$((\$genes_total-\$genes_dn_only)) ;
	 genes_ann=\$((\$genes_gg-\$genes_gg_only)) ;
	 rm -rf $output ;
	 echo "Genes found:" >> $output ;
         echo "-------------" >> $output ;
	 echo "\$genes_ann\t\$ann_text" >> $output ;
	 echo "\$genes_gg\t\$gg_text" >> $output ;
	 echo "\$genes_total\t\$dn_text" >> $output ;
      """
      }
   }
}

get_bp_info = {
   output.dir=stats_dir
   produce("size_info.txt"){
      from("clusters.txt",
	   "SuperDuper.fasta",
	   "genome_merged.flattened.gft",
	   "ref_annotations_combined.gtf"){

      exec """
         ann_text="Annotation" ; 
	 gg_text=`echo "\$ann_text + Genome guided assembly"` ;
	 dn_text=`echo "\$gg_text + De novo assembly"` ;
	 genes_dn_only=`cut -f2 $input1 | grep "^REL-" | sort -u  | wc -l` ;
	 genes_gg_only=`cut -f2 $input1 | grep "^MSTRG" | sort -u  | wc -l` ;
   	 genes_total=`cut -f2 $input1 | sort -u  | wc -l` ;
	 rm -rf $output ;
	 bp_all=`grep -v "^>" $input2 | wc -c` ;
	 bp_all=\$((\$bp_all-\$genes_total)) ;
	 bp_gg=`awk 'END { print s } { s += \$5 - \$4 + 1 }' $input3` ;
	 bp_ann=`$stringtie --merge $stringtie_merge_options $input4 | $gtf2flatgtf /dev/stdin /dev/stdout | awk 'END { print s } { s += \$5 - \$4 + 1 }' - ` ; 
	 echo "SuperTranscriptome Size (bp):" >> $output ;
	 echo "-----------------------------" >> $output ;
	 echo "\$bp_ann\t\$ann_text" >> $output ;
	 echo "\$bp_gg\t\$gg_text" >> $output ;
	 echo \$bp_all"\t\$dn_text" >> $output ;
      """
      }
   }
}

get_mapping_info = {
   output.dir=stats_dir
   produce("mapping_info.txt"){
      from("mapped2genome.sum","*.summary"){
      exec """
	 pairs=`cut -f1 -d " " $input1 | head -n1` ;
	 reads=\$((\$pairs*2)) ;
	 not_aligning_g=`grep "aligned 0 times\$" $input1 | cut -f1 -d "(" | awk 'END { print s } { s += \$1 }' - ` ;
	 not_aligning_st=`cat $inputs.summary | grep "aligned 0 times\$" | cut -f1 -d "(" | awk 'END { print s } { s += \$1 }' - ` ;
	 rm -rf $output ;
	 echo "Read:" >> $output ;
	 echo "-----" >> $output ;
	 echo "\$reads\tTotal" >> $output ;
	 echo \$((\$reads-\$not_aligning_g))"\tMap to genome" >> $output ;
	 echo \$((\$reads-\$not_aligning_st))"\tMap to superTranscriptome" >> $output ;
      """
      }
   }
}

get_stats = segment { get_gene_info + get_bp_info + get_mapping_info }

