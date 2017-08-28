/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 
 *********************************************************/

superTranscriptome_dir="superTranscriptome"

run_lace = {
   output.dir=superTranscriptome_dir
   from("clusters.txt","sequences.fasta") 
   produce ("SuperDuper.fasta"){
   exec """
      $lace --cores $threads --tidy
            --outputDir $output.dir
            $input2 $input1 ;
      rm -rf $output.dir/SuperFiles
      ""","lace"
   }
}


