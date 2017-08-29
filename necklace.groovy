/***********************************************************
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 
 *********************************************************/

VERSION="0.80"

codeBase = file(bpipe.Config.config.script).parentFile.absolutePath
load codeBase+"/tools.groovy"
load args[0]

load codeBase+"/bpipe_stages/genome_guided_assembly.groovy"
load codeBase+"/bpipe_stages/build_genome_superTranscriptome.groovy"
load codeBase+"/bpipe_stages/build_relatives_superTranscriptome.groovy"
load codeBase+"/bpipe_stages/de_novo_assembly.groovy"
load codeBase+"/bpipe_stages/cluster.groovy"
load codeBase+"/bpipe_stages/run_lace.groovy"
load codeBase+"/bpipe_stages/map_reads.groovy"
load codeBase+"/bpipe_stages/get_counts.groovy"
load codeBase+"/bpipe_stages/get_stats.groovy"


/******************* Here are the pipeline stages **********************/

run_check = {
    doc "check that the data files exist"
    if (output_dir) {
        output.dir=output_dir
    }
    produce("checks") {
        exec """
            echo "Running necklace version $VERSION" ;
            echo "Checking for the data files..." ;
	    for i in $genome $annotation $reads_R1 $reads_R2 $prot_related_species ; 
                 do ls $i 2>/dev/null || { echo "CAN'T FIND ${i}..." ;
		 echo "PLEASE FIX PATH... STOPPING NOW" ; exit 1  ; } ; 
	    done ;
            echo "All looking good" ;
            echo "running  necklace version $VERSION.. checks passed" > $output
        ""","checks"
    }
}

nthreads=8

run { run_check + //single thread 
    [ build_genome_guided_assembly + build_genome_superTranscriptome, 
    de_novo_assemble.using(threads: nthreads-2), 
    build_relatives_superTranscriptome ] + 
    cluster_files +
    run_lace.using(threads: nthreads) + 
    map_reads + 
    get_counts.using(threads: nthreads) + //
    get_stats } //single thread



