/**
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 15/09/2019
 **
 **/ 

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <iterator>
#include <sstream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

//const static int max_id_length=80;

typedef struct Interval {
  vector<string> name;
  int start;
  int end;
} interval;

/** function to return all the names of overlapping interval **/
vector<string> get_overlaps(const vector<interval> & target, const interval & query){
  vector<string> result;
  //loop over the targer intervals and check if any overlap
  for(int i=0; i<target.size(); i++){
    interval t=target.at(i);
    int total_length=(t.end-t.start)+(query.end-query.start);
    int union_length=max(t.end,query.end)-min(t.start,query.start);
    if(total_length > union_length) result.insert(result.end(),t.name.begin(),t.name.end());
  }
  return result;
}

/** function which takes a set of intervals and returns the gaps
    in coverage (excluding the start and end) **/
vector<interval> get_regions(const vector<interval> & alignments){
  vector<bool> bases;
  //loop over the alignments for each contig
  for(int m=0 ; m<alignments.size(); m++){
    interval current_match=alignments.at(m);
    //resize the vector if the match extends past the end of the bases vector
    if(bases.size()<=current_match.end)
      bases.resize(current_match.end+1, 0);
    for(int b=current_match.start; b<=current_match.end; b++)
      bases.at(b)=true;
  }
  //now loop over the bases to work out if there are multiple intervals
  vector<interval> uinterval;
  interval region;
  for(int b=0; b<bases.size(); b++){
    //check for the start of an alignment region
    if(bases.at(b)==true && ( b==0 || bases.at(b-1)==false ))
      region.start=b;
     //and the end
    if(bases.at(b)==true && (b==(bases.size()-1) || bases.at(b+1)==false)){
      region.end=b;
      //get the names of all genes that overlap the region
      region.name=get_overlaps(alignments,region);
      uinterval.push_back(region);
    } 
  }
  return(uinterval);
};

//function to insert a new alignment into the map that checks if alignment to same
//gene already exists and extends it
void insert_alignment( vector<interval> & alignments, interval new_match){
  //check for existing alignment
  //loop over alignments
  int i=0;
  while(i<alignments.size()){
    if(alignments.at(i).name[0]==new_match.name[0]){
      alignments.at(i).start=min(alignments.at(i).start,new_match.start);
      alignments.at(i).end=max(alignments.at(i).end,new_match.end);
      return;
    }
    i++;
  }
  //if no existing match found
  alignments.push_back(new_match);
}


int main(int argc, char **argv){
  
  bool DO_CUTTING=true;
  char * cluster_output_file=NULL;
  char * fasta_output_file=NULL;
  float min_coverage=0;

  int params=1,c;
  while((c =  getopt(argc, argv, "nc:f:m:")) != EOF){
    switch(c){
    case 'n': { //don't cut
      DO_CUTTING=false;
      cout << "Setting to no cutting"<<endl;
      params+=1;
      break; }
    case 'c': {//set the cluster file
      cout << "Setting output cluster file to " << optarg << endl;
      cluster_output_file=optarg;
      params+=2;
      break; }
    case 'f': {//set the cluster file
      cout << "Setting output fasta file to " << optarg << endl;
      fasta_output_file=optarg;
      params+=2;
      break; } 
    case 'm': {//set the cluster file
      cout << "Setting minimum coverage to " << optarg << endl;
      min_coverage=atof(optarg);
      params+=2;
      break; }
    case '?': //other unknown options
      cerr << "Unknown option.. stopping" << endl;
      exit(1);
      break;
    }
  }
  if((argc-params)!=2){
    cerr << "Wrong number of arguments" << endl;
    cerr << "Usage: " << endl;
    cerr << "   chimera_breaker [options] <blat result> <input_fasta_file>" << endl;
    cerr << "   -n    no cutting" << endl;
    cerr << "   -c    cluster file output name" << endl;
    cerr << "   -f    fasta file output name" << endl;
    cerr << "   -m    minimum coverage to count an alignemtn (alignment length / length of smaller sequence)" << endl;
    exit(1);
  }
  
  ifstream file;
  char * filename=argv[params];
  //Open the blat table
  file.open(filename);
  if(!(file.good())){
    cerr << "Unable to open file " << filename << endl;
    exit(1);
  } 
  /**********  Read all the blat alignments ****************/
  cerr << "Reading the input alignment file, "<<filename << endl;
  string line;
  map<string,vector<interval> > contig_alignments;
  map<string,vector<interval> > gene_alignments;

  for(int i=0; i<5 ; i++) getline(file,line) ; //skip the first 5 lines
  while(getline(file,line) ){
    istringstream line_stream(line);
    vector<string> columns;
    string temp;
    while(getline(line_stream, temp,'\t'))
      columns.push_back(temp);
    if(columns.size()!=21){ 
      cerr << columns.size() ; 
      cerr << "Warning: unexpected number of columns in psl file, "
	   << filename << endl; continue ; 
    }
    if(atof(columns[0].c_str()) / min(atof(columns[10].c_str()),atof(columns[14].c_str())) >= min_coverage){
      interval this_match; //temporary variable for alignment
      //add the match of contig onto gene
      this_match.name.push_back(columns[13]);
      this_match.start=atoi(columns[11].c_str());
      this_match.end=atoi(columns[12].c_str());
      insert_alignment(gene_alignments[columns[9]],this_match);
      //gene_alignments[columns[9]].push_back(this_match);
      //add the match of the gene onto the contig
      interval this_match2; //temporary variable for alignment
      this_match2.name.push_back(columns[9]);
      this_match2.start=atoi(columns[15].c_str());
      this_match2.end=atoi(columns[16].c_str());
      insert_alignment(contig_alignments[columns[13]],this_match2);
    //contig_alignments[columns[13]].push_back(this_match);
    }
  }
  file.close();

  /************* Read the fasta file ***************************/
  //Now open the fasta file and get the IDs and the sequence
  char * fname=argv[params+1];
  file.open(fname);
  if(!(file.good())){
    cerr << "Unable to open file " << fname << endl;
    exit(1);
  }
  cerr << "Reading the input fasta file, " << fname << endl;
  string id;
  string sequence;
  map<string,string> contigs;
  while ( getline (file,line) ){
    int start=line.find(">")+1;
    if(start==1){ //if this is the ID line...
      //push back the pervious entry
      if(id!="") contigs[id]=sequence;
      //reset the variables
      sequence=""; 
      int end=line.find_first_of("\t\n ")-1;
      id=line.substr(start,end);
    } else { sequence+=line; }
  } //push back the final sequence
  contigs[id]=sequence ;
  file.close();

  /******* Loop through the contigs and decide if they need to be cut ****/
  cerr << "Checking for chimeric contigs to break" << endl;
  int ccount=0;
  map<string,string> new_contigs;
  map<string,vector<string > > gene_contig_mapping; //holds the clusters
  vector<vector<string > > genes_to_merge;
  ostringstream gene_id;
  
  for(map<string,string>::iterator citr=contigs.begin(); citr!=contigs.end(); citr++){
    vector<interval> align_regions;
    if(DO_CUTTING) align_regions=get_regions(contig_alignments[citr->first]);
    else { //if we don't want to break the chimeras, just make a single interval will
      // the names of all genes that match.
      interval new_interval; 
      for(int i=0; i<contig_alignments[citr->first].size() ; i++){
	new_interval.name.push_back(contig_alignments[citr->first].at(i).name.at(0));
      }
      align_regions.push_back(new_interval);
    }
    interval contig_range;
    //cut the contig if a chimera was detected and update the id
    int last_cut_position=-1;
    for(int i=1; i < align_regions.size() ; i++){
       ccount++;
       int start=last_cut_position+1;
       int cut_position=(align_regions[i].start+align_regions[i-1].end)/2;
       stringstream new_id;
       //new_id << align_regions.at(i-1).name
       new_id << citr->first << "." << last_cut_position+1 << "-" << cut_position;
       new_contigs[new_id.str()]=citr->second.substr(last_cut_position+1,cut_position-last_cut_position-1);
       last_cut_position=cut_position;

       //loop over the aligned genes and add this contig to the cluster list
       vector<string> genes = align_regions.at(i-1).name;
       genes_to_merge.push_back(genes);
       for(int n=0; n<genes.size(); n++)
	 gene_contig_mapping[genes.at(n)].push_back(new_id.str());
    }
    stringstream new_id; //output the final region (which might be whole contig for case of no breaks)
    //if(align_regions.size()>0) //if the contig aligns to a gene
    //  new_id << align_regions.at(align_regions.size()-1).name;
    new_id << citr->first;
    if(align_regions.size()>1) //if the contig was broken
      new_id << "." << last_cut_position+1 << "-" << citr->second.size()-1;
    new_contigs[new_id.str()]=citr->second.substr(last_cut_position+1);
    if(align_regions.size()>0){
      vector<string> genes = align_regions.at(align_regions.size()-1).name;
      genes_to_merge.push_back(genes);
      for(int n=0; n<genes.size(); n++)
	gene_contig_mapping[genes.at(n)].push_back(new_id.str());
    }
  }

  /**** work out the clusters *****/
  ofstream ofile;
  if(cluster_output_file){
    ofile.open(cluster_output_file);
    if(!(ofile.good())){
      cerr << "Unable to create file " << cluster_output_file << endl;
      exit(1);
    }
    cerr << "Creating the annotation cluster file, " << cluster_output_file << endl;
    //start by clustering the annotated genes/transcripts
    vector<vector<string > >::iterator sitr = genes_to_merge.begin(); 
    
    // Loop over the vectors
    for(;genes_to_merge.size()>0 && sitr != genes_to_merge.end(); sitr++){
      bool finished_cluster=false;
      while(!finished_cluster){ //we need check the same cluster over after adding in new annotations.
	finished_cluster=true;
	//clean up the list by removing redundancy
	sort( sitr->begin(), sitr->end() );
	sitr->erase( unique( sitr->begin(), sitr->end() ), sitr->end() );
	//for each vector loop over the genes to merge
	vector<string> these_genes = *sitr; //copy the list (can't use direct as we're modifying it
	for(vector<string>::iterator gitr=these_genes.begin() ; gitr!=these_genes.end(); gitr++){
	  //loop again over the vector to see if the gene can be found in another group
	  vector< vector< string > >::iterator sitr2=sitr+1;
	  for(; sitr2 < genes_to_merge.end(); sitr2++){
	    if(find(sitr2->begin(), sitr2->end(), *gitr) != sitr2->end()){
	      //merge the groups
	      sitr->insert(sitr->end(),sitr2->begin(),sitr2->end());
	      //remove duplicate
	      sitr2->clear();
	      finished_cluster=false;
	    }
	  }
	}
      }
      //make a single string for the gene name
      ostringstream annotation_oss;
      copy(sitr->begin(),sitr->end(),ostream_iterator<string>(annotation_oss,":"));
      //now get all the transcripts that match the genes
      vector<string> trans_ids;
      vector<string>::iterator gene_itr=sitr->begin();
      for(; gene_itr!=sitr->end(); gene_itr++){
	vector<string> this_genes_trans=gene_contig_mapping[*gene_itr];
	trans_ids.insert(trans_ids.end(),this_genes_trans.begin(),this_genes_trans.end());
      }
      //clean up the transcript id list
      sort( trans_ids.begin(), trans_ids.end() );
      trans_ids.erase( unique( trans_ids.begin(), trans_ids.end() ), trans_ids.end() );
      vector<string>::iterator trans_itr=trans_ids.begin();
      for(;trans_itr!=trans_ids.end();trans_itr++){
	ofile << *trans_itr << "\t" << annotation_oss.str() << endl;
      }
    }
    ofile.close();
  }
  
  /***** output the new contigs ****/
  if(fasta_output_file){
    ofile.open(fasta_output_file);
    if(!(ofile.good())){
      cerr << "Unable to create file " << fasta_output_file << endl;
      exit(1);
    }
    cerr << "Creating the chimera broken fasta file, " << fasta_output_file << endl;
    for(map<string,string>::iterator citr=new_contigs.begin(); citr!=new_contigs.end(); citr++){  
      ofile << ">" << citr->first << endl;
      ofile << citr->second << endl;
    }
    ofile.close();
  }
  
  cerr << contigs.size() << " contigs processed" << endl;
  cerr << gene_alignments.size() << " known genes/transcripts annotated " << endl;
  cerr << ccount << " cuts made " << endl;
  
  cerr << "All done" << endl;
}
