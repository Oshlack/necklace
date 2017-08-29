/**
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 26 October 2016
 **
 ** when finished can use the result like so:
 ** gffread <result> -g /group/bioi1/shared/genomes/hg38/fasta/hg38.fa -w temp.fasta
 **/ 

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <map>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#include "misc.h"

using namespace std;

vector<g_interval> merge_intervals( vector<g_interval> interval){
  //loop over each interval
  vector<g_interval> result;
  if(interval.size()<2) return interval;
  for(int i=0; i<(interval.size()-1); i++){
    bool merged=false;
    for(int j=(i+1); (j<interval.size()) && !merged; j++){
      //check if the two intervals overlap
      int x1=interval.at(i).start;
      int y1=interval.at(i).end;
      int x2=interval.at(j).start;
      int y2=interval.at(j).end;
      if((interval.at(i).chrom==interval.at(j).chrom) &&
	 !(y1<x2) && !(y2<x1) &&
	 (interval.at(i).strand==interval.at(j).strand || 
	  interval.at(i).strand=="." || 
	  interval.at(j).strand==".")
	 ){
	interval.at(j).start=min(x1,x2);
	interval.at(j).end=max(y1,y2);
	merged=true;
      }
    }
    if(!merged) result.push_back(interval.at(i));
  }
  result.push_back(interval.at(interval.size()-1));
  return result;
}

// the real stuff starts here.
int main(int argc, char **argv){

  if(argc!=3){
    cerr << "Wrong number of arguments" << endl;
    cerr << "Usage: " << endl;
    cerr << "   gtf2flatgtf <in_gtf> <out_gtf_flattened>" << endl; // <out_chain>" << endl;
    exit(1);
  }

  ifstream file;
  //Open the exons file
  file.open(argv[1]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  // reading the gtf file
  string line;
  map<string,vector<g_interval> > exons;
  map<string,vector<string> > reference_ids; 
  while(getline(file,line) ){
    istringstream line_stream(line);
    vector<string> columns;
    string temp;
    while(getline(line_stream, temp,'\t'))  
      columns.push_back(temp);
    if(columns.size()!=9){ continue; }//cout << columns.size() ; cout << "Warning: unexpected number of columns in gtf file, "<<argv[1]<<endl;}    
    if(columns[2]=="exon"){
      g_interval new_exon;
      new_exon.chrom=columns[0];
      new_exon.start=atoi(columns[3].c_str());
      new_exon.end=atoi(columns[4].c_str());
      new_exon.strand=columns[6];
      //find the gene_id in the meta data.
      size_t gene_id_start=columns[8].find('"',columns[8].find("gene_id"));
      size_t gene_id_end=temp.find('"',gene_id_start+1);
      string symbol=columns[8].substr(gene_id_start+1,gene_id_end-gene_id_start-1);
      exons[symbol].push_back(new_exon);
      //find the ref_gene_id in the meta data if it exists
      size_t ref_gene_id_start=columns[8].find('"',columns[8].find("ref_gene_id"));
      if(ref_gene_id_start!=string::npos){
	size_t ref_gene_id_end=temp.find('"',ref_gene_id_start+1);
	string ref_symbol=columns[8].substr(ref_gene_id_start+1,ref_gene_id_end-ref_gene_id_start-1);
	reference_ids[symbol].push_back(ref_symbol);
      }
    }
  }
  file.close();

  const static int MAX_INTRON = 1000000;

  //neaten up the map to reference IDs
  map<string,string> ref_ids_string;
  map<string,vector<string> >::iterator idItr=reference_ids.begin();
  map<string,vector<string> >::iterator idItrEnd=reference_ids.end();
  for(;idItr!=idItrEnd; idItr++){
    sort_vector(idItr->second);
    vector<string> ref_ids=idItr->second;
    stringstream label;
    for(int r=0; r < ref_ids.size(); r++)
      label << ref_ids[r] << ":" ;
    label << idItr->first ;
    ref_ids_string[idItr->first]=label.str();
  }

  //loop over exons
  map<string,vector<g_interval> >::iterator exonItr=exons.begin();
  map<string,vector<g_interval> >::iterator exonItrEnd=exons.end();
  //Open the output file
  ofstream flat_file;
  static int chain_id=1;
  flat_file.open(argv[2]);
  for( ; exonItr!=exonItrEnd ; exonItr++){ //looping over the genes
    vector<g_interval> exonInterval = merge_intervals(exonItr->second); //merging exons
    sort(exonInterval.begin(),exonInterval.end(),g_interval_compare);
    int last_start, last_end=0;
    vector<int> genome_starts,genome_ends,block_sizes;
    string chrom=exonInterval.at(0).chrom;
    string strand=exonInterval.at(0).strand;
    int st_length=0;
    for(int i=0; i<exonInterval.size(); i++){ //loop over the exons in the gene
      string this_chrom=exonInterval.at(i).chrom;
      int this_start=exonInterval.at(i).start;
      int this_end=exonInterval.at(i).end;
      string this_strand=exonInterval.at(i).strand;
      if(last_end==0 
	 || ( chrom==this_chrom && abs(this_start-last_end)<MAX_INTRON && strand==this_strand)){
	genome_starts.push_back(this_start);
	genome_ends.push_back(this_end);
	st_length+=this_end-this_start+1;
	block_sizes.push_back(this_end-this_start+1);
	last_start=this_start; last_end=this_end; 
      }
    }
    for(int i=0; i<genome_starts.size(); i++){ 
      flat_file << chrom << "\t" ;
      flat_file << "gtf2flatgtf\texon\t";
      flat_file << genome_starts.at(i) << "\t" ;
      flat_file << genome_ends.at(i) << "\t" ;
      flat_file << ".\t"<< strand << "\t.\t";
      string gene_id=ref_ids_string[exonItr->first]; //use the reference ids if found
      if(gene_id=="") gene_id=exonItr->first; //else the assembly id
      flat_file << "gene_id \""<< gene_id << "\"; " ; 
      flat_file << "transcript_id \""<< gene_id << "\"" << endl;
    } }
  flat_file.close();

}
