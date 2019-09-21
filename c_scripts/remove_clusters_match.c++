/**
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 9th May 2017
 **
 **/ 

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <vector>
#include <stdlib.h>

using namespace std;

int main(int argc, char **argv){
  if(argc!=3){
    cerr << "Wrong number of arguments" << endl;
    cerr << "Usage: " << endl;
    cerr << "   remove_clusters_with_match <ref.clusters> <denovo.clusters>" << endl;
    exit(1);
  }

  vector<string>contigs_in_anno;
  ifstream file;
  //Open the cluster table
  file.open(argv[1]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  cerr << "Reading the input file, "<< argv[1] << endl;
  string line; 
  string contig_id;
  while(getline(file,line) ){
    istringstream line_stream(line);
    getline(line_stream, contig_id,'\t');
    contigs_in_anno.push_back(contig_id);
  }
  file.close();

  //open the next file
  map<string,string> denovo_clust;
  file.open(argv[2]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[2] << endl;
    exit(1);
  }
  cerr << "Reading the input file, "<< argv[2] << endl;
  while(getline(file,line) ){
    istringstream line_stream(line);
    string transcript_id;
    string gene_id;
    getline(line_stream, transcript_id,'\t');
    getline(line_stream, gene_id,'\t');
    denovo_clust[transcript_id]=gene_id;
  }
  file.close();

  //loop through the denovo list and work out which clusters have the previous contigs
  cerr << "Get rel-clusters with a reference match" << endl;
  vector<string> denovo_clust_in_anno;
  for(vector<string>::iterator it = contigs_in_anno.begin(); it != contigs_in_anno.end(); ++it){
    map<string,string>::iterator dc_it = denovo_clust.find( *it );
    if(dc_it != denovo_clust.end())
      denovo_clust_in_anno.push_back(dc_it->second);
  }
  
  //loop a final time to remove the clusters
  cerr << "Remove rel-clusters with a reference match and output" << endl;
  for(map<string,string>::iterator it = denovo_clust.begin(); it != denovo_clust.end(); ++it){
    if(find(denovo_clust_in_anno.begin(), denovo_clust_in_anno.end(), it->second) == denovo_clust_in_anno.end()){
      cout << it->first << "\t" << it->second << endl;      
    }
  }

  cerr << "All done" << endl;
}
