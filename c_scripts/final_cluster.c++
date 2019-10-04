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

const int max_id_length=80;

map< string, string> get_alignments( string filename, string prefix){
  map<string,string> clusters;
  ifstream file;
  //Open the cluster table
  file.open(filename.c_str());
  if(!(file.good())){
    cerr << "Unable to open file " << filename << endl;
    exit(1);
  }
  // Read the alignments
  cerr << "Reading the input file, "<<filename << endl;
  string line;
  while(getline(file,line) ){
    istringstream line_stream(line);
    string transcript_id;
    string gene_id;
    string new_gene_id;
    getline(line_stream, transcript_id,'\t');
    if(!getline(line_stream, gene_id,'\t'))
      new_gene_id=transcript_id;
    else
      new_gene_id=prefix+gene_id;
    if(new_gene_id.size()>max_id_length)
      new_gene_id=new_gene_id.substr(0,max_id_length);
    clusters[transcript_id]=new_gene_id;
  }
  file.close();
  return clusters;
}

int main(int argc, char **argv){
  if(argc!=4){
    cerr << "Wrong number of arguments" << endl;
    cerr << "Usage: " << endl;
    cerr << "   final_cluster <genome_superT.names> <ref.clusters> <denovo.clusters>" << endl;
    exit(1);
  }

  //order matters here as final clustering will override earlier ones
  map<string,string> anno_clust=get_alignments(argv[1],"");
  map<string,string> denovo_anno_clust=get_alignments(argv[2],"");
  map<string,string> denovo_clust=get_alignments(argv[3],"REL-");
  
  //merge the clusters prioritizing the annotation over related species clusters.

  cerr << "Merge clusters" << endl;
  //make the final clustering
  map<string,string> clusters;
  clusters.insert(denovo_anno_clust.begin(),denovo_anno_clust.end());
  clusters.insert(anno_clust.begin(),anno_clust.end());
  clusters.insert(denovo_clust.begin(),denovo_clust.end());

  //Output the clusters
  cerr << "Outputting final clusters" << endl;
  for(map<string,string>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    cout << it->first << "\t" << it->second << endl;

  cerr << "All done" << endl;
}
