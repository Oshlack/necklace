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
//#include <vector>
//#include <algorithm>
#include <stdlib.h>

using namespace std;

void get_alignments( string filename, string prefix, map<string,string> & clusters){
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
    getline(line_stream, transcript_id,'\t');
    if(!getline(line_stream, gene_id,'\t')){
      clusters[transcript_id]=transcript_id;
    } else {
      clusters[transcript_id]=prefix+gene_id;
    }
  }
  file.close();
}

int main(int argc, char **argv){
  if(argc!=4){
    cerr << "Wrong number of arguments" << endl;
    cerr << "Usage: " << endl;
    cerr << "   final_cluster <genome_superT.names> <ref.clusters> <denovo.clusters>" << endl;
    exit(1);
  }

  //order matters here as final clustering will override earlier ones
  map<string,string> clusters;
  get_alignments(argv[3],"REL-",clusters);
  get_alignments(argv[1],"",clusters);
  get_alignments(argv[2],"",clusters);
  
  //Output the clusters
  for(map<string,string>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    cout << it->first << "\t" << it->second << endl;

  cerr << "All done" << endl;
}
