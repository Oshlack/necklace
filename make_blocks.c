/**
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 22 August 2017
 **
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

//Very simple script which takes a list of gene lengths
//And a list of splice sites and produces a gtf with
//"splice blocks"
int main(int argc, char **argv){

  if(argc<3){
    cerr << "Wrong number of arguments" << endl;
    cerr << "Usage: " << endl;
    cerr << "   make_blocks <gene.sizes> <splice.sites1> <splice.sites2> ..." << endl; 
    exit(1);
  }

  ifstream file;
  string line;
  map<string,vector<int> > blocks;

  //Open and read gene sizes files
  file.open(argv[1]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  while(getline(file,line) ){
    istringstream line_stream(line);
    string gene; int size ;
    line_stream >> gene ;
    line_stream >> size ;
    blocks[gene].push_back(1);
    blocks[gene].push_back(size+1);
  }
  file.close();  

  //Open and read the splice site file
  for(int f=2; f < argc; f++){
    file.open(argv[f]);
    if(!(file.good())){
      cerr << "Unable to open file " << argv[f] << endl;
      exit(1);
    }
    while(getline(file,line) ){
      istringstream line_stream(line);
      string gene; int start ; int end ;
      line_stream >> gene ;
      line_stream >> start ;
      line_stream >> end ;
      blocks[gene].push_back(start+2);
      blocks[gene].push_back(end+1);
    }
    file.close();
  }

  //Now output the gtf to stdout
  map<string,vector<int> >::iterator blockItr=blocks.begin();
  map<string,vector<int> >::iterator blockItrEnd=blocks.end();
  for(;blockItr!=blockItrEnd; blockItr++){
    sort_vector(blockItr->second);
    vector<int> boundaries=blockItr->second;    
    for(int b=0; b<(boundaries.size()-1) ; b++){ 
      cout << blockItr->first << "\t" ;
      cout << "superTranscript_blocks\texon\t";
      cout << boundaries.at(b) << "\t" ;
      cout << boundaries.at(b+1)-1 << "\t" ;
      cout << ".\t+\t.\t";
      cout << "gene_id \""<< blockItr->first << "\"; " ; 
      cout << "transcript_id \""<< blockItr->first << "\"" << endl;
    } 
  }
}
