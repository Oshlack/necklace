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
#include <vector>
#include <algorithm>
#include <stdlib.h>

using namespace std;

const static int max_id_length=100;

map<string,vector<string> > get_alignments( string filename, string prefix1="", string prefix2=""){
  ifstream file;
  //Open the blat table
  file.open(filename.c_str());
  if(!(file.good())){
    cerr << "Unable to open file " << filename << endl;
    exit(1);
  } 
  // Read the alignments
  cerr << "Reading the input file, "<<filename << endl;
  string line;
  map<string,vector<string> > alignments;
  for(int i=0; i<5 ; i++) getline(file,line) ; //skip the first 5 lines
  while(getline(file,line) ){
    istringstream line_stream(line);
    vector<string> columns;
    string temp;
    while(getline(line_stream, temp,'\t'))
      columns.push_back(temp);
    if(columns.size()!=21){ cerr << columns.size() ; cerr << "Warning: unexpected number of columns in psl file, "
				<< filename << endl; continue ; }
    alignments[prefix1 + columns[9]].push_back(prefix2 + columns[13]);
  }
  file.close();

  //remove duplicate entries
  for(map<string,vector<string> >::iterator it=alignments.begin(); it!=alignments.end(); ++it){
    sort(it->second.begin(), it->second.end()); 
    vector<string>::iterator last = unique(it->second.begin(), it->second.end());
    it->second.resize(distance(it->second.begin(),last));
  }
  return alignments;
}

int main(int argc, char **argv){
  if(argc!=5){
    cerr << "Wrong number of arguments" << endl;
    cerr << "Usage: " << endl;
    cerr << "   cluster <genomeST_denovo.psl> <relST_denovo.psl> <relST_genomeST.psl> <genome_superT.names>" << endl;
    exit(1);
  }

  map<string,vector<string> > genomeST_denovo = get_alignments(argv[1],"","");
  map<string,vector<string> > relST_denovo = get_alignments(argv[2],"","REL-");
  map<string,vector<string> > relST_genomeST = get_alignments(argv[3],"","REL-");

  //read the gene IDs from the genome and use these are cluster names to start with
  map<string,vector<string> > clusters;
  ifstream file;
  file.open(argv[4]);
  if(!(file.good())){
    cerr << "Unable to open file " << argv[4] << endl;
    exit(1);
  }
  string line;
  while(getline(file,line) )
    clusters[line].push_back(line);
  file.close();

  //loop over all the contigs which align to a known genome superTranscript
  for(map<string,vector<string> >::iterator it=genomeST_denovo.begin(); it!=genomeST_denovo.end(); ++it){
    //does the contig also align to the superTranscript of the related species?
    map<string,vector<string> >::iterator it2=relST_denovo.find(it->first);
    //if it is found, delete it from the related species list (this species takes priority)
    if(it2!=relST_denovo.end()) relST_denovo.erase(it2);
    //Now add the contig ID to the cluster list for those than only match 1 genome gene.
    //De novo contigs that match two or more genes have a high chance of being false chimeras.
    if(it->second.size()==1) clusters[it->second[0]].push_back(it->first);
  }  

  //loop over the alignments between our species's genome and the related
  //we will remove these from the related species cluster list later
  vector<string> relST_with_genomeST_match;
  for(map<string,vector<string> >::iterator it=relST_genomeST.begin(); it!=relST_genomeST.end(); ++it){
    relST_with_genomeST_match.insert(relST_with_genomeST_match.end(), 
				     it->second.begin(), it->second.end() );
  }
  sort(relST_with_genomeST_match.begin(), relST_with_genomeST_match.end());
  vector<string>::iterator last = unique(relST_with_genomeST_match.begin(),relST_with_genomeST_match.end());
  relST_with_genomeST_match.resize(distance(relST_with_genomeST_match.begin(),last));

  //now loop over the de novo contigs alignments to the related species (only the novel stuff should be left)
  for(map<string,vector<string> >::iterator it=relST_denovo.begin(); it!=relST_denovo.end(); ++it){
    //Again, remove all de novo contigs that match two or more gene (these are unreliable).
    if(it->second.size()==1){
      //Is this gene actually missing from the genome?
      //Check whether the related gene aligns to the genome.
      vector<string>::iterator found_position = find(relST_with_genomeST_match.begin(),relST_with_genomeST_match.end(),it->second[0]);
      if(found_position==relST_with_genomeST_match.end()) //not found
	clusters[it->second[0]].push_back(it->first);
    }
  }

  //Now we are finished making the clusters. Just need to output them
  for(map<string,vector<string> >::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    string cluster_id = it->first;
    //check if the cluster_id it too long (lace will make a filename that's too long in this case)
    if(cluster_id.size()>max_id_length){
      cluster_id=cluster_id.substr(0,max_id_length);
    }
    vector<string> seqs = it->second;
    for(int s=0; s < seqs.size(); s++)
      cout << seqs.at(s) << "\t" << cluster_id << endl;
  }

  cerr << "All done" << endl;
}
