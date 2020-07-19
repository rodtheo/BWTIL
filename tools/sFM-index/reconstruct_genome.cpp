#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <sdsl/bit_vectors.hpp>

// Heng Li - KSEQ LIB for parsing FASTA files
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "kstring.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace sdsl;

int main(int argc, char *argv[]){
  // const int NUM_MISMATCHES = 10000;
  string inFileName = "outM.txt";
  ifstream inFile;
  inFile.open(inFileName.c_str());

  inFile.seekg (0, inFile.end);
  int length = inFile.tellg();
  inFile.seekg (0, inFile.beg);
  length = ceil(length/2 + 1);

  int numberOfMis;
  std::vector<char> M( length );

  cout << "LENGTH OF VECTOR IS = " << length << endl;
  int count = 0;
  while (!inFile.eof()){
    inFile >> M[count];
    // std::cout << M[count] << " ";
    count = count + 1;
  }

  sd_vector<> B;
  load_from_file(B, "sdbB.sdsl");
  cout << "Bit vector B (1s marking the mismatches positions) has size of "<< B.size() << " bits" << endl;

  sd_vector<> I;
  load_from_file(I, "sdbI.sdsl");
  cout << "Bit vector I indicating the first character of each phrase, has size of "<< I.size() << " bits" << endl;

  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(argv[1], "r");
  seq = kseq_init(fp);
  char *s;
  while ((l = kseq_read(seq)) >= 0){
    s = seq->seq.s;
    // cout << "seq: \n" <<  s << endl;
  }
  kseq_destroy(seq);
  gzclose(fp);

  cout << "seq: \n" <<  s << endl;


  return 0;
}
