/*
 *  This file is part of BWTIL.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>
 *
 *   BWTIL is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   BWTIL is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

//============================================================================
// Name        : sFM-index.cpp
// Author      : Nicola Prezza
// Version     : 1.0
// Copyright   : GNU General Public License (http://www.gnu.org/copyleft/gpl.html)
// Description : test class for the succinct FM index
//============================================================================

#include <iostream>
#include <queue>
#include <string>

#include "../../data_structures/succinctFMIndex.h"
#include "../../data_structures/RMaxQBlockDecomp.h"

using namespace bwtil;
using namespace std;

static inline int aln_score(const int m, const int o) {
	// return m*p->mm_score + o*p->gapo_score + e*p->gape_score;
  return m*3 + o*11;
}

struct query_bases_t {
  uchar b;
  ulint L;
  ulint U;
  ulint maxAr;

  query_bases_t(){}
  query_bases_t(uchar b, ulint L, ulint U, ulint maxAr){}
};

struct aln_entry_t {
	// bwtint_t L;
	// bwtint_t U; // (L,U): SA interval of [i,n-1]
  ulint i;
  ulint L;
  ulint U;
  ulint next_p;
	uint32_t num_mm, num_gapo, num_snps;
	uint32_t score, state, aln_length;
  uchar c;
  uchar b;
	//uint8_t score; // aln score so far
	//uint8_t i; // position in the read
	//uint8_t num_mm;
	//uint8_t num_gapo;
	//uint8_t num_gape;
	//uint8_t state;
	//int n_seed_mm;
	//uint8_t last_diff_pos;
	//uint8_t padding;

	// edit transcript
	//uint8_t aln_length;
	// char aln_path[ALN_PATH_ALLOC];

  aln_entry_t(){} // DEFAULT CONSTRUCTOR

  aln_entry_t(ulint i, ulint L, ulint U, uint32_t num_mm, uint32_t num_gapo, uint32_t next_p, uint32_t state, uint32_t num_snps, uint32_t aln_length, uchar& uu, uchar& uuu): i(i), L(L), U(U), num_mm(num_mm), num_gapo(num_gapo), next_p(next_p), state(state), num_snps(num_snps), aln_length(aln_length), c(uu), b(uuu)
  {
    // c = *uu;
    // c = strcpy(c, uu);
  }

};

struct CompareAlnScores {
  bool operator()(aln_entry_t const& p1, aln_entry_t const& p2){
    // return true if alignment p1 has lower score
    // than alignment p2
    return aln_score(p1.num_mm, p1.num_gapo) < aln_score(p2.num_mm, p2.num_gapo);
  }
};

// Function which return string by concatenating it.
string repeat(string s, int n)
{
    // Copying given string to temparory string.
    string s1 = s;

    for (int i=1; i<n;i++)
        s += s1; // Concatinating strings

    return s;
}

int main(int argc,char** argv) {

	if(argc != 4 and argc != 3 and argc != 5){
		cout << "*** succinct FM-index data structure : a wavelet-tree based uncompressed FM index ***\n";
		cout << "Usage: sFM-index option file [pattern]\n";
		cout << "where:\n";
		cout <<	"- option = build|search|lz. \n";
		cout << "- file = path of the text file (if build mode) or .sfm sFM-index file (if search mode). \n";
		cout << "- pattern = must be specified in search mode. It is the pattern to be searched in the index.\n";
		exit(0);
	}

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

	int build=0,search=1,lz=2;

	int mode;

	if(string(argv[1]).compare("build")==0)
		mode=build;
	else if(string(argv[1]).compare("search")==0)
		mode=search;
	else if(string(argv[1]).compare("lz")==0)
		mode=lz;
	else{
		cout << "Unrecognized option "<<argv[1]<<endl;
		exit(0);
	}

	string in(argv[2]);
	string out(in);
	out.append(".sfm");

	string pattern;

	if(mode==search){
		pattern = string(argv[3]);
	}

  int k_seed;
  int debug = 0;
  string inS;


  if(mode==lz){
    k_seed = atoi(argv[3]);
    inS = argv[4];
  }

    auto t1 = high_resolution_clock::now();

	succinctFMIndex SFMI;
  // RL
  succinctFMIndex SFMI_S;

	if(mode==build){

		cout << "Building succinct FM-index of file "<< in << endl;
		SFMI = succinctFMIndex(in,true);

		cout << "\nStoring succinct FM-index in "<< out << endl;
		SFMI.saveToFile(out);

		cout << "Done.\n";

	}

	if(mode==search){

		cout << "Loading succinct FM-index from file "<< in <<endl;
		SFMI = succinctFMIndex::loadFromFile(in);
		cout << "Done." << endl;

		cout << "\nSearching pattern \""<< pattern << "\""<<endl;

		vector<ulint> occ = SFMI.getOccurrencies( pattern );

		cout << "The pattern occurs " << occ.size() << " times in the text at the following positions : \n";

		for(uint i=0;i<occ.size();i++)
			cout << occ.at(i) << " ";

		cout << "\n\nDone.\n";

	}

	if (mode==lz){

    ofstream outfile;
    outfile.open("out.txt");

		cout << "Loading succinct FM-index from file " << in << endl;
		SFMI = succinctFMIndex::loadFromFile(in);
		cout << "Done." << endl;
		IndexedBWT* idxBWT = SFMI.get_idxBWTPtr();
		RMaxQBlockDecomp RMQ(idxBWT);

    // RLZ - start
    SFMI_S = succinctFMIndex::loadFromFile(inS);
    IndexedBWT* idxBWTS = SFMI_S.get_idxBWTPtr();

		// vector<factor> enc;

		ulint i = 0;
		uint j = 0;
		ulint p = 1;
		uchar c;
		ulint cidx;
		pair<ulint, ulint> interval;
		ulint max_Ar_idx, max_Ar;
		uint iprime;
    // LZ-parser - take the Text length
		// ulint n = SFMI.textLength();
    // RLZ-parser - take the Sample text length
    ulint n = SFMI_S.textLength();
    ulint nR = SFMI.textLength();
    // RLZ-parser - END snippet

    cout << "TEXT LENGTH = " << nR << endl;
		ulint sp, ep;
    ulint next_p;
		int perc, last_perc= -1;
		// cout << "Letter " << c << endl;
		priority_queue<aln_entry_t, vector<aln_entry_t>, CompareAlnScores> heap;

    int best_score, num_best, max_entries, total_entries, last_num_entries, alns_get_num_entries;

    int diff_left, best_diff;
    uint32_t matches, max_diff;

    if (debug){
      cout << "###################################################" << endl;
      cout << "INITIALIZING ALGORITHM" << endl;
      cout << "###################################################\n\n" << endl;
    }
    string s = "##";
    int nhashtags;
    int allow_diff;

    int ii, old_p;
    p = 1;
    // for(ii=0; ii <= n; ii++){
    //   old_p = p;
    //   c = idxBWTS->at( p );
    //   p = idxBWTS->LF( p );
    //   cout << "CHAR " << c << " int(CHAR) = " << int(c) << " WITH P = " << old_p <<", NEXT P = " << p << endl;
    // }

    // int ref_len = 13;
    // int k_seed = 3;
		while (i < (n-1)){

      if (debug){
        cout << "###################################################" << endl;
        cout << " FOR i = " << i << endl;
      }
      nhashtags = 1;
			// cout << (int)(n/10.) << endl;
			// perc = (100*j)/n;
			// if (perc > last_perc and (perc%5)) {
			// 	cout << " " << perc << "% done." << endl;
			// 	last_perc=perc;
			// }

			j = i + 1;



			c = idxBWTS->at( p );
			// p = idxBWT->LF( p );

			interval = SFMI.arrayC(c);
			// cout << "LETTER =" << c << ", INTERVAL FIRST = " << interval.first << ", SA=" << idxBWT->convertToTextCoordinate(interval.first) << "| SECOND = " << interval.second << ", SA=" << idxBWT->convertToTextCoordinate(interval.second-1) << endl;
			// cout << "RMQ[sp,ep] = " << RMQ.query(interval.first, interval.second-1) << endl;

			sp = interval.first;
			ep = interval.second;

			max_Ar = RMQ.query(sp, ep-1);
			iprime = 0;

			heap = priority_queue<aln_entry_t, vector<aln_entry_t>, CompareAlnScores>();

			heap.push(aln_entry_t(i, sp, ep, 0, 0, p, 0, j, 0, c, c));

      best_score = aln_score(2, 2);
      best_diff = 2;
      max_diff = 1;
      num_best = 0;
      max_entries = 0;

      total_entries = 0;
      last_num_entries = 1;


      alns_get_num_entries = 0;

			// cout << "Max_Ar = " << max_Ar << endl;
			// cout << "i= " << i << " j= " << j << " c= " << c << " p= " << p << " sp=" << sp << " ep=" << ep << endl;

			// while ((j + 1 < n) && (sp <= ep) && (max_Ar >= (n-(i+1)))){
      aln_entry_t a;

      // allow_diff = 0;
      while (!heap.empty()){

        a = heap.top();
        heap.pop();

        allow_diff = a.state;

        diff_left = max_diff - a.num_mm;
        sp = a.L;
        ep = a.U;
        matches = a.aln_length;

        c = a.c;
        p = idxBWTS->LF( a.next_p );

        max_Ar = RMQ.query(sp, ep-1);

        if (debug){
          cout << "Pop i = " << a.i << " c=" << a.c << " and c =" << c << "int(c) =" << int(a.c) << " , actual base = "<< a.b << ", max Ar=" << max_Ar << " p =" << a.next_p << " next c=" << idxBWT->at( a.next_p ) << ", LF(p)=" << p << "  j=" << a.num_snps <<", j new=" << a.i+a.aln_length << " aln len="<< matches << " iprime =" << a.num_gapo << " L=" << a.L << " U=" << a.U << " num_mm = " << a.num_mm << " n-(a.i+1) =" << n-(a.i+1) << endl;

          // cout << "DUPLICATED Pop i + allow_diff = " << a.i + a.state << " allow_diff = " << a.state << " c=" << a.c << " and c =" << c << " , actual base = "<< a.b << ", max Ar=" << max_Ar << " p =" << a.next_p << " next c=" << idxBWT->at( a.next_p ) << ", LF(p)=" << p << "  j=" << a.num_snps <<", j new=" << a.i+a.aln_length+allow_diff << " aln len="<< matches << " iprime =" << a.num_gapo << " L=" << a.L << " U=" << a.U << " num_mm = " << a.num_mm << endl;

          cout << "Allow diff = " << allow_diff << endl;
        }

        if (diff_left < 0){
          if (debug){
            cout << " INSIDE DIFF LEFT !!!" << endl;
          }
	  // j = j - 1;
	  // a.aln_length = a.aln_length - 1;
          continue;}

        if (sp > (ep-1)){
          if (debug){
            cout << " INSIDE SP > (EP-1) !!!" << endl;
          }
          continue;}

        // LZ
        // if (max_Ar < n-(a.i+1)){
        // RLZ
        if (max_Ar < 0){
          if (debug){
          cout << " INSIDE THE MAX AR !!!" << endl;
        }
          continue;
        }

        // LZ
        // if (a.num_snps >= (n-1)){
        // RLZ
        if (a.num_snps >= (n)){
          if (debug){
          cout << " INSIDE THE NUM SNPS !!! " << a.num_snps << endl;
        }
          // a.num_snps = j + a.num_mm;

          // if (allow_diff >= k_seed){
          //   j = a.num_snps;
          // }
          // j = j + 1;
          continue;
        }

        if ((int(a.c) == 36) || (int(a.c) == 0)){
          if (debug){
            cout << " INSIDE $ BASE DETECTION !!! " << a.c << endl;
          }
          continue;
        }

        // in RLZ we comment all the above code, up to the empty line
        // if (a.next_p == 1) {
        //   if (debug){
        //     cout << " CYCLING THROUGH BWT !! STOP BEFORE GETTING TO FIRST CHAR !!" << endl;
        //   }
        //   break;
        // }

        // cout << "## Pop i = " << a.i << " c=" << a.c << " and c =" << c << " , actual base = "<< a.b << " max Ar=" << max_Ar << " p =" << a.next_p << " next c=" << idxBWT->at( a.next_p ) << ", LF(p)=" << p << "  j=" << a.num_snps << " aln len="<< matches << " iprime =" << a.num_gapo << " L=" << a.L << " U=" << a.U << endl;
        nhashtags ++;
        if (debug){
          // cout << repeat(s, a.aln_length+2) << " ALLOWING INNER MATCHES" << endl;
          cout << "## " << " ALLOWING INNER MATCHES" << endl;
        }
        // LZ - start
				// iprime = n - max_Ar;
        // LZ - end

        // RLZ - begin
        iprime = nR - max_Ar;
        // RLZ - end

        a.num_gapo = iprime;
        j = a.num_snps;
        // c = idxBWT->at( p );
        c = idxBWTS->at( p );

        if ((int(c) == 0) || (int(c) == 36)){
          continue;
        }

        if (debug){
        // cout << repeat(s, a.aln_length+2) << " IN c=" << c << " with iprime=" << iprime << " and i=" << a.i+allow_diff << ", p=" << p << endl;
        cout << "### " << " IN c=" << c << " with iprime=" << iprime << " and i=" << a.i+allow_diff << ", p=" << p << endl;
        }

        next_p =  a.next_p;
        // p = a.next_p;

        uchar alphabet[] = {'A', 'C', 'G', 'T'};
        struct query_bases_t qbases[4];

        for(int i=0; i < 4; i++){
          interval = idxBWT->exact_match( alphabet[i], a.L, a.U);
  				sp = interval.first;
  				ep = interval.second;

          qbases[i].b = alphabet[i];
          qbases[i].L = sp;
          qbases[i].U = ep;
          qbases[i].maxAr =RMQ.query(sp, ep-1);
          // cout << "letter = " << alphabet[i] << " with L = " << sp << " and U = " << ep << " and max Ar=" <<   qbases[i].maxAr << endl;
        }

				// j = j + 1;
        int is_mm;
        // int allow_diff = 1;
        // allow_diff = allow_diff - a.num_mm;
        // cout << "PROXIMO ALLOW DIFF = " << allow_diff << endl;

        if (allow_diff > k_seed){

          if (debug){
            // cout << repeat(s, a.aln_length+2) << "&&& BEGIN MISMATCH STYLE ALLOWED !!!" << endl;
            // cout << repeat(s, a.aln_length+2) << " max Ar = " << max_Ar << " n-(i+1) = " << n-(a.i+allow_diff+1) << " allow_diff =" << allow_diff << endl;

            cout << "## "<< "&&& BEGIN MISMATCH STYLE ALLOWED !!!" << endl;
            cout << "## "<< " max Ar = " << max_Ar << " n-(i+1) = " << n-(a.i+allow_diff+1) << " allow_diff =" << allow_diff << endl;

          }

            for(int i=0; i < 4; i++){

              if (debug) {
            			cout << "BEFORE CONDITION BEGIN" << endl;
            			// cout << repeat(s, a.aln_length+2) << " Pushed base " << qbases[i].b << " with " << a.num_mm << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << "Allow_diff = " << allow_diff << endl;
                  cout << "## " << " Pushed base " << qbases[i].b << " with " << a.num_mm << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << "Allow_diff = " << allow_diff << endl;
            			cout << "BEFORE CONDITION END" << endl;
            		  }

              // LZ
              // if ((qbases[i].L < qbases[i].U) && (max_Ar >= n-(a.i+1)) ){
              // RLZ
              if ( (qbases[i].L < qbases[i].U) ){

                // cout << "AAAAQUI " << n-(a.i+1)<< endl;
                is_mm = 0;

                if (c != qbases[i].b){
                  is_mm = 1;
                }

                if (is_mm == 0){

                  if (debug){
                    // cout << repeat(s, a.aln_length+2) << " Pushed base " << qbases[i].b << " with " << a.num_mm << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << "Allow_diff = " << allow_diff << endl;
                    cout << "## " << " Pushed base " << qbases[i].b << " with " << a.num_mm << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << "Allow_diff = " <<  allow_diff << " max(n-(i+1), n-(ref_len+1)) = " << n-(a.i+1) << endl;
                  }

                  // LZ
                  // if (qbases[i].maxAr >= n-(a.i+1))
                  // RLZ
                  if (qbases[i].maxAr >= 0)
                    heap.push(aln_entry_t(a.i, qbases[i].L, qbases[i].U, a.num_mm, (nR-qbases[i].maxAr), p, a.state+1, j+1, a.aln_length+1, c, qbases[i].b));

                } else {

                  if (debug){
                    // cout << repeat(s, a.aln_length+2) << " Pushed base " << qbases[i].b << " with " << a.num_mm+1 << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << endl;
                    cout << "## " << " Pushed base " << qbases[i].b << " with " << a.num_mm+1 << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << " max(n-(i+1), n-(ref_len+1)) = " << std::max(n-(a.i+1), n-(13+1)) << endl;
                  }

                  // LZ
                  // if (qbases[i].maxAr >= n-(a.i+1))
                  // RLZ
                  if (qbases[i].maxAr >= 0)
                    heap.push(aln_entry_t(a.i, qbases[i].L, qbases[i].U, a.num_mm+1, (nR-qbases[i].maxAr), p, 0, j+1, a.aln_length+1, c, qbases[i].b ));
                }

              }
            }

        if (debug){
          // cout << repeat(s, a.aln_length+2) << "&&& END MISMATCH STYLE ALLOWED !!!" << endl;
          cout << "## " << "&&& END MISMATCH STYLE ALLOWED !!!" << endl;
        }
        } else {
          interval = idxBWT->exact_match( c, a.L, a.U);
          sp = interval.first;
          ep = interval.second;
          ulint maxArm = RMQ.query(sp, ep-1);
          // maxAr = RMQ.query(sp, ep-1);
          if (debug){
            // cout << repeat(s, a.aln_length+2) << "BEFORE Pushed MATCHED base " << c << " with " << a.num_mm << " mis" << " L=" << sp << " , U=" << ep<< " , maxAr=" << maxArm << " max_Ar_ori=" << max_Ar << " n-(i+1) = " << n-(a.i+allow_diff+1) << endl;
            // cout << repeat(s, a.aln_length+2) << "@@ TESTING FOR n-(i+1) = " <<  n-(a.i+allow_diff+1) << " BUT " << allow_diff << endl;

            cout << "## " << "BEFORE Pushed MATCHED base " << c << " with " << a.num_mm << " mis" << " L=" << sp << " , U=" << ep<< " , maxAr=" << maxArm << " max_Ar_ori=" << max_Ar << " n-(i+1) = " << n-(a.i+allow_diff+1) << endl;
            cout << "## " << "@@ TESTING FOR n-(i+1) = " <<  n-(a.i+allow_diff+1) << " BUT " << allow_diff << endl;
          }
          // LZ
          // if ( (sp < ep) && (maxArm >= n-(a.i+1)) ){
          // RLZ
          if ( (sp < ep) && (maxArm >= 0) ){

          // if ( (sp < ep) && (maxArm >= (a.num_gapo+1)) ){
            if (debug){
              // cout << repeat(s, a.aln_length+2) << " Pushed MATCHED base " << c << " with " << a.num_mm << " mis" << " L=" << sp << " , U=" << ep<< " , maxAr=" << maxArm << ", allow_diff=" << a.state+1 << " iprime = " << a.num_gapo << " n-(a.i+1) = " << n << endl;
              cout << "## " << " Pushed MATCHED base " << c << " with " << a.num_mm << " mis" << " L=" << sp << " , U=" << ep<< " , maxAr=" << maxArm << ", allow_diff=" << a.state+1 << " iprime = " << a.num_gapo << " n-(a.i+1) = " << n << endl;
            }

            heap.push(aln_entry_t(a.i, sp, ep, a.num_mm, (nR-maxArm), p, a.state+1, j+1, a.aln_length+1, c, c ));
            // allow_diff ++;
            i ++;
          }
        }


        if (debug){
          // cout << repeat(s, a.aln_length) << " Heap SIZE = " << heap.size() << endl;
            cout << "## " << " Heap SIZE = " << heap.size() << endl;
        }

				// perc = (100*j)/n;
				// if (perc > last_perc and (perc%5)) {
				// 	cout << " " << perc << "% done." << endl;
				// 	last_perc=perc;
				// }
        //
				// c = idxBWT->at( p );
				// p = idxBWT->LF( p );
        //
				// interval = idxBWT->exact_match( c, sp, ep );
				// sp = interval.first;
				// ep = interval.second;
        //
				// // cout << "#### sp = " << sp << ", ep = " << ep << endl;
				// max_Ar = RMQ.query(sp, ep - 1);
				// // cout << "Max_Ar = " << max_Ar << endl;
				// // cout << "i= " << i << " j= " << j << " c= " << c << " p= " << p << " sp=" << sp << " ep=" << ep << endl;
        // aln = a;
			}

      // j = j;
      if (debug){
        cout << "J final=" << j << endl;
      }

			if (a.num_gapo - 1 == -1) {
				outfile << j << "\t" << c << "\t-\t-" <<endl;
        cout << j << "\t" << c << "\t-\t-" <<endl;
			} else {
				// cout << "(" << iprime-(j-i) << "," << j-i-1 << ", " << c << ")" << endl;
        if (debug){
          cout << "c=" << a.b << " num_gapo (iprime) = " << a.num_gapo << ", num_snps = " << a.num_snps << " a.i = " << a.i << endl;
        }

        // if ( (a.num_gapo-(a.num_snps-a.i)) == 0) {
        //   cout << "BIRL" << endl;
        //   cout << "(" << (a.num_gapo-(a.num_snps-a.i)) << "," << a.aln_length + 1 << ")" <<  endl;
        // }
        // else {
          if (a.num_mm > 0){
                outfile << j << "\t" << ((a.num_gapo)-(a.num_snps-a.i)) << "\t" << a.aln_length + 1 << "\tM" <<  endl;

                cout << j << "\t" << ((a.num_gapo)-(a.aln_length)) << "\t" << a.aln_length + 1 << "\tM" <<  endl;
          } else {
            outfile << j<< "\t" << ((a.num_gapo)-(a.num_snps-a.i)) << "\t" << a.aln_length + 1 << "\tN" << endl;

            cout << j<< "\t" << ((a.num_gapo)-(a.aln_length)) << "\t" << a.aln_length + 1 << "\tN" << endl;
          }
        // }
			}

      a.state = 0;
			i = j;

		}

    outfile.close();
	}

	printRSSstat();
	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();

	if(total>=3600){

		uint h = total/3600;
		uint m = (total%3600)/60;
		uint s = (total%3600)%60;

		cout << "Total time: " << h << "h " << m << "m " << s << "s" << endl;

	}else if (total>=60){

		uint m = total/60;
		uint s = total%60;

		cout << "Total time: " << m << "m " << s << "s" << endl;

	}else{

		cout << "Total time: " << total << "s" << endl;

	}
}
