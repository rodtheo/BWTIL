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
	uint32_t num_mm:8, num_gapo:8, num_snps:8;
	uint32_t score:8, state:8, aln_length:8;
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

	if(argc != 4 and argc != 3){
		cout << "*** succinct FM-index data structure : a wavelet-tree based uncompressed FM index ***\n";
		cout << "Usage: sFM-index option file [pattern]\n";
		cout << "where:\n";
		cout <<	"- option = build|search. \n";
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

    auto t1 = high_resolution_clock::now();

	succinctFMIndex SFMI;

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

		cout << "Loading succinct FM-index from file " << in << endl;
		SFMI = succinctFMIndex::loadFromFile(in);
		cout << "Done." << endl;
		IndexedBWT* idxBWT = SFMI.get_idxBWTPtr();
		RMaxQBlockDecomp RMQ(idxBWT);

		// vector<factor> enc;

		ulint i = 0;
		uint j = 0;
		ulint p = 1;
		uchar c;
		ulint cidx;
		pair<ulint, ulint> interval;
		ulint max_Ar_idx, max_Ar;
		uint iprime;
		ulint n = SFMI.textLength();
    cout << "TEXT LENGTH = " << n << endl;
		ulint sp, ep;
    ulint next_p;
		int perc, last_perc= -1;
		// cout << "Letter " << c << endl;
		priority_queue<aln_entry_t, vector<aln_entry_t>, CompareAlnScores> heap;

    int best_score, num_best, max_entries, total_entries, last_num_entries, alns_get_num_entries;

    int diff_left, best_diff;
    uint32_t matches, max_diff;

    cout << "###################################################" << endl;
    cout << "INITIALIZING ALGORITHM" << endl;
    cout << "###################################################\n\n" << endl;

    string s = "##";
    int nhashtags;

		while (i < (n-1)){

      cout << "###################################################" << endl;
      cout << " FOR i = " << i << endl;
      nhashtags = 1;
			// cout << (int)(n/10.) << endl;
			// perc = (100*j)/n;
			// if (perc > last_perc and (perc%5)) {
			// 	cout << " " << perc << "% done." << endl;
			// 	last_perc=perc;
			// }

			j = i + 1;



			c = idxBWT->at( p );
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
      while (!heap.empty()){

        a = heap.top();
        heap.pop();

        diff_left = max_diff - a.num_mm;
        sp = a.L;
        ep = a.U;
        matches = a.aln_length;

        c = a.c;
        p = idxBWT->LF( a.next_p );

        max_Ar = RMQ.query(sp, ep-1);

        cout << "Pop i = " << a.i << " c=" << a.c << " and c =" << c << " , actual base = "<< a.b << ", max Ar=" << max_Ar << " p =" << a.next_p << " next c=" << idxBWT->at( a.next_p ) << ", LF(p)=" << p << "  j=" << a.num_snps <<", j new=" << a.i+a.aln_length << " aln len="<< matches << " iprime =" << a.num_gapo << " L=" << a.L << " U=" << a.U << endl;

        if (diff_left < 0)
          continue;

        if (sp > (ep-1))
          continue;

        if (max_Ar < n-(i+1))
          continue;

        if (a.num_snps >= (n-1))
          continue;

        // cout << "## Pop i = " << a.i << " c=" << a.c << " and c =" << c << " , actual base = "<< a.b << " max Ar=" << max_Ar << " p =" << a.next_p << " next c=" << idxBWT->at( a.next_p ) << ", LF(p)=" << p << "  j=" << a.num_snps << " aln len="<< matches << " iprime =" << a.num_gapo << " L=" << a.L << " U=" << a.U << endl;
        nhashtags ++;
        cout << repeat(s, a.aln_length+2) << " ALLOWING INNER MATCHES" << endl;

				iprime = n - max_Ar;
        a.num_gapo = iprime;
        j = a.num_snps;
        // c = idxBWT->at( p );
        c = idxBWT->at( p );
        cout << repeat(s, a.aln_length+2) << " IN c=" << c << " with iprime=" << iprime << " and i=" << a.i << ", p=" << p << endl;
        next_p =  a.next_p;
        // p = a.next_p;

        uchar alphabet[] = {'G', 'O', 'L'};
        struct query_bases_t qbases[3];

        for(int i=0; i < 3; i++){
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
        int allow_diff = 1;

        if (allow_diff){
            cout << repeat(s, a.aln_length+2) << " max Ar = " << max_Ar << " n-(i+1) = " << n-(a.i+1) << endl;
            for(int i=0; i < 3; i++){
              if ((qbases[i].L < qbases[i].U) && (max_Ar >= n-(a.i+1))){
                // cout << "AAAAQUI" << endl;
                is_mm = 0;
                if (c != qbases[i].b){
                  is_mm = 1;
                }

                if (is_mm == 0){
                  cout << repeat(s, a.aln_length+2) << " Pushed base " << qbases[i].b << " with " << a.num_mm << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << endl;
                  if (qbases[i].maxAr >= n-(a.i+1))
                    heap.push(aln_entry_t(a.i, qbases[i].L, qbases[i].U, a.num_mm, (n-qbases[i].maxAr), p, a.state, j+1, a.aln_length+1, c, qbases[i].b));
                } else {
                  cout << repeat(s, a.aln_length+2) << " Pushed base " << qbases[i].b << " with " << a.num_mm+1 << " mis" << " L=" << qbases[i].L << " , U=" << qbases[i].U << " , maxAr=" << qbases[i].maxAr << endl;
                  if (qbases[i].maxAr >= n-(a.i+1))
                    heap.push(aln_entry_t(a.i, qbases[i].L, qbases[i].U, a.num_mm+1, (n-qbases[i].maxAr), p, a.state, j+1, a.aln_length+1, c, qbases[i].b ));
                }

              }
            }

        }

        cout << repeat(s, a.aln_length) << " Heap SIZE = " << heap.size() << endl;

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
      cout << "J final=" << j << endl;

			if (iprime - 1 == -1) {
				cout << "(" << c << ")"<< endl;
			} else {
				// cout << "(" << iprime-(j-i) << "," << j-i-1 << ", " << c << ")" << endl;
        cout << "(" << (a.num_gapo-(a.num_snps-a.i))-1 << "," << a.aln_length + 1 << ")" <<  endl;
			}

			i = j;

		}


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
