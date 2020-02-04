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

#include "../../data_structures/succinctFMIndex.h"
#include "../../data_structures/RMaxQBlockDecomp.h"

using namespace bwtil;
using namespace std;


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
		ulint sp, ep;
		int perc, last_perc= -1;
		// cout << "Letter " << c << endl;
		while (i < n){
			// cout << (int)(n/10.) << endl;
			perc = (100*j)/n;
			if (perc > last_perc and (perc%5)) {
				cout << " " << perc << "% done." << endl;
				last_perc=perc;
			}

			j = i + 1;



			c = idxBWT->at( p );
			p = idxBWT->LF( p );

			interval = SFMI.arrayC(c);
			// cout << "LETTER =" << c << ", INTERVAL FIRST = " << interval.first << ", SA=" << idxBWT->convertToTextCoordinate(interval.first) << "| SECOND = " << interval.second << ", SA=" << idxBWT->convertToTextCoordinate(interval.second-1) << endl;
			// cout << "RMQ[sp,ep] = " << RMQ.query(interval.first, interval.second-1) << endl;

			sp = interval.first;
			ep = interval.second;

			max_Ar = RMQ.query(sp, ep-1);
			iprime = 0;

			// cout << "Max_Ar = " << max_Ar << endl;
			// cout << "i= " << i << " j= " << j << " c= " << c << " p= " << p << " sp=" << sp << " ep=" << ep << endl;

			while ((j + 1 < n) && (sp <= ep) && (max_Ar >= (n-(i+1)))){

				iprime = n - max_Ar;
				j = j + 1;

				perc = (100*j)/n;
				if (perc > last_perc and (perc%5)) {
					cout << " " << perc << "% done." << endl;
					last_perc=perc;
				}

				c = idxBWT->at( p );
				p = idxBWT->LF( p );

				interval = idxBWT->exact_match( c, sp, ep );
				sp = interval.first;
				ep = interval.second;

				// cout << "#### sp = " << sp << ", ep = " << ep << endl;
				max_Ar = RMQ.query(sp, ep - 1);
				// cout << "Max_Ar = " << max_Ar << endl;
				// cout << "i= " << i << " j= " << j << " c= " << c << " p= " << p << " sp=" << sp << " ep=" << ep << endl;


			}

			// if (iprime - (j-i) == -1) {
			// 	cout << "(" << c << ")"<< endl;
			// } else {
			// 	cout << "(" << iprime-(j-i) << "," << j-i-1 << ", " << c << ")" << endl;
			// }

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
