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

/*
 * PartialSums.h
 *
 *  Created on: Jun 22, 2014
 *  Author: nicola
 *
 *  Description: maintains sigma partial sums; when incrementing counter i, all counters j>=i are incremented.
 *  Complexity: the structure is a packed B-tree. increment/access time = O( log_d(sigma), where d=w/log w and w is the word size (w=64) )
 *
 */

#ifndef PARTIALSUMS_H_
#define PARTIALSUMS_H_

#include "../common/common.h"

namespace bwtil {

class PartialSums {

public:

	PartialSums(){nr_of_nodes=0;};

	//the structure mantains one counter for each symbol in {0,1,...,sigma-1}.
	PartialSums(uint sigma, ulint n){//size of the alphabet and maximum number to be stored in a counter

		nr_of_nodes=0;

		if(n==0){//empty counters
			empty=true;
			return;
		}

		base_counter = 0;

		empty=false;
		this->sigma=sigma;

		log2n = ceil(log2(n+1));

		uint w=64;
		d = w/log2n;

		ones = (ulint)0;
		for(uint i=0;i<d;i++)
			ones = ones | (((ulint)1)<<(i*log2n));

		nr_of_leafs = sigma/d + (sigma%d==0?0:1);//leafs

		uint nr_of_nodes_in_level = nr_of_leafs;
		nr_of_nodes = nr_of_nodes_in_level;

		//round nr_of_nodes_in_level to the next power of d+1

		uint nodes_pow = 1;
		while(nodes_pow<nr_of_nodes_in_level)
			nodes_pow *= (d+1);

		nr_of_nodes_in_level = nodes_pow;
		uint height=1;

		while(nr_of_nodes_in_level>1){

			height++;
			nr_of_nodes_in_level = nr_of_nodes_in_level/(d+1);
			nr_of_nodes += nr_of_nodes_in_level;

		}

		/*cout<<"log2n="<<(uint)log2n<<endl;
		cout << "d="<<(uint)d<<endl;
		cout << "nr_of_leafs="<<(uint)nr_of_leafs<<endl;
		cout << "height="<<height<<endl;
		cout << "nr_of_nodes="<<(uint)nr_of_nodes<<endl;*/

		nodes = vector<ulint>(nr_of_nodes);
		for(uint i=0;i<nr_of_nodes;i++)//reset all counters
			nodes[i]=0;

	}

	string toString(){

	#ifdef DEBUG
		if(empty){
			cout << "ERROR (PartialSums): toString() called on empty counter\n";
			exit(0);
		}
	#endif

		stringstream ss;

		 for(uint i =0;i<sigma;i++)
			ss << getCount(i) << " ";

		 return ss.str();

	}

	void increment(symbol s){

	#ifdef DEBUG
		if(empty){
			cout << "ERROR (PartialSums): increment() called on empty counter\n";
			exit(0);
		}
	#endif

	#ifdef DEBUG
		if(s>=sigma){
			cout << "ERROR (PartialSums): symbol " << s << " not in alphabet.\n";
			exit(0);
		}
	#endif

		uint current_node = (nr_of_nodes - nr_of_leafs) + (s/d);//offset leafs + leaf number
		uint offset_in_node = s%d;//number of the counter inside the node

		incrementFrom(&nodes, current_node, offset_in_node);//increment counters in the leaf

		while(current_node>0){//repeat while current node is not the root

			offset_in_node = childNumber(current_node);
			current_node = parent(current_node);

			incrementFrom(&nodes, current_node, offset_in_node);

		}

	}

	//number of symbols <s inserted. Eg: getCount(2) returns number of occurrencies of 0 and 1. getCount(0) returns always 0.
	ulint getCount(symbol s){

	#ifdef DEBUG
		if(empty){
			cout << "ERROR (PartialSums): getCount() called on empty counter\n";
			exit(0);
		}
	#endif

		if(s==0)
			return base_counter;

		s--;

	#ifdef DEBUG
		if(s>=sigma){
			cout << "ERROR (PartialSums): symbol " << s << " not in alphabet.\n";
			exit(0);
		}
	#endif

		uint current_node = (nr_of_nodes - nr_of_leafs) + (s/d);//offset leafs + leaf number

		uint offset_in_node = s%d;//number of the counter inside the node

		ulint count = getCounterNumber(nodes[current_node],offset_in_node);

		while(current_node>0){//repeat while current node is not the root

			offset_in_node = childNumber(current_node);
			current_node = parent(current_node);

			if(offset_in_node>0)
				count += getCounterNumber(nodes[current_node],offset_in_node-1);

		}

		return count + base_counter;

	}

	void setBaseCounter(){base_counter=1;};

	uint bitSize(){//return size in bits

		return CHAR_BIT*(sizeof(this) + nr_of_nodes*sizeof(ulint));

	}

private:

	//increment counters in node by 1 starting from counter number i (from left)
	inline void incrementFrom(vector<ulint> * nodes, ulint node_nr, uint i){

	#ifdef DEBUG
		if(empty){
			cout << "ERROR (PartialSums): incrementFrom() called on empty counter\n";
			exit(0);
		}
	#endif

	#ifdef DEBUG
		if(i>d){
			cout << "ERROR (PartialSums): incrementing counter i>d (" << i << ">" << d << ")\n";
			exit(0);
		}
	#endif

		ulint MASK = ~((ulint)0);

		if(i>0) MASK = ((ulint)1<<((d-i)*log2n))-1;

		//printWord(MASK&ones);

		nodes->at(node_nr) += (MASK&ones);

	}

	//get value of the counter number 0<=i<d in the node n
	inline ulint getCounterNumber(ulint n, uint i){

	#ifdef DEBUG
		if(empty){
			cout << "ERROR (PartialSums): getCounterNumber() called on empty counter\n";
			exit(0);
		}
	#endif

	#ifdef DEBUG
		if(i>=d){
			cout << "ERROR (PartialSums): get counter i>=d (" << i << ">=" << d << ")\n";
			exit(0);
		}
	#endif

		ulint MASK = (((ulint)1)<<log2n)-1;

		return (n>>((d-i-1)*log2n))&MASK;

	}

	//return child number i of node n
	inline uint16_t child(uint16_t n, uint8_t i){//return child number 0<=i<=d of node n

	#ifdef DEBUG
		if(empty){
			cout << "ERROR (PartialSums): child() called on empty counter\n";
			exit(0);
		}
	#endif

	#ifdef DEBUG
		if(i>d){
			cout << "ERROR (PartialSums): child number i>d (" << i << ">" << d << ")\n";
			exit(0);
		}
	#endif

		return (n*(d+1))+i+1;

	}

	//return parent of node n
	inline uint16_t parent(uint16_t n){//return parent of node n

	#ifdef DEBUG
		if(empty){
			cout << "ERROR (PartialSums): parent() called on empty counter\n";
			exit(0);
		}
	#endif

	#ifdef DEBUG
		if(n>=nr_of_nodes){
			cout << "ERROR (PartialSums): parent of inexistent node " << n <<"\n";
			exit(0);
		}
	#endif

		return (n-(((n-1)%(d+1))+1))/(d+1);

	}

	//return which children number is n in his parent
	inline uint8_t childNumber(uint16_t n){

	#ifdef DEBUG
		if(n>=nr_of_nodes){
			cout << "ERROR (PartialSums): parent of inexistent node " << n << "\n";
			exit(0);
		}
	#endif

		return (n-1)%(d+1);

	}

	uint8_t sigma;
	uint8_t nr_of_leafs;
	uint8_t log2n;//number of bits of each counter
	uint8_t d;//number of counters per node = 64/log2(n)

	bool base_counter;//added to each count (1 only in the last context in the text, to count for one terminator)
	bool empty;

	uint16_t nr_of_nodes;

	ulint ones;//1^d in base 2^log2n
	vector<ulint> nodes;//each node is a 64-bits word and stores d counters

};


} /* namespace bwtil */
#endif /* PARTIALSUMS_H_ */
