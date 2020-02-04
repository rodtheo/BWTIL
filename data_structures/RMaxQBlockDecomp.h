#ifndef RMAQBLOCKDECOMP_H_
#define RMAQBLOCKDECOMP_H_

#include <iostream>

#include "IndexedBWT.h"

using namespace bwtil;
using namespace std;

namespace bwtil {

  class RMaxQBlockDecomp {

    public:

      RMaxQBlockDecomp(IndexedBWT* idxBWT){

        this->idxBWT = idxBWT;
        n = idxBWT->length();
        b = ceil(log2(n));
        if (b<1) b=1;

        block = vector<ulint>((int)(n/b), 0);

        ulint blk_idx = -1;
        uint blk_sz = b;

        cout << "Block size = " << block.size() << endl;
        ulint sa_val;

        for (int i=0; i < n; i++){
          if ((i % b) == 0){
            blk_idx ++;
          }
          sa_val = idxBWT->convertToTextCoordinate(i);
          if (block[blk_idx] < sa_val){
            block[blk_idx] = sa_val;
          }
        }
      }

      ulint query(ulint l, ulint r){
        ulint res = 0;

        int blk_idx_l = (int)(l/b);
        int blk_idx_r = (int)(r/b);

        ulint sa_val;
        while( (l<r) && ((l%b) != 0) && (l != 0)){
          sa_val = idxBWT->convertToTextCoordinate(l);
          if (sa_val > res){
            res = sa_val;
          }
          l++;
        }

        while((l + b) <= r) {
          if (block[(int)(l/b)] > res)
            res = block[(int)(l/b)];
          l = l + b;
        }

        while((l <= r)){
          sa_val = idxBWT->convertToTextCoordinate(l);
          if (sa_val > res)
            res = sa_val;
          l ++;
        }

        return res;
      }



      IndexedBWT* idxBWT;
      uint n;
      uint b; // size of a block = log2 n
      vector<ulint> block;

  };
} // namespace data_structures
#endif
