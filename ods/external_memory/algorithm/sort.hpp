#pragma once
#include "external_memory/emvector.hpp"

// UNDONE(): This file should implement the sort optimized for the external
// memory model.
//
namespace EM::Algorithm {
  using EM::Vector::Vector;

  template <typename T, typename Compare>
  bool IsSorted(Vector<T>& v, Compare cmp) {
    bool ret = true;
    for (uint64_t i=1; i<v.size(); i++) {
      ret = ret 
          * ( cmp(v[i-1],v[i])
              + (!cmp(v[i], v[i-1]))
            );
    }
    return ret;
  }

  template <typename T, typename Compare>
  void Sort(Vector<T>& v, Compare cmp) {
    Assert(GetNextPowerOfTwo(v.size()) == v.size(), v.size());
    // UNDONE: use the efficient oblivious butterfly sort described in
    // https://arxiv.org/pdf/2008.00332.pdf
    // We are currently just doing bitonic sort, which doesn't has neither the complexity
    // nor the optimal external_memory transfers.
    uint64_t n = v.size();
    for (uint64_t k = 2; k <= n; k *= 2) {
      for (uint64_t j = k/2; j > 0; j /= 2) {
        for (uint64_t i = 0; i < n; ++i) {
          uint64_t l = i ^ j;
          if (l > i) {
            bool llti = cmp(v[l], v[i]);
            bool swap = 
              (((i & k) == 0) * llti) +
              (((i & k) != 0) * (!llti));
            CXCHG(swap, v[i], v[l]);
          }
        }
      }
    }
  }

  template <typename T, uint64_t Z=512>
  void MergeSplit(Vector<TaggedT<T> >& curr,  Vector<TaggedT<T> >& next, uint64_t r0, uint64_t r1, uint64_t w0, uint64_t w1, uint64_t mask) {
    // What we want to do is sort M[r0..r0+Z; r1..r1+Z] using a compare function that sorts and splits in half.
    //

    Vector<TaggedT<T> > toSort(2*Z, TaggedT<T>::DUMMY());
    uint64_t extendedMask = (mask|(mask>>1));
    uint64_t selectBit = mask ^ extendedMask;
    uint64_t reals[2] = {0};
    const uint64_t v0[2] = {r0,r1};
    const uint64_t v1[2] = {w0,w1};
    const uint64_t o0[2] = {0, Z};

    // DEBUG:
    uint64_t localMask = 0;
    //.

    // Pass 1) Put everything in the same array:
    //
    for (uint64_t dir=0; dir<=1; dir++) {
      for (uint64_t i=0; i<Z; i++) {
        toSort[Z*dir+i] = curr[v0[dir]+i];

        // DEBUG ONLY:
        if (!toSort[Z*dir+i].isDummy && localMask == 0) {
          localMask = toSort[Z*dir+i].tag&mask;
        }
        // Assert(!toSort[Z*dir+i].isDummy || ((toSort[Z*dir+i].tag&mask)==localMask));
        //.
        bool isRealDir = 
            (!toSort[Z*dir+i].isDummy)
          * ((toSort[Z*dir+i].tag & selectBit) == dir);
        uint64_t successor = reals[dir]+1;
        CMOV(isRealDir, reals[dir], successor);
      }
    }

    // Pass 2) Update empty blocks offset, they will fall in the end of each array:
    //
    uint64_t maskedLargeTag = mask | ((mask|(mask>>1)) - 1);
    for (uint64_t i=0; i<2*Z; i++) {      
      bool setEmptyBit;
      setEmptyBit = toSort[i].isDummy;
      for (uint64_t dir=0; dir<=1; dir++) {
        bool dirRemaining = setEmptyBit 
                          * (reals[dir]>0);
        CMOV(dirRemaining, toSort[i].tag, maskedLargeTag | (dir ? selectBit : 0));
        CMOV(dirRemaining, reals[dir], reals[dir]-1);
        setEmptyBit = setEmptyBit * (!dirRemaining);
      }
      // Assert(!setEmptyBit);
    }
    // Assert (reals[0] == 0);
    // Assert (reals[1] == 0);

    // Pass 3) Oblivious sort the small array
    //
    auto cmp = [](const TaggedT<T>& b1, const TaggedT<T>& b2) {
      return 
          (b1.tag < b2.tag)
        + ( (b1.tag == b2.tag)
          * (b1.isDummy < b2.isDummy)
          );
    };
    Sort(toSort, cmp);

    // Pass 4) Write back:
    //
    for (uint64_t dir=0; dir<=1; dir++) {
      for (uint64_t i=0; i<Z; i++) {
        next[v1[dir]+i] = toSort[Z*dir+i];
      }
    }
  }


  template <typename T, uint64_t Z=512>
  void BucketObliviousSort_Internal(Vector<TaggedT<T> >& b, uint64_t keyBits=64) {
    Assert((2*b.size())%Z == 0);
    uint64_t N = b.size();
    uint64_t B = GetNextPowerOfTwo((2*b.size())/Z);
    uint64_t lB = GetLogBaseTwo(B);
    
    Vector<TaggedT<T> > tv(B*Z);
    Vector<TaggedT<T> > tvNext(B*Z);

    
    // 1) First divide into B groups:
    uint64_t i0 = 0;
    uint64_t Z0 = N/B;
    uint64_t N0 = N - Z0*B; // number of blocks that will have one extra value than the others

    for (uint64_t i=0; i<B*Z; i+=Z) {
      uint64_t j=0;
      for (; j<Z0 + (N0 > 0 ? 1 : 0); j++) {
        tv[i+j] = b[i0++];
      }
      if (N0 > 0) N0--;
      for (; j<Z; j++) {
        tv[i+j] = TaggedT<T>::DUMMY();
      }
    }
    Assert(N0 == 0, N0);
    Assert(i0 == N, i0, N);

    uint64_t maskShift = (keyBits+lB-1)/lB;
    uint64_t mask = ((1ULL<<maskShift)-1ULL) << (keyBits);

    // Sort by key MSB's iteratively:
    //
    for (uint64_t i=0; i<lB-1; i++) {
      Vector<TaggedT<T> >& curr = (i%2 == 0) ? tv : tvNext;
      Vector<TaggedT<T> >& next = (i%2 == 0) ? tvNext : tv;

      for (uint64_t j=0; j<B/2; j++) { // 2|B
        uint64_t jp = (j / (1<<i)) * (1<<i);
        MergeSplit(curr, next, (j+jp)*Z, (j+jp+(1<<i))*Z,(2*j)*Z, (2*j+1)*Z, mask);
      }

      mask |= (mask>>maskShift);
    }
  }

  template <typename T>
  void BucketObliviousRandomPermutation(Vector<T>& b) {
    uint64_t intermidiateSize = GetNextPowerOfTwo(2*b.size());
    Vector<TaggedT<T> > tv(intermidiateSize, TaggedT<T>::DUMMY());
    
    for (uint64_t i=0; i<intermidiateSize; i++) {
      tv[i].tag = UniformRandom(b.size()-1);
      if (i < b.size()) {
        tv[i].isDummy = false;
        tv[i].v = b[i];
      } else {
        tv[i].tag = UniformRandom(b.size(), intermidiateSize);
      }
    }

    BucketObliviousSort_Internal(tv, GetLogBaseTwo(intermidiateSize));

    for (uint64_t i=0; i<b.size(); i++) {
      b[i] = tv[i].v;
    }
  }

  template <typename T>
  void BucketObliviousSort(Vector<TaggedT<T> >& b, uint64_t keyBits=64) {
    BucketObliviousRandomPermutation(b);
    BucketObliviousSort_Internal(b, keyBits);
  }
}