#include "oram/pathoram/oram.hpp"
#include "otree/otree.hpp"
#include <gtest/gtest.h>
using namespace _ORAM;
using namespace std;

using K = _OBST::K;
using V = _OBST::V;

using ORAMClient_t = typename _ORAM::PathORAM::ORAMClient::ORAMClient<_OBST::Node,ORAM__Z,false,4>;
using OramClient_t = typename _OBST::OramClient::OramClient<ORAMClient_t>;
using OBST = typename _OBST::OBST::OBST<OramClient_t>;

TEST(TestOTreeInit, BasicAssertions) {
// int main() {
  for (uint64_t sz=10; sz<26; sz++) {
    const uint64_t size = (1<<sz);
    const K maxK = size * 100;
    std::cout << size << std::endl;
    EM::Vector::Vector<std::pair<K, V> > emVec(size/2, std::pair<K,V>(0,0));
    for (uint64_t i=0; i<size/2; i++) {
      emVec[i].first = i*100 + (random() % 100);
      emVec[i].second = random() % 100;
    }

    {
      const clock_t t0 = clock();
      OBST t(size,emVec);
      const clock_t t1 = clock();
      double s = (t1 - t0) / (double)CLOCKS_PER_SEC;
      cout << endl;
      cout << "[Report]" << endl;
      cout << "max depth:" << t.maxDepth << endl;
      cout << "max nodes:" << t.maxNodes << endl;
      cout << "oram.N:" << t.oram.N_ << endl;
      cout << "oram.L:" << t.oram.L_ << endl;
      cout << "oram.Z:" << ORAM__Z << endl;
      cout << "s initialization time:" << s << endl;
      cout << "us/node:" << s / (size) * 1000000 << endl;
      cout << "[/Report]" << endl;
      // ASSERT_GE(t.maxNodes, size);
    }
  }
}