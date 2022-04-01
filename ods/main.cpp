#include "otree/otree.hpp"
#include "oram/pathoram/oram.hpp"
#include <cassert>
#include <iostream>

int main() {
  using ORAMClient = _ORAM::PathORAM::ORAMClient::ORAMClient<_OBST::Node,ORAM__Z,false,false,4>;
  using OramClient = _OBST::OramClient::OramClient<ORAMClient>;
  _OBST::OBST::OBST<OramClient> x(4);
  x.Insert(0, 32);
  x.Insert(1, 10);
  x.Insert(2, 10);
  _OBST::V a;
  std::cout << x.Get(0, a) << std::endl;
  std::cout << a << std::endl;
  return 0;
}
