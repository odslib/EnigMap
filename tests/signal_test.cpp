#include "doram.hpp"
#include "oram.hpp"
#include "otree.hpp"
#include "recdoram.hpp"
#include "recoram.hpp"
#include "signal.hpp"
#include "utils.hpp"
#include <chrono>
#include <gtest/gtest.h>
#include <map>
using namespace _OBST;
using namespace _ORAM;
using namespace std;

template <typename T> class TestSignal : public ::testing::Test {
protected:
  uint64_t num_users_ = 100000;
  std::unique_ptr<Signal<T>> signal_;

  void SetUp() override {
    signal_ = std::make_unique<Signal<T>>(num_users_);
    signal_->RegisterUser(213);
    signal_->RegisterUser(410);
    signal_->RegisterUser(829);
  }
};
// Note that this should work with every client
// using MyTypes = ::testing::Types<RecDoramClient, RecOramClient,
// DoramPosClient, OramPosClient>;
using MyTypes = ::testing::Types<RecDoramClient>;
TYPED_TEST_SUITE(TestSignal, MyTypes);
TYPED_TEST(TestSignal, Correctness) {
  EXPECT_EQ(1, this->signal_->QueryUser(213));
  EXPECT_EQ(1, this->signal_->QueryUser(410));
  EXPECT_EQ(1, this->signal_->QueryUser(829));
  EXPECT_EQ(0, this->signal_->QueryUser(122));
}

/*
template <typename T>
class TestSignalPerf10000 : public ::testing::Test {
 protected:
  uint64_t num_users_ = 10000;
  std::unique_ptr<Signal<T>> signal_;

  void SetUp() override {
    signal_ = std::make_unique<Signal<T>>(num_users_);
    signal_->RegisterUser(213);
    signal_->RegisterUser(410);
    signal_->RegisterUser(829);
  }
};
TYPED_TEST_SUITE(TestSignalPerf10000, MyTypes);
TYPED_TEST(TestSignalPerf10000, Keke) {
  using clock = std::chrono::system_clock;
  using ms = std::chrono::duration<double, std::milli>;
  const auto before = clock::now();

  std::vector<uint64_t> cands = {213, 410, 829, 251, 451, 751};
  srand(time(0));
  for (uint64_t i = 0; i < 1000; ++i) {
    uint64_t ind = rand() % cands.size();
    uint64_t uid = cands[ind];
    bool res = this->signal_->QueryUser(uid);
    EXPECT_EQ(res, uid == 213 || uid == 410 || uid == 829);
  }

  const ms duration = clock::now() - before;
  std::cout << "It took " << duration.count() << " ms.\n";
}
*/

template <typename T> class TestSignalPerf100000 : public ::testing::Test {
protected:
  uint64_t num_users_ = 100000;
  std::unique_ptr<Signal<T>> signal_;

  void SetUp() override {
    signal_ = std::make_unique<Signal<T>>(num_users_);
    signal_->RegisterUser(213);
    signal_->RegisterUser(410);
    signal_->RegisterUser(829);
  }
};
TYPED_TEST_SUITE(TestSignalPerf100000, MyTypes);
TYPED_TEST(TestSignalPerf100000, Keke) {
  using clock = std::chrono::system_clock;
  using ms = std::chrono::duration<double, std::milli>;
  const auto before = clock::now();

  std::vector<uint64_t> cands = {213, 410, 829, 251, 451, 751};
  srand(time(0));
  for (uint64_t i = 0; i < 1000; ++i) {
    uint64_t ind = rand() % cands.size();
    uint64_t uid = cands[ind];
    bool res = this->signal_->QueryUser(uid);
    EXPECT_EQ(res, uid == 213 || uid == 410 || uid == 829);
  }

  const ms duration = clock::now() - before;
  std::cout << "It took " << duration.count() << " ms.\n";
}

template <typename T> class TestSignalPerf1000000 : public ::testing::Test {
protected:
  uint64_t num_users_ = 1000000;
  std::unique_ptr<Signal<T>> signal_;

  void SetUp() override {
    signal_ = std::make_unique<Signal<T>>(num_users_);
    signal_->RegisterUser(213);
    signal_->RegisterUser(410);
    signal_->RegisterUser(829);
  }
};
TYPED_TEST_SUITE(TestSignalPerf1000000, MyTypes);
TYPED_TEST(TestSignalPerf1000000, Keke) {
  using clock = std::chrono::system_clock;
  using ms = std::chrono::duration<double, std::milli>;
  const auto before = clock::now();

  std::vector<uint64_t> cands = {213, 410, 829, 251, 451, 751};
  srand(time(0));
  for (uint64_t i = 0; i < 1000; ++i) {
    uint64_t ind = rand() % cands.size();
    uint64_t uid = cands[ind];
    bool res = this->signal_->QueryUser(uid);
    EXPECT_EQ(res, uid == 213 || uid == 410 || uid == 829);
  }

  const ms duration = clock::now() - before;
  std::cout << "It took " << duration.count() << " ms.\n";
}

template <typename T> class TestSignalPerf10000000 : public ::testing::Test {
protected:
  uint64_t num_users_ = 10000000;
  std::unique_ptr<Signal<T>> signal_;

  void SetUp() override {
    signal_ = std::make_unique<Signal<T>>(num_users_);
    signal_->RegisterUser(213);
    signal_->RegisterUser(410);
    signal_->RegisterUser(829);
  }
};
TYPED_TEST_SUITE(TestSignalPerf10000000, MyTypes);
TYPED_TEST(TestSignalPerf10000000, Keke) {
  using clock = std::chrono::system_clock;
  using ms = std::chrono::duration<double, std::milli>;
  const auto before = clock::now();

  std::vector<uint64_t> cands = {213, 410, 829, 251, 451, 751};
  srand(time(0));
  for (uint64_t i = 0; i < 1000; ++i) {
    uint64_t ind = rand() % cands.size();
    uint64_t uid = cands[ind];
    bool res = this->signal_->QueryUser(uid);
    EXPECT_EQ(res, uid == 213 || uid == 410 || uid == 829);
  }

  const ms duration = clock::now() - before;
  std::cout << "It took " << duration.count() << " ms.\n";
}

template <typename T> class TestSignalPerf100000000 : public ::testing::Test {
protected:
  uint64_t num_users_ = 100000000;
  std::unique_ptr<Signal<T>> signal_;

  void SetUp() override {
    signal_ = std::make_unique<Signal<T>>(num_users_);
    signal_->RegisterUser(213);
    signal_->RegisterUser(410);
    signal_->RegisterUser(829);
  }
};
TYPED_TEST_SUITE(TestSignalPerf100000000, MyTypes);
TYPED_TEST(TestSignalPerf100000000, Keke) {
  using clock = std::chrono::system_clock;
  using ms = std::chrono::duration<double, std::milli>;
  const auto before = clock::now();

  std::vector<uint64_t> cands = {213, 410, 829, 251, 451, 751};
  srand(time(0));
  for (uint64_t i = 0; i < 1000; ++i) {
    uint64_t ind = rand() % cands.size();
    uint64_t uid = cands[ind];
    bool res = this->signal_->QueryUser(uid);
    EXPECT_EQ(res, uid == 213 || uid == 410 || uid == 829);
  }

  const ms duration = clock::now() - before;
  std::cout << "It took " << duration.count() << " ms.\n";
}

template <typename T> class TestSignalPerf1000000000 : public ::testing::Test {
protected:
  uint64_t num_users_ = 1000000000;
  std::unique_ptr<Signal<T>> signal_;

  void SetUp() override {
    signal_ = std::make_unique<Signal<T>>(num_users_);
    signal_->RegisterUser(213);
    signal_->RegisterUser(410);
    signal_->RegisterUser(829);
  }
};
TYPED_TEST_SUITE(TestSignalPerf1000000000, MyTypes);
TYPED_TEST(TestSignalPerf1000000000, Keke) {
  using clock = std::chrono::system_clock;
  using ms = std::chrono::duration<double, std::milli>;
  const auto before = clock::now();

  std::vector<uint64_t> cands = {213, 410, 829, 251, 451, 751};
  srand(time(0));
  for (uint64_t i = 0; i < 1000; ++i) {
    uint64_t ind = rand() % cands.size();
    uint64_t uid = cands[ind];
    bool res = this->signal_->QueryUser(uid);
    EXPECT_EQ(res, uid == 213 || uid == 410 || uid == 829);
  }

  const ms duration = clock::now() - before;
  std::cout << "It took " << duration.count() << " ms.\n";
}
