#include "testutils.hpp"
#include "external_memory/algorithm/shift.hpp"

TEST(ShiftByConstant, correctness) {
    for (uint64_t N=1; N<64; N++) {
        for (uint64_t shift=0; shift<= N; shift++) {
            std::cerr << "N=" << N << " shift=" << shift << std::endl;
            std::vector<uint32_t> orig;
            for (uint64_t i=0; i<N; i++) {
                orig.push_back(i);
            }
            maybeShiftByK(false, orig, shift);
            for (uint64_t i=0; i<N; i++) {
                ASSERT_EQ(orig[i], i);
            }
            maybeShiftByK(true, orig, shift);
            for (uint64_t i=0; i<N; i++) {
                std::cerr << orig[i] << " ";
            }
            std::cerr << std::endl;
            for (uint64_t i=0; i<N-shift; i++) {
                ASSERT_EQ(orig[i], i+shift);
            }
            for (uint64_t i=N-shift; i<N; i++) {
                ASSERT_EQ(orig[i], i-(N-shift));
            }
        }
    }
}

TEST(ShiftByHiddenConstant, correctness) {
    for (uint64_t N=1; N<64; N++) {
        for (uint64_t shift=0; shift<= N; shift++) {
            std::cerr << "N=" << N << " shift=" << shift << std::endl;
            std::vector<uint32_t> orig;
            for (uint64_t i=0; i<N; i++) {
                orig.push_back(i);
            }

            cyclicShift(orig, shift,20);
            for (uint64_t i=0; i<N; i++) {
                std::cerr << orig[i] << " ";
            }
            std::cerr << std::endl;
            for (uint64_t i=0; i<N-shift; i++) {
                ASSERT_EQ(orig[i], i+shift);
            }
            for (uint64_t i=N-shift; i<N; i++) {
                ASSERT_EQ(orig[i], i-(N-shift));
            }
        }
    }
}

TEST(Knapsack, correctness) {
    for (uint64_t instance=0; instance<1000; instance++) {
        uint64_t capacity;
        uint64_t N=2 + (rand()%1000);
        uint64_t W=1 + (rand() % 1234);
        uint64_t V=1+(rand() % (2<<(rand() % 10)));
        uint64_t C=1+(rand() % (2<<(rand() % 10)));
        std::vector<uint64_t> weights;
        std::vector<uint64_t> values;
        for (uint64_t i=0; i<N; i++) {
            weights.push_back(rand() % W);
            values.push_back(rand() % V);
        }

        uint64_t k1 = knapsack_val<false>(weights, values, C, W);
        uint64_t k2 = knapsack_val<true>(weights, values, C, W);
        ASSERT_EQ(k1, k2);
    }
}