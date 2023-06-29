#pragma once
#include <cmath>
#include <vector>

#include "block_for_sort.hpp"

namespace EM::Algorithm {
static double addLogs(double logA, double logB) {
  double bigger = std::max(logA, logB);
  double smaller = std::min(logA, logB);
  if (bigger == -INFINITY) {
    return bigger;
  }
  if (bigger - smaller > 100) {
    return bigger;
  }
  return bigger + log2(1.0 + pow(2, (smaller - bigger)));
}

enum BOUND_TYPE { EXACT, UPPER, LOWER };

static double binary_entropy(double p) {
  Assert(0 < p);
  Assert(p < 1);
  return -p * log2(p) - (1 - p) * log2(1 - p);
}

template <BOUND_TYPE type = EXACT>
static double logCombin(uint64_t k, uint64_t n) {
  Assert(k <= n);
  if (k * 2 > n) {
    k = n - k;
  }

  double res = (lgamma(n + 1) - lgamma(n - k + 1) - lgamma(k + 1)) / log(2);
  return res;
}

template <BOUND_TYPE type = EXACT>
static double binomLogPmf(uint64_t k, uint64_t n, double p) {
  double logCkn = logCombin<type>(k, n);

  if (p > 1e-5) {
    return logCkn + k * log2(p) + (n - k) * log2(1 - p);
  }
  return logCkn + k * log2(p) -
         (n - k) / log(2) * (p + p * p / 2 + p * p * p / 3);
}

template <BOUND_TYPE type = EXACT>
static double binomLogSf(uint64_t k, uint64_t n, double p) {
  double sf = -INFINITY;
  double pmf = binomLogPmf<type>(k, n, p);
  double eps = pmf - 40;
  while (pmf > eps && k < n) {
    ++k;
    pmf += log2((double)(n - k + 1) / k) + log2(p / (1 - p));
    sf = addLogs(sf, pmf);
  }
  return sf;
}

template <BOUND_TYPE type = EXACT>
static double binomLogCdf(uint64_t k, uint64_t n, double p) {
  double sf = -INFINITY;
  double pmf = binomLogPmf<type>(k, n, p);
  double eps = pmf - 40;
  while (pmf > eps) {
    sf = addLogs(sf, pmf);
    if (k == 0) {
      break;
    }
    pmf -= log2((double)(n - k + 1) / k) + log2(p / (1 - p));
    --k;
  }
  return sf;
}

template <BOUND_TYPE type = EXACT>
static double hypergeomLogPmf(uint64_t k, uint64_t M, uint64_t n, uint64_t N) {
  Assert(N <= M);
  if (k > N || k > n) {
    return -INFINITY;
  }

  double logCkn = logCombin<type>(k, n);
  double logCN_kM_n = logCombin<type>(N - k, M - n);
  static size_t cachedM = 0;
  static size_t cachedN = 0;
  static size_t cachedLogMN = 0;
  double logCMN;
  if (M == cachedM && N == cachedN) {
    logCMN = cachedLogMN;
  } else {
    if constexpr (type == EXACT) {
      logCMN = logCombin<EXACT>(N, M);
    } else if constexpr (type == LOWER) {
      logCMN = logCombin<UPPER>(N, M);
    } else {
      logCMN = logCombin<LOWER>(
          N, M);  // could change to exact here to get a tigher bound
    }
    cachedM = M;
    cachedN = N;
    cachedLogMN = logCMN;
  }
  return logCkn + logCN_kM_n - logCMN;
}

template <BOUND_TYPE type = EXACT>
static double hypergeomLogSf(uint64_t k, uint64_t M, uint64_t n, uint64_t N) {
  double sf = -INFINITY;
  double pmf = hypergeomLogPmf<type>(k, M, n, N);
  double eps = pmf - 40;
  while (pmf > eps && k < n) {
    ++k;
    pmf += log2((double)(n - k + 1) / k) -
           log2((double)(M - n - N + k) / (N - k + 1));
    sf = addLogs(sf, pmf);
  }
  return sf;
}

static size_t getBucketCount(size_t Z, size_t N) {
  return GetNextPowerOfTwo((N + Z - 1) / Z);
}

static double failureProbBucketSort(size_t Z, size_t N, size_t numTotalBucket) {
  size_t NTotal = numTotalBucket * Z;
  double failProb = -INFINITY;
  size_t numRealPerBucket = (N - 1) / numTotalBucket + 1;
  for (size_t numBucket = numTotalBucket; numBucket >= 2; numBucket /= 2) {
    size_t n = numRealPerBucket * numBucket;  // num real on first level
    if (Z >= n) {
      break;
    }
    double p = 1.0 / numBucket;
    double bucketFailProb = binomLogSf(Z, n, p);
    double layerFailProb = log2(numTotalBucket) + bucketFailProb;
    failProb = addLogs(failProb, layerFailProb);
  }
  return failProb;
}

static double failureProbBucketSort(size_t Z, size_t N) {
  size_t numTotalBucket = getBucketCount(Z, N);
  return failureProbBucketSort(Z, N, numTotalBucket);
}

/** Z: bucket size (ine element)
N: total number of real elements
M: number of elements that can fit into memory
return IO cost, unit is amortized cost to read/write one element
*/
static double IOCostBucketSort(size_t Z, size_t N, size_t M) {
  size_t numTotalBucket = getBucketCount(Z, N);
  M -= Z * 2;  // since we copy two buckets to temp every time
  size_t numLayerButterfly = ceil(log2(numTotalBucket) / floor(log2(M / Z)));
  if (numLayerButterfly < 2) {
    numLayerButterfly = 2;
  }
  size_t numLayerMerge = 2;
  size_t cost =
      numTotalBucket * Z * (numLayerButterfly - 1) + N * numLayerMerge;
  return cost * 2;
}

/**
return computation cost, unit is number of obli swaps
*/
static double ComputeCostBucketSort(size_t Z, size_t N) {
  size_t numTotalBucket = getBucketCount(Z, N);
  double logZ = log2(Z);
  size_t numLevelButterfly = GetLogBaseTwo(numTotalBucket);
  size_t cost = numTotalBucket * (numLevelButterfly * 0.5 * Z * logZ +
                                  0.25 * Z * logZ * (logZ + 1));
  return cost;
}

template <typename T, class Check>
static T lowerBound(T left, T right, const Check& satisfy, T prec = 1) {
  while (true) {
    T mid = (left + right) / 2;
    if (satisfy(mid)) {
      right = mid;
      if (right - left <= prec) {
        return mid;
      }
    } else {
      left = mid;
      if (right - left <= prec) {
        return right;
      }
    }
  }
}

template <typename T, class Target>
static T convexMinimize(T left, T right, const Target& target, T prec) {
  Assert(prec * 2 < right - left) T curr = (left + right) / 2;
  T step = (right - left) / 4;
  while (true) {
    if (curr + prec > right || curr - prec < left) {
      return curr;
    }
    if (target(curr + prec) < target(curr - prec)) {
      curr += step;
    } else {
      curr -= step;
    }
    step /= 2;
    if (step <= prec) {
      break;
    }
  }
  return curr;
}

static double IOCostOQSort(size_t N, size_t M, size_t layer,
                           double recursive_ratio, double eps) {
  return N * (1 + recursive_ratio) *
         (2 * layer * (1 + eps) + 3 + recursive_ratio);
}

static double ComputeCostOQSort(size_t N, size_t M, size_t layer,
                                double recursive_ratio, double eps) {
  size_t numTotalPart = divRoundUp(N * (1 + eps), M);
  M = divRoundUp(N * (1 + eps), numTotalPart);
  double logM = log2(M);
  double avgRealPerM = M / (1 + eps);
  double logAvgRealPerM = log2(M);
  size_t way = ceil(pow(numTotalPart, 1.0 / layer));
  double logWay = log2(way);
  size_t cost =
      (1 + recursive_ratio) * numTotalPart *
      (layer * 0.5 * M * (logM - logWay / 2) * logWay + 0.5 * M * logM +
       0.25 * avgRealPerM * logAvgRealPerM * (logAvgRealPerM + 1));
  return cost;
}

static double TotalCostOQSort(size_t N, size_t M, size_t layer,
                              double recursive_ratio, double eps) {
  double ioCost = IOCostOQSort(N, M, layer, recursive_ratio, eps);
  double computeCost = ComputeCostOQSort(N, M, layer, recursive_ratio, eps);
  return ioCost * 10 + computeCost;
}

static double ComputeCostDistriOSort(size_t Z, size_t N, size_t layer,
                                     double recursive_ratio) {
  size_t numTotalBucket = getBucketCount(Z, N);
  size_t total = numTotalBucket * Z;
  size_t M = total >> layer;
  size_t avgReal = N >> layer;
  size_t cost = (1 + recursive_ratio) *
                (total * layer * 0.5 * log2(Z) + 0.5 * total * log2(M) +
                 0.25 * N * log2(avgReal) * (log2(avgReal) + 1));
  return cost;
}

static double IOCostDistriOSort(size_t Z, size_t N, size_t M,
                                double recursive_ratio) {
  size_t numTotalBucket = getBucketCount(Z, N);
  M -= Z * 2;  // since we copy two buckets to temp every time
  size_t numLayerButterfly = ceil(log2(numTotalBucket) / floor(log2(M / Z)));
  if (numLayerButterfly < 2) {
    numLayerButterfly = 2;
  }
  size_t cost = 2 * numTotalBucket * Z * (numLayerButterfly - 1) +
                N * (3 + recursive_ratio);
  return cost * (1 + recursive_ratio);
}

static double TotalCostDistriOSort(size_t Z, size_t N, size_t M, size_t layer,
                                   double recursive_ratio) {
  double ioCost = IOCostDistriOSort(Z, N, M, recursive_ratio);
  double computeCost = ComputeCostDistriOSort(Z, N, layer, recursive_ratio);
  return ioCost * 10 + computeCost;
}

struct OQSortParams {
  double alpha;
  double slack_sampling;
  double eps;
  size_t p_max;
  size_t M;
  size_t layer;
};

struct DistriOSortParams {
  double alpha;
  double slack_sampling;
  double eps;
  size_t Z;
  size_t layer;
};

static double logP_t1_data_between_pivots(size_t total, double sampleratio,
                                          double p, size_t t1) {
  // P[max(element between any two pivots) = t1]
  // <= P[there exist two pivots sandwiching t1 elements]
  // <= P[there exists a consecutive chunk of t1-1 elements where
  // total*sampleratio/p elements are sampled]
  // <= 2^prob
  double prob = binomLogPmf<UPPER>(ceil(total * sampleratio / p) - 1, t1 - 1,
                                   sampleratio) +
                log2(total) + 2 * log2(sampleratio);
  return std::min(0.0, prob);
}

static double logP_more_than_t1_data_between_pivots(size_t total,
                                                    double sampleratio,
                                                    double p, size_t t1) {
  // P[max(element between any two pivots) >= t1]
  // <= P[exists a consecutive chunk of t1 elements at most total*sampleratio/p
  // elements are sampled]
  // <= 2^prob
  double prob =
      binomLogCdf<UPPER>(ceil(total * sampleratio / p), t1, sampleratio) +
      log2(total);
  return std::min(0.0, prob);
}

static double logP_Overflow_when_t1_data_between_pivots(size_t total, double p,
                                                        double eps, size_t t1,
                                                        size_t B) {
  // P[overflow | (max(element between any two pivots) = t1)]
  // <= P[overflow | (element between any two pivots = t1)]  (i.e. replacing
  // fillers with real elements only increases failure probability)
  // <= 2^prob
  double prob = hypergeomLogSf<UPPER>(int(B * (1 + eps)), total, t1, B * p) +
                log2(total / B);
  return std::min(0.0, prob);
}

static bool withinFailureProbOQPartition(size_t N, size_t M, size_t layer,
                                         double alpha, double eps,
                                         size_t maxFinalPartSize,
                                         double target = -60) {
  size_t logTotal = ceil(log2(N));
  size_t numTotalPart = divRoundUp(N * (1 + eps), maxFinalPartSize);
  // M = divRoundUp(N * (1+eps), numTotalPart);
  size_t way = ceil(pow(numTotalPart, 1.0 / layer));
  numTotalPart =
      divRoundUp(numTotalPart, pow(way, layer - 1)) * pow(way, layer - 1);

  // printf("numTotalPart=%ld, layer=%ld, way %ld\n", numTotalPart, layer, way);
  double logWay = log2(way);
  size_t p = numTotalPart;
  double q = -INFINITY;
  for (size_t lyr = 0; lyr < layer; ++lyr) {
    Assert(p > 1);
    // num real element per bucket on the first layer
    // N / (number of batches * bucket per batch)
    size_t B = divRoundUp(N, divRoundUp(N * (1 + eps), M) * std::min(way, p));
    // printf("B = %ld\n", B);

    auto satisfy_start = [&](int64_t negt1) {
      return logP_Overflow_when_t1_data_between_pivots(N, p, eps, -negt1, B) <
             target - layer - 3;
    };
    size_t t1_start = -lowerBound(-(int64_t)(N * (1 + eps) / p),
                                  -(int64_t)(N / p), satisfy_start);
    // printf("t1_start = %ld\n", t1_start);

    auto satisfy_end = [&](int64_t t1) {
      return logP_more_than_t1_data_between_pivots(N, alpha, p, t1) <
             target - layer - 3;
    };
    size_t t1_end =
        lowerBound((int64_t)(N / p), (int64_t)(N * (1 + eps) / p), satisfy_end);
    // printf("t1_end = %ld\n", t1_end);
    if (t1_start >= t1_end) {
      t1_start = t1_end = (t1_start + t1_end) / 2;
    }

    // add failure prob on two ends
    double qLayer =
        logP_Overflow_when_t1_data_between_pivots(N, p, eps, t1_start - 1, B);
    qLayer = addLogs(
        qLayer, logP_more_than_t1_data_between_pivots(N, alpha, p, t1_end));

    size_t step = std::max(1UL, (t1_end - t1_start) / 100);
    for (size_t t1 = t1_start; t1 < t1_end; t1 += step) {
      double qi = logP_t1_data_between_pivots(N, alpha, p, t1);
      double qj =
          logP_Overflow_when_t1_data_between_pivots(N, p, eps, t1 + step, B);
      qLayer = addLogs(qLayer, qi + qj + log2(std::min(step, t1_end - t1)));
    }
    // qLayer += log2(log2(std::min(way, p)));
    q = addLogs(q, qLayer);
    if (q > target) {
      return false;
    }
    double ub = addLogs(q, qLayer + log2(layer - lyr - 1));
    if (ub < target) {
      return true;
    }
    p /= way;
  }
  return q < target;
}

static DistriOSortParams bestDistriOSortParams(size_t N, size_t M,
                                               double target = -60) {
  DistriOSortParams params = {};
  double bestCost = INFINITY;
  for (double alpha = 0.02; alpha <= 0.1; alpha += 0.01) {
    size_t k = 128;
    size_t Z;
    auto samplingSatisify = [&](double slack) {
      return binomLogSf<UPPER>(slack * M * alpha, M, alpha) +
                 log2(divRoundUp(N, M)) <
             target - 4;
    };
    static const double slack_sampling_prec = 1e-4;
    double slack_sampling =
        lowerBound(1.0, 1.5, samplingSatisify, slack_sampling_prec);
    for (;; k *= 2) {
      Z = (N - 1) / getBucketCount(k, N) / 2 * 2;
      size_t numBucket = getBucketCount(Z, N);
      size_t numBucketFit = 1UL << GetLogBaseTwo(M / Z);
      double eps = double(Z * numBucket) / N - 1;
      size_t layer = GetLogBaseTwo(numBucket / numBucketFit);
      if (withinFailureProbOQPartition(N, Z * 2, layer, alpha, eps, M,
                                       target)) {
        break;
      }
    }
    // printf("Z init = %ld\n", Z);
    auto satisfy = [=](size_t _Z) {
      _Z += _Z & 1;  // only use even bucket sizes
      size_t numBucket = getBucketCount(_Z, N);
      size_t numBucketFit = 1UL << GetLogBaseTwo(M / _Z);
      double eps = double(_Z * numBucket) / N - 1;
      size_t layer = GetLogBaseTwo(numBucket / numBucketFit);
      return withinFailureProbOQPartition(N, _Z * 2, layer, alpha, eps, M,
                                          target - 0.05);
    };  // -0.05 due to failure prob in sampling
    Z = lowerBound(Z / 2, Z, satisfy);
    Z += Z & 1;
    size_t numBucket = getBucketCount(Z, N);
    size_t numBucketFit = 1UL << GetLogBaseTwo(M / Z);
    size_t layer = GetLogBaseTwo(numBucket / numBucketFit);
    // printf("alpha=%f, slack=%f\n", alpha, slack_sampling);
    double cost = TotalCostDistriOSort(Z, N, M, layer, slack_sampling * alpha);
    double nextCost = 0;
    while (true) {
      size_t nextZ_min = N / (getBucketCount(Z, N) / 2) + 1;
      size_t nextZ_max = 2 * Z;
      size_t nextZ = lowerBound(nextZ_min, nextZ_max, satisfy);
      nextZ += nextZ & 1;
      // printf("nextZ = %ld\n", nextZ);
      if (nextZ >= M / 32 || nextZ >= N / 32) {
        break;
      }
      size_t numBucket = getBucketCount(nextZ, N);
      size_t numBucketFit = 1UL << GetLogBaseTwo(M / nextZ);
      size_t layer = GetLogBaseTwo(numBucket / numBucketFit);
      nextCost =
          TotalCostDistriOSort(nextZ, N, M, layer, slack_sampling * alpha);
      // printf("nextCost=%f\n", nextCost);

      if (nextCost > cost) {
        break;
      }
      Z = nextZ;
      cost = nextCost;
    }
    if (cost < bestCost) {
      bestCost = cost;
      numBucket = getBucketCount(Z, N);
      numBucketFit = 1UL << GetLogBaseTwo(M / Z);
      double eps = double(Z * numBucket) / N - 1;
      params.layer = GetLogBaseTwo(numBucket / numBucketFit);
      params.Z = Z;
      params.eps = eps;
      params.alpha = alpha;
      params.slack_sampling = slack_sampling;
    } else {
      break;
    }
  }
  return params;
}

static OQSortParams bestOQSortParams(size_t N, size_t M, double target = -60) {
  static const double eps_prec = 1e-4;
  static const double slack_sampling_prec = 1e-4;
  static const double alpha_prec = 1e-4;
  static const double logM_prec = 0.25;
  const size_t layer_max = std::max(1UL, (size_t)pow(log(N) / log(M), 2));
  OQSortParams bestParams = {};
  Assert(N > M);
  double bestCost = (double)INFINITY;
  for (size_t layer = 1; layer <= layer_max; ++layer) {
    OQSortParams params = {};
    auto MTarget = [&](int64_t M) {
      auto alphaTarget = [&](double alpha) {
        // printf("alpha = %f\n", alpha);
        auto samplingSatisify = [&](double slack) {
          return binomLogSf<UPPER>(slack * M * alpha, M, alpha) +
                     log2(divRoundUp(N, M)) <
                 target - 4;
        };
        double slack_sampling =
            lowerBound(1.0, 1.5, samplingSatisify, slack_sampling_prec);

        auto epsSatisfy = [&](double _eps) {
          return withinFailureProbOQPartition(N, M, layer, alpha, _eps, M,
                                              target - 0.05);
        };
        double eps = lowerBound(0.0, 1.0, epsSatisfy, eps_prec);
        if (eps >= 1.0 - 2 * eps_prec) {
          eps = lowerBound(0.0, 10.0, epsSatisfy, eps_prec * 10);
          ;
        }
        params.eps = eps;
        params.slack_sampling = slack_sampling;
        // printf("alpha = %f, eps = %f\n", alpha, eps);
        double cost = TotalCostOQSort(N, M, layer, alpha * slack_sampling, eps);
        // printf("cost = %f\n", cost);
        return cost;
      };
      static const double alpha_min = 0.01;
      static const double alpha_max = 0.1;
      double best_alpha = alpha_min;
      double min_cost = INFINITY;
      for (double alpha = alpha_min; alpha <= alpha_max; alpha += 0.01) {
        double cost = alphaTarget(alpha);
        if (cost < min_cost) {
          best_alpha = alpha;
          min_cost = cost;
        } else {
          break;
        }
      }
      params.alpha = best_alpha;
      return alphaTarget(best_alpha);
    };

    size_t bestM = M;

    double cost = MTarget(bestM);

    params.M = bestM;

    if (cost < bestCost) {
      bestCost = cost;
      params.layer = layer;
      bestParams = params;
    }
  }
  bestParams.p_max = ceil(
      pow((1 + bestParams.eps) * N / M,
          1.0 / (bestParams.layer - 0.2)));  // -0.2 to avoid rounding issue
  return bestParams;
}

class ButterflyWaySolver {
 private:
  // options of mergesplit ways
  std::vector<uint64_t> options;
  // amortized cost on each bucket in option-way mergesplit
  std::vector<double> unitCosts;
  // amortized cost to read&write a bucket from/to ext mem
  double ioCost;
  // maximum number of ways that can be done in internal memory
  uint64_t maxWayInternal;

  double solverHelper(std::vector<uint64_t>& optimalWays, size_t optionCount,
                      uint64_t targetProduct, uint64_t remainWayInternal) {
    optimalWays.clear();
    if (targetProduct <= 1) {
      return 0;
    }
    if (!optionCount) {
      return INFINITY;
    }
    std::vector<uint64_t> childrenOptimalWays;
    double minCost = INFINITY;
    for (size_t optionIdx = 0; optionIdx < optionCount; ++optionIdx) {
      uint64_t option = options[optionIdx];
      if (option <= remainWayInternal) {
        // for efficiency only allow to use options <= current option
        double costChildren = solverHelper(childrenOptimalWays, optionIdx + 1,
                                           divRoundUp(targetProduct, option),
                                           remainWayInternal / option);
        uint64_t product = option;
        for (uint64_t way : childrenOptimalWays) {
          product *= way;
        }
        double cost = costChildren * option + unitCosts[optionIdx] * product;
        if (cost < minCost) {
          minCost = cost;
          optimalWays.resize(childrenOptimalWays.size() + 1);
          std::copy(childrenOptimalWays.begin(), childrenOptimalWays.end(),
                    optimalWays.begin());
          optimalWays.back() = option;
        }
        if (option == targetProduct || option == maxWayInternal) {
          return minCost;  // no need to try next option
        }
      }
    }
    if (options[optionCount - 1] > remainWayInternal) {
      // start new layer, reset all options and remainWayInternal
      double costWithoutIO = solverHelper(childrenOptimalWays, options.size(),
                                          targetProduct, maxWayInternal);
      uint64_t product = 1;
      for (uint64_t way : childrenOptimalWays) {
        product *= way;
      }
      double cost = costWithoutIO + product * ioCost;
      if (cost < minCost) {
        minCost = cost;
        optimalWays.resize(childrenOptimalWays.size());
        std::copy(childrenOptimalWays.begin(), childrenOptimalWays.end(),
                  optimalWays.begin());
      }
    }
    return minCost;
  }

 public:
  ButterflyWaySolver(const std::vector<uint64_t>& options,
                     const std::vector<double>& unitCosts, double ioCost,
                     uint64_t maxWayInternal)
      : options(options),
        unitCosts(unitCosts),
        ioCost(ioCost),
        maxWayInternal(maxWayInternal) {
    Assert(unitCosts.size() == options.size());
    std::sort(this->options.begin(), this->options.end());
  }

  double solve(std::vector<std::vector<uint64_t>>& output,
               uint64_t targetProduct) {
    std::vector<uint64_t> optimalWays;
    double cost = solverHelper(optimalWays, options.size(), targetProduct,
                               maxWayInternal);
    output.clear();
    uint64_t accWays = 1;
    auto begin = optimalWays.begin();
    auto it = begin;
    for (; it != optimalWays.end(); ++it) {
      if (accWays * (*it) > maxWayInternal) {
        if (begin != it) {
          output.emplace_back(begin, it);
        }
        begin = it;
        accWays = *it;
      } else {
        accWays *= (*it);
      }
    }
    output.emplace_back(begin, it);
    return cost;
  }
};

struct KWayButterflyParams {
  std::vector<std::vector<uint64_t>> ways;
  size_t Z;
  size_t totalBucket;
};

static KWayButterflyParams bestKWayButterflyParams(size_t N,
                                                   size_t M = 0x7000000 / 136,
                                                   double target = -60) {
  double minCost = INFINITY;
  KWayButterflyParams optimalParams;
  const std::vector<uint64_t> choices = {2, 3, 4, 5, 6, 7, 8};
  const std::vector<double> bitonicCost = {
      0.183596, 0.444128, 0.985086, 2.275917, 5.363025, 12.508203, 29.262215};
  const std::vector<std::vector<double>> costs = {
      {0.067040, 0.092965, 0.100130, 0.098910, 0.092490, 0.096466, 0.098223},
      {0.108723, 0.161996, 0.174614, 0.186288, 0.195849, 0.207293, 0.212523},
      {0.234463, 0.346844, 0.379626, 0.408319, 0.433684, 0.459960, 0.471675},
      {0.509743, 0.784653, 0.856935, 0.910044, 0.948897, 0.994666, 1.013502},
      {1.139603, 1.725107, 1.854331, 2.058958, 2.064610, 2.118721, 2.149046},
      {2.839750, 3.790226, 3.941216, 4.139045, 4.287033, 4.438938, 4.513650},
      {6.552275, 7.858482, 8.353425, 8.662583, 8.997599, 9.329178, 9.564495}};

#ifdef DISK_IO
  static const double ioCostPerElement = 0.00096857;
#else
  static const double ioCostPerElement = 0.00061467;
#endif
  for (int i = 0; i < 7; ++i) {
    size_t Z = 256UL << i;
    if (Z >= N / 2 || Z >= M / 16) {
      if (i == 0) {
        printf("Input size / memory size too small\n");
        abort();
      }
      break;
    }
    auto satisfy = [Z = Z, N = N, target = target](size_t bucketCount) {
      return failureProbBucketSort(Z, N, bucketCount) < target;
    };
    size_t minBucketCount = lowerBound(N / Z + 1, 3 * N / Z, satisfy);
    size_t maxWayInternal = (M - Z * 8) / Z;
    ButterflyWaySolver solver(choices, costs[i], ioCostPerElement * Z,
                              maxWayInternal);
    KWayButterflyParams params;
    double cost = solver.solve(params.ways, minBucketCount);
    uint64_t actualBucketCount = 1;
    for (auto& vec : params.ways) {
      for (auto way : vec) {
        actualBucketCount *= way;
      }
    }
    cost += bitonicCost[i] * actualBucketCount;
    if (cost < minCost) {
      minCost = cost;
      optimalParams = params;
      optimalParams.Z = Z;
      optimalParams.totalBucket = actualBucketCount;
    }
  }
  return optimalParams;
}
}  // namespace EM::Algorithm