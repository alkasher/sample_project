#pragma once 
#include <cmath>
#include <unordered_set>
#include <cstdint> // needed for uint*_t types
#include <iostream>
#include <time.h>
#include  <random>
#include <algorithm>
#include  <iterator>
#include "tms-nets-3.0.1/include/tms-nets.hpp"
using namespace tms;

namespace std {
	template <>
	struct hash<std::vector<int>> {
		size_t operator()(const vector<int>& v) const {
			std::hash<int> hasher;
			size_t seed = 0;
			for (int i : v) {
				seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};
}


void findFactorImpl(int N, int K, std::vector<std::vector<int>>& temp,
	unsigned f_min) {
	if (K <= 1) {
		temp[temp.size() - 1].push_back(N);
		temp.push_back({});

	}
	else {
		for (int f = f_min; f * K <= N; ++f) {
			temp[temp.size() - 1].push_back(f);
			findFactorImpl(N - f, K - 1, temp, f);
		}
	}
}

std::vector<std::vector<int>> findFactor(unsigned N, unsigned K) {
	std::vector<std::vector<int>> temp;
	temp.push_back({});
	findFactorImpl(N, K, temp, 0);
	return temp;
}


void permute(std::vector<int>& nums, std::vector<int>& currPerm, std::vector<bool>& used, std::unordered_set<std::vector<int>>& result) {
	if (currPerm.size() == nums.size()) {
		result.insert(currPerm);
		return;
	}

	for (int i = 0; i < nums.size(); ++i) {
		if (used[i]) {
			continue;
		}
		if (i > 0 && nums[i] == nums[i - 1] && !used[i - 1]) {
			continue; // Skip duplicates
		}
		used[i] = true;
		currPerm.push_back(nums[i]);
		permute(nums, currPerm, used, result);
		currPerm.pop_back();
		used[i] = false;
	}
}

void test_net(std::vector<std::pair<Point, Point>> result, std::vector<Point> vect, int b, int t) {
	int s = result.size();
	std::vector<int> count(s);
	for (int n = 0; n < vect.size(); ++n) {
		for (int j = 0; j < result.size(); ++j) {
			bool flag = true;
			for (int i = 0; i < vect[n].size(); ++i) {
				if (vect[n][i] < result[j].first[i] or vect[n][i] >= result[j].second[i])
					flag = false;
			}
			if (flag == true) { count[j]++; }
		}
	}
	int k = count[0];
	for (int i = 0; i < count.size(); ++i) {
		if (k != count[i]) {
			std::cout << "\nt can't be used\n";
		}
	}
	std::cout << "\nt-m-s net\n";
	return;
}



void generateIntervals(int s, int d, int b, std::vector<Point> vect, int t) {
	//const std::unordered_set<std::vector<int>>& set;
	std::vector<std::vector<int>> temp = findFactor(d, s);
	temp.pop_back();
	for (auto nums : temp) {
		std::vector<int> currPerm;
		std::vector<bool> used(nums.size(), false);
		std::unordered_set<std::vector<int>> res;
		std::sort(nums.begin(), nums.end());
		permute(nums, currPerm, used, res);

		for (auto d : res) {
			// ¬ектор дл€ хранени€ текущих значений a[i] дл€ каждого измерени€
			std::vector<int> a(s, 0);

			// ћаксимальное значение a[i] дл€ каждого измерени€
			std::vector<int> max_a(s);
			for (int i = 0; i < s; ++i) {
				max_a[i] = std::pow(b, d[i]) - 1;
			}
			std::vector<std::pair<Point, Point>> result;

			while (true) {
				// ¬ывод текущего интервала
				std::pair<Point, Point> tmp;
				//std::cout << "[";
				for (int i = 0; i < s; ++i) {
					Real start = static_cast<Real>(a[i]) / std::pow(b, d[i]);
					Real end = static_cast<Real>(a[i] + 1) / std::pow(b, d[i]);
					tmp.first.push_back(start);
					tmp.second.push_back(end);
					//std::cout << "[" << start << ", " << end << ")";
					//if (i < s - 1) std::cout << " x ";
				}
				result.push_back(tmp);
				//	std::cout << "]\n";

					// ќбновление значений a[i]
				int k = s - 1;
				while (k >= 0 && a[k] == max_a[k]) {
					a[k] = 0;
					--k;
				}
				if (k < 0) { break; }
				++a[k];
			}
			test_net(result, vect, b, t);

		}
	}

}