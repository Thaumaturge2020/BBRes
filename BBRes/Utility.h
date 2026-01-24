#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <tuple>
#include <queue>
#include <set>

#define NDEBUG
#define NEW_BRANCH
#define FLOW_PRUNE
#define COLOR_PRUNE
#define COLOR_BASED_PRUNE
// #define MDC_BRANCH
// #define WODC_PIVOT
#define COLOR_AHEAD_PRUNE

#include <cassert>

using ui = unsigned int; // vertex type
using ept = unsigned int; // edge pointer type; unsigned int can be used to process upto two billion undirected edges

#define pb push_back
#define mp std::make_pair

#define mmax(a,b) ((a)>(b)?(a):(b))
#define mmin(a,b) ((a)<(b)?(a):(b))

class Utility {
public:
	static FILE *open_file(const char *file_name, const char *mode) {
		FILE *f = fopen(file_name, mode);
		if(f == nullptr) {
			printf("Can not open file: %s\n", file_name);
			exit(1);
		}

		return f;
	}

	static std::string integer_to_string(long long number) {
		std::vector<long long> sequence;
		sequence.pb(number);

		char buf[30];
		std::string res;
		for(ui i = sequence.size();i > 0;i --) {
			if(i == sequence.size()) sprintf(buf, "%llu", sequence[i-1]);
			else sprintf(buf, ",%03llu", sequence[i-1]);
			res += std::string(buf);
		}
		return res;
	}
};

#endif
