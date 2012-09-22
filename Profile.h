// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#ifndef PROFILE_H_
#define PROFILE_H_

//#define ENABLE_PROFILE

#ifdef ENABLE_PROFILE
#include <cstdio>
#include <map>
#include <windows.h>
static struct Profile{
	std::map<const char *, std::pair<int, LONGLONG> > pset; 
	Profile(){
	}
	~Profile(){
		FILE *fp = std::fopen("d:/profile.txt", "w");
		LARGE_INTEGER freq;
		::QueryPerformanceFrequency(&freq);
		double dfreq = static_cast<double>(freq.QuadPart);
		std::fprintf(fp, "fname\tcount\ttotal time(msec)\ttime per call(msec)\n======\n");
		for (std::map<const char *, std::pair<int, LONGLONG> >::iterator iter = pset.begin(); iter != pset.end(); ++iter){
			double ttime = static_cast<double>(iter->second.second) * 1000.0 / dfreq;
			std::fprintf(fp, "%s\t%d\t%f\t%f\n", iter->first, iter->second.first, ttime, ttime / static_cast<double>(iter->second.first));
		}
		fclose(fp);
	}
	struct Item{
		const char *name;
		LARGE_INTEGER c1, c2;
		Item(const char *name_) : name(name_){
			::QueryPerformanceCounter(&c1);
		}
		~Item(){
			::QueryPerformanceCounter(&c2);
			std::pair<int, LONGLONG> &p = profile.pset[name];
			++p.first;
			p.second += c2.QuadPart - c1.QuadPart;
		}
	};
}profile;
#define PROF(X) Profile::Item item(X)
#else
#define PROF(X)
#endif


#endif // PROFILE_H_
