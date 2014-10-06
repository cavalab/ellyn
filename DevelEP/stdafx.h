// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
// DEVELEP

#pragma once
#include <stdio.h>
#include <numeric>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <boost/random.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <boost/regex.hpp>
#include "time.h"
#include <math.h>
#if defined(_WIN32)
	#include "targetver.h"
	#include <tchar.h>
	#include <boost/shared_ptr.hpp>
	#include <direct.h>
#else
	#include <memory>
	#include <sys/stat.h>
	#include <unistd.h>
#endif


//#include <stdlib.h>
#include <cstdlib>
//#include "ppl.h"
#include <omp.h>

#include <exception>
/*#include "pop.h"
#include "pareto.h"
#include "p_archive.h"
#include "data.h"
#include "general_fns.h"*/
// TODO: reference additional headers your program requires here
