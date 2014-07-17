#pragma once
#ifndef FITNESS_H
#define FITNESS_H

//#include "pop.h"
//#include "params.h"
//#include "data.h"
//
////#include "evaluator.h"
//
//#include "EvalEqnStr.h"
//#include <unordered_map>
//#include "Line2Eqn.h"
#include "state.h"
#include "FitnessEstimator.h"

class CException
{
public:
	char* message;
	CException( char* m ) { message = m; };
	void Report(){cout << "error calculating fitness\n";};

};

//void getEqnForm(std::string& eqn,std::string& eqn_form);
float getCorr(vector<float>& output,vector<float>& target,float meanout,float meantarget,int off,float& target_std);
//int getComplexity(string& eqn);
float std_dev(vector<float>& target,float& meantarget);
void Fitness(vector<ind>& pop,params& p,data& d,state& s,FitnessEstimator& FE);


#endif