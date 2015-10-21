#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "data.h"
#include "rnd.h"
#include "state.h"
#include "Line2Eqn.h"
#include "EvalEqnStr.h"
#include <unordered_map>
#include <bitset>
#include <math.h>  
#include "matrix.h"
#include <Eigen/Dense>
#include "general_fns.h"
using Eigen::MatrixXf;
using Eigen::ArrayXf;
//#include "runEllenGP.h"

#if defined(_WIN32)
	#include <regex>
#else
	#include <boost/regex.hpp>
#endif
#include "FitnessEstimator.h"
//#include "Fitness.h"

class CException
{
public:
	char* message;
	CException( char* m ) { message = m; };
	void Report(){cout << "error calculating fitness\n";};

};

//void getEqnForm(std::string& eqn,std::string& eqn_form);
//float getCorr(vector<float>& output,vector<float>& target,float meanout,float meantarget,int off);
//int getComplexity(string& eqn);

//void getEqnForm(std::string& eqn,std::string& eqn_form)
//{
////replace numbers with the letter c
//#if defined(_WIN32)
//	std::regex e ("(([0-9]+)(\\.)([0-9]+))|([0-9]+)");
//	std::basic_string<char> tmp = "c";
//	std::regex_replace (std::back_inserter(eqn_form), eqn.begin(), eqn.end(), e,tmp);
//#else
//	boost::regex e ("(([0-9]+)(\\.)([0-9]+))|([0-9]+)");
//	eqn_form = boost::regex_replace(eqn,e,"c");
//#endif 
//	//(\d+\.\d+)|(\d+)
//	//eqn_form = eqn;
//	//eqn_form=std::tr1::regex_replace(eqn,e,tmp.c_str(),std::tr1::regex_constants::match_default);
//	//std::string result;
//
//
//	//std::regex_replace(std::back_inserter(eqn_form),eqn.begin(),eqn.end(),e,"c",std::regex_constants::match_default);
//	//std::regex_replace(std::back_inserter(eqn_form),eqn.begin(),eqn.end(),e,"c",std::tr1::regex_constants::match_default
//    //std::cout << result;
//	//std::cout << eqn << "\t" << eqn_form <<"\n";
//}

const int bitlen( int n )  
{  //return bitstring length needed to represent n
    // log(n)/log(2) is log2.  
    return ceil(log( double(n) ) / log( float(2.0) ));  
}

void CalcOutput(ind& me,params& p,vector<vector<float>>& vals,vector<float>& dattovar,vector<float>& target,state& s);
void CalcClassOutput(ind& me,params& p,vector<vector<float>>& vals,vector<float>& dattovar,vector<float>& target,state& s);
void Calc_M3GP_Output(ind& me,params& p,vector<vector<float>>& vals,vector<float>& dattovar,vector<float>& target,state& s);
float getCorr(vector<float>& output,vector<float>& target,float meanout,float meantarget,int off,float& target_std)
{
	float v1,v2;
	float var_target=0;
	float var_ind = 0;
	float q=0;
	float ndata = float(output.size());
	float corr;
	
	//calculate correlation coefficient
	for (unsigned int c = 0; c<output.size(); ++c)
	{

		v1 = target.at(c+off)-meantarget;
		v2 = output.at(c)-meanout;

		q += v1*v2;
		var_target+=pow(v1,2);
		var_ind+=pow(v2,2);
	}
	q = q/(ndata-1); //unbiased esimator
	var_target=var_target/(ndata-1); //unbiased esimator
	var_ind =var_ind/(ndata-1); //unbiased esimator
	if(abs(var_target)<0.0000001 || abs(var_ind)<0.0000001)
		corr = 0;
	else{
		corr = pow(q,2)/(var_target*var_ind);
		//corr = pow(q,2)/(sqrt(var_target)*var_ind); //for normalizing the error by the standard deviation of the target
	}
	target_std = sqrt(var_target);
	return corr;
}
float VAF_loud(vector<float>& output,vector<float>& target,float meantarget,int off,state& s)
{
	float v1,v2;
	float var_target=0;
	float var_diff = 0;
	float q=0;
	float ndata = float(output.size());
	float vaf=0;
//	float diff;
	float diffmean=0;
	s.out << "output.size() = " << output.size() << "\n";
	s.out << "target.size() = " << target.size() << "\n";

	for (unsigned int c = 0; c<output.size(); ++c)
		diffmean += target.at(c+off)-output.at(c);
	diffmean = diffmean/output.size();
	s.out << "diffmean = " << diffmean << "\n";
	//calculate correlation coefficient
	s.out << "var_diff: ";
	for (unsigned int c = 0; c<output.size(); ++c)
	{		
		v1 = target.at(c+off)-meantarget;
		v2 = (target.at(c+off)-output.at(c))-diffmean;
		var_target+=pow(v1,2);
		var_diff+=pow(v2,2);
		s.out << var_diff << "\t";
	}
	//q = q/(ndata-1); //unbiased estimator
	var_target=var_target/(ndata-1); //unbiased estimator
	var_diff =var_diff/(ndata-1); //unbiased estimator
	s.out << "\nvar_diff = " << var_diff << "\n";
	s.out << "var_target = " << var_target << "\n";
	s.out << "var_diff / var_target = " << var_diff/var_target << "\n";
	if(var_target<0.0000001)
		return 0;
	else{
		float tmp = (1-var_diff/var_target)*100;
		s.out << "tmp = " << tmp << "\n";
		vaf = std::max(float(0),tmp);
		s.out << "vaf = " << vaf << "\n";
	}
	return vaf;
}
float VAF(vector<float>& output,vector<float>& target,float meantarget,int off)
{
	if (*min_element(output.begin(),output.end())==*max_element(output.begin(),output.end()))
		return 0;

	float v1,v2;
	float var_target=0;
	float var_diff = 0;
	float q=0;
	float ndata = float(output.size());
	float vaf=0;
//	float diff;
	float diffmean=0;

	for (unsigned int c = 0; c<output.size(); ++c)
		diffmean += target.at(c+off)-output.at(c);
	diffmean = diffmean/output.size();
	//calculate correlation coefficient
	for (unsigned int c = 0; c<output.size(); ++c)
	{		
		v1 = target.at(c+off)-meantarget;
		v2 = (target.at(c+off)-output.at(c))-diffmean;
		var_target+=pow(v1,2);
		var_diff+=pow(v2,2);
	}
	//q = q/(ndata-1); //unbiased estimator
	var_target=var_target/(ndata-1); //unbiased estimator
	var_diff =var_diff/(ndata-1); //unbiased estimator
	if(var_target<0.0000001)
		return 0;
	else{
		float tmp = (1-var_diff/var_target)*100;
		vaf = std::max(float(0),tmp);
	}
	return vaf;
}
float std_dev(vector<float>& target,float& meantarget)
{
	float s=0;
	//calculate correlation coefficient
	for (unsigned int c = 0; c<target.size(); ++c)
	{
		s+=pow(abs(target[c]-meantarget),2);
	}
	s = s/(target.size()-1);
	return sqrt(s);
}

int recComplexity(vector<node>& line,int i, int &j)
{
	while (!line[i].on || line[i].intron) 
		--i;
	j = i;
	int added = 0;
	if (line[i].arity() == 0){
		//cout << line[i].c ;
		return line[i].c;
	}
	else{
		int comp = 0; 
		while (added<line[i].arity()){
			++added; // can't guarantee that updated i-added will traverse the last value added. need better way to do this...				
			comp += line[i].c * recComplexity(line,j-1,j);
		}
		comp += line[i].c;
		return comp;
	}
}
int recComplexity2(vector<node>& line, int i, int& j, char type)
{
	
	while (!line[i].on || line[i].intron || line[i].return_type != type){
		if (line[i].return_type != type && !line[i].intron && line[i].on){
			if (type == 'f' && line[i].arity_float > 0) {
				i -= line[i].arity_float;
			}
			else if (type == 'b' && line[i].arity_bool > 0) {
				i -= line[i].arity_bool;
			}
			else
				--i;
		}
		else
			--i;
		assert(i > -1);
	}
	int k = i;

	j = i;
	
	if (line[i].arity() == 0) 
		return line[i].c;
	else {
		
		//int added = 0;
		int comp = 0;
		//int k = i;		
		
		for (int f = 0; f < line[i].arity_float; ++f) {		
			comp += line[i].c * recComplexity2(line, j - 1, j,'f');
		}
		
		if (k != i)
			cout << "wth\n";
		for (int b = 0; b < line[i].arity_bool; ++b) {				
			comp += line[i].c * recComplexity2(line, k - 1, k, 'b');
		}
		comp += line[i].c;
		return comp;
	}
}
void eval_complexity(const node& n, vector<int>& c_float, vector<int>& c_bool)
{
	//float n1, n2;
	//bool b1, b2;
	if (c_float.size() >= n.arity_float && c_bool.size() >= n.arity_bool) {
		int c_f = 0;
		int c_b = 0;

		if (n.type == 'n' || n.type == 'v')
			c_float.push_back(n.c);
		else {
			for (size_t i = 0; i < n.arity_float; ++i) {
				c_f += c_float.back();
				c_float.pop_back();
			}
			for (size_t i = 0; i < n.arity_bool; ++i) {
				c_b += c_bool.back();
				c_bool.pop_back();
			}
			if (n.return_type == 'f')
				c_float.push_back(n.c*(c_f + c_b + 1));
			else
				c_bool.push_back(n.c*(c_f + c_b + 1));
		}
	}

}
int getComplexity(ind& me, params& p)
{	
	int complexity = 0;

	vector<int> c_float;
	vector<int> c_bool;

	for (size_t i = 0; i < me.line.size(); ++i) 
		if (me.line[i].on) eval_complexity(me.line[i], c_float, c_bool);

	if (c_float.empty())
		complexity = p.max_fit;
	else if (p.classification && p.class_m3gp) 
		complexity = accumulate(c_float.begin(), c_float.end(), 0);
	else
		complexity = c_float.back();

	return complexity;
}

//int getComplexity(string& eqn)
//{
//	int complexity=0;
//	char c;
//	for(int m=0;m<eqn.size();++m){
//		c=eqn[m];
//		
//		if(c=='/')
//			complexity=complexity+2;
//		else if (c=='s'){
//			
//			if(m+2<eqn.size()){
//				if ( eqn[m+1]=='i' && eqn[m+2] == 'n'){
//					complexity=complexity+3;
//					m=m+2;
//				}
//			}
//			if (m+3<eqn.size()){ //sqrt
//				if ( eqn[m+1]=='q' && eqn[m+2] == 'r' && eqn[m+3]=='t'){
//					complexity=complexity+3;
//					m=m+2;
//				}
//			}
//			
//		}
//		else if (c=='c'){
//			if(m+2<eqn.size()){
//				if ( eqn[m+1]=='o' && eqn[m+2] == 's'){
//					complexity=complexity+3;
//					m=m+2;
//				}
//			}
//		}
//		else if (c=='e'){
//			if(m+2<eqn.size()){
//				if ( eqn[m+1]=='x' && eqn[m+2] == 'p'){
//					complexity=complexity+4;
//					m=m+2;
//				}
//			}
//		}
//		else if (c=='l'){
//			if(m+2<eqn.size()){
//				if ( eqn[m+1]=='o' && eqn[m+2] == 'g'){
//					complexity=complexity+4;
//					m=m+2;
//				}
//			}
//		}
//		else if (isalpha(c) && (m+1)<eqn.size()){
//			bool pass=true;
//			while ((m+1)<eqn.size() && pass){
//				if (isalpha(eqn[m+1]) || isdigit(eqn[m+1])) ++m; 
//				else pass=0;
//			}
//			++complexity;
//		}
//		else
//			++complexity;
//	}
//
//	return complexity;
//}
int getEffSize(vector<node>& line)
{
	int eff_size=0;
	for(int m=0;m<line.size();++m){
		if(line.at(m).on)
			++eff_size;
	}
	return eff_size;
}
void eval(node& n,vector<float>& stack_float,vector<bool>& stack_bool)
{
	float n1,n2;
	bool b1,b2; 
	if (stack_float.size()>=n.arity_float && stack_bool.size()>=n.arity_bool){
		n.intron=false;
		switch(n.type) 
		{
		case 'n':
			stack_float.push_back(n.value);
			break;
		case 'v':
			stack_float.push_back(*n.valpt);
			break;
		case '+':
			n1 = stack_float.back(); stack_float.pop_back();
			n2 = stack_float.back(); stack_float.pop_back();
			stack_float.push_back(n2+n1);
			break;
		case '-':
			n1 = stack_float.back(); stack_float.pop_back();
			n2 = stack_float.back(); stack_float.pop_back();
			stack_float.push_back(n2-n1);
			break;
		case '*':
			n1 = stack_float.back(); stack_float.pop_back();
			n2 = stack_float.back(); stack_float.pop_back();
			stack_float.push_back(n2*n1);
			break;
		case '/':
			n1 = stack_float.back(); stack_float.pop_back();
			n2 = stack_float.back(); stack_float.pop_back();
			if(abs(n1)<0.000001)
				stack_float.push_back(1);
			else
				stack_float.push_back(n2/n1);
			break;
		case 's':
			n1 = stack_float.back(); stack_float.pop_back();
			stack_float.push_back(sin(n1));
			break;
		case 'c':
			n1 = stack_float.back(); stack_float.pop_back();
			stack_float.push_back(cos(n1));
			break;
		case 'e':
			n1 = stack_float.back(); stack_float.pop_back();
			stack_float.push_back(exp(n1));
			break;
		case 'l':
			n1 = stack_float.back(); stack_float.pop_back();
			if (abs(n1)<0.000001)
				stack_float.push_back(0);
			else 
				// safe log of absolute value of n1
				stack_float.push_back(log(abs(n1)));
				// unsafe log of real value
				//stack_float.push_back(log(n1));
			break;
		case 'q':
			n1 = stack_float.back(); stack_float.pop_back();
			// safe sqrt of absolute value of n1
			stack_float.push_back(sqrt(abs(n1)));
			break;
		case '=': // equals			
			n1 = stack_float.back(); stack_float.pop_back();
			n2 = stack_float.back(); stack_float.pop_back();
			stack_bool.push_back(n1==n2);		
			break;
		case '!': // does not equal			
			b1 = stack_bool.back(); stack_bool.pop_back();
			stack_bool.push_back(!b1);
			break;
		case '>': //greater than
			n1 = stack_float.back(); stack_float.pop_back();
			n2 = stack_float.back(); stack_float.pop_back();
			stack_bool.push_back(n2 > n1);
			break;
		case '<': //less than			
			n1 = stack_float.back(); stack_float.pop_back();
			n2 = stack_float.back(); stack_float.pop_back();
			stack_bool.push_back(n2 < n1);		
			break;
		case '}': //greater than or equal
			n1 = stack_float.back(); stack_float.pop_back();
			n2 = stack_float.back(); stack_float.pop_back();
			stack_bool.push_back(n2 >= n1);
			break;
		case '{': //less than or equal
			n1 = stack_float.back(); stack_float.pop_back();
			n2 = stack_float.back(); stack_float.pop_back();
			stack_bool.push_back(n2 <= n1);
			break;
		case 'i': // if (arity 2). if stack_bool true, leave top element of floating stack. otherwise, pop it.
			b1 = stack_bool.back(); stack_bool.pop_back();
			if (!b1)
				stack_float.pop_back();
			break;
		case 't': // if-then-else (arity 3). if stack_bool true, leave 2nd element of floating stack and pop first. otherwise, pop 2nd and leave first.
			b1 = stack_bool.back(); stack_bool.pop_back();
			//n1 = stack_float.back(); stack_float.pop_back();
			//n2 = stack_float.back(); stack_float.pop_back();
			if (b1)
				stack_float.pop_back();
			else{
				swap(stack_float[stack_float.size()-2],stack_float.back());
				stack_float.pop_back();
			}
			break;
		case '&':
			b1 = stack_bool.back(); stack_bool.pop_back();
			b2 = stack_bool.back(); stack_bool.pop_back();
			stack_bool.push_back(b1 && b2);
			break;
		case '|':
			b1 = stack_bool.back(); stack_bool.pop_back();
			b2 = stack_bool.back(); stack_bool.pop_back();
			stack_bool.push_back(b1 || b2);
			break;
		default:
			cout << "eval error\n";
			break;
		}
	}
	else
		n.intron= n.intron && true; // only set it to intron if it isn't used in any of the execution

}

void FitnessEstimate(vector<ind>& pop,params& p,Data& d,state& s,FitnessEstimator& FE);
void StandardFitness(ind& me,params& p,Data& d,state& s,FitnessEstimator& FE, unordered_map<string,float*>& datatable, vector<float>& dattovar);
void LexicaseFitness(ind& me,params& p,Data& d,state& s,FitnessEstimator& FE);

void Fitness(vector<ind>& pop,params& p,Data& d,state& s,FitnessEstimator& FE)
{
	////
	//set up data table for conversion of symbolic variables
	unordered_map <string,float*> datatable;
	vector<float> dattovar(d.label.size());

	for (unsigned int i=0;i<d.label.size(); ++i)
			datatable.insert(pair<string,float*>(d.label[i],&dattovar[i]));

	//#pragma omp parallel for private(e)
	for(int count = 0; count<pop.size(); ++count)
	{
		if (p.print_protected_operators){
			pop.at(count).eqn_matlab = Line2Eqn(pop.at(count).line,pop.at(count).eqn_form,p,true);
			pop.at(count).eqn = Line2Eqn(pop.at(count).line,pop.at(count).eqn_form,p,false);
		}
		else
			pop.at(count).eqn = Line2Eqn(pop.at(count).line,pop.at(count).eqn_form,p,true);
		//getEqnForm(pop.at(count).eqn,pop.at(count).eqn_form);
		
		pop.at(count).eff_size = getEffSize(pop.at(count).line);

		//if(p.sel!=3){
		StandardFitness(pop.at(count),p,d,s,FE,datatable,dattovar);		

		if (p.classification && p.class_m3gp)
			pop[count].dim = pop[count].M.cols();
		//} // if p.sel!=3
		//else //LEXICASE FITNESS 
		//{
		//	LexicaseFitness(pop.at(count),p,d,s,FE);
		//}//LEXICASE FITNESS
		if (pop[count].eqn.compare("unwriteable")==0)
			pop[count].complexity=p.max_fit;
		else
			pop.at(count).complexity= getComplexity(pop.at(count),p);
		/*if (p.estimate_generality && pop.at(count).genty != abs(pop[count].fitness-pop[count].fitness_v)/pop[count].fitness && pop.at(count).genty != p.max_fit) 
			std::cerr << "genty error, line 300 Fitness.cpp\n";*/
	}//for(int count = 0; count<pop.size(); ++count)
	s.numevals[omp_get_thread_num()]=s.numevals[omp_get_thread_num()]+pop.size();
	//cout << "\nFitness Time: ";
}
void StandardFitness(ind& me,params& p,Data& d,state& s,FitnessEstimator& FE,unordered_map<string,float*>& datatable, vector<float>& dattovar)
{
    
	me.abserror = 0;
	me.abserror_v = 0;
	
	// set data table and pointers to data in program nodes
	for(int m=0;m<me.line.size();++m){
		if(me.line.at(m).type=='v')
			{// set pointer to dattovar 
				//float* set = datatable.at(static_pointer_cast<n_sym>(me.line.at(m)).varname);
				float* set = datatable.at(me.line.at(m).varname);
				if(set==NULL)
					cout<<"hmm";
				/*static_pointer_cast<n_sym>(me.line.at(m)).setpt(set);*/
				me.line.at(m).setpt(set);
				/*if (static_pointer_cast<n_sym>(me.line.at(m)).valpt==NULL)
					cout<<"wth";*/
			}
	}
	//cout << "Equation" << count << ": f=" << me.eqn << "\n";
	bool pass=true;
	if(!me.eqn.compare("unwriteable")==0){
							
// calculate error 
		if(!p.EstimateFitness){
			if(p.classification && p.class_m3gp)  Calc_M3GP_Output(me,p,d.vals,dattovar,d.target,s);
			else if (p.classification)  CalcClassOutput(me,p,d.vals,dattovar,d.target,s);
			else CalcOutput(me,p,d.vals,dattovar,d.target,s);
		} // if not estimate fitness 
		else{ 
// use fitness estimator subset of d.vals ===========================================================	
			vector<vector<float>> FEvals;
			vector<float> FEtarget;
			setFEvals(FEvals,FEtarget,FE,d);
			if(p.classification && p.class_m3gp)  Calc_M3GP_Output(me,p,FEvals,dattovar,FEtarget,s);
			else if(p.classification) CalcClassOutput(me,p,FEvals,dattovar,FEtarget,s);
			else CalcOutput(me,p,FEvals,dattovar,FEtarget,s);
		} // if estimate fitness
		if (p.estimate_generality || p.PS_sel==2){
				if (me.fitness == p.max_fit || me.fitness_v== p.max_fit || boost::math::isnan(me.fitness_v) || boost::math::isinf(me.fitness_v))
					me.genty = p.max_fit;
				else{
					if (p.G_sel==1) // MAE
						me.genty = abs(me.abserror-me.abserror_v)/me.abserror;
					else if (p.G_sel==2) // R2
						me.genty = abs(me.corr-me.corr_v)/me.corr;
					else if (p.G_sel==3) // MAE R2 combo
						me.genty = abs(me.abserror/me.corr-me.abserror_v/me.corr_v)/(me.abserror/me.corr);
					else if (p.G_sel==3) // VAF
						me.genty = abs(me.VAF-me.VAF_v)/me.VAF;
				}
				if ( boost::math::isnan(me.genty) || boost::math::isinf(me.genty))
					me.genty=p.max_fit;
		}
	} // if not unwriteable equation
	else{ // bad equation, assign maximum fitness
		me.abserror=p.max_fit;
		me.corr = p.min_fit;
		me.VAF = p.min_fit;
		me.fitness = p.max_fit;
		me.fitness_v = p.max_fit;
		me.genty = p.max_fit;
		if (p.train){
			me.abserror_v=p.max_fit;
			me.corr_v = p.min_fit;
			me.VAF_v = p.min_fit;
		}
	}
			
			 
}
void CalcOutput(ind& me,params& p,vector<vector<float>>& vals,vector<float>& dattovar,vector<float>& target,state& s)
{		
	vector<float> stack_float;
	vector<bool> stack_bool;
	me.reset_introns();
	me.output.resize(0); 
	me.output_v.resize(0);
	if (p.sel==3) 
		me.error.resize(0);
	float SStot=0;
	float SSreg=0;
	float SSres=0;
	float q = 0;
	float var_target = 0;
	float var_ind = 0;
	bool pass = true;
	float meanout=0;
	float meantarget=0;
	float meanout_v=0;
	float meantarget_v=0;
	float target_std=1000;
	float target_std_v=1000;
	int ptevals=0;
	unsigned int ndata_t,ndata_v; // training and validation data sizes
	if (p.train){
		ndata_t = vals.size()*p.train_pct;
		ndata_v = vals.size()-ndata_t;
	}
	else{
		ndata_t = vals.size();
		ndata_v=0;
	}
	me.output.reserve(ndata_t);
	me.output_v.reserve(ndata_v);
	//if (p.eHC_on && p.eHC_slim){ me.stack_float.resize(0); me.stack_float.reserve((me.eff_size)*vals.size());}
	if (p.eHC_on && p.eHC_slim)
	{ 
		me.stack_float.resize((me.eff_size)*vals.size());
		me.stack_floatlen.resize(0);
		me.stack_floatlen.reserve(me.eff_size);
	}
	//stack_float.reserve(me.eff_size/2);
	 
	for(unsigned int sim=0;sim<vals.size();++sim)
		{
			//if (p.eHC_slim) me.stack_float.push_back(vector<float>());
			int k_eff=0;

			for (unsigned int j=0; j<p.allvars.size()-p.AR_n;++j) //wgl: add time delay of output variable here 
				dattovar.at(j) = vals[sim][j]; // can we replace this with a pointer to the data so we don't have to copy the whole thing every time?
			if (p.AR){ // auto-regressive output variables
				int ARstart = p.allvars.size()-p.AR_n; 
				for (unsigned int h=0; h<p.AR_n; ++h){
					if (sim<ndata_t){
						if (me.output.size()>h) // add distinction for training / validation data
							if (p.AR_lookahead) dattovar[ARstart+h] = target[sim-1-h];
							else dattovar[ARstart+h] = me.output[sim-1-h];
						else dattovar[ARstart+h] = 0;
					}
					else{
						if (me.output_v.size()>h) // add distinction for training / validation data
							if (p.AR_lookahead) dattovar[ARstart+h] = target[sim-1-h];
							else dattovar[ARstart+h] = me.output_v[sim-1-h-ndata_t];
						else dattovar[ARstart+h] = 0;
					}
				}
			}
			for(int k=0;k<me.line.size();++k){
				/*if(me.line.at(k).type=='v'){
					if (static_pointer_cast<n_sym>(me.line.at(k)).valpt==NULL)
						cout<<"WTF";
				}*/
				if (me.line.at(k).on){
					//me.line.at(k).eval(stack_float);
					eval(me.line.at(k),stack_float,stack_bool);
					++ptevals;
					if (p.eHC_on && p.eHC_slim) // stack tracing
					{
						if(!stack_float.empty()) me.stack_float[k_eff*vals.size() + sim] = stack_float.back();
						else me.stack_float[k_eff*vals.size() + sim] = 0;
						if (sim==0) me.stack_floatlen.push_back(stack_float.size());
						/*if(!stack_float.empty()) me.stack_float.push_back(stack_float.back());
						else me.stack_float.push_back(0);*/
					}
					++k_eff;
					/*else
						cout << "hm";*/
				}
			}
			/*if (stack_float.size()>1)
				cout << "non-tree\n";*/

			//if(!(!p.classification && stack_float.empty()) && !(p.classification && stack_bool.empty())){
			  if(!stack_float.empty()){	

				if ((p.train && sim<ndata_t) || (!p.train)){
						me.output.push_back(stack_float.back());						
						
						if (p.weight_error)
							me.abserror += p.error_weight[sim]*abs(target.at(sim)-me.output.at(sim));
						else
							me.abserror += abs(target.at(sim)-me.output.at(sim));

						if (p.sel==3) // lexicase error vector
							me.error.push_back(abs(target.at(sim)-me.output.at(sim)));						
						
						meantarget += target.at(sim);
						meanout += me.output[sim];
					}
					else //validation set
					{
						me.output_v.push_back(stack_float.back());

						if(p.weight_error)
							me.abserror_v += p.error_weight[sim]*abs(target[sim]-me.output_v[sim-ndata_t]);
						else
							me.abserror_v += abs(target[sim]-me.output_v[sim-ndata_t]);
						
						meantarget_v += target.at(sim);
						meanout_v += me.output_v[sim-ndata_t];
					}
			} //stack_float empty check
			else{
				pass=false;
				break;
				}
			stack_float.resize(0);
			stack_bool.resize(0);
		}
		if (pass){
			assert(me.output.size()==ndata_t);
			// mean absolute error
			me.abserror = me.abserror/ndata_t;
			meantarget = meantarget/ndata_t;
			meanout = meanout/ndata_t;
			//calculate correlation coefficient
			me.corr = getCorr(me.output,target,meanout,meantarget,0,target_std);
			me.VAF = VAF(me.output,target,meantarget,0);
			if (me.corr < 0.001 && me.VAF/100 > .7){
				s.out << "Error in VAF calculation (probably):\n";
				s.out << me.eqn; 
				s.out << "meantarget: " << meantarget << "\n";
				s.out << "corr: " << me.corr << "\n";
				s.out << "VAF: " << me.VAF << "\n";
				s.out << "output: \t";

				for (int i = 0; i<me.output.size(); ++i){
					s.out << me.output[i] << "\t";
				}
				s.out << "\ntarget: \t";
				for (int i = 0; i<me.output.size(); ++i){
					s.out << target[i] << "\t";
				}
				s.out << "\n recalculating with VAF_loud...\n";
				me.VAF = VAF_loud(me.output,target,meantarget,0,s);
				s.out << "\n setting VAF = 0 and continuing...\n";
				me.VAF=0;
			}
			if (p.train)
			{
				q = 0;
				var_target = 0;
				var_ind = 0;
					// mean absolute error
				me.abserror_v = me.abserror_v/ndata_v;
				meantarget_v = meantarget_v/ndata_v;
				meanout_v = meanout_v/ndata_v;
				//calculate correlation coefficient
				me.corr_v = getCorr(me.output_v,target,meanout_v,meantarget_v,ndata_t,target_std_v);
				me.VAF_v = VAF(me.output_v,target,meantarget_v,ndata_t);
			}
		}

		if (me.eqn.compare("1")==0 && me.corr > 0.0001)
				cout << "caught\n";

			if (!pass){
				me.corr = 0;
				me.VAF = 0;
				me.abserror = p.max_fit;
			}
						

			if(me.corr < p.min_fit)
				me.corr=p.min_fit;
			if(me.VAF < p.min_fit)
				me.VAF=p.min_fit;

		    if(me.output.empty())
				me.fitness=p.max_fit;
			else if ( boost::math::isnan(me.abserror) || boost::math::isinf(me.abserror) || boost::math::isnan(me.corr) || boost::math::isinf(me.corr))
				me.fitness=p.max_fit;
			else{
				if (p.fit_type==1)
					me.fitness = me.abserror;
				else if (p.fit_type==2)
					me.fitness = 1-me.corr;
				else if (p.fit_type==3)
					me.fitness = me.abserror/me.corr;
				else if (p.fit_type==4)
					me.fitness = 1-me.VAF/100;
				if (p.norm_error)
					me.fitness = me.fitness/target_std;
			}

		
			if(me.fitness>p.max_fit)
				me.fitness=p.max_fit;
			else if(me.fitness<p.min_fit)
				(me.fitness=p.min_fit);

			if(p.train){ //assign validation fitness
				if (!pass){
					me.corr_v = 0;
					me.VAF_v = 0;
				}
						

				if(me.corr_v < p.min_fit)
					me.corr_v=p.min_fit;
				if(me.VAF_v < p.min_fit)
					me.VAF_v=p.min_fit;

				if(me.output_v.empty())
					me.fitness_v=p.max_fit;
				/*else if (*std::max_element(me.output_v.begin(),me.output_v.end())==*std::min_element(me.output_v.begin(),me.output_v.end()))
					me.fitness_v=p.max_fit;*/
				else if ( boost::math::isnan(me.abserror_v) || boost::math::isinf(me.abserror_v) || boost::math::isnan(me.corr_v) || boost::math::isinf(me.corr_v))
					me.fitness_v=p.max_fit;
				else{
					if (p.fit_type==1)
						me.fitness_v = me.abserror_v;
					else if (p.fit_type==2)
						me.fitness_v = 1-me.corr_v;
					else if (p.fit_type==3)
						me.fitness_v = me.abserror_v/me.corr_v;
					else if (p.fit_type==4)
						me.fitness_v = 1-me.VAF_v/100;
					if (p.norm_error)
						me.fitness_v = me.fitness_v/target_std_v;

				}

		
				if(me.fitness_v>p.max_fit)
					me.fitness_v=p.max_fit;
				else if(me.fitness_v<p.min_fit)
					(me.fitness_v=p.min_fit);
			}
			else{ // if not training, assign copy of regular fitness to the validation fitness variables
				me.corr_v=me.corr;
				me.abserror_v=me.abserror;
				me.VAF_v = me.VAF;
				me.fitness_v=me.fitness;
			}
			//if (p.estimate_generality || p.PS_sel==2){
			//	if (me.fitness == p.max_fit || me.fitness_v== p.max_fit)
			//		me.genty = p.max_fit;
			//	else{
			//		if (p.G_sel==1) // MAE
			//			me.genty = abs(me.abserror-me.abserror_v)/me.abserror;
			//		else if (p.G_sel==2) // R2
			//			me.genty = abs(me.corr-me.corr_v)/me.corr;
			//		else if (p.G_sel==3) // MAE R2 combo
			//			me.genty = abs(me.abserror/me.corr-me.abserror_v/me.corr_v)/(me.abserror/me.corr);
			//		else if (p.G_sel==3) // VAF
			//			me.genty = abs(me.VAF-me.VAF_v)/me.VAF;
			//	}
			//}
			int tmp = omp_get_thread_num();
		s.ptevals[omp_get_thread_num()]=s.ptevals[omp_get_thread_num()]+ptevals;
}
void Calc_M3GP_Output(ind& me,params& p,vector<vector<float>>& vals,vector<float>& dattovar,vector<float>& target,state& s)
{
	vector<float> stack_float;
	vector<bool> stack_bool;
	
	me.reset_introns();
	me.output.resize(0); 
	me.output_v.resize(0);
	if (p.sel==3){ 
		if (p.lex_class) 
			me.error.assign(p.number_of_classes,0);
		else
			me.error.resize(0);
	}
	else if (p.sel==4 && p.PS_sel>=4) // each class is an objective
		me.error.assign(p.number_of_classes,0);

	float var_target = 0;
	float var_ind = 0;
	bool pass = true;
	int ptevals=0;
	unsigned int ndata_t,ndata_v; // training and validation data sizes
	if (p.train){
		ndata_t = vals.size()*p.train_pct;
		ndata_v = vals.size()-ndata_t;
	}
	else{
		ndata_t = vals.size();
		ndata_v=0;
	}
	me.output.reserve(ndata_t);
	me.output_v.reserve(ndata_v);
	//if (p.eHC_on && p.eHC_slim){ me.stack_float.resize(0); me.stack_float.reserve((me.eff_size)*vals.size());}
	if (p.eHC_on && p.eHC_slim)
	{ 
		me.stack_float.resize((me.eff_size)*vals.size());
		me.stack_floatlen.resize(0);
		me.stack_floatlen.reserve(me.eff_size);
	}
	//stack_float.reserve(me.eff_size/2);

	//vector<vector<float>> Z; // n_t x d output for m3gp
	MatrixXf Z(ndata_t,1);
	//vector<vector<float>> Z_v; // n_v x d output for m3gp
	MatrixXf Z_v(ndata_v,1);
	vector<MatrixXf> K(p.number_of_classes); // output subsets by class label
	
	// Calculate Output
	for(unsigned int sim=0;sim<vals.size();++sim)
		{
			//if (p.eHC_slim) me.stack_float.push_back(vector<float>());
			int k_eff=0;

			for (unsigned int j=0; j<p.allvars.size();++j) 
				dattovar.at(j)= vals[sim][j];
			
			for(int k=0;k<me.line.size();++k){
				if (me.line.at(k).on){
					eval(me.line.at(k),stack_float,stack_bool);
					++ptevals;
					if (p.eHC_on && p.eHC_slim) // stack tracing
					{
						if(!stack_float.empty()) me.stack_float[k_eff*vals.size() + sim] = stack_float.back();
						else me.stack_float[k_eff*vals.size() + sim] = 0;
						if (sim==0) me.stack_floatlen.push_back(stack_float.size());
						/*if(!stack_float.empty()) me.stack_float.push_back(stack_float.back());
						else me.stack_float.push_back(0);*/
					}
					++k_eff;
					/*else
						cout << "hm";*/
				}
			}
			/*if (stack_float.size()>1)
				cout << "non-tree\n";*/

			//if(!(!p.classification && stack_float.empty()) && !(p.classification && stack_bool.empty())){
			if(!stack_float.empty()){	
				Z.resize(ndata_t,stack_float.size());
				if ((p.train && sim<ndata_t) || (!p.train)){
																	
							//Z.push_back(vector<float>());
							for (int z =0;z<stack_float.size();++z){
								//Z[sim].push_back(stack_float[z]);	
								Z(sim,z) = stack_float[z];
								//K[target[sim]].push_back(stack_float[z]);	
							}
							//VectorXf Kadd;					
							
							/*cout << K[target[sim]], 
								         Z.row(sim);*/
							K[target[sim]].resize(K[target[sim]].rows(),Z.cols());

					//		cout << "K size: " << K[target[sim]].rows() << " x " <<  K[target[sim]].cols() << endl;
							
					//		cout << "K new size: " << K[target[sim]].rows() << " x " <<  K[target[sim]].cols() << endl;
					//		cout << "Z: " << Z.row(sim) << endl;
							MatrixXf Kadd(K[target[sim]].rows()+1,Z.cols());
							Kadd << K[target[sim]], 
								    Z.row(sim); 
					//		cout << "Kadd: " << Kadd <<endl;
							//K[target[sim]].bottomRows(1) = Z.row(sim);
							K[target[sim]].swap(Kadd);
							//K[target[sim]].set( (MatrixXd(K[target[sim]].rows()+1,Z) << mat, vec.transpose()).finished() );
					//		cout << "K: " << K[target[sim]] << endl;
							//K[target[sim]].set( (MatrixXf(
							//K[target[sim]] << K[target[sim]], 
								              //Z.row(sim);
			// for multiclass classification, error should be based solely on whether or not the right class was assigned						
						
					}
					else //validation set
					{
						Z_v.resize(ndata_v,stack_float.size());
						//Z_v.push_back(vector<float>());
						for (int z =0;z<stack_float.size();++z)
							//Z_v[sim-ndata_t].push_back(stack_float[z]);	
							Z_v(sim-ndata_t,z) = stack_float[z];

					}
			}
			else{
				pass=false;
				break;
				}
			//stack_float.clear();
			stack_float.resize(0);
			stack_bool.resize(0);
		}
	// Calculate Covariance and Centroids of model output data by class 
	/////////////////////////////////////////////////////////////////////
	VectorXf D(p.number_of_classes);
	VectorXf D_v(p.number_of_classes);
	if (pass){
		bool pass2=true; // pass check for invertible covariance matrix
		me.M.resize(p.number_of_classes,Z.cols());
		
		std::vector<MatrixXf> Cinv(p.number_of_classes, MatrixXf(Z.cols(),Z.cols()));
		for (int i = 0; i<p.number_of_classes; ++i){
			me.C.push_back(MatrixXf(Z.cols(),Z.cols()));
			//me.M.push_back(vector<float>());
			//cout << K[i].colwise().mean() << endl;
			me.M.row(i) << K[i].colwise().mean();
			//cout << "M centroids for class " << i << ":\n" << me.M.row(i) << endl;
			/*s.out << "/////////////////////////////\n";
			s.out << "K(" << i << ")\n";
			s.out << K[i] << "\n";*/

			cov(K[i],me.C[i]);
			Eigen::FullPivLU<MatrixXf> check(me.C[i]);
			if (check.isInvertible()){
				//s.out << "C: \n";
				//s.out << me.C[i] << "\n";
				Cinv[i] = me.C[i].inverse();
				//s.out << "C^-1:\n";
				//s.out << Cinv[i] << "\n";
			}
			else{
				//s.out << "C: \n";
				//s.out << me.C[i] << "\n";
				//Cinv[i] = me.C[i].inverse();
				Cinv[i] = MatrixXf::Identity(Cinv[i].rows(),Cinv[i].cols());
				//s.out << "C^-1:\n";
				//s.out << Cinv[i] << "\n";
				//pass=false;
				//pass2=false;
			}
		 
		
		
		}
		//vector<float> D; // mahalanobis distance on training set
		//vector<float> D_v; // mahalanobis distance on validation set
		
		// Calculate Fitness
		////////////////////////////////////////////////////////////////////////
		if (pass2){
			for (int sim=0;sim<vals.size();++sim){
				if ((p.train && sim<ndata_t) || (!p.train)){
					//D.resize(0);
					//training set mahalanobis distance
			
					MahalanobisDistance(Z.row(sim),Cinv,me.M,D,s);
			

					/*s.out << "Z:\n";
					s.out << Z.row(sim) << "\n"; 
					s.out << "D:\n";
					s.out << D << "\n";*/
					//assign class based on minimum distance 
					//vector<float>::iterator it = min_element(D.begin(),D.end());
					MatrixXf::Index min_i;
					float min = D.minCoeff(&min_i);
					//float ans = min_i;
					me.output.push_back(float(min_i));

					// assign error
					if (target[sim]!=me.output[sim]){
						++me.abserror;
						if (p.sel==3){ // lexicase error vector
							if (p.lex_class) 
								++me.error[target[sim]];
							else
								me.error.push_back(1);
						}
						else if (p.sel==4 && p.PS_sel>=4) // each class is an objective
							++me.error[target[sim]];

					}
					else if (p.sel==3 && !p.lex_class)
						me.error.push_back(0);
				//assign error
				}
				else{	
					//D_v.resize(0);
					//validation set mahalanobis distance
					MahalanobisDistance(Z_v.row(sim-ndata_t),Cinv,me.M,D_v,s);
			
					//assign class based on minimum distance in D_v[sim]
					//vector<float>::iterator it = min_element(D_v.begin(),D_v.end());
					//me.output_v.push_back(float(it - D_v.begin()));
					MatrixXf::Index min_i;
					float min = D_v.minCoeff(&min_i);
					//float ans = (min_i);
					me.output_v.push_back(float(min_i));

					if (target[sim]!=me.output_v[sim-ndata_t])
							++me.abserror_v;													
						
		
				}
			}
		}
	}
		/////////////////////////////////////////////////////////////////////////////////
		if (pass){
			assert(me.output.size()==ndata_t);
			// mean absolute error
			me.abserror = me.abserror/ndata_t;
			
			if (p.train)//mean absolute error
				me.abserror_v = me.abserror_v/ndata_v;
			
		}

			if (!pass){
				me.abserror = p.max_fit;
			}


		    if(me.output.empty())
				me.fitness=p.max_fit;
			else if ( boost::math::isnan(me.abserror) || boost::math::isinf(me.abserror) )
				me.fitness=p.max_fit;
			else{
				if (p.fit_type!=1){
					s.out << "warning: fit_type not set to error. using error anyway (because classification is being output)\n";
					p.fit_type=1;
				}
				me.fitness = me.abserror;
									
				/*if (p.norm_error)
					me.fitness = me.fitness/target_std;*/
			}

			for (unsigned z=0;z<D.size();++z){
				if(!boost::math::isfinite(D(z)))
					me.fitness = p.max_fit;
			}
			if(me.fitness>p.max_fit)
				me.fitness=p.max_fit;
			else if(me.fitness<p.min_fit)
				(me.fitness=p.min_fit);

			if(p.train){ //assign validation fitness
				
				for (unsigned z=0;z<D_v.size();++z){
				if(!boost::math::isfinite(D_v(z)))
					me.fitness = p.max_fit;
				}

				if(me.output_v.empty())
					me.fitness_v=p.max_fit;
				/*else if (*std::max_element(me.output_v.begin(),me.output_v.end())==*std::min_element(me.output_v.begin(),me.output_v.end()))
					me.fitness_v=p.max_fit;*/
				else if ( boost::math::isnan(me.abserror_v) || boost::math::isinf(me.abserror_v))
					me.fitness_v=p.max_fit;
				else{
					if (p.fit_type!=1){
						s.out << "WARNING: fit_type not set to error. using error anyway (because classification is being output)\n";
						p.fit_type=1;
					}
					me.fitness_v = me.abserror_v;
					/*if (p.norm_error)
						me.fitness_v = me.fitness_v/target_std_v;*/

				}

		
				if(me.fitness_v>p.max_fit)
					me.fitness_v=p.max_fit;
				else if(me.fitness_v<p.min_fit)
					(me.fitness_v=p.min_fit);
			}
			else{ // if not training, assign copy of regular fitness to the validation fitness variables
				me.abserror_v=me.abserror;
				me.fitness_v=me.fitness;
			}
			//if (p.estimate_generality || p.PS_sel==2){
			//	if (me.fitness == p.max_fit || me.fitness_v== p.max_fit)
			//		me.genty = p.max_fit;
			//	else{
			//		if (p.G_sel==1) // MAE
			//			me.genty = abs(me.abserror-me.abserror_v)/me.abserror;
			//		else if (p.G_sel==2) // R2
			//			me.genty = abs(me.corr-me.corr_v)/me.corr;
			//		else if (p.G_sel==3) // MAE R2 combo
			//			me.genty = abs(me.abserror/me.corr-me.abserror_v/me.corr_v)/(me.abserror/me.corr);
			//		else if (p.G_sel==3) // VAF
			//			me.genty = abs(me.VAF-me.VAF_v)/me.VAF;
			//	}
			//}
			int tmp = omp_get_thread_num();
		s.ptevals[omp_get_thread_num()]=s.ptevals[omp_get_thread_num()]+ptevals;
}

void CalcClassOutput(ind& me,params& p,vector<vector<float>>& vals,vector<float>& dattovar,vector<float>& target,state& s)
{		
	vector<float> stack_float;
	vector<bool> stack_bool;
	vector<vector<float>> Z; // n x d output for m3gp
	vector<vector<float>> K; // output subsets by class label
	me.reset_introns();
	me.output.resize(0); 
	me.output_v.resize(0);
	me.abserror = 0;
	me.abserror_v = 0;

	if (p.sel==3){ 
		if (p.lex_class) 
			me.error.assign(p.number_of_classes,0);
		else
			me.error.resize(0);
	}
	else if (p.sel==4 && p.PS_sel>=4) // each class is an objective
		me.error.assign(p.number_of_classes,0);

	float SStot=0;
	float SSreg=0;
	float SSres=0;
	float q = 0;
	float var_target = 0;
	float var_ind = 0;
	bool pass = true;
	float meanout=0;
	float meantarget=0;
	float meanout_v=0;
	float meantarget_v=0;
	float target_std=1000;
	float target_std_v=1000;
	int ptevals=0;
	unsigned int ndata_t,ndata_v; // training and validation data sizes
	if (p.train){
		ndata_t = vals.size()*p.train_pct;
		ndata_v = vals.size()-ndata_t;
	}
	else{
		ndata_t = vals.size();
		ndata_v=0;
	}
	me.output.reserve(ndata_t);
	me.output_v.reserve(ndata_v);
	//if (p.eHC_on && p.eHC_slim){ me.stack_float.resize(0); me.stack_float.reserve((me.eff_size)*vals.size());}
	if (p.eHC_on && p.eHC_slim)
	{ 
		me.stack_float.resize((me.eff_size)*vals.size());
		me.stack_floatlen.resize(0);
		me.stack_floatlen.reserve(me.eff_size);
	}
	//stack_float.reserve(me.eff_size/2);
	 
	for(unsigned int sim=0;sim<vals.size();++sim)
		{
			//if (p.eHC_slim) me.stack_float.push_back(vector<float>());
			int k_eff=0;

			for (unsigned int j=0; j<p.allvars.size();++j) //wgl: add time delay of output variable here 
				dattovar.at(j)= vals[sim][j];
			
			for(int k=0;k<me.line.size();++k){
				if (me.line.at(k).on){
					eval(me.line.at(k),stack_float,stack_bool);
					++ptevals;
					if (p.eHC_on && p.eHC_slim) // stack tracing
					{
						if(!stack_float.empty()) me.stack_float[k_eff*vals.size() + sim] = stack_float.back();
						else me.stack_float[k_eff*vals.size() + sim] = 0;
						if (sim==0) me.stack_floatlen.push_back(stack_float.size());
						/*if(!stack_float.empty()) me.stack_float.push_back(stack_float.back());
						else me.stack_float.push_back(0);*/
					}
					++k_eff;
					/*else
						cout << "hm";*/
				}
			}
			/*if (stack_float.size()>1)
				cout << "non-tree\n";*/

			//if(!(!p.classification && stack_float.empty()) && !(p.classification && stack_bool.empty())){
			if((p.class_bool && !stack_bool.empty()) || (!p.class_bool && !stack_float.empty())){	

				if ((p.train && sim<ndata_t) || (!p.train)){
						 if (p.class_bool){
							// use stack_bool 
							std::bitset<4> bitout;//(p.number_of_classes,0);
							for (int i = 0;i<bitlen(p.number_of_classes); ++i){
								if (stack_bool.size()>i) // copy from the back of the stack towards the front
									bitout.set(i,stack_bool[stack_bool.size()-1-i]); 
							}
							// interpret output as the integer represented by the bitstring produced by stack_bool
							me.output.push_back(float(bitout.to_ulong()));
						}
						else{ 
							// interpret output as the index of the largest float in stack_float
							vector<float>::iterator it = max_element(stack_float.begin(),stack_float.end());
							me.output.push_back(float(it - stack_float.begin()));
						}
							
						// error is based solely on whether or not the right class was assigned
						if (target[sim]!=me.output[sim]){
							++me.abserror;
							if (p.sel==3){ // lexicase error vector
								if (p.lex_class) 
									++me.error[target[sim]];
								else
									me.error.push_back(1);
							}
							else if (p.sel==4 && p.PS_sel>=4) // each class is an objective
								++me.error[target[sim]];
						}
						else if (p.sel==3 && !p.lex_class)
							me.error.push_back(0);						
					}
					else //validation set
					{
						if (p.class_bool){	
							//bool stack							
							std::bitset<4> bitout; //(p.number_of_classes,0);
							for (int i = 0;i<bitlen(p.number_of_classes); ++i){
								if (stack_bool.size()>i) // copy from back of stack forward
									bitout.set(i,stack_bool[stack_bool.size()-1-i]); 
							}
							me.output_v.push_back(float(bitout.to_ulong()));
						}
						else{
							// class is index of largest element in floating point stack
							vector<float>::iterator it = max_element(stack_float.begin(),stack_float.end());
							me.output_v.push_back(float(it - stack_float.begin()));	
						}

					 // error based solely on whether or not the right class was assigned
						if (target[sim]!=me.output_v[sim-ndata_t])
							++me.abserror_v;													
					}
			} // if stack empty
			else{
				pass=false;
				break;
				}
			stack_float.resize(0);
			stack_bool.resize(0);
		} // end for

		if (pass){
			assert(me.output.size()==ndata_t);
			// mean absolute error
			me.abserror = me.abserror/ndata_t;
						
			if (p.train)
			{
				// mean absolute error
				me.abserror_v = me.abserror_v/ndata_v;
				meantarget_v = meantarget_v/ndata_v;
				meanout_v = meanout_v/ndata_v;
			}
		}

		if (!pass){
			me.abserror = p.max_fit;
		}


		if(me.output.empty())
			me.fitness=p.max_fit;
		else if ( boost::math::isnan(me.abserror) || boost::math::isinf(me.abserror) )
			me.fitness=p.max_fit;
		else{
			if (p.fit_type!=1){
				s.out << "warning: fit_type not set to error. using error anyway (because classification is being output)\n";
				p.fit_type=1;
			}
			me.fitness = me.abserror;								
		}
		
		if(me.fitness>p.max_fit)
			me.fitness=p.max_fit;
		else if(me.fitness<p.min_fit)
			(me.fitness=p.min_fit);

		if(p.train){ //assign validation fitness				
			if(me.output_v.empty())
				me.fitness_v=p.max_fit;
			else if ( boost::math::isnan(me.abserror_v) || boost::math::isinf(me.abserror_v))
				me.fitness_v=p.max_fit;
			else{
				if (p.fit_type!=1){
					s.out << "WARNING: fit_type not set to error. using error anyway (because classification is being output)\n";
					p.fit_type=1;
				}
				me.fitness_v = me.abserror_v;
			}

		
			if(me.fitness_v>p.max_fit)
				me.fitness_v=p.max_fit;
			else if(me.fitness_v<p.min_fit)
				(me.fitness_v=p.min_fit);
		}
		else{ // if not training, assign copy of regular fitness to the validation fitness variables
			me.abserror_v=me.abserror;
			me.fitness_v=me.fitness;
		}
		//if (p.estimate_generality || p.PS_sel==2){
		//	if (me.fitness == p.max_fit || me.fitness_v== p.max_fit)
		//		me.genty = p.max_fit;
		//	else{
		//		if (p.G_sel==1) // MAE
		//			me.genty = abs(me.abserror-me.abserror_v)/me.abserror;
		//		else if (p.G_sel==2) // R2
		//			me.genty = abs(me.corr-me.corr_v)/me.corr;
		//		else if (p.G_sel==3) // MAE R2 combo
		//			me.genty = abs(me.abserror/me.corr-me.abserror_v/me.corr_v)/(me.abserror/me.corr);
		//		else if (p.G_sel==3) // VAF
		//			me.genty = abs(me.VAF-me.VAF_v)/me.VAF;
		//	}
		//}
		int tmp = omp_get_thread_num();
	s.ptevals[omp_get_thread_num()]=s.ptevals[omp_get_thread_num()]+ptevals;
}
bool CalcSlimOutput(ind& me,params& p,vector<vector<float>>& vals,vector<float>& dattovar,vector<float>& target,state& s,int linestart, float orig_fit)
{
	vector<float> stack_float;// = init_stack;
	vector<bool> stack_bool;
	me.reset_introns();
	me.output.clear();
	me.output_v.clear();
	me.error.clear();
	float SStot=0;
	float SSreg=0;
	float SSres=0;
	float q = 0;
	float var_target = 0;
	float var_ind = 0;
	bool pass = true;
	float meanout=0;
	float meantarget=0;
	float meanout_v=0;
	float meantarget_v=0;
	float target_std=1000;
	float target_std_v=1000;
	int ptevals=0;
	int ndata_t,ndata_v; // training and validation data sizes
	if (p.train){
		ndata_t = vals.size()*p.train_pct;
		ndata_v = vals.size()-ndata_t;
	}
	else{
		ndata_t = vals.size();
		ndata_v=0;
	}
	// initialize stack for new individual from old stack
	
	/*for(int i=0;i<vals.size();++i){
		me.stack_float.push_back(vector<float>());
		for(int j=0;j<outstart;++j)
			me.stack_float[i].push_back(init_stack[i].at(j));
	}*/
	int outstart = int(me.stack_float.size()/vals.size());
	me.stack_float.resize(me.eff_size*vals.size());
	int k_eff = outstart;
	int k_fill = outstart;
	me.stack_floatlen.reserve(me.eff_size);
	//if(!me.stack_floatlen.empty()) stack_float.reserve(*std::max(me.stack_floatlen.begin(),me.stack_floatlen.end()));
	for(unsigned int sim=0;sim<vals.size();++sim)
		{
			//initialize stack_float with correct number of elements
			if (me.stack_floatlen.size()>0 && outstart>0){
				for (unsigned int i=me.stack_floatlen.at(outstart-1);i>0;--i) 
					stack_float.push_back(me.stack_float[(outstart-i)*vals.size()+sim]);
			}
			
			/*if (outstart>0)
				stack_float.push_back(init_stack[sim][outstart-1]);*/

			//me.stack_float.push_back(vector<float>());

			for (unsigned int j=0; j<p.allvars.size();++j)
				dattovar.at(j)= vals[sim][j];
			
			k_eff=outstart;
			for(int k=linestart;k<me.line.size();++k){
				/*if(me.line.at(k).type=='v'){
					if (static_pointer_cast<n_sym>(me.line.at(k)).valpt==NULL)
						cout<<"WTF";
				}*/
				if (me.line.at(k).on){
					//me.line.at(k).eval(stack_float);
					eval(me.line.at(k),stack_float, stack_bool);
					++ptevals;
					if(!stack_float.empty()) me.stack_float[k_eff*vals.size() + sim] = stack_float.back(); 
					else me.stack_float[k_eff*vals.size() + sim] = 0;
					/*if(!stack_float.empty()) me.stack_float.push_back(stack_float.back());
					else me.stack_float.push_back(0);*/
					if (sim==0) me.stack_floatlen.push_back(stack_float.size());
					++k_eff;
				}
					
			}

			if(!stack_float.empty()){
				if (p.train){
					if(sim<ndata_t){
						me.output.push_back(stack_float.back());
						me.abserror += abs(target.at(sim)-me.output.at(sim));
						if(me.abserror/vals.size() > orig_fit)
							return 0;
						meantarget += target.at(sim);
						meanout += me.output[sim];

					}
					else
					{
						me.output_v.push_back(stack_float.back());
						me.abserror_v += abs(target.at(sim)-me.output_v.at(sim-ndata_t));
						meantarget_v += target.at(sim);
						meanout_v += me.output_v[sim-ndata_t];

					}

				}
				else {
					me.output.push_back(stack_float.back());
					me.abserror += abs(target.at(sim)-me.output.at(sim));
					if(me.abserror/vals.size() > orig_fit)
							return 0;
					meantarget += target.at(sim);
					meanout += me.output[sim];
				}
			}
			else{
				return 0;
				}
			//stack_float.clear();
			stack_float.resize(0);
		}
		if (!pass) return 0;
		else
		{
			// mean absolute error
			me.abserror = me.abserror/ndata_t;
			meantarget = meantarget/ndata_t;
			meanout = meanout/ndata_t;
			//calculate correlation coefficient
			me.corr = getCorr(me.output,target,meanout,meantarget,0,target_std);
			me.VAF = VAF(me.output,target,meantarget,0);

			if (p.train)
			{
				q = 0;
				var_target = 0;
				var_ind = 0;
					// mean absolute error
				me.abserror_v = me.abserror_v/ndata_v;
				meantarget_v = meantarget_v/ndata_v;
				meanout_v = meanout_v/ndata_v;
				//calculate correlation coefficient
				me.corr_v = getCorr(me.output_v,target,meanout_v,meantarget_v,ndata_t,target_std_v);
				me.VAF_v = VAF(me.output_v,target,meantarget_v,ndata_t);
			}
		}

		if (me.eqn.compare("1")==0 && me.corr > 0.0001)
				cout << "caught\n";

			if (!pass){
				me.corr = 0;
				me.VAF = 0;
			}

			if(me.corr < p.min_fit)
				me.corr=p.min_fit;
			if(me.VAF < p.min_fit)
				me.VAF=p.min_fit;

		    if(me.output.empty())
				me.fitness=p.max_fit;
			else if ( boost::math::isnan(me.abserror) || boost::math::isinf(me.abserror) || boost::math::isnan(me.corr) || boost::math::isinf(me.corr))
			{
				me.fitness=p.max_fit;
				me.abserror=p.max_fit;
				me.corr = p.min_fit;
				me.VAF = p.min_fit;
			}
			else{
				if (p.fit_type==1)
					me.fitness = me.abserror;
				else if (p.fit_type==2)
					me.fitness = 1-me.corr;
				else if (p.fit_type==3)
					me.fitness = me.abserror/me.corr;
				else if (p.fit_type==4)
					me.fitness = 1-me.VAF/100;
				if (p.norm_error)
					me.fitness = me.fitness/target_std;
			}

		
			if(me.fitness>p.max_fit)
				me.fitness=p.max_fit;
			else if(me.fitness<p.min_fit)
				(me.fitness=p.min_fit);

			if(p.train){ //assign validation fitness
				if (!pass){
					me.corr_v = 0;
					me.VAF_v = 0;
				}

				if(me.corr_v < p.min_fit)
					me.corr_v=p.min_fit;
				if(me.VAF_v < p.min_fit)
					me.VAF_v=p.min_fit;

				if(me.output_v.empty())
					me.fitness_v=p.max_fit;
				/*else if (*std::max_element(me.output_v.begin(),me.output_v.end())==*std::min_element(me.output_v.begin(),me.output_v.end()))
					me.fitness_v=p.max_fit;*/
				else if ( boost::math::isnan(me.abserror_v) || boost::math::isinf(me.abserror_v) || boost::math::isnan(me.corr_v) || boost::math::isinf(me.corr_v))
					me.fitness_v=p.max_fit;
				else{
					if (p.fit_type==1)
						me.fitness_v = me.abserror_v;
					else if (p.fit_type==2)
						me.fitness_v = 1-me.corr_v;
					else if (p.fit_type==3)
						me.fitness_v = me.abserror_v/me.corr_v;
					else if (p.fit_type==4)
						me.fitness_v = 1-me.VAF/100;
					if (p.norm_error)
						me.fitness_v = me.fitness_v/target_std_v;

				}

		
				if(me.fitness_v>p.max_fit)
					me.fitness_v=p.max_fit;
				else if(me.fitness_v<p.min_fit)
					(me.fitness_v=p.min_fit);
			}
			else{ // if not training, assign copy of regular fitness to the validation fitness variables
				me.corr_v=me.corr;
				me.VAF_v = me.VAF;
				me.abserror_v=me.abserror;
				me.fitness_v=me.fitness;
			}
			/*if (p.estimate_generality || p.PS_sel==2){
				if (me.fitness == p.max_fit || me.fitness_v== p.max_fit)
					me.genty = p.max_fit;
				else
					me.genty = abs(me.fitness-me.fitness_v)/me.fitness;
			}*/
			int tmp = omp_get_thread_num();
		s.ptevals[omp_get_thread_num()]=s.ptevals[omp_get_thread_num()]+ptevals;
		return true;
}
bool getSlimFit(ind& me,params& p,Data& d,state& s,FitnessEstimator& FE,int linestart, float orig_fit)
	
{
    //set up data table for conversion of symbolic variables
	unordered_map <string,float*> datatable;
	vector<float> dattovar(d.label.size());

	for (unsigned int i=0;i<d.label.size(); ++i)
			datatable.insert(pair<string,float*>(d.label[i],&dattovar[i]));

	vector<vector<float>> FEvals;
	vector<float> FEtarget;
	setFEvals(FEvals,FEtarget,FE,d);
		
	me.abserror = 0;
	me.abserror_v = 0;
	
	// set data table
	for(int m=0;m<me.line.size();++m){
		if(me.line.at(m).type=='v')
			{// set pointer to dattovar 
				//float* set = datatable.at(static_pointer_cast<n_sym>(me.line.at(m)).varname);
				float* set = datatable.at(me.line.at(m).varname);
				if(set==NULL)
					cout<<"hmm";
				/*static_pointer_cast<n_sym>(me.line.at(m)).setpt(set);*/
				me.line.at(m).setpt(set);
				/*if (static_pointer_cast<n_sym>(me.line.at(m)).valpt==NULL)
					cout<<"wth";*/
			}
	}
	//cout << "Equation" << count << ": f=" << me.eqn << "\n";
			//bool pass=true;
			if(me.eqn.compare("unwriteable")==0)
				return 0;
			else
			{
			    // calculate error 
				if(!p.EstimateFitness){
					return CalcSlimOutput(me,p,d.vals,dattovar,d.target,s,linestart,orig_fit);
				} // if not estimate fitness 	
				else{ 
					return CalcSlimOutput(me,p,FEvals,dattovar,FEtarget,s,linestart,orig_fit);
				} // if estimate fitness
			} // if not unwriteable equation
			 
}
bool SlimFitness(ind& me,params& p,Data& d,state& s,FitnessEstimator& FE, int linestart,  float orig_fit)
{

	me.eqn = Line2Eqn(me.line,me.eqn_form,p,true);
	//getEqnForm(me.eqn,me.eqn_form);
	
	me.eff_size = getEffSize(me.line);

	bool pass = getSlimFit(me,p,d,s,FE,linestart,orig_fit);
	
	me.complexity= getComplexity(me,p);
	s.numevals[omp_get_thread_num()]=s.numevals[omp_get_thread_num()]+1;
	return pass;
}
//void LexicaseFitness(ind& me,params& p,Data& d,state& s,FitnessEstimator& FE)
//{
//	//set up data table for conversion of symbolic variables
//	unordered_map <string,float*> datatable;
//	vector<float> dattovar(d.label.size());
//
//	for (unsigned int i=0;i<d.label.size(); ++i)
//			datatable.insert(pair<string,float*>(d.label[i],&dattovar[i]));
//
//	int ndata_t,ndata_v; // training and validation data sizes
//	
//	vector<vector<float>> FEvals;
//	vector<float> FEtarget;
//	setFEvals(FEvals,FEtarget,FE,d);
//	
//	if (p.train){
//		if(!p.EstimateFitness){
//			ndata_t = d.vals.size()*p.train_pct;
//			ndata_v = d.vals.size()-ndata_t;
//		}
//		else{
//			ndata_t = FEvals.size()*p.train_pct;
//			ndata_v = FEvals.size()-ndata_t;
//		}
//	}
//	else{
//		if(!p.EstimateFitness){
//			ndata_t = d.vals.size();
//			ndata_v=0;
//		}
//		else{
//			ndata_t = FEvals.size();
//			ndata_v=0;
//		}
//
//	}
//
//	int ptevals=0;
//	//get equation and equation form
//			me.eqn = Line2Eqn(me.line,me.eqn_form,p);
//			//getEqnForm(me.eqn,me.eqn_form);
//		//Get effective size and complexity
//			me.eff_size=0;
//			for(int m=0;m<me.line.size();++m){
//
//				if(me.line.at(m).on)
//					++me.eff_size;
//			}
//			// Get Complexity
//			me.complexity= getComplexity(me.line);
//			
//			// set data pointers
//			for(int m=0;m<me.line.size();++m){
//				if(me.line.at(m).type=='v')
//					{// set pointer to dattovar 
//						/*float* set = datatable.at(static_pointer_cast<n_sym>(me.line.at(m)).varname);*/
//						float* set = datatable.at(me.line.at(m).varname);
//						if(set==NULL)
//							cout<<"hmm";
//						/*static_pointer_cast<n_sym>(me.line.at(m)).setpt(set);*/
//						me.line.at(m).setpt(set);
//						/*if (static_pointer_cast<n_sym>(me.line.at(m)).valpt==NULL)
//							cout<<"wth";*/
//					}
//			}
//
//			me.abserror = 0;
//			me.abserror_v = 0;
//			float meanout=0;
//			float meanoutlex;
//			float meantarget=0;
//			float meantargetlex;
//			float meanout_v=0;
//			float meantarget_v=0;
//			float target_std_lex, target_std_lex_v;
//			// set pointer to dattovar in symbolic functions
//
//			//cout << "Equation" << count << ": f=" << me.eqn << "\n";
//			me.fitlex.resize(p.numcases);
//			bool pass=true;
//			if(!me.eqn.compare("unwriteable")==0){
//				
//				vector<float> stack_float;
//				vector<bool> stack_bool;
//				me.output.clear();
//				me.output_v.clear();
//				
//				vector<float> tarlex; // target data arranged for comparison with output
//				vector<float> tarlex_v; // target data arranged for comparison with output_v
//
//				vector<vector<float>> outlex;
//				vector<float> errorlex;
//				vector<float> corrlex;
//				vector<float> VAFlex;
//				//vector<vector<float>> outlex_v;
//				//vector<float> errorlex_v;
//				//vector<float> corrlex_v;
//
//		// loop through numcases
//				for (int lex=0; lex<d.lexvals.size(); ++lex)
//				{
//					
//					meanoutlex=0;
//					meantargetlex=0;
//					outlex.push_back(vector<float>());
//					errorlex.push_back(0);
//					corrlex.push_back(0);
//					VAFlex.push_back(0);
//
//					if (p.train){
//						ndata_t = d.lexvals[lex].size()*p.train_pct;
//						ndata_v = d.lexvals[lex].size()-ndata_t;
//					}
//					else{
//						ndata_t = d.lexvals[lex].size();
//						ndata_v=0;
//					}	
//
//	// calculate error 
//					for(unsigned int sim=0;sim<d.lexvals[lex].size();++sim)
//					{
//						//outlex_v.push_back(vector<float>());
//
//						for (unsigned int j=0; j<p.allvars.size();++j)
//							dattovar.at(j)= d.lexvals[lex][sim][j];
//
//						for(int k=0;k<me.line.size();++k){
//							if (me.line.at(k).on){
//								//me.line.at(k)->eval(stack_float);
//								eval(me.line.at(k),stack_float,stack_bool);
//								++ptevals;}
//						}
//
//						if(!stack_float.empty()){
//							if (p.train){
//								if(sim<ndata_t){
//									me.output.push_back(stack_float.back());
//									tarlex.push_back(d.targetlex[lex][sim]);
//									outlex[lex].push_back(stack_float.back());
//
//									me.abserror += abs(d.targetlex[lex].at(sim)-me.output.back());
//									errorlex[lex]+=abs(d.targetlex[lex].at(sim)-me.output.back());
//
//									meantarget += d.targetlex[lex].at(sim);
//									meantargetlex += d.targetlex[lex].at(sim);
//									meanout += me.output.back();
//									meanoutlex+= me.output.back();
//
//								}
//								else
//								{
//									me.output_v.push_back(stack_float.back());
//									tarlex_v.push_back(d.targetlex[lex][sim]);
//									//outlex_v[lex].push_back(stack_float.back());
//									me.abserror_v += abs(d.targetlex[lex].at(sim)-me.output_v.back());
//									//errorlex_v[lex]+=abs(d.targetlex[lex]lex[lex].at(sim)-me.output.at(sim-ndata_t));
//									meantarget_v += d.targetlex[lex][sim];
//									//meantargetlex += d.targetlex[lex]lex[lex].at(sim);
//									meanout_v += me.output_v.back();
//									//meanoutlex_v += me.output_v[sim-ndata_t];
//
//								}
//							}
//							else {
//								me.output.push_back(stack_float.back());
//								tarlex.push_back(d.targetlex[lex][sim]);
//								outlex[lex].push_back(stack_float.back());
//								me.abserror += abs(d.targetlex[lex].at(sim)-me.output.back());
//								errorlex[lex]+=abs(d.targetlex[lex].at(sim)-me.output.back());
//								meantarget += d.targetlex[lex].at(sim);
//								meantargetlex += d.targetlex[lex].at(sim);
//								meanout += me.output.back();
//								meanoutlex+= me.output.back();
//							}
//						}
//						else{
//							pass=false;
//							break;
//							}
//						stack_float.clear();
//					}//for(unsigned int sim=0;sim<d.lexvals[lex].size();++sim)
//
//					if (pass){
//					
//						// mean absolute error
//						errorlex[lex] = errorlex[lex]/outlex[lex].size();
//						meantargetlex = meantargetlex/outlex[lex].size();
//						meanoutlex= meanoutlex/outlex[lex].size();
//
//						//calculate correlation coefficient
//						string tmp = me.eqn;
//						corrlex[lex] = getCorr(outlex[lex],d.targetlex[lex],meanoutlex,meantargetlex,0,target_std_lex);
//						VAFlex[lex] = VAF(outlex[lex],d.targetlex[lex],meantargetlex,0);
//
//					}
//					else{
//						errorlex[lex]=p.max_fit;
//						corrlex[lex] = p.min_fit;
//						VAFlex[lex] = p.min_fit;
//					}
//			
//					if(corrlex[lex] < p.min_fit)
//						corrlex[lex]=p.min_fit;
//
//					if(VAFlex[lex] < p.min_fit)
//						VAFlex[lex]=p.min_fit;
//
//					if(me.output.empty())
//						me.fitlex[lex]=p.max_fit;
//					/*else if (*std::max_element(me.output.begin(),me.output.end())==*std::min_element(me.output.begin(),me.output.end()))
//						me.fitness=p.max_fit;*/
//					else if ( boost::math::isnan(errorlex[lex]) || boost::math::isinf(errorlex[lex]) || boost::math::isnan(corrlex[lex]) || boost::math::isinf(corrlex[lex]))
//						me.fitlex[lex]=p.max_fit;
//					else{
//						if (p.fit_type==1)
//							me.fitlex[lex] = errorlex[lex];
//						else if (p.fit_type==2)
//							me.fitlex[lex] = 1-corrlex[lex];
//						else if (p.fit_type==3)
//							me.fitlex[lex] = errorlex[lex]/corrlex[lex];
//						else if (p.fit_type==4)
//							me.fitlex[lex] = 1-VAFlex[lex]/100;
//						if (p.norm_error)
//							me.fitlex[lex] = me.fitlex[lex]/target_std_lex;
//					}
//
//		
//					if(me.fitlex[lex]>p.max_fit)
//						me.fitlex[lex]=p.max_fit;
//					else if(me.fitlex[lex]<p.min_fit)
//						(me.fitlex[lex]=p.min_fit);
//
//				} //for (int lex=0; lex<d.lexvals.size(); lex++)
//	// fill in overall fitness information
//				if (pass){
//					
//						// mean absolute error
//						me.abserror = me.abserror/me.output.size();
//						meantarget = meantarget/me.output.size();
//						meanout = meanout/me.output.size();
//						//meanoutlex= meanoutlex/me.output.size();
//
//						//calculate correlation coefficient
//						me.corr = getCorr(me.output,tarlex,meanout,meantarget,0,target_std_lex);
//						me.VAF = VAF(me.output,tarlex,meantarget,0);
//					
//						if (p.train)
//						{
//							// mean absolute error
//							me.abserror_v = me.abserror_v/me.output_v.size();
//							meantarget_v = meantarget_v/tarlex_v.size();
//							meanout_v = meanout_v/me.output_v.size();
//							//calculate correlation coefficient
//							me.corr_v = getCorr(me.output_v,tarlex_v,meanout_v,meantarget_v,0,target_std_lex_v);
//							me.VAF_v = VAF(me.output_v,tarlex_v,meantarget_v,0);
//						}
//
//					}
//					else{
//						me.abserror=p.max_fit;
//						me.corr = p.min_fit;
//						me.VAF = p.min_fit;
//
//						if (p.train){
//							me.abserror_v=p.max_fit;
//							me.corr_v = p.min_fit;
//							me.VAF_v = p.min_fit;
//						}
//					}						
//
//				if(me.corr < p.min_fit)
//					me.corr=p.min_fit;
//				if(me.VAF < p.min_fit)
//					me.VAF = p.min_fit;
//
//				if(me.output.empty())
//					me.fitness=p.max_fit;
//				/*else if (*std::max_element(me.output.begin(),me.output.end())==*std::min_element(me.output.begin(),me.output.end()))
//					me.fitness=p.max_fit;*/
//				else if ( boost::math::isnan(me.abserror) || boost::math::isinf(me.abserror) || boost::math::isnan(me.corr) || boost::math::isinf(me.corr))
//					me.fitness=p.max_fit;
//				else{
//					if (p.fit_type==1)
//						me.fitness = me.abserror;
//					else if (p.fit_type==2)
//						me.fitness = 1-me.corr;
//					else if (p.fit_type==3)
//						me.fitness = me.abserror/me.corr;
//					else if (p.fit_type==4)
//						me.fitness = 1-me.VAF/100;
//					if (p.norm_error)
//						me.fitness = me.fitness/target_std_lex;
//
//				}
//
//		
//				if(me.fitness>p.max_fit)
//					me.fitness=p.max_fit;
//				else if(me.fitness<p.min_fit)
//					(me.fitness=p.min_fit);
//
//				if(p.train){
//					if (!pass){
//						me.corr_v = 0;
//						me.VAF_v = 0;
//					}
//
//					if(me.corr_v < p.min_fit)
//						me.corr_v=p.min_fit;
//					if(me.VAF_v < p.min_fit)
//						me.VAF_v=p.min_fit;
//
//					if(me.output_v.empty())
//						me.fitness_v=p.max_fit;
//					/*else if (*std::max_element(me.output_v.begin(),me.output_v.end())==*std::min_element(me.output_v.begin(),me.output_v.end()))
//						me.fitness_v=p.max_fit;*/
//					else if ( boost::math::isnan(me.abserror_v) || boost::math::isinf(me.abserror_v) || boost::math::isnan(me.corr_v) || boost::math::isinf(me.corr_v))
//					{
//						me.fitness_v=p.max_fit;
//						me.abserror=p.max_fit;
//						me.corr = p.min_fit;
//						me.VAF = p.min_fit;
//					}
//					else{
//						if (p.fit_type==1)
//							me.fitness_v = me.abserror_v;
//						else if (p.fit_type==2)
//							me.fitness_v = 1-me.corr_v;
//						else if (p.fit_type==3)
//							me.fitness_v = me.abserror_v/me.corr_v;
//						else if (p.fit_type==4)
//							me.fitness = 1-me.VAF/100;
//						if (p.norm_error)
//							me.fitness_v = me.fitness_v/target_std_lex_v;
//					}
//
//		
//					if(me.fitness_v>p.max_fit)
//						me.fitness_v=p.max_fit;
//					else if(me.fitness_v<p.min_fit)
//						(me.fitness_v=p.min_fit);
//				}
//				else{ // if not training, assign copy of regular fitness to validation variales
//					me.corr_v=me.corr;
//					me.VAF_v = me.VAF;
//					me.abserror_v=me.abserror;
//					me.fitness_v=me.fitness;
//				}
//			}//if(!me.eqn.compare("unwriteable")==0)
//	else
//		{
//			for (int i=0;i<p.numcases;++i)
//				me.fitlex[i]=p.max_fit;
//			me.fitness = p.max_fit;
//			me.fitness_v = p.max_fit;
//
//		}
//}
