#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "data.h"
#include "rnd.h"
#include "state.h"
#include "Line2Eqn.h"
#include "EvalEqnStr.h"
#include <unordered_map>
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
void CalcOutput(ind& me,params& p,vector<vector<float>>& vals,vector<float>& dattovar,vector<float>& target,state& s);
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
float VAF(vector<float>& output,vector<float>& target,float meantarget,int off)
{
	float v1,v2;
	float var_target=0;
	float var_diff = 0;
	float q=0;
	float ndata = float(output.size());
	float vaf;
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
	//q = q/(ndata-1); //unbiased esimator
	var_target=var_target/(ndata-1); //unbiased esimator
	var_diff =var_diff/(ndata-1); //unbiased esimator
	if(abs(var_target)<0.0000001 || abs(var_diff)<0.0000001)
		vaf = 0;
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
int getComplexity(string& eqn)
{
	int complexity=0;
	char c;
	for(int m=0;m<eqn.size();++m){
		c=eqn[m];
		
		if(c=='/')
			complexity=complexity+2;
		else if (c=='s'){
			if(m+2<eqn.size()){
				if ( eqn[m+1]=='i' && eqn[m+2] == 'n'){
					complexity=complexity+3;
					m=m+2;
				}
			}
		}
		else if (c=='c'){
			if(m+2<eqn.size()){
				if ( eqn[m+1]=='o' && eqn[m+2] == 's'){
					complexity=complexity+3;
					m=m+2;
				}
			}
		}
		else if (c=='e'){
			if(m+2<eqn.size()){
				if ( eqn[m+1]=='x' && eqn[m+2] == 'p'){
					complexity=complexity+4;
					m=m+2;
				}
			}
		}
		else if (c=='l'){
			if(m+2<eqn.size()){
				if ( eqn[m+1]=='o' && eqn[m+2] == 'g'){
					complexity=complexity+4;
					m=m+2;
				}
			}
		}
		else if (isalpha(c) && (m+1)<eqn.size()){
			bool pass=true;
			while ((m+1)<eqn.size() && pass){
				if (isalpha(eqn[m+1])) ++m; 
				else pass=0;
			}
			++complexity;
		}
		else
			++complexity;
	}

	return complexity;
}
int getEffSize(vector<node>& line)
{
	int eff_size=0;
	for(int m=0;m<line.size();++m){
		if(line.at(m).on)
			++eff_size;
	}
	return eff_size;
}
void eval(node& n,vector<float>& outstack)
{
	switch(n.type) 
	{
	case 'n':
		outstack.push_back(n.value);
		break;
	case 'v':
		/*if (n.valpt==NULL)
			cout<<"problem";
		else*/
		outstack.push_back(*n.valpt);
		break;
	case '+':
		if(outstack.size()>=2){
				float n1 = outstack.back(); outstack.pop_back();
				float n2 = outstack.back(); outstack.pop_back();

				outstack.push_back(n2+n1);
		}
		break;
	case '-':
		if(outstack.size()>=2){
				float n1 = outstack.back(); outstack.pop_back();
				float n2 = outstack.back(); outstack.pop_back();

				outstack.push_back(n2-n1);
		}
		break;
	case '*':
		if(outstack.size()>=2){
				float n1 = outstack.back(); outstack.pop_back();
				float n2 = outstack.back(); outstack.pop_back();

				outstack.push_back(n2*n1);
		}
		break;
	case '/':
		if(outstack.size()>=2){
				float n1 = outstack.back(); outstack.pop_back();
				float n2 = outstack.back(); outstack.pop_back();
				if(abs(n1)<0.0001)
					outstack.push_back(0);
				else
					outstack.push_back(n2/n1);
		}
		break;
	case 's':
		if(outstack.size()>=1){
				float n1 = outstack.back(); outstack.pop_back();
				outstack.push_back(sin(n1));
		}
		break;
	case 'c':
		if(outstack.size()>=1){
				float n1 = outstack.back(); outstack.pop_back();
				outstack.push_back(cos(n1));
		}
		break;
	case 'e':
		if(outstack.size()>=1){
				float n1 = outstack.back(); outstack.pop_back();
				outstack.push_back(exp(n1));
		}
		break;
	case 'l':
		if(outstack.size()>=1){
				float n1 = outstack.back(); outstack.pop_back();
				if (abs(n1)<0.0001)
					outstack.push_back(0);
				else
					outstack.push_back(log(n1));
		}
		break;
	}
}
void FitnessEstimate(vector<ind>& pop,params& p,data& d,state& s,FitnessEstimator& FE);
void StandardFitness(ind& me,params& p,data& d,state& s,FitnessEstimator& FE, unordered_map<string,float*>& datatable, vector<float>& dattovar);
void LexicaseFitness(ind& me,params& p,data& d,state& s,FitnessEstimator& FE);

void Fitness(vector<ind>& pop,params& p,data& d,state& s,FitnessEstimator& FE)
{
	
	//set up data table for conversion of symbolic variables
	unordered_map <string,float*> datatable;
	vector<float> dattovar(d.label.size());

	for (unsigned int i=0;i<d.label.size(); ++i)
			datatable.insert(pair<string,float*>(d.label[i],&dattovar[i]));

	//#pragma omp parallel for private(e)
	for(int count = 0; count<pop.size(); ++count)
	{
		pop.at(count).eqn = Line2Eqn(pop.at(count).line,pop.at(count).eqn_form);
		//getEqnForm(pop.at(count).eqn,pop.at(count).eqn_form);
		pop.at(count).complexity= getComplexity(pop.at(count).eqn_form);
		pop.at(count).eff_size = getEffSize(pop.at(count).line);

		if(p.sel!=3){
			StandardFitness(pop.at(count),p,d,s,FE,datatable,dattovar);		
		} // if p.sel!=3
		else //LEXICASE FITNESS 
		{
			LexicaseFitness(pop.at(count),p,d,s,FE);
		}//LEXICASE FITNESS

		/*if (p.estimate_generality && pop.at(count).genty != abs(pop[count].fitness-pop[count].fitness_v)/pop[count].fitness && pop.at(count).genty != p.max_fit) 
			std::cerr << "genty error, line 300 Fitness.cpp\n";*/
	}//for(int count = 0; count<pop.size(); ++count)
	s.numevals[omp_get_thread_num()]=s.numevals[omp_get_thread_num()]+pop.size();
	//cout << "\nFitness Time: ";
}
void StandardFitness(ind& me,params& p,data& d,state& s,FitnessEstimator& FE,unordered_map<string,float*>& datatable, vector<float>& dattovar)
{
    

	
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
	bool pass=true;
	if(!me.eqn.compare("unwriteable")==0){
							
// calculate error 
		if(!p.EstimateFitness){
			CalcOutput(me,p,d.vals,dattovar,d.target,s);
		} // if not estimate fitness 
		else{ 
// use fitness estimator subset of d.vals ===========================================================	
			vector<vector<float>> FEvals;
			vector<float> FEtarget;
			setFEvals(FEvals,FEtarget,FE,d);
			CalcOutput(me,p,FEvals,dattovar,FEtarget,s);
		} // if estimate fitness
		if (p.estimate_generality || p.PS_sel==2){
				if (me.fitness == p.max_fit || me.fitness_v== p.max_fit)
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
	vector<float> outstack;
	me.output.resize(0); 
	me.output_v.resize(0); 
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
		ndata_t = vals.size()/2;
		ndata_v = vals.size()-ndata_t;
	}
	else{
		ndata_t = vals.size();
		ndata_v=0;
	}
	me.output.reserve(ndata_t);
	me.output_v.reserve(ndata_v);
	//if (p.eHC_on && p.eHC_slim){ me.outstack.resize(0); me.outstack.reserve((me.eff_size)*vals.size());}
	if (p.eHC_on && p.eHC_slim)
	{ 
		me.outstack.resize((me.eff_size)*vals.size());
		me.outstacklen.resize(0);
		me.outstacklen.reserve(me.eff_size);
	}
	//outstack.reserve(me.eff_size/2);
	 
	for(unsigned int sim=0;sim<vals.size();++sim)
		{
			//if (p.eHC_slim) me.outstack.push_back(vector<float>());
			int k_eff=0;

			for (unsigned int j=0; j<p.allvars.size();++j)
				dattovar.at(j)= vals[sim][j];

			for(int k=0;k<me.line.size();++k){
				/*if(me.line.at(k).type=='v'){
					if (static_pointer_cast<n_sym>(me.line.at(k)).valpt==NULL)
						cout<<"WTF";
				}*/
				if (me.line.at(k).on){
					//me.line.at(k).eval(outstack);
					eval(me.line.at(k),outstack);
					++ptevals;
				if (p.eHC_on && p.eHC_slim) // stack tracing
				{
					if(!outstack.empty()) me.outstack[k_eff*vals.size() + sim] = outstack.back();
					else me.outstack[k_eff*vals.size() + sim] = 0;
					if (sim==0) me.outstacklen.push_back(outstack.size());
					/*if(!outstack.empty()) me.outstack.push_back(outstack.back());
					else me.outstack.push_back(0);*/
				}
				++k_eff;
					/*else
						cout << "hm";*/
				}
			}

			if(!outstack.empty()){
				

				if (p.train){
					if(sim<ndata_t){
						me.output.push_back(outstack.back());
						me.abserror += abs(target.at(sim)-me.output.at(sim));
						meantarget += target.at(sim);
						meanout += me.output[sim];

					}
					else
					{
						me.output_v.push_back(outstack.back());
						me.abserror_v += abs(target.at(sim)-me.output_v.at(sim-ndata_t));
						meantarget_v += target.at(sim);
						meanout_v += me.output_v[sim-ndata_t];

					}

				}
				else {
					me.output.push_back(outstack.back());
					me.abserror += abs(target.at(sim)-me.output.at(sim));
					meantarget += target.at(sim);
					meanout += me.output[sim];
				}
			}
			else{
				pass=false;
				break;
				}
			//outstack.clear();
			outstack.resize(0);
		}
		if (pass){
			
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

			if (!pass)
				me.corr = 0;
						

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
				if (!pass)
					me.corr_v = 0;
						

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

		s.ptevals[omp_get_thread_num()]=s.ptevals[omp_get_thread_num()]+ptevals;
}

bool CalcSlimOutput(ind& me,params& p,vector<vector<float>>& vals,vector<float>& dattovar,vector<float>& target,state& s,int linestart, float orig_fit)
{
	vector<float> outstack;// = init_stack;
	me.output.clear();
	me.output_v.clear();
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
		ndata_t = vals.size()/2;
		ndata_v = vals.size()-ndata_t;
	}
	else{
		ndata_t = vals.size();
		ndata_v=0;
	}
	// initialize stack for new individual from old stack
	
	/*for(int i=0;i<vals.size();++i){
		me.outstack.push_back(vector<float>());
		for(int j=0;j<outstart;++j)
			me.outstack[i].push_back(init_stack[i].at(j));
	}*/
	int outstart = int(me.outstack.size()/vals.size());
	me.outstack.resize(me.eff_size*vals.size());
	int k_eff = outstart;
	int k_fill = outstart;
	me.outstacklen.reserve(me.eff_size);
	//if(!me.outstacklen.empty()) outstack.reserve(*std::max(me.outstacklen.begin(),me.outstacklen.end()));
	for(unsigned int sim=0;sim<vals.size();++sim)
		{
			//initialize outstack with correct number of elements
			if (me.outstacklen.size()>0 && outstart>0){
				for (unsigned int i=me.outstacklen.at(outstart-1);i>0;--i) 
					outstack.push_back(me.outstack[(outstart-i)*vals.size()+sim]);
			}
			
			/*if (outstart>0)
				outstack.push_back(init_stack[sim][outstart-1]);*/

			//me.outstack.push_back(vector<float>());

			for (unsigned int j=0; j<p.allvars.size();++j)
				dattovar.at(j)= vals[sim][j];
			
			k_eff=outstart;
			for(int k=linestart;k<me.line.size();++k){
				/*if(me.line.at(k).type=='v'){
					if (static_pointer_cast<n_sym>(me.line.at(k)).valpt==NULL)
						cout<<"WTF";
				}*/
				if (me.line.at(k).on){
					//me.line.at(k).eval(outstack);
					eval(me.line.at(k),outstack);
					++ptevals;
					if(!outstack.empty()) me.outstack[k_eff*vals.size() + sim] = outstack.back(); 
					else me.outstack[k_eff*vals.size() + sim] = 0;
					/*if(!outstack.empty()) me.outstack.push_back(outstack.back());
					else me.outstack.push_back(0);*/
					if (sim==0) me.outstacklen.push_back(outstack.size());
					++k_eff;
				}
					
			}

			if(!outstack.empty()){
				if (p.train){
					if(sim<ndata_t){
						me.output.push_back(outstack.back());
						me.abserror += abs(target.at(sim)-me.output.at(sim));
						if(me.abserror/vals.size() > orig_fit)
							return 0;
						meantarget += target.at(sim);
						meanout += me.output[sim];

					}
					else
					{
						me.output_v.push_back(outstack.back());
						me.abserror_v += abs(target.at(sim)-me.output_v.at(sim-ndata_t));
						meantarget_v += target.at(sim);
						meanout_v += me.output_v[sim-ndata_t];

					}

				}
				else {
					me.output.push_back(outstack.back());
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
			//outstack.clear();
			outstack.resize(0);
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

			if (!pass)
				me.corr = 0;
						

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
			}
			else{
				if (p.fit_type==1)
					me.fitness = me.abserror;
				else if (p.fit_type==2)
					me.fitness = 1-me.corr;
				else if (p.fit_type==3)
					me.fitness = me.abserror/me.corr;
				else if (p.fit_type==4)
					me.fitness = 1-me.VAF;
				if (p.norm_error)
					me.fitness = me.fitness/target_std;
			}

		
			if(me.fitness>p.max_fit)
				me.fitness=p.max_fit;
			else if(me.fitness<p.min_fit)
				(me.fitness=p.min_fit);

			if(p.train){ //assign validation fitness
				if (!pass)
					me.corr_v = 0;
						

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
						me.fitness = 1-me.VAF;
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
				me.fitness_v=me.fitness;
			}
			/*if (p.estimate_generality || p.PS_sel==2){
				if (me.fitness == p.max_fit || me.fitness_v== p.max_fit)
					me.genty = p.max_fit;
				else
					me.genty = abs(me.fitness-me.fitness_v)/me.fitness;
			}*/
		s.ptevals[omp_get_thread_num()]=s.ptevals[omp_get_thread_num()]+ptevals;
		return true;
}
bool getSlimFit(ind& me,params& p,data& d,state& s,FitnessEstimator& FE,int linestart, float orig_fit)
	
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
bool SlimFitness(ind& me,params& p,data& d,state& s,FitnessEstimator& FE, int linestart,  float orig_fit)
{

	me.eqn = Line2Eqn(me.line,me.eqn_form);
	//getEqnForm(me.eqn,me.eqn_form);
	me.complexity= getComplexity(me.eqn_form);
	me.eff_size = getEffSize(me.line);

	bool pass = getSlimFit(me,p,d,s,FE,linestart,orig_fit);
		
	s.numevals[omp_get_thread_num()]=s.numevals[omp_get_thread_num()]+1;
	return pass;
}
void LexicaseFitness(ind& me,params& p,data& d,state& s,FitnessEstimator& FE)
{
	//set up data table for conversion of symbolic variables
	unordered_map <string,float*> datatable;
	vector<float> dattovar(d.label.size());

	for (unsigned int i=0;i<d.label.size(); ++i)
			datatable.insert(pair<string,float*>(d.label[i],&dattovar[i]));

	int ndata_t,ndata_v; // training and validation data sizes
	
	vector<vector<float>> FEvals;
	vector<float> FEtarget;
	setFEvals(FEvals,FEtarget,FE,d);
	
	if (p.train){
		if(!p.EstimateFitness){
			ndata_t = d.vals.size()/2;
			ndata_v = d.vals.size()-ndata_t;
		}
		else{
			ndata_t = FEvals.size()/2;
			ndata_v = FEvals.size()-ndata_t;
		}
	}
	else{
		if(!p.EstimateFitness){
			ndata_t = d.vals.size();
			ndata_v=0;
		}
		else{
			ndata_t = FEvals.size();
			ndata_v=0;
		}

	}

	int ptevals=0;
	//get equation and equation form
			me.eqn = Line2Eqn(me.line,me.eqn_form);
			//getEqnForm(me.eqn,me.eqn_form);
		//Get effective size and complexity
			me.eff_size=0;
			for(int m=0;m<me.line.size();++m){

				if(me.line.at(m).on)
					++me.eff_size;
			}
			// Get Complexity
			me.complexity= getComplexity(me.eqn_form);
			
			// set data pointers
			for(int m=0;m<me.line.size();++m){
				if(me.line.at(m).type=='v')
					{// set pointer to dattovar 
						/*float* set = datatable.at(static_pointer_cast<n_sym>(me.line.at(m)).varname);*/
						float* set = datatable.at(me.line.at(m).varname);
						if(set==NULL)
							cout<<"hmm";
						/*static_pointer_cast<n_sym>(me.line.at(m)).setpt(set);*/
						me.line.at(m).setpt(set);
						/*if (static_pointer_cast<n_sym>(me.line.at(m)).valpt==NULL)
							cout<<"wth";*/
					}
			}

			me.abserror = 0;
			me.abserror_v = 0;
			float meanout=0;
			float meanoutlex;
			float meantarget=0;
			float meantargetlex;
			float meanout_v=0;
			float meantarget_v=0;
			float target_std_lex, target_std_lex_v;
			// set pointer to dattovar in symbolic functions

			//cout << "Equation" << count << ": f=" << me.eqn << "\n";
			me.fitlex.resize(p.numcases);
			bool pass=true;
			if(!me.eqn.compare("unwriteable")==0){
				
				vector<float> outstack;
				me.output.clear();
				me.output_v.clear();
				
				vector<float> tarlex; // target data arranged for comparison with output
				vector<float> tarlex_v; // target data arranged for comparison with output_v

				vector<vector<float>> outlex;
				vector<float> errorlex;
				vector<float> corrlex;
				vector<float> VAFlex;
				//vector<vector<float>> outlex_v;
				//vector<float> errorlex_v;
				//vector<float> corrlex_v;

		// loop through numcases
				for (int lex=0; lex<d.lexvals.size(); ++lex)
				{
					
					meanoutlex=0;
					meantargetlex=0;
					outlex.push_back(vector<float>());
					errorlex.push_back(0);
					corrlex.push_back(0);
					VAFlex.push_back(0);

					if (p.train){
						ndata_t = d.lexvals[lex].size()/2;
						ndata_v = d.lexvals[lex].size()-ndata_t;
					}
					else{
						ndata_t = d.lexvals[lex].size();
						ndata_v=0;
					}	

	// calculate error 
					for(unsigned int sim=0;sim<d.lexvals[lex].size();++sim)
					{
						//outlex_v.push_back(vector<float>());

						for (unsigned int j=0; j<p.allvars.size();++j)
							dattovar.at(j)= d.lexvals[lex][sim][j];

						for(int k=0;k<me.line.size();++k){
							if (me.line.at(k).on){
								//me.line.at(k)->eval(outstack);
								eval(me.line.at(k),outstack);
								++ptevals;}
						}

						if(!outstack.empty()){
							if (p.train){
								if(sim<ndata_t){
									me.output.push_back(outstack.back());
									tarlex.push_back(d.targetlex[lex][sim]);
									outlex[lex].push_back(outstack.back());

									me.abserror += abs(d.targetlex[lex].at(sim)-me.output.back());
									errorlex[lex]+=abs(d.targetlex[lex].at(sim)-me.output.back());

									meantarget += d.targetlex[lex].at(sim);
									meantargetlex += d.targetlex[lex].at(sim);
									meanout += me.output.back();
									meanoutlex+= me.output.back();

								}
								else
								{
									me.output_v.push_back(outstack.back());
									tarlex_v.push_back(d.targetlex[lex][sim]);
									//outlex_v[lex].push_back(outstack.back());
									me.abserror_v += abs(d.targetlex[lex].at(sim)-me.output_v.back());
									//errorlex_v[lex]+=abs(d.targetlex[lex]lex[lex].at(sim)-me.output.at(sim-ndata_t));
									meantarget_v += d.targetlex[lex][sim];
									//meantargetlex += d.targetlex[lex]lex[lex].at(sim);
									meanout_v += me.output_v.back();
									//meanoutlex_v += me.output_v[sim-ndata_t];

								}
							}
							else {
								me.output.push_back(outstack.back());
								tarlex.push_back(d.targetlex[lex][sim]);
								outlex[lex].push_back(outstack.back());
								me.abserror += abs(d.targetlex[lex].at(sim)-me.output.back());
								errorlex[lex]+=abs(d.targetlex[lex].at(sim)-me.output.back());
								meantarget += d.targetlex[lex].at(sim);
								meantargetlex += d.targetlex[lex].at(sim);
								meanout += me.output.back();
								meanoutlex+= me.output.back();
							}
						}
						else{
							pass=false;
							break;
							}
						outstack.clear();
					}//for(unsigned int sim=0;sim<d.lexvals[lex].size();++sim)

					if (pass){
					
						// mean absolute error
						errorlex[lex] = errorlex[lex]/outlex[lex].size();
						meantargetlex = meantargetlex/outlex[lex].size();
						meanoutlex= meanoutlex/outlex[lex].size();

						//calculate correlation coefficient
						string tmp = me.eqn;
						corrlex[lex] = getCorr(outlex[lex],d.targetlex[lex],meanoutlex,meantargetlex,0,target_std_lex);
						VAFlex[lex] = VAF(outlex[lex],d.targetlex[lex],meantargetlex,0);

					}
					else{
						errorlex[lex]=p.max_fit;
						corrlex[lex] = p.min_fit;
						VAFlex[lex] = p.min_fit;
					}
			
					if(corrlex[lex] < p.min_fit)
						corrlex[lex]=p.min_fit;

					if(VAFlex[lex] < p.min_fit)
						VAFlex[lex]=p.min_fit;

					if(me.output.empty())
						me.fitlex[lex]=p.max_fit;
					/*else if (*std::max_element(me.output.begin(),me.output.end())==*std::min_element(me.output.begin(),me.output.end()))
						me.fitness=p.max_fit;*/
					else if ( boost::math::isnan(errorlex[lex]) || boost::math::isinf(errorlex[lex]) || boost::math::isnan(corrlex[lex]) || boost::math::isinf(corrlex[lex]))
						me.fitlex[lex]=p.max_fit;
					else{
						if (p.fit_type==1)
							me.fitlex[lex] = errorlex[lex];
						else if (p.fit_type==2)
							me.fitlex[lex] = 1-corrlex[lex];
						else if (p.fit_type==3)
							me.fitlex[lex] = errorlex[lex]/corrlex[lex];
						else if (p.fit_type==4)
							me.fitlex[lex] = 1-VAFlex[lex];
						if (p.norm_error)
							me.fitlex[lex] = me.fitlex[lex]/target_std_lex;
					}

		
					if(me.fitlex[lex]>p.max_fit)
						me.fitlex[lex]=p.max_fit;
					else if(me.fitlex[lex]<p.min_fit)
						(me.fitlex[lex]=p.min_fit);

				} //for (int lex=0; lex<d.lexvals.size(); lex++)
	// fill in overall fitness information
				if (pass){
					
						// mean absolute error
						me.abserror = me.abserror/me.output.size();
						meantarget = meantarget/me.output.size();
						meanout = meanout/me.output.size();
						//meanoutlex= meanoutlex/me.output.size();

						//calculate correlation coefficient
						me.corr = getCorr(me.output,tarlex,meanout,meantarget,0,target_std_lex);
						me.VAF = VAF(me.output,tarlex,meantarget,0);
					
						if (p.train)
						{
							// mean absolute error
							me.abserror_v = me.abserror_v/me.output_v.size();
							meantarget_v = meantarget_v/tarlex_v.size();
							meanout_v = meanout_v/me.output_v.size();
							//calculate correlation coefficient
							me.corr_v = getCorr(me.output_v,tarlex_v,meanout_v,meantarget_v,0,target_std_lex_v);
							me.VAF_v = VAF(me.output_v,tarlex_v,meantarget_v,0);
						}

					}
					else{
						me.abserror=p.max_fit;
						me.corr = p.min_fit;
						me.VAF = p.min_fit;

						if (p.train){
							me.abserror_v=p.max_fit;
							me.corr_v = p.min_fit;
							me.VAF_v = p.min_fit;
						}
					}						

				if(me.corr < p.min_fit)
					me.corr=p.min_fit;
				if(me.VAF < p.min_fit)
					me.VAF = p.min_fit;

				if(me.output.empty())
					me.fitness=p.max_fit;
				/*else if (*std::max_element(me.output.begin(),me.output.end())==*std::min_element(me.output.begin(),me.output.end()))
					me.fitness=p.max_fit;*/
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
						me.fitness = 1-me.VAF;
					if (p.norm_error)
						me.fitness = me.fitness/target_std_lex;

				}

		
				if(me.fitness>p.max_fit)
					me.fitness=p.max_fit;
				else if(me.fitness<p.min_fit)
					(me.fitness=p.min_fit);

				if(p.train){
					if (!pass)
						me.corr_v = 0;
						

					if(me.corr_v < p.min_fit)
						me.corr_v=p.min_fit;
					if(me.VAF_v < p.min_fit)
						me.VAF_v=p.min_fit;

					if(me.output_v.empty())
						me.fitness_v=p.max_fit;
					/*else if (*std::max_element(me.output_v.begin(),me.output_v.end())==*std::min_element(me.output_v.begin(),me.output_v.end()))
						me.fitness_v=p.max_fit;*/
					else if ( boost::math::isnan(me.abserror_v) || boost::math::isinf(me.abserror_v) || boost::math::isnan(me.corr_v) || boost::math::isinf(me.corr_v))
					{
						me.fitness_v=p.max_fit;
						me.abserror=p.max_fit;
						me.corr = p.min_fit;
					}
					else{
						if (p.fit_type==1)
							me.fitness_v = me.abserror_v;
						else if (p.fit_type==2)
							me.fitness_v = 1-me.corr_v;
						else if (p.fit_type==3)
							me.fitness_v = me.abserror_v/me.corr_v;
						else if (p.fit_type==4)
							me.fitness = 1-me.VAF;
						if (p.norm_error)
							me.fitness_v = me.fitness_v/target_std_lex_v;
					}

		
					if(me.fitness_v>p.max_fit)
						me.fitness_v=p.max_fit;
					else if(me.fitness_v<p.min_fit)
						(me.fitness_v=p.min_fit);
				}
				else{ // if not training, assign copy of regular fitness to validation variales
					me.corr_v=me.corr;
					me.abserror_v=me.abserror;
					me.fitness_v=me.fitness;
				}
			}//if(!me.eqn.compare("unwriteable")==0)
	else
		{
			for (int i=0;i<p.numcases;++i)
				me.fitlex[i]=p.max_fit;
			me.fitness = p.max_fit;
			me.fitness_v = p.max_fit;

		}
}
