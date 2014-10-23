#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "data.h"
#include "rnd.h"
#include "state.h"
#include "Line2Eqn.h"
#include "EvalEqnStr.h"
#include <unordered_map>
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

void getEqnForm(std::string& eqn,std::string& eqn_form)
{
//replace numbers with the letter c
#if defined(_WIN32)
	std::regex e ("(([0-9]+)(\.)([0-9]+))|([0-9]+)");
	std::basic_string<char> tmp = "c";
	std::regex_replace (std::back_inserter(eqn_form), eqn.begin(), eqn.end(), e,tmp);
#else
	boost::regex e ("(([0-9]+)(\.)([0-9]+))|([0-9]+)");
	eqn_form = boost::regex_replace(eqn,e,"c");
#endif 
	//(\d+\.\d+)|(\d+)
	//eqn_form = eqn;
	//eqn_form=std::tr1::regex_replace(eqn,e,tmp.c_str(),std::tr1::regex_constants::match_default);
	//std::string result;


	//std::regex_replace(std::back_inserter(eqn_form),eqn.begin(),eqn.end(),e,"c",std::regex_constants::match_default);
	//std::regex_replace(std::back_inserter(eqn_form),eqn.begin(),eqn.end(),e,"c",std::tr1::regex_constants::match_default
    //std::cout << result;
	//std::cout << eqn << "\t" << eqn_form <<"\n";
}
float getCorr(vector<float>& output,vector<float>& target,float meanout,float meantarget,int off,float& target_std)
{
	float v1,v2;
	float var_target=0;
	float var_ind = 0;
	float q=0;
	float ndata = float(output.size());
	float corr;
	
	//calculate correlation coefficient
	for (unsigned int c = 0; c<output.size(); c++)
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
	float diff;
	float diffmean=0;

	for (unsigned int c = 0; c<output.size(); c++)
		diffmean += target.at(c+off)-output.at(c);
	diffmean = diffmean/output.size();
	//calculate correlation coefficient
	for (unsigned int c = 0; c<output.size(); c++)
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
	for (unsigned int c = 0; c<target.size(); c++)
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
	for(int m=0;m<eqn.size();m++){
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
				if (isalpha(eqn[m+1])) m++; 
				else pass=0;
			}
			complexity++;
		}
		else
			complexity++;
	}

	return complexity;
}

void Fitness(vector<ind>& pop,params& p,data& d,state& s,FitnessEstimator& FE)
{
	//set up data table for conversion of symbolic variables
	unordered_map <string,float*> datatable;
	vector<float> dattovar(d.label.size());

	for (unsigned int i=0;i<d.label.size(); i++)
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
	//#pragma omp parallel for private(e)
	for(int count = 0; count<pop.size(); count++)
	{
		if(p.sel!=3){
							
			pop.at(count).abserror = 0;
			pop.at(count).abserror_v = 0;
			float meanout=0;
			float meantarget=0;
			float meanout_v=0;
			float meantarget_v=0;
			float target_std=1000;
			float target_std_v=1000;
			//float sumout=0;
			//float meantarget=0; 
	//get equation and equation form
			pop.at(count).eqn = Line2Eqn(pop.at(count).line);
			getEqnForm(pop.at(count).eqn,pop.at(count).eqn_form);
			//if (pop.at(count).eqn.compare("1")==0)
				//cout << "caught one\n";
			// set pointer to dattovar in symbolic functions

	//effective size
			pop.at(count).eff_size=0;
			for(int m=0;m<pop.at(count).line.size();m++){

				if(pop.at(count).line.at(m)->on)
					pop.at(count).eff_size++;
			}
	// Get Complexity
			pop.at(count).complexity= getComplexity(pop.at(count).eqn_form);
		
		// set data table
			for(int m=0;m<pop.at(count).line.size();m++){
				if(pop.at(count).line.at(m)->type=='v')
					{// set pointer to dattovar 
						float* set = datatable.at(static_pointer_cast<n_sym>(pop.at(count).line.at(m))->varname);
						if(set==NULL)
							cout<<"hmm";
						static_pointer_cast<n_sym>(pop.at(count).line.at(m))->setpt(set);
						/*if (static_pointer_cast<n_sym>(pop.at(count).line.at(m))->valpt==NULL)
							cout<<"wth";*/
					}
			}
			//cout << "Equation" << count << ": f=" << pop.at(count).eqn << "\n";
			bool pass=true;
			if(!pop.at(count).eqn.compare("unwriteable")==0){
				vector<float> outstack;
				pop.at(count).output.clear();
				pop.at(count).output_v.clear();
				float SStot=0;
				float SSreg=0;
				float SSres=0;
				float q = 0;
				float var_target = 0;
				float var_ind = 0;
				
	// calculate error 
				if(!p.EstimateFitness){
					for(unsigned int sim=0;sim<d.vals.size();sim++)
					{
						for (unsigned int j=0; j<p.allvars.size();j++)
							dattovar.at(j)= d.vals[sim][j];

						for(int k=0;k<pop.at(count).line.size();k++){
							/*if(pop.at(count).line.at(k)->type=='v'){
								if (static_pointer_cast<n_sym>(pop.at(count).line.at(k))->valpt==NULL)
									cout<<"WTF";
							}*/
							if (pop.at(count).line.at(k)->on){
								pop.at(count).line.at(k)->eval(outstack);
								ptevals++;}
						}

						if(!outstack.empty()){
							if (p.train){
								if(sim<ndata_t){
									pop.at(count).output.push_back(outstack.back());
									pop.at(count).abserror += abs(d.target.at(sim)-pop.at(count).output.at(sim));
									meantarget += d.target.at(sim);
									meanout += pop.at(count).output[sim];

								}
								else
								{
									pop.at(count).output_v.push_back(outstack.back());
									pop.at(count).abserror_v += abs(d.target.at(sim)-pop.at(count).output_v.at(sim-ndata_t));
									meantarget_v += d.target.at(sim);
									meanout_v += pop.at(count).output_v[sim-ndata_t];

								}

							}
							else {
								pop.at(count).output.push_back(outstack.back());
								pop.at(count).abserror += abs(d.target.at(sim)-pop.at(count).output.at(sim));
								meantarget += d.target.at(sim);
								meanout += pop.at(count).output[sim];
							}
						}
						else{
							pass=false;
							break;
							}
						outstack.clear();
					}
					if (pass){
						
						// mean absolute error
						pop.at(count).abserror = pop.at(count).abserror/ndata_t;
						meantarget = meantarget/ndata_t;
						meanout = meanout/ndata_t;
						//calculate correlation coefficient
						pop.at(count).corr = getCorr(pop.at(count).output,d.target,meanout,meantarget,0,target_std);
						pop.at(count).VAF = VAF(pop.at(count).output,d.target,meantarget,0);

						if (p.train)
						{
							q = 0;
							var_target = 0;
							var_ind = 0;
								// mean absolute error
							pop.at(count).abserror_v = pop.at(count).abserror_v/ndata_v;
							meantarget_v = meantarget_v/ndata_v;
							meanout_v = meanout_v/ndata_v;
							//calculate correlation coefficient
							pop.at(count).corr_v = getCorr(pop.at(count).output_v,d.target,meanout_v,meantarget_v,ndata_t,target_std_v);
							pop.at(count).VAF_v = VAF(pop.at(count).output_v,d.target,meantarget_v,ndata_t);
						}

					}
				} // if not estimate fitness 
	
				else{ 
	// use fitness estimator subset of d.vals ===========================================================		
					/*pop.at(count).output.clear();
					pop.at(count).output_v.clear();
					int tmp1= pop.at(count).output.size() ;
					int tmp2= pop.at(count).output_v.size();*/
					for(unsigned int sim=0;sim<FEvals.size();sim++)
					{
						for (unsigned int j=0; j<p.allvars.size();j++)
							dattovar.at(j)= FEvals[sim][j];

						for(int k=0;k<pop.at(count).line.size();k++){
							/*if(pop.at(count).line.at(k)->type=='v'){
								if (static_pointer_cast<n_sym>(pop.at(count).line.at(k))->valpt==NULL)
									cout<<"WTF";
							}*/
							if (pop.at(count).line.at(k)->on){
								pop.at(count).line.at(k)->eval(outstack);
								ptevals++;}
						}

						if(!outstack.empty()){
							if (p.train){
								if(sim<ndata_t){
									if (pop.at(count).output.size() != sim)
										cout << "hold up\n";

									pop.at(count).output.push_back(outstack.back());
									pop.at(count).abserror += abs(FEtarget.at(sim)-pop.at(count).output.at(sim));
									meantarget += FEtarget.at(sim);
									meanout += pop.at(count).output[sim];
									}
									else
									{
										pop.at(count).output_v.push_back(outstack.back());
										pop.at(count).abserror_v += abs(FEtarget.at(sim)-pop.at(count).output_v.at(sim-ndata_t));
										meantarget_v += FEtarget.at(sim);
										meanout_v += pop.at(count).output_v[sim-ndata_t];
								}
							} // if p.train
							else {
								pop.at(count).output.push_back(outstack.back());
								pop.at(count).abserror += abs(FEtarget.at(sim)-pop.at(count).output.at(sim));
								meantarget += FEtarget.at(sim);
								meanout += pop.at(count).output[sim];
							}
						}
						else{
							pass=false;
							break;
							}
						outstack.clear();
					}
					if (pass){
						// mean absolute error
						pop.at(count).abserror = pop.at(count).abserror/ndata_t;
						meantarget = meantarget/ndata_t;
						meanout = meanout/ndata_t;
						//calculate correlation coefficient
						pop.at(count).corr = getCorr(pop.at(count).output,FEtarget,meanout,meantarget,0,target_std);
						pop.at(count).VAF = VAF(pop.at(count).output,FEtarget,meantarget,0);

						if (p.train)
						{
							q = 0;
							var_target = 0;
							var_ind = 0;
								// mean absolute error
							pop.at(count).abserror_v = pop.at(count).abserror_v/ndata_v;
							meantarget_v = meantarget_v/ndata_v;
							meanout_v = meanout_v/ndata_v;
							//calculate correlation coefficient
							
							pop.at(count).corr_v = getCorr(pop.at(count).output_v,FEtarget,meanout_v,meantarget_v,ndata_t,target_std_v);
							pop.at(count).VAF_v = VAF(pop.at(count).output_v,FEtarget,meantarget_v,0);
						}

					} // if pass
				} // if estimate fitness
			} // if not unwriteable equation
			else{ // bad equation, assign maximum fitness
				pop.at(count).abserror=p.max_fit;
				pop.at(count).corr = p.min_fit;
				pop.at(count).VAF = p.min_fit;
				if (p.train){
					pop.at(count).abserror_v=p.max_fit;
					pop.at(count).corr_v = p.min_fit;
					pop.at(count).VAF_v = p.min_fit;
				}
			}
			
			if (pop.at(count).eqn.compare("1")==0 && pop.at(count).corr > 0.0001)
				cout << "caught\n";

			if (!pass)
				pop.at(count).corr = 0;
						

			if(pop.at(count).corr < p.min_fit)
				pop.at(count).corr=p.min_fit;
			if(pop.at(count).VAF < p.min_fit)
				pop.at(count).VAF=p.min_fit;

		    if(pop.at(count).output.empty())
				pop.at(count).fitness=p.max_fit;
			/*else if (*std::max_element(pop.at(count).output.begin(),pop.at(count).output.end())==*std::min_element(pop.at(count).output.begin(),pop.at(count).output.end()))
				pop.at(count).fitness=p.max_fit;*/
			else if ( boost::math::isnan(pop.at(count).abserror) || boost::math::isinf(pop.at(count).abserror) || boost::math::isnan(pop.at(count).corr) || boost::math::isinf(pop.at(count).corr))
				pop.at(count).fitness=p.max_fit;
			else{
				if (p.fit_type==1)
					pop.at(count).fitness = pop.at(count).abserror;
				else if (p.fit_type==2)
					pop.at(count).fitness = 1-pop.at(count).corr;
				else if (p.fit_type==3)
					pop.at(count).fitness = pop.at(count).abserror/pop.at(count).corr;
				else if (p.fit_type==4)
					pop.at(count).fitness = 1-pop.at(count).VAF;
				if (p.norm_error)
					pop.at(count).fitness = pop.at(count).fitness/target_std;
			}

		
			if(pop.at(count).fitness>p.max_fit)
				pop.at(count).fitness=p.max_fit;
			else if(pop.at(count).fitness<p.min_fit)
				(pop.at(count).fitness=p.min_fit);

			if(p.train){ //assign validation fitness
				if (!pass)
					pop.at(count).corr_v = 0;
						

				if(pop.at(count).corr_v < p.min_fit)
					pop.at(count).corr_v=p.min_fit;
				if(pop.at(count).VAF_v < p.min_fit)
					pop.at(count).VAF_v=p.min_fit;

				if(pop.at(count).output_v.empty())
					pop.at(count).fitness_v=p.max_fit;
				/*else if (*std::max_element(pop.at(count).output_v.begin(),pop.at(count).output_v.end())==*std::min_element(pop.at(count).output_v.begin(),pop.at(count).output_v.end()))
					pop.at(count).fitness_v=p.max_fit;*/
				else if ( boost::math::isnan(pop.at(count).abserror_v) || boost::math::isinf(pop.at(count).abserror_v) || boost::math::isnan(pop.at(count).corr_v) || boost::math::isinf(pop.at(count).corr_v))
					pop.at(count).fitness_v=p.max_fit;
				else{
					if (p.fit_type==1)
						pop.at(count).fitness_v = pop.at(count).abserror_v;
					else if (p.fit_type==2)
						pop.at(count).fitness_v = 1-pop.at(count).corr_v;
					else if (p.fit_type==3)
						pop.at(count).fitness_v = pop.at(count).abserror_v/pop.at(count).corr_v;
					else if (p.fit_type==4)
						pop.at(count).fitness = 1-pop.at(count).VAF;
					if (p.norm_error)
						pop.at(count).fitness_v = pop.at(count).fitness_v/target_std_v;

				}

		
				if(pop.at(count).fitness_v>p.max_fit)
					pop.at(count).fitness_v=p.max_fit;
				else if(pop.at(count).fitness_v<p.min_fit)
					(pop.at(count).fitness_v=p.min_fit);
			}
			else{ // if not training, assign copy of regular fitness to the validation fitness variables
				pop.at(count).corr_v=pop.at(count).corr;
				pop.at(count).abserror_v=pop.at(count).abserror;
				pop.at(count).fitness_v=pop.at(count).fitness;
			} 
		
		} // if p.sel!=3
		else //LEXICASE FITNESS ==================================================================
		{
		//get equation and equation form
			pop.at(count).eqn = Line2Eqn(pop.at(count).line);
			getEqnForm(pop.at(count).eqn,pop.at(count).eqn_form);
		//Get effective size and complexity
			pop.at(count).eff_size=0;
			for(int m=0;m<pop.at(count).line.size();m++){

				if(pop.at(count).line.at(m)->on)
					pop.at(count).eff_size++;
			}
			// Get Complexity
			pop.at(count).complexity= getComplexity(pop.at(count).eqn_form);
			
			// set data pointers
			for(int m=0;m<pop.at(count).line.size();m++){
				if(pop.at(count).line.at(m)->type=='v')
					{// set pointer to dattovar 
						float* set = datatable.at(static_pointer_cast<n_sym>(pop.at(count).line.at(m))->varname);
						if(set==NULL)
							cout<<"hmm";
						static_pointer_cast<n_sym>(pop.at(count).line.at(m))->setpt(set);
						/*if (static_pointer_cast<n_sym>(pop.at(count).line.at(m))->valpt==NULL)
							cout<<"wth";*/
					}
			}

			pop.at(count).abserror = 0;
			pop.at(count).abserror_v = 0;
			float meanout=0;
			float meanoutlex;
			float meantarget=0;
			float meantargetlex;
			float meanout_v=0;
			float meantarget_v=0;
			float target_std_lex, target_std_lex_v;
			// set pointer to dattovar in symbolic functions

			//cout << "Equation" << count << ": f=" << pop.at(count).eqn << "\n";
			pop.at(count).fitlex.resize(p.numcases);
			bool pass=true;
			if(!pop.at(count).eqn.compare("unwriteable")==0){
				
				vector<float> outstack;
				pop.at(count).output.clear();
				pop.at(count).output_v.clear();
				
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
				for (int lex=0; lex<d.lexvals.size(); lex++)
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
					for(unsigned int sim=0;sim<d.lexvals[lex].size();sim++)
					{
						//outlex_v.push_back(vector<float>());

						for (unsigned int j=0; j<p.allvars.size();j++)
							dattovar.at(j)= d.lexvals[lex][sim][j];

						for(int k=0;k<pop.at(count).line.size();k++){
							if (pop.at(count).line.at(k)->on){
								pop.at(count).line.at(k)->eval(outstack);
								ptevals++;}
						}

						if(!outstack.empty()){
							if (p.train){
								if(sim<ndata_t){
									pop.at(count).output.push_back(outstack.back());
									tarlex.push_back(d.targetlex[lex][sim]);
									outlex[lex].push_back(outstack.back());

									pop.at(count).abserror += abs(d.targetlex[lex].at(sim)-pop.at(count).output.back());
									errorlex[lex]+=abs(d.targetlex[lex].at(sim)-pop.at(count).output.back());

									meantarget += d.targetlex[lex].at(sim);
									meantargetlex += d.targetlex[lex].at(sim);
									meanout += pop.at(count).output.back();
									meanoutlex+= pop.at(count).output.back();

								}
								else
								{
									pop.at(count).output_v.push_back(outstack.back());
									tarlex_v.push_back(d.targetlex[lex][sim]);
									//outlex_v[lex].push_back(outstack.back());
									pop.at(count).abserror_v += abs(d.targetlex[lex].at(sim)-pop.at(count).output_v.back());
									//errorlex_v[lex]+=abs(d.targetlex[lex]lex[lex].at(sim)-pop.at(count).output.at(sim-ndata_t));
									meantarget_v += d.targetlex[lex][sim];
									//meantargetlex += d.targetlex[lex]lex[lex].at(sim);
									meanout_v += pop.at(count).output_v.back();
									//meanoutlex_v += pop.at(count).output_v[sim-ndata_t];

								}
							}
							else {
								pop.at(count).output.push_back(outstack.back());
								tarlex.push_back(d.targetlex[lex][sim]);
								outlex[lex].push_back(outstack.back());
								pop.at(count).abserror += abs(d.targetlex[lex].at(sim)-pop.at(count).output.back());
								errorlex[lex]+=abs(d.targetlex[lex].at(sim)-pop.at(count).output.back());
								meantarget += d.targetlex[lex].at(sim);
								meantargetlex += d.targetlex[lex].at(sim);
								meanout += pop.at(count).output.back();
								meanoutlex+= pop.at(count).output.back();
							}
						}
						else{
							pass=false;
							break;
							}
						outstack.clear();
					}//for(unsigned int sim=0;sim<d.lexvals[lex].size();sim++)

					if (pass){
					
						// mean absolute error
						errorlex[lex] = errorlex[lex]/outlex[lex].size();
						meantargetlex = meantargetlex/outlex[lex].size();
						meanoutlex= meanoutlex/outlex[lex].size();

						//calculate correlation coefficient
						string tmp = pop.at(count).eqn;
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

					if(pop.at(count).output.empty())
						pop.at(count).fitlex[lex]=p.max_fit;
					/*else if (*std::max_element(pop.at(count).output.begin(),pop.at(count).output.end())==*std::min_element(pop.at(count).output.begin(),pop.at(count).output.end()))
						pop.at(count).fitness=p.max_fit;*/
					else if ( boost::math::isnan(errorlex[lex]) || boost::math::isinf(errorlex[lex]) || boost::math::isnan(corrlex[lex]) || boost::math::isinf(corrlex[lex]))
						pop.at(count).fitlex[lex]=p.max_fit;
					else{
						if (p.fit_type==1)
							pop.at(count).fitlex[lex] = errorlex[lex];
						else if (p.fit_type==2)
							pop.at(count).fitlex[lex] = 1-corrlex[lex];
						else if (p.fit_type==3)
							pop.at(count).fitlex[lex] = errorlex[lex]/corrlex[lex];
						else if (p.fit_type==4)
							pop.at(count).fitlex[lex] = 1-VAFlex[lex];
						if (p.norm_error)
							pop.at(count).fitlex[lex] = pop.at(count).fitlex[lex]/target_std_lex;
					}

		
					if(pop.at(count).fitlex[lex]>p.max_fit)
						pop.at(count).fitlex[lex]=p.max_fit;
					else if(pop.at(count).fitlex[lex]<p.min_fit)
						(pop.at(count).fitlex[lex]=p.min_fit);

				} //for (int lex=0; lex<d.lexvals.size(); lex++)
	// fill in overall fitness information
				if (pass){
					
						// mean absolute error
						pop.at(count).abserror = pop.at(count).abserror/pop.at(count).output.size();
						meantarget = meantarget/pop.at(count).output.size();
						meanout = meanout/pop.at(count).output.size();
						//meanoutlex= meanoutlex/pop.at(count).output.size();

						//calculate correlation coefficient
						pop.at(count).corr = getCorr(pop.at(count).output,tarlex,meanout,meantarget,0,target_std_lex);
						pop.at(count).VAF = VAF(pop.at(count).output,tarlex,meantarget,0);
					
						if (p.train)
						{
							// mean absolute error
							pop.at(count).abserror_v = pop.at(count).abserror_v/pop.at(count).output_v.size();
							meantarget_v = meantarget_v/tarlex_v.size();
							meanout_v = meanout_v/pop.at(count).output_v.size();
							//calculate correlation coefficient
							pop.at(count).corr_v = getCorr(pop.at(count).output_v,tarlex_v,meanout_v,meantarget_v,0,target_std_lex_v);
							pop.at(count).VAF_v = VAF(pop.at(count).output_v,tarlex_v,meantarget_v,0);
						}

					}
					else{
						pop.at(count).abserror=p.max_fit;
						pop.at(count).corr = p.min_fit;
						pop.at(count).VAF = p.min_fit;

						if (p.train){
							pop.at(count).abserror_v=p.max_fit;
							pop.at(count).corr_v = p.min_fit;
							pop.at(count).VAF_v = p.min_fit;
						}
					}						

				if(pop.at(count).corr < p.min_fit)
					pop.at(count).corr=p.min_fit;
				if(pop.at(count).VAF < p.min_fit)
					pop.at(count).VAF = p.min_fit;

				if(pop.at(count).output.empty())
					pop.at(count).fitness=p.max_fit;
				/*else if (*std::max_element(pop.at(count).output.begin(),pop.at(count).output.end())==*std::min_element(pop.at(count).output.begin(),pop.at(count).output.end()))
					pop.at(count).fitness=p.max_fit;*/
				else if ( boost::math::isnan(pop.at(count).abserror) || boost::math::isinf(pop.at(count).abserror) || boost::math::isnan(pop.at(count).corr) || boost::math::isinf(pop.at(count).corr))
					pop.at(count).fitness=p.max_fit;
				else{
					if (p.fit_type==1)
						pop.at(count).fitness = pop.at(count).abserror;
					else if (p.fit_type==2)
						pop.at(count).fitness = 1-pop.at(count).corr;
					else if (p.fit_type==3)
						pop.at(count).fitness = pop.at(count).abserror/pop.at(count).corr;
					else if (p.fit_type==4)
						pop.at(count).fitness = 1-pop.at(count).VAF;
					if (p.norm_error)
						pop.at(count).fitness = pop.at(count).fitness/target_std_lex;

				}

		
				if(pop.at(count).fitness>p.max_fit)
					pop.at(count).fitness=p.max_fit;
				else if(pop.at(count).fitness<p.min_fit)
					(pop.at(count).fitness=p.min_fit);

				if(p.train){
					if (!pass)
						pop.at(count).corr_v = 0;
						

					if(pop.at(count).corr_v < p.min_fit)
						pop.at(count).corr_v=p.min_fit;
					if(pop.at(count).VAF_v < p.min_fit)
						pop.at(count).VAF_v=p.min_fit;

					if(pop.at(count).output_v.empty())
						pop.at(count).fitness_v=p.max_fit;
					/*else if (*std::max_element(pop.at(count).output_v.begin(),pop.at(count).output_v.end())==*std::min_element(pop.at(count).output_v.begin(),pop.at(count).output_v.end()))
						pop.at(count).fitness_v=p.max_fit;*/
					else if ( boost::math::isnan(pop.at(count).abserror_v) || boost::math::isinf(pop.at(count).abserror_v) || boost::math::isnan(pop.at(count).corr_v) || boost::math::isinf(pop.at(count).corr_v))
						pop.at(count).fitness_v=p.max_fit;
					else{
						if (p.fit_type==1)
							pop.at(count).fitness_v = pop.at(count).abserror_v;
						else if (p.fit_type==2)
							pop.at(count).fitness_v = 1-pop.at(count).corr_v;
						else if (p.fit_type==3)
							pop.at(count).fitness_v = pop.at(count).abserror_v/pop.at(count).corr_v;
						else if (p.fit_type==4)
							pop.at(count).fitness = 1-pop.at(count).VAF;
						if (p.norm_error)
							pop.at(count).fitness_v = pop.at(count).fitness_v/target_std_lex_v;
					}

		
					if(pop.at(count).fitness_v>p.max_fit)
						pop.at(count).fitness_v=p.max_fit;
					else if(pop.at(count).fitness_v<p.min_fit)
						(pop.at(count).fitness_v=p.min_fit);
				}
				else{ // if not training, assign copy of regular fitness to validation variales
					pop.at(count).corr_v=pop.at(count).corr;
					pop.at(count).abserror_v=pop.at(count).abserror;
					pop.at(count).fitness_v=pop.at(count).fitness;
				}
			}//if(!pop.at(count).eqn.compare("unwriteable")==0)
	else
		{
			for (int i=0;i<p.numcases;i++)
				pop.at(count).fitlex[i]=p.max_fit;
			pop.at(count).fitness = p.max_fit;
			pop.at(count).fitness_v = p.max_fit;

		}
	}//LEXICASE FITNESS
	}//for(int count = 0; count<pop.size(); count++)
	s.numevals[omp_get_thread_num()]=s.numevals[omp_get_thread_num()]+pop.size();
	s.ptevals[omp_get_thread_num()]=s.ptevals[omp_get_thread_num()]+ptevals;
	//cout << "\nFitness Time: ";
}
