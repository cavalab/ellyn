#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "data.h"
#include <cstdlib>
#include <math.h>
//#include "evaluator.h"
#include "state.h"
#include "EvalEqnStr.h"
#include <unordered_map>
#include "Line2Eqn.h"
//#include "exprtk.hpp"

class CException
{
public:
	char* message;
	CException( char* m ) { message = m; };
	void Report(){cout << "error calculating fitness\n";};

};
void getEqnForm(std::string& eqn,std::string& eqn_form);
float getCorr(vector<float>& output,vector<float>& target,float meanout,float meantarget,int off);

void Fitness(vector<ind>& pop,params& p,data& d,state& s)
{
	//boost::progress_timer timer;
	//int count;
	/*evaluator e;
	e.init(p,d);*/
	//int wtf = pop.size();

	//set up data table for conversion of symbolic variables
	unordered_map <string,float*> datatable;
	vector<float> dattovar(d.label.size());

	for (unsigned int i=0;i<d.label.size(); i++)
			datatable.insert(pair<string,float*>(d.label[i],&dattovar[i]));

	int ndata_t,ndata_v; // training and validation data sizes
	if (p.train){
		ndata_t = d.vals.size()/2;
		ndata_v = d.vals.size()-ndata_t;
	}
	else{
		ndata_t = d.vals.size();
		ndata_v=0;
	}


	//#pragma omp parallel for private(e)
	for(int count = 0; count<pop.size(); count++)
	{
		/*for (unsigned int y = 0;y<pop.at(count).line.size();y++)
		{
			if (pop.at(count).line.at(y)->type == 'v')
			{}

		}*/
		if(!p.sel==3){
							
			pop.at(count).abserror = 0;
			pop.at(count).abserror_v = 0;
			float meanout=0;
			float meantarget=0;
			float meanout_v=0;
			float meantarget_v=0;
			//float sumout=0;
			//float meantarget=0; 
			//get equation and equation form
			pop.at(count).eqn = Line2Eqn(pop.at(count).line);
			getEqnForm(pop.at(count).eqn,pop.at(count).eqn_form);
			
			// set pointer to dattovar in symbolic functions

			//GET EFFECTIVE SIZE
			pop.at(count).eff_size=0;
			for(int m=0;m<pop.at(count).line.size();m++){

				if(pop.at(count).line.at(m)->on)
					pop.at(count).eff_size++;
			}
			// Get Complexity
			pop.at(count).complexity=0;
			for(int m=0;m<pop.at(count).line.size();m++){
				if(pop.at(count).line.at(m)->on)
				{
					pop.at(count).complexity++;
					if (pop.at(count).line.at(m)->type=='/')
						pop.at(count).complexity++;
					else if (pop.at(count).line.at(m)->type=='s' || pop.at(count).line.at(m)->type=='c')
						pop.at(count).eff_size=pop.at(count).eff_size+2;
					else if (pop.at(count).line.at(m)->type=='e' || pop.at(count).line.at(m)->type=='l')
						pop.at(count).eff_size=pop.at(count).eff_size+3;
				}
			}
			
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
				for(unsigned int sim=0;sim<d.vals.size();sim++)
				{
					for (unsigned int j=0; j<p.allvars.size();j++)
						dattovar.at(j)= d.vals[sim][j];

					for(int k=0;k<pop.at(count).line.size();k++){
						/*if(pop.at(count).line.at(k)->type=='v'){
							if (static_pointer_cast<n_sym>(pop.at(count).line.at(k))->valpt==NULL)
								cout<<"WTF";
						}*/
						if (pop.at(count).line.at(k)->on)
							pop.at(count).line.at(k)->eval(outstack);
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
					pop.at(count).corr = getCorr(pop.at(count).output,d.target,meanout,meantarget,0);
					

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
						pop.at(count).corr_v = getCorr(pop.at(count).output_v,d.target,meanout_v,meantarget_v,ndata_t);
					}

				}
			}
			else{
				pop.at(count).abserror=p.max_fit;
				pop.at(count).corr = p.min_fit;
				if (p.train){
					pop.at(count).abserror_v=p.max_fit;
					pop.at(count).corr_v = p.min_fit;
				}
				}
			
			if (!pass)
				pop.at(count).corr = 0;
						

			if(pop.at(count).corr < p.min_fit)
				pop.at(count).corr=p.min_fit;

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
				}

		
				if(pop.at(count).fitness_v>p.max_fit)
					pop.at(count).fitness_v=p.max_fit;
				else if(pop.at(count).fitness_v<p.min_fit)
					(pop.at(count).fitness_v=p.min_fit);
			}
			else{
				pop.at(count).corr_v=pop.at(count).corr;
				pop.at(count).abserror_v=pop.at(count).abserror;
				pop.at(count).fitness_v=pop.at(count).fitness;
			} 
		
	}
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
			pop.at(count).complexity=0;
			for(int m=0;m<pop.at(count).line.size();m++){
				if(pop.at(count).line.at(m)->on)
				{
					pop.at(count).complexity++;
					if (pop.at(count).line.at(m)->type=='/')
						pop.at(count).complexity++;
					else if (pop.at(count).line.at(m)->type=='s' || pop.at(count).line.at(m)->type=='c')
						pop.at(count).complexity=pop.at(count).complexity+2;
					else if (pop.at(count).line.at(m)->type=='e' || pop.at(count).line.at(m)->type=='l')
						pop.at(count).complexity=pop.at(count).complexity+3;
				}
			}
			
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
				
			// set pointer to dattovar in symbolic functions

			//cout << "Equation" << count << ": f=" << pop.at(count).eqn << "\n";
			pop.at(count).fitlex.resize(p.numcases);
			bool pass=true;
			if(!pop.at(count).eqn.compare("unwriteable")==0){
				
				vector<float> outstack;
				pop.at(count).output.clear();
				pop.at(count).output_v.clear();
				
				vector<vector<float>> outlex;
				vector<float> errorlex;
				vector<float> corrlex;

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
							if (pop.at(count).line.at(k)->on)
								pop.at(count).line.at(k)->eval(outstack);
						}

						if(!outstack.empty()){
							if (p.train){
								if(sim<ndata_t){
									pop.at(count).output.push_back(outstack.back());
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
									//outlex_v[lex].push_back(outstack.back());
									pop.at(count).abserror_v += abs(d.targetlex[lex].at(sim)-pop.at(count).output_v.back());
									//errorlex_v[lex]+=abs(d.targetlex[lex]lex[lex].at(sim)-pop.at(count).output.at(sim-ndata_t));
									//meantarget_v += d.targetlex[lex].at(sim);
									//meantargetlex += d.targetlex[lex]lex[lex].at(sim);
									meanout_v += pop.at(count).output_v.back();
									//meanoutlex_v += pop.at(count).output_v[sim-ndata_t];

								}
							}
							else {
								pop.at(count).output.push_back(outstack.back());
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
						corrlex[lex] = getCorr(outlex[lex],d.targetlex[lex],meanoutlex,meantargetlex,0);

					}
					else{
						errorlex[lex]=p.max_fit;
						corrlex[lex] = p.min_fit;
					}
			
					if(corrlex[lex] < p.min_fit)
						corrlex[lex]=p.min_fit;

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
						pop.at(count).corr = getCorr(pop.at(count).output,d.target,meanout,meantarget,0);
					
						if (p.train)
						{
							// mean absolute error
							pop.at(count).abserror_v = pop.at(count).abserror_v/pop.at(count).output_v.size();
							meantarget_v = meantarget_v/pop.at(count).output_v.size();
							meanout_v = meanout_v/pop.at(count).output_v.size();
							//calculate correlation coefficient
							pop.at(count).corr_v = getCorr(pop.at(count).output_v,d.target,meanout_v,meantarget_v,pop.at(count).output.size());
						}

					}
					else{
						pop.at(count).abserror=p.max_fit;
						pop.at(count).corr = p.min_fit;
						if (p.train){
							pop.at(count).abserror_v=p.max_fit;
							pop.at(count).corr_v = p.min_fit;
						}
					}						

				if(pop.at(count).corr < p.min_fit)
					pop.at(count).corr=p.min_fit;

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
					}

		
					if(pop.at(count).fitness_v>p.max_fit)
						pop.at(count).fitness_v=p.max_fit;
					else if(pop.at(count).fitness_v<p.min_fit)
						(pop.at(count).fitness_v=p.min_fit);
				}
				else{
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
	//cout << "\nFitness Time: ";
}

void getEqnForm(std::string& eqn,std::string& eqn_form)
{
//replace numbers with the letter c
	//boost::regex re("\d|[\d\.\d]");
	//(([0-9]+)(\.)([0-9]+))|([0-9]+)
	boost::regex re("(([0-9]+)(\.)([0-9]+))|([0-9]+)");
	//(\d+\.\d+)|(\d+)
	eqn_form=boost::regex_replace(eqn,re,"c");
	//std::cout << eqn << "\t" << eqn_form <<"\n";
}
float getCorr(vector<float>& output,vector<float>& target,float meanout,float meantarget,int off)
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
	if(var_target==0 || var_ind==0)
		corr = 0;
	else
		corr = pow(q,2)/(var_target*var_ind);

	return corr;
}
