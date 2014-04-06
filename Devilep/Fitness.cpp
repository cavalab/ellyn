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

void Fitness(vector<ind>& pop,params& p,data& d,state& s)
{
	//boost::progress_timer timer;
	//int count;
	/*evaluator e;
	e.init(p,d);*/
	//int wtf = pop.size();
	
	//#pragma omp parallel for private(e)
	for(int count = 0; count<pop.size(); count++)
	{
			try 
			{
			pop.at(count).abserror = 0;
			float meanout=0;
			float meantarget=0;
			//float sumout=0;
			//float meantarget=0; 
			//get equation and equation form
			pop.at(count).eqn = Line2Eqn(pop.at(count).line);
			getEqnForm(pop.at(count).eqn,pop.at(count).eqn_form);
			//GET EFFECTIVE SIZE
			pop.at(count).eff_size=0;
			for(int m=0;m<pop.at(count).line.size();m++){
				if(pop.at(count).line.at(m)->on)
					pop.at(count).eff_size++;
			}
			//cout << "Equation" << count << ": f=" << pop.at(count).eqn << "\n";
			if(!pop.at(count).eqn.compare("unwriteable")==0){
				vector<float> outstack;

				for(unsigned int sim=0;sim<d.vals.size();sim++)
				{
					for (unsigned int j=0; j<p.allvars.size();j++)
						d.dattovar.at(j)= d.vals[sim][j];

					for(int k=0;k<pop.at(count).line.size();k++){
						if (pop.at(count).line.at(k)->on)
							pop.at(count).line.at(k)->eval(outstack);
					}

					if(!outstack.empty()){
						pop.at(count).output.push_back(outstack.front());
						pop.at(count).abserror += abs(d.target.at(sim)-pop.at(count).output.at(sim));
						meantarget += d.target.at(sim)/d.vals.size();
						meanout += pop.at(count).output[sim]/d.vals.size();
					}

					outstack.clear();
				}
			}
			else
				pop.at(count).abserror=p.max_fit;

			pop.at(count).corr = 0;
						

			if(pop.at(count).corr < p.min_fit)
				pop.at(count).corr=p.min_fit;

		    if(pop.at(count).output.empty())
				pop.at(count).fitness=p.max_fit;
			else if (*std::max_element(pop.at(count).output.begin(),pop.at(count).output.end())==*std::min_element(pop.at(count).output.begin(),pop.at(count).output.end()))
				pop.at(count).fitness=p.max_fit;
			else if ( boost::math::isnan(pop.at(count).abserror) || boost::math::isinf(pop.at(count).abserror) || boost::math::isnan(pop.at(count).corr) || boost::math::isinf(pop.at(count).corr))
				pop.at(count).fitness=p.max_fit;
			else
				pop.at(count).fitness = pop.at(count).abserror;
			//else
				//pop.at(count).fitness = pop.at(count).abserror/pop.at(count).corr;
		
			if(pop.at(count).fitness>p.max_fit)
				pop.at(count).fitness=p.max_fit;
			else if(pop.at(count).fitness<p.min_fit)
				(pop.at(count).fitness=p.min_fit);

			}
			catch(CException ex)
			{
				ex.Report();
				pop.at(count).fitness=p.max_fit;
			}
		/*}
		else
			pop.at(count).fitness=p.max_fit;*/
		
		//cout << count << "\t" << pop.at(count).fitness << endl; 
		
	}
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
