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
//#include "exprtk.hpp"

class CException
{
public:
	char* message;
	CException( char* m ) { message = m; };
	void Report(){cout << "error calculating fitness\n";};

};

void Fitness(vector<ind>& pop,params& p,data& d,state& s)
{
	//boost::progress_timer timer;
	//int count;
	/*evaluator e;
	e.init(p,d);*/
	//int wtf = pop.size();
	vector<float> dattovar(p.allvars.size());
	unordered_map <string,float*> datatable;

	for (unsigned int i=0;i<d.label.size(); i++)
			datatable.insert(pair<string,float*>(d.label[i],&dattovar[i]));

	//#pragma omp parallel for private(e)
	for(int count = 0; count<pop.size(); count++)
	{
		///*if(count==0)
		//	cout<<"number of threads: " << omp_get_num_threads() << endl;*/
		//e.parser.compile(pop.at(count).eqn,e.expression);
		//bool pass=1;
		//if (!e.parser.compile(pop.at(count).eqn,e.expression))
		//{
		//	  // A compilation error has occured. Attempt to
		//	  // print all errors to the stdout.
		//	 // typedef exprtk::parser_error::type error_t;
		//	 s.out << "Expression " + pop.at(count).eqn + " did not compile \n";
		//	  //s.out << "Error: " +  string(e.parser.error().c_str()) + "\tExpression: " + pop.at(count).eqn + "\n";
		//	  //

		//	  //for (std::size_t i = 0; i < e.parser.error_count(); ++i)
		//	  //{
		//		 //// Include the specific nature of each error
		//		 //// and its position in the expression string.
		//		 // 
		//		 //error_t error = e.parser.get_error(i);

		//		 //s.out << "Error: " + to_string(static_cast<long long>(i)) + "\nPosition:  " + to_string(static_cast<long long>(error.token.position)) + "\nType: " + exprtk::parser_error::to_str(error.mode).c_str() + "\nMessage: " + error.diagnostic.c_str() +	"\nExpression: " + pop.at(count).eqn + "\n";
		//	  //}
		//	  pass=false;
		//   }
		/*if( pass )
		{*/
			try 
			{
			pop.at(count).abserror = 0;
			float meanout=0;
			float meantarget=0;
			//float sumout=0;
			//float meantarget=0; 
			
				for(unsigned int sim=0;sim<d.vals.size();sim++)
				{
					for (unsigned int j=0; j<p.allvars.size();j++)
						dattovar.at(j)= d.vals[sim][j];
			
				
						//pop.at(count).output.push_back(e.expression.value());
						pop.at(count).output.push_back(EvalEqnStr(pop.at(count).eqn,datatable));
				
			
					pop.at(count).abserror += abs(d.target.at(sim)-pop.at(count).output.at(sim));
					meantarget += d.target.at(sim)/d.vals.size();
					meanout += pop.at(count).output[sim]/d.vals.size();
				}
		
			////runge kutta
			//float x0=0;
			//float delt=.01;
			//float y1,y2,y3,y4;
			//for(unsigned int sim=0;sim<d.vals.size();sim++)
			//{
			//	/*y1 = [0 t(n); 0 0]*x(:,n)*delt;
			//	y2 = [0 t(n); 0 0]*(x(:,n)+.5*y1)*delt;
			//	y3 = [0 t(n); 0 0]*(x(:,n)+.5*y2)*delt;
			//	y4 = [0 t(n); 0 0]*(x(:,n)+y3)*delt;
	  //  
			//	x(:,n+1) = x(:,n)+1/6*(y1+2*y2+2*y3+y4);*/

			///*
			//}
			//correlation coefficient
			pop.at(count).corr = 0;
			/*for(unsigned int sim=0;sim<pop.at(count).output.size();sim++)
			{
				pop.at(count).corr += ((pop[count].output[sim]-meanout)*(d.target[sim]-meantarget))/pop.at(count).output.size();
			}
			pop.at(count).corr=pop.at(count).corr*pop.at(count).corr;*/

			if(pop.at(count).corr < p.min_fit)
				pop.at(count).corr=p.min_fit;
		
			pop.at(count).fitness = pop.at(count).abserror;

			if (*std::max_element(pop.at(count).output.begin(),pop.at(count).output.end())==*std::min_element(pop.at(count).output.begin(),pop.at(count).output.end()))
				pop.at(count).fitness=p.max_fit;

			if ( boost::math::isnan(pop.at(count).abserror) || boost::math::isinf(pop.at(count).abserror) || boost::math::isnan(pop.at(count).corr) || boost::math::isinf(pop.at(count).corr))
				pop.at(count).fitness=p.max_fit;
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
	dattovar.clear();
	datatable.clear();
	s.numevals[omp_get_thread_num()]=s.numevals[omp_get_thread_num()]+pop.size();
	//cout << "\nFitness Time: ";
		
}
