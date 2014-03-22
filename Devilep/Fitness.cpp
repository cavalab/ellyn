#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "data.h"
#include <cstdlib>
#include <math.h>
#include "evaluator.h"
#include "state.h"


void Fitness(vector<ind>& pop,params& p,data& d,state& s)
{
	//boost::progress_timer timer;
	//int count;
	evaluator e;
	e.init(p,d);
	int wtf = pop.size();

	//#pragma omp parallel for private(e)
	for(unsigned int count = 0; count<pop.size(); count++)
	{
		///*if(count==0)
		//	cout<<"number of threads: " << omp_get_num_threads() << endl;*/

		e.parser.compile(pop.at(count).eqn,e.expression);
		pop.at(count).abserror = 0;
		float meanout=0;
		float meantarget=0;
		//float sumout=0;
		//float meantarget=0;
		for(unsigned int sim=0;sim<d.vals.size();sim++)
		{
			for (unsigned int j=0; j<p.allvars.size();j++)
				d.dattovar.at(j)= d.vals[sim][j];
			
			pop.at(count).output.push_back(e.expression.value());
			
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
		else
			//pop.at(count).fitness = pop.at(count).abserror/pop.at(count).corr;
		

		if(pop.at(count).fitness>p.max_fit)
			pop.at(count).fitness=p.max_fit;
		else if(pop.at(count).fitness<p.min_fit)
			(pop.at(count).fitness=p.min_fit);
		else{}
		//cout << count << "\t" << pop.at(count).fitness << endl; 
	}
	
	s.numevals[omp_get_thread_num()]=s.numevals[omp_get_thread_num()]+pop.size();
		
}

	//for(count = 0; count<pop.size(); count++)
	//{
	//	d.parser.compile(pop.at(count).eqn,pop.at(count).expression);
	// /*   if (!d.parser.compile(pop.at(count).eqn,pop.at(count).expression))
	//    {
	//	   printf("Error: %s\tExpression: %s\n",
	// 			 d.parser.error().c_str(),
 //				 pop.at(count).eqn.c_str());
 // 
 //		  for (std::size_t i = 0; i < d.parser.error_count(); ++i)
	//	  {
	//		 error_t error = d.parser.get_error(i);
	//		 printf("Error: %02d Position: %02d Type: [%s] Msg: %s Expr: %s\n",
	//				static_cast<int>(i),
	//				static_cast<int>(error.token.position),
	//				exprtk::parser_error::to_str(error.mode).c_str(),
	//				error.diagnostic.c_str(),
	//				pop.at(count).eqn.c_str());
	//	  }
	//   }*/
	//}