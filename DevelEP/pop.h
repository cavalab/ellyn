// header file for ind struct
#pragma once
#ifndef POP_H
#define POP_H

//#include "params.h"
//#include "data.h"

//#include "RPN_class.h"
#include "op_node.h"
#include "rnd.h"
#include "strdist.h"
//#include "general_fns.h"
//#include "pareto.h"

struct ind {
	/* 
	===================================================================
	BE SURE TO ADD ANY NEW VARIABLES TO THE SWAP FUNCTION FOR COPYING!!
	===================================================================
	*/
	//unsigned int id;
	/*vector <std::shared_ptr<node> > line;*/
	
	vector <node> line;
	vector<float> output;
	vector<float> output_v;
	vector<float> fitlex; // fitnesses for lexicase selection
	vector<float> f; // vector of objectives for multi-objective implementations (PS_sel)
	std::vector<unsigned int> outstacklen;
	std::vector<float> outstack; // linearized outstack
	std::vector<int> dominated; //for spea2 strength
	std::string eqn;
	std::string eqn_form; // equation form for string distance comparison to other forms	
	float abserror;
	float abserror_v;
	float corr;
	float corr_v;
	float VAF;
	float VAF_v;
	float fitness;
	float fitness_v;	
	float FEvar; //variance in fitness estimates (for sorting purposes)
	float GEvar; //variance in generality estimates (for sorting purposes)
	float genty; //generality
	float spea_fit;
	float parentfitness;
	int eff_size;
	int age;
	int rank;
	int complexity;
	char origin; // x: crossover, m: mutation, i: initialization

	/* 
	===================================================================
	BE SURE TO ADD ANY NEW VARIABLES TO THE SWAP FUNCTION FOR COPYING!!
	===================================================================
	*/

	ind()
	{
		abserror = 0;
		corr = 0;
		age = 1;
		genty = 1;
	}
	/*ind(const ind& x)
	{
		*this = x;
	}*/
	~ind() {
		//cout << "ind destructor\n";
		//if(!line.empty())
		//{
		//	for(vector<node*>::iterator it= line.begin();it!=line.end();it++)
		//		delete(*it);
		//	line.clear();
		//	//cout << "ind destructor deleted line nodes\n";
		//}
		//
	}
	ind & operator = (ind s) // over-ride copy construction with swap 
    {
      s.swap (*this); // Non-throwing swap
      return *this;
    }
 
    void swap (ind &s) 
	{	
		
		line.swap(s.line);						// vectors
		output.swap(s.output);
		output_v.swap(s.output_v);
		fitlex.swap(s.fitlex);
		f.swap(s.f);
		outstacklen.swap(s.outstacklen);
		outstack.swap(s.outstack);		
		dominated.swap(dominated);
		eqn.swap(s.eqn);						// strings
		eqn_form.swap(s.eqn_form);

		using std::swap;
		
		swap(this->abserror,s.abserror);		// floats
		swap(this->abserror_v,s.abserror_v);
		swap(this->corr,s.corr);
		swap(this->corr_v,s.corr_v);
		swap(this->VAF,s.VAF);
		swap(this->VAF_v,s.VAF_v);
		swap(this->fitness,s.fitness);
		swap(this->fitness_v,s.fitness_v);
		swap(this->FEvar,s.FEvar);
		swap(this->GEvar,s.GEvar);
		swap(this->genty,s.genty);
		swap(this->spea_fit,s.spea_fit);
		swap(this->parentfitness,s.parentfitness);
		swap(this->eff_size,s.eff_size);		// ints
		swap(this->age,s.age);
		swap(this->rank,s.rank);
		swap(this->complexity,s.complexity);
		
		swap(this->origin,s.origin);			// chars
		

	}//throw (); // Also see the non-throwing swap idiom
	////swap optimization
	//void swap(ind&) throw();
	//void init(string& nom_mod)
	//{
	//	eqn = nom_mod;
	//	ptr.push_back(1);
	//	ptr.push_back(nom_mod.size()-2);
	//	//nominal_model=nom_mod;
	//	//expression.register_symbol_table(d.symbol_table);
	//}

	void clrPhen()
	{
		abserror = 0;
		abserror_v=0;
		corr = 0;
		corr_v=0;
		fitness=0;
		fitness_v=0;
		eqn = "";
		eqn_form="";
		output.clear();
		output_v.clear();
		genty = 1;
		//outstack.clear();
		// nominal model must be encased in set of parenthesis. the pointer points to that which is encased.
		//ptr[0]= 1;
		//ptr[1] = nom_mod.size()-2;
	}
//private:
//	string& nominal_model;
};
struct sub_ind
{
	float fitness;
	float abserror;
	float corr;
	float VAF;

	float abserror_v;
	float corr_v;
	float VAF_v;
	string eqn;
	int age;
	int complexity;
	sub_ind(){}
	void init(ind& x){
		fitness = x.fitness; 
		abserror = x.abserror; 
		abserror_v = x.abserror_v;
		corr = x.corr; 
		corr_v = x.corr_v; 
		VAF = x.VAF;
		VAF_v = x.VAF_v;
		eqn = x.eqn; 
		age=x.age; 
		complexity = x.complexity;
	}
	~sub_ind(){}
};
//swap optimization
inline void swap(ind& lhs, ind& rhs) { lhs.swap(rhs); }
namespace std { template<> inline void swap<struct ind>(ind& lhs, ind& rhs)	{lhs.swap(rhs);	}}
////using std::swap;
struct SortFit{ bool operator() (const ind& i,const ind& j) { return (i.fitness<j.fitness);} };
struct SortFit2{ bool operator() (const sub_ind& i,const sub_ind& j) { return (i.fitness<j.fitness);} };
struct SortRank{ bool operator() (const ind& i,const ind& j) { return (i.rank<j.rank);} };
struct SortGenty{ bool operator() (const ind& i,const ind& j) { return (i.genty<j.genty);} };
struct revSortRank{	bool operator() (ind& i,ind& j) { return (i.rank>j.rank);} };

struct SortEqnSize{	bool operator() (const ind& i,const ind& j) { return (i.eqn.size()<j.eqn.size());} };
struct SortFEVar{	bool operator() (const ind& i,const ind& j) { return (i.FEvar>j.FEvar);} };
struct SortGEVar{	bool operator() (const ind& i,const ind& j) { return (i.GEvar>j.GEvar);} };
struct SortComplexity{bool operator() (const ind& i,const ind& j) { return (i.complexity<j.complexity);}};
struct SortFit_v{	bool operator() (const ind& i,const ind& j) { return (i.fitness_v<j.fitness_v);}};
struct SortSize{	bool operator() (const ind& i,const ind& j) { return (i.line.size()<j.line.size());}};
struct SortAge{	bool operator() (const ind& i,const ind& j) { return (i.age < j.age );}};

struct sameEqn{	bool operator() (ind& i,ind& j) { return i.eqn==j.eqn;} };

struct sameEqnSize{	bool operator() (ind& i,ind& j) { return i.eqn.size()==j.eqn.size();} };

struct sameSizeFit{	bool operator() (ind& i,ind& j) { return (i.fitness==j.fitness && i.eqn.size()==j.eqn.size());} };

struct sameFit{	bool operator() (ind& i,ind& j) { return (i.fitness==j.fitness);} };
struct sameFit2{	bool operator() (sub_ind& i,sub_ind& j) { return (i.fitness==j.fitness);} };
struct sameComplexity{bool operator() (const ind& i,const ind& j) { return (i.complexity==j.complexity);} };

struct sameFitComplexity{bool operator() (const ind& i,const ind& j) { return (i.fitness==j.fitness && i.complexity==j.complexity);} };


struct tribe{ 

	vector <ind> pop; // population
	float best;
	float worst;

	tribe(int size,float& max_fit,float& min_fit)	
	{
		pop.resize(size);
		best=max_fit;
		worst=min_fit;
		/*for(unsigned int i = 1;i<pop.size();++i)
			pop.at(i).init(nom_mod);*/
		maxf=max_fit;
		minf=min_fit;
	}
	~tribe() {}

	float bestFit() // returns best fitness value
	{

		/*#pragma omp parallel
		{
		   float localbest = maxf;

		   #pragma omp for schedule(static)
		   for(int i = 0; i < pop.size(); ++i)
			   localbest = min(localbest, pop.at(i).fitness);

		   #pragma omp critical
		   {
			  best = min(localbest, best);
		   }
		}*/
		best = maxf;
		for(int i = 0; i < pop.size(); ++i)
			   best = min(best, pop.at(i).fitness);
		return best;
	}
	float bestFit_v() // returns best fitness value
	{

		/*#pragma omp parallel
		{
		   float localbest = maxf;

		   #pragma omp for schedule(static)
		   for(int i = 0; i < pop.size(); ++i)
			   localbest = min(localbest, pop.at(i).fitness);

		   #pragma omp critical
		   {
			  best = min(localbest, best);
		   }
		}*/
		best = maxf;
		for(int i = 0; i < pop.size(); ++i)
			   best = min(best, pop.at(i).fitness_v);
		return best;
	}
	float worstFit() //worst fitness
	{
		worst = minf;
		/*#pragma omp parallel
		{
		   float localworst = minf;

		   #pragma omp for schedule(static)
		   for(int i = 0; i < pop.size(); ++i)
			   localworst = max(localworst, pop.at(i).fitness);

		   #pragma omp critical
		   {
			  worst = max(localworst, worst);
		   }
		}*/
		 for(int i = 0; i < pop.size(); ++i)
			 worst = max(worst, pop.at(i).fitness);
		return worst;
	}	
	float medFit() //median fitness
	{
		vector<float> fitness(pop.size());
		for(int i =0; i < pop.size(); i++)
			fitness.at(i) = pop.at(i).fitness;
		sort(fitness.begin(),fitness.end());
		if (pop.size() % 2==0) //even
			return fitness.at((int)floor((float)pop.size()/2));
		else
			return (fitness.at(pop.size()/2)+fitness.at(pop.size()/2-1))/2;
	} 
	float medFit_v() //median fitness
	{
		vector<float> fitness(pop.size());
		for(int i =0; i < pop.size(); i++)
			fitness.at(i) = pop.at(i).fitness_v;
		sort(fitness.begin(),fitness.end());
		if (pop.size() % 2==0) //even
			return fitness.at((int)floor((float)pop.size()/2));
		else
			return (fitness.at(pop.size()/2)+fitness.at(pop.size()/2-1))/2;

	}
	float meanFit() // mean fitness
	{
		float answer=0;
		//#pragma omp parallel for reduction(+ : answer)
		for(int i=0; i<pop.size(); ++i)
		{
			answer+=pop.at(i).fitness;
		}
		return (float)answer/pop.size();
	}
			
	float meanSize() // mean line length
	{
		float answer=0;
		//#pragma omp parallel for reduction(+ : answer)
		for(int i=0; i<pop.size(); ++i)
		{
			answer+=pop.at(i).line.size();
		}
		return (float)answer/pop.size();
	}
	float meanEffSize()
	{
		float answer=0;
		//#pragma omp parallel for reduction(+ : answer)
		for(int i=0; i<pop.size(); ++i)
		{
			answer+=pop.at(i).eff_size;
		}
		return (float)answer/pop.size();
	}
	int medSize() // median line length
	{		
		//vector<ind> tmppop = pop;
		sort(pop.begin(),pop.end(),SortSize());
		int index = (int)floor((float)pop.size()/2);
		return int(pop.at(index).line.size()); 
	}

	void topTen(vector <sub_ind>& eqns) //returns address to vector of equation strings
	{
		vector<sub_ind> tmppop(pop.size());
		for (int i = 0;i<pop.size();++i) tmppop[i].init(pop.at(i));
		//vector<ind> tmppop = pop;
		sort(tmppop.begin(),tmppop.end(),SortFit2());
		unique(tmppop.begin(),tmppop.end(),sameFit2());
		
		for (int i=0;i<10;++i)
			eqns.push_back(tmppop.at(i));

		/*vector <float> fitnesses;
		int i=0;
		bool pass=true;
		while(eqns.size()<10 && i<pop.size())
		{
			fitnesses.push_back(pop.at(i).fitness);
			for(unsigned int j=0;j<fitnesses.size()-1;++j)		 
			{
				if(fitnesses.at(j)==fitnesses.back())
				{
					fitnesses.pop_back();
					pass=0;
					break;
				}
			}

			if (pass)
				eqns.push_back(pop.at(i));
			else
				pass=1;
			++i;
		}*/
	}
	void getbestsubind(sub_ind& bestind)
	{
		vector<sub_ind> subpop(pop.size());
		for (int i = 0;i<pop.size();++i) subpop[i].init(pop.at(i));
		//vector<ind> tmppop = pop;
		sort(subpop.begin(),subpop.end(),SortFit2());
		bestind = subpop.front();
	}// address of best individual
	void getbestind(ind& bestind)
	{
		//vector<ind> tmppop = pop;
		sort(pop.begin(),pop.end(),SortFit());
		bestind = pop.front();
	}// address of best individual
	void sortpop()
	{
		sort(pop.begin(),pop.end(),SortFit());
	}
	void hom(vector<Randclass>& r, float& tot_hom, float& on_hom, float& off_hom)
	{
		tot_hom = 0; on_hom=0; off_hom=0;
		//float sum_strdist=0;
		int c1, c2, s_tot,s_on,s_off;
		float tot_tmp=0,on_tmp=0,off_tmp=0;
		int samplesize=200;
		std::string tmp1, tmp2, tmp1on, tmp2on, tmp1off, tmp2off;
		//std::string tmp2;
		for (int i=0; i<samplesize; ++i)
		{
			//reset temporary strings
			tmp1.resize(0); tmp2.resize(0); 
			tmp1on.resize(0); tmp2on.resize(0); 
			tmp1off.resize(0); tmp2off.resize(0);
			tot_tmp = 0; on_tmp = 0; off_tmp = 0;
			c1 = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
			c2 = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);

			for (int j=pop.at(c1).line.size()-1; j>=0;--j){
				if (pop.at(c1).line.at(j).type=='v') tmp1 += pop.at(c1).line.at(j).varname;
				else tmp1 += pop.at(c1).line.at(j).type;

				if(pop.at(c1).line.at(j).on){
					if (pop.at(c1).line.at(j).type=='v') tmp1on += pop.at(c1).line.at(j).varname;
					else tmp1on += pop.at(c1).line.at(j).type;
				}
				/*else
					tmp1on += ' ';*/

				if(!pop.at(c1).line.at(j).on){
					if (pop.at(c1).line.at(j).type=='v') tmp1off += pop.at(c1).line.at(j).varname;
					else tmp1off += pop.at(c1).line.at(j).type;
				}
				/*else
					tmp1off += ' ';*/
			}

			for (int j=pop.at(c2).line.size()-1; j>=0;--j){
				if (pop.at(c2).line.at(j).type=='v') tmp2 += pop.at(c2).line.at(j).varname;
				else tmp2 += pop.at(c2).line.at(j).type;

				if(pop.at(c2).line.at(j).on){
					if (pop.at(c2).line.at(j).type=='v') tmp2on += pop.at(c2).line.at(j).varname;
					else tmp2on += pop.at(c2).line.at(j).type;
				}
				/*else
					tmp2on += ' ';*/

				if(!pop.at(c2).line.at(j).on){
					if (pop.at(c2).line.at(j).type=='v') tmp2off += pop.at(c2).line.at(j).varname;
					else tmp2off += pop.at(c2).line.at(j).type;
				}
				/*else
					tmp2off += ' ';*/
			}

			s_tot = strdist(tmp1,tmp2);
			s_on = strdist(tmp1on,tmp2on);
			//s_off = s_tot-s_on;
			if (tmp1off.size()>0 && tmp2off.size()>0) s_off = strdist(tmp1off,tmp2off);
			else s_off = 13785;

			tot_tmp = float(s_tot)/float(std::max(tmp1.size(),tmp2.size()));
			on_tmp  = float(s_on)/float(std::max(tmp1on.size(),tmp2on.size()));
			if (s_off!= 13785) off_tmp = float(s_off)/float(std::max(tmp1off.size(),tmp2off.size()));
			else off_tmp = 1;

			tot_hom += tot_tmp;
			on_hom += on_tmp;
			off_hom += off_tmp;
		}

		tot_hom = 1-tot_hom/samplesize;
		on_hom = 1-on_hom/samplesize;
		off_hom = 1-off_hom/samplesize;
		
	}
	/*float on_hom(vector<Randclass>& r){
		float sum_strdist=0;
		int c1, c2;
		int samplesize = 100;
		std::string tmp1;
		std::string tmp2;
		for (int i=0; i<samplesize; ++i)
		{
			c1 = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
			c2 = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);

			for(int j=pop.at(c1).line.size()-1; j>=0;--j){
				if(pop.at(c1).line.at(j).on)
					tmp1 += pop.at(c1).line.at(j).type;
				else
					tmp1 += ' ';
			}
			for (int j=pop.at(c2).line.size()-1; j>=0;--j){
				if(pop.at(c2).line.at(j).on)
					tmp2 += pop.at(c2).line.at(j).type;
				else
					tmp2 += ' ';
			}
			sum_strdist += strdist(tmp1,tmp2)/float(std::max(tmp1.size(),tmp2.size()));
		}

		return 1-sum_strdist/samplesize;
	}
	float off_hom(vector<Randclass>& r){
		float sum_strdist=0;
		int c1, c2;
		int samplesize=100;
		std::string tmp1;
		std::string tmp2;
		for (int i=0; i<samplesize; ++i)
		{
			c1 = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);
			c2 = r[omp_get_thread_num()].rnd_int(0,pop.size()-1);

			for(int j=pop.at(c1).line.size(); j>0;--j){
				if(!pop.at(c1).line.at(j-1).on)
					tmp1 += pop.at(c1).line.at(j-1).type;
				else
					tmp1 += ' ';
			}
			for (int j=pop.at(c2).line.size(); j>0;--j){
				if(!pop.at(c2).line.at(j-1).on)
					tmp2 += pop.at(c2).line.at(j-1).type;
				else
					tmp2 += ' ';
			}
			sum_strdist += strdist(tmp1,tmp2)/float(std::max(tmp1.size(),tmp2.size()));
		}

		return 1-sum_strdist/samplesize;
	}*/
	/*
private:
	
	bool fitlow (ind& i,ind& j) { return (i.fitness<j.fitness); }
	bool eqncomp(ind& i,ind& j) { return (i.eqn_form.compare(j.eqn_form)==0); }
	bool fitcomp (ind& i,ind& j) { return (i.fitness==j.fitness); }
	*/

private:
		float maxf;
		float minf;


};


#endif
