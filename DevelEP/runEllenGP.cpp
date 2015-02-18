#include "stdafx.h"
#include <string>
// mine
#include "pop.h"
#include "params.h"
#include "rnd.h"
#include "data.h"
#include "state.h"
#include "logger.h"

#include "InitPop.h"
#include "FitnessEstimator.h"
#include "Fitness.h"
#include "Generation.h"
#include "instructionset.h"
//#include "omp.h"
#include "Generationfns.h"
#include "strdist.h"
#include <time.h>
#include <cstring>
#include "p_archive.h"
#include "Eqn2Line.h"
#include "general_fns.h"

//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

using namespace std;

// global parameters structure


//void load_params(params &p, std::ifstream& is);
//void load_data(data &d, std::ifstream& is,params&);
//void load_lexdata(data &d, std::ifstream& fs,params& p);
//bool stopcondition(tribe& T,params& p,data& d,state& s,FitnessEstimator& FE);
//void printstats(tribe& T,int &i,state& s,params& p,paretoarchive& A);
//void printbestind(tribe& T,params& p,state& s,string& logname);
//void printpop(vector<ind>& pop,params& p,state& s,string& logname,int type);
//void shuffle_data(data& d, params& p, vector<Randclass>& r,state& s);

bool check_genty(vector<ind>& pop,params& p){
	
	for(int count = 0; count<pop.size(); ++count)
	{
		float tmp = abs(pop[count].fitness-pop[count].fitness_v)/pop[count].fitness;
		if (pop.at(count).genty != tmp && pop.at(count).genty != p.max_fit) 
			return 0;
	}
	return 1;
}
void printGenome(tribe& T,int gen,string& logname,data& d,params& p)
{
	// print genome to file in csv format
	// create data file
	// print format: print in head-to-tail fashion with grid locations
	// x y gene
	// gene encoding: +:1,-:2,*:3,/:4,sin:5,cos:6,exp:7,log:8,constant:9,vars:10-length(vars)

	 std::ofstream gout;
	 string filename = logname.substr(0,logname.size()-4)+"_g.csv."+std::to_string(static_cast<long long>(gen));
	 gout.open(filename);
	 gout << "x,y,fit,age,g,e\n";
	 string out;
	 string tmp;
	 int k;
	 float y;
	 float max_fit=0;
	 float min_fit=10000;
	 float max_age=0;
	 float min_age=100000;
	 float max_line=0;
	 float min_line = 1000;
	 for (size_t i =0; i<T.pop.size(); ++i){
		if (log(1+T.pop[i].fitness)>max_fit) max_fit=log(1+T.pop[i].fitness);
		if (log(1+T.pop[i].fitness)<min_fit) min_fit=log(1+T.pop[i].fitness);
		if (T.pop[i].age>max_age) max_age=T.pop[i].age;
		if (T.pop[i].age<min_age) min_age=T.pop[i].age;
		for (size_t j = 0; j<T.pop[i].line.size(); ++j){
			if (T.pop[i].line.size()>max_line) max_line=T.pop[i].line.size();
			if (T.pop[i].line.size()<min_line) min_line=T.pop[i].line.size();
		}
	}
	 T.sortpop();
	// T.sortpop_age();
	 for (size_t i=0; i<T.pop.size(); ++i){
		 y=0;
		 for (int j=T.pop[i].line.size()-1; j>=0; --j){
			 switch(T.pop[i].line[j].type){
			 case '+':
				 out="1";
				 break;
			 case '-':
				 out="2";
				 break;
			 case '*':
				 out="3";
				 break;
			 case '/':
				 out="4";
				 break;
			 case 's':
				 out="5";
				 break;
			 case 'c':
				 out="6";
				 break;
			 case 'e':
				 out="7";
				 break;
			 case 'l':
				 out="8";
				 break;
			 case 'n':
				 out="9";
				 break;
			 case 'v':
				 k=0;
				 while (T.pop[i].line[j].varname.compare(d.label[k])!=0)
					 ++k;
				 tmp=std::to_string(static_cast<long long>(k+10));
				 out = tmp;
				 break;
			 }
			 gout << float(float(i)/float(T.pop.size())) << "," << y/p.max_len << "," << (log(1+T.pop[i].fitness)-min_fit)/(max_fit-min_fit) << "," << (float(T.pop[i].age)-min_age)/(max_age-min_age) << "," << out << "," << T.pop[i].line[j].on << "\n";
			++y;
		}
	 }
	 gout.close();
}
void printstats(tribe& T,int &i,state& s,params& p,paretoarchive& A)
{
	
	//boost::progress_timer timer;
s.out << "--- Generation " << i << "---------------------------------------------------------------" << "\n";
s.out << "population size: " << T.pop.size() << "\n";
s.out << "Number of evals: " << s.genevals.back() << "\n";
s.out << "Best Fitness: " << T.bestFit() <<"\n";
s.out << "Best Fitness (v): " << T.bestFit_v() <<"\n";
s.out << "Median Fitness: " << T.medFit_v()<<"\n";
s.out << "Median Fitness (v): " << T.medFit()<<"\n";
s.out << "Mean Size: " << T.meanSize() << "\n";
s.out << "Mean Eff Size: " << T.meanEffSize() << "\n";
s.out << "Pareto Front Equations: " << A.optimal_size << "\n";
if(p.pHC_on)
	s.out << "Parameter updates: " << float(s.setPHCupdates())/float(p.popsize)*100 << "%%\n";
if(p.eHC_on)
	s.out << "Epigenetic updates: " << float(s.setEHCupdates())/float(p.popsize)*100 << "%%\n";
s.out << "Beneficial Genetics: " << s.getGoodCrossPct() << "%%\n";
s.out << "Neutral Genetics: " << s.getNeutCrossPct() << "%%\n";
s.out << "Bad Genetics: " << s.getBadCrossPct() << "%\n";
//s.clearCross();
//float totalshares = 0;
//float c1=0;
//for (int i = 0; i<T.pop.size();++i)
//{
//	for (int j=0;j<T.pop.at(i).line.size();++j){
//		totalshares+=float(T.pop.at(i).line.at(j).use_count()); c1++;}
//}
//s.out << "Average shared pointer use count: " << totalshares/(c1) << "\n";
s.out << "MAE \t R^2 \t Fitness \t Equation \n";
vector <sub_ind> besteqns;
T.topTen(besteqns);
//for(unsigned int j=0;j<min(10,int(A.pop.size()));++j)
//	s.out <<A.pop.at(j).abserror_v << "\t" << A.pop.at(j).corr_v << "\t" << A.pop.at(j).eqn <<"\n";
for(unsigned int j=0;j<besteqns.size();++j){
	s.out <<besteqns.at(j).abserror << "\t" << besteqns.at(j).corr << "\t" << besteqns.at(j).fitness << "\t" << besteqns.at(j).eqn <<"\n";
	/*if(boost::math::isnan(besteqns.at(j).abserror))
	{
		cout << "equation with NaN error: " + besteqns.at(j).eqn + "\n";	
	}*/
}

s.out << "-------------------------------------------------------------------------------" << "\n";
}
//void printstatsP(tribe T,int &i,state s,params& p,paretoarchive A)
//{
//	
//	//boost::progress_timer timer;
//s.out << "--- Generation " << i << "---------------------------------------------------------------" << "\n";
//s.out << "Number of evals: " << s.genevals.back() << "\n";
//s.out << "Best Fitness: " << T.bestFit() <<"\n";
//s.out << "Best Fitness (v): " << T.bestFit_v() <<"\n";
//s.out << "Median Fitness: " << T.medFit_v()<<"\n";
//s.out << "Median Fitness: " << T.medFit()<<"\n";
//s.out << "Mean Size: " << T.meanSize() << "\n";
//s.out << "Mean Eff Size: " << T.meanEffSize() << "\n";
//s.out << "Pareto Front Equations: " << A.optimal_size << "\n";
//if(p.pHC_on)
//	s.out << "Parameter updates: " << s.setpHCupdates() << "\n";
//if(p.eHC_on)
//	s.out << "Epigenetic updates: " << s.geteHCupdates() << "\n";
//s.out << "Beneficial Genetics: " << s.getGoodCrossPct() << "%%\n";
//s.out << "Neutral Genetics: " << s.getNeutCrossPct() << "%%\n";
//s.out << "Bad Genetics: " << s.getBadCrossPct() << "%\n";
//
////float totalshares = 0;
////float c1=0;
////for (int i = 0; i<T.pop.size();++i)
////{
////	for (int j=0;j<T.pop.at(i).line.size();++j){
////		totalshares+=float(T.pop.at(i).line.at(j).use_count()); c1++;}
////}
////s.out << "Average shared pointer use count: " << totalshares/(c1) << "\n";
//s.out << "MAE \t R^2 \t Fitness \t Equation \n";
//vector <ind> besteqns;
//T.topTen(besteqns);
////for(unsigned int j=0;j<min(10,int(A.pop.size()));++j)
////	s.out <<A.pop.at(j).abserror_v << "\t" << A.pop.at(j).corr_v << "\t" << A.pop.at(j).eqn <<"\n";
//for(unsigned int j=0;j<besteqns.size();++j)
//	s.out <<besteqns.at(j).abserror << "\t" << besteqns.at(j).corr << "\t" << besteqns.at(j).fitness << "\t" << besteqns.at(j).eqn <<"\n";
//
//s.out << "-------------------------------------------------------------------------------" << "\n";
//}

void printbestind(tribe& T,params& p,state& s,string& logname)
{
	s.out << "saving best ind... \n";
	ind best;
	T.getbestind(best);
	string bestname = logname.substr(0,logname.size()-4)+".best_ind";
	std::ofstream fout;
	fout.open(bestname);
	//boost::progress_timer timer;
	fout << "--- Best Individual ---------------------------------------------------------------" << "\n";
	fout << "Corresponding Logfile: " + logname + "\n";
	fout << "Total Evaluations: " << s.totalevals() << "\n";
	fout << "f = " + best.eqn + "\n";
	fout << "gline: ";
	for(unsigned int i =0;i<best.line.size();++i)
	{
		if (best.line[i].type=='n')
			fout << best.line[i].value << "\t";
		else if (best.line[i].type=='v')
			fout << best.line[i].varname << "\t";
		else
			fout << best.line[i].type << "\t";
	}
	fout << endl;
	fout << "eline: ";
	for(unsigned int i =0;i<best.line.size();++i){
		if (best.line.at(i).on)
			fout <<"1\t";
		else
			fout <<"0\t";
	}
	fout << endl;
	fout << "size: " << best.line.size() << "\n";
	fout << "eff size: " << best.eff_size << "\n";
	fout << "training abs error: " << best.abserror<< "\n";
	fout << "training correlation: " << best.corr<< "\n";
	fout << "training fitness: " << best.fitness<< "\n";
	fout << "training VAF: " << best.VAF <<"\n";
	fout << "validation abs error: " << best.abserror_v<< "\n";
	fout << "validation correlation: " << best.corr_v<< "\n";
	fout << "validation fitness: " << best.fitness_v<< "\n";
	fout << "validation VAF: " << best.VAF_v <<"\n";
	fout << "parent fitness: " << best.parentfitness << "\n";
	fout << "origin: " << best.origin << "\n";
	fout << "age: " << best.age << "\n";
	fout << "eqn form: " << best.eqn_form << "\n";
	/*fout << "output: ";
	for(unsigned int i =0;i<best.output.size();++i)
	{
		fout << best.output.at(i) << " ";
	}
	fout<<"\n";*/
}
void initdatafile(std::ofstream& dfout,string & logname)
{
	string dataname = logname.substr(0,logname.size()-4)+".data";
	dfout.open(dataname,std::ofstream::out | std::ofstream::app);
	//dfout.open(dataname,std::ofstream::app);
	//dfout << "pt_evals \t best_eqn \t best_fit \t best_fit_v \t med_fit \t med_fit_v \t best_MAE \t best_MAE_v \t best_R2 \t best_R2_v \t best_VAF \t best_VAF_v \t size \t eff_size \t pHC_pct \t eHC_pct \t good_g_pct \t neut_g_pct \t bad_g_pct \t tot_hom \t on_hom \t off_hom\n";
	dfout << "gen,pt_evals,best_eqn,best_fit,best_fit_v,med_fit,med_fit_v,best_MAE,best_MAE_v,best_R2,best_R2_v,size,eff_size,pHC_pct,eHC_pct,good_g_pct,neut_g_pct,bad_g_pct,tot_hom,on_hom,off_hom\n";
	//fout.close(dataname);
}
void printdatafile(tribe& T,state& s,params& p, vector<Randclass>& r,std::ofstream& dfout, int gen)
{
	//string dataname = logname.substr(0,logname.size()-4)+".data";
	//std::ofstream fout;
	//dfout.open(dataname,std::ofstream::app);

	sub_ind best_ind;
	T.getbestsubind(best_ind);

	/*dfout << s.totalptevals() << "\t" << best_ind.eqn << "\t" << T.bestFit() << "\t" << T.bestFit_v() << "\t" << T.medFit() << "\t" << T.medFit_v() << "\t" << best_ind.abserror << "\t" << best_ind.abserror_v << "\t" << best_ind.corr << "\t" << best_ind.corr_v << "\t" << T.meanSize() << "\t" << T.meanEffSize() << "\t" << s.current_pHC_updates/float(p.popsize)*100.0 << "\t" << s.current_eHC_updates/float(p.popsize)*100.0 << "\t" <<  s.good_cross_pct << "\t" << s.neut_cross_pct << "\t" << s.bad_cross_pct;*/
	dfout << gen << "," << s.totalptevals() << "," << best_ind.eqn << "," << T.bestFit() << "," << T.bestFit_v() << "," << T.medFit() << "," << T.medFit_v() << "," << best_ind.abserror << "," << best_ind.abserror_v << "," << best_ind.corr << "," << best_ind.corr_v << "," << best_ind.VAF << "," << best_ind.VAF_v << T.meanSize() << "," << T.meanEffSize() << "," << s.current_pHC_updates/float(p.popsize)*100.0 << "," << s.current_eHC_updates/float(p.popsize)*100.0 << "," <<  s.good_cross_pct << "," << s.neut_cross_pct << "," << s.bad_cross_pct;
	if (p.print_homology){
		float tot_hom, on_hom, off_hom;
		T.hom(r,tot_hom,on_hom,off_hom);
		dfout << "," << tot_hom << "," << on_hom << "," << off_hom;
	}
	dfout <<"\n";

	//s.clearCross();
}
void printpop(vector<ind>& pop,params& p,state& s,string& logname,int type)
{
	string bestname;
	if (type==1){
		//s.out << "saving pareto archive... \n";
		bestname = logname.substr(0,logname.size()-4)+".archive";
		sort(pop.begin(),pop.end(),SortComplexity());
		stable_sort(pop.begin(),pop.end(),SortRank());
	}
	else if (type == 2){
		string gen = to_string(static_cast<long long>(s.genevals.size()));
		s.out << "saving pop... \n";
		bestname = logname.substr(0,logname.size()-4) + "gen" + gen + ".pop";
		sort(pop.begin(),pop.end(),SortFit());
	}
	else if (type == 3){
		string gen = to_string(static_cast<long long>(s.genevals.size()));
		s.out << "saving initial pop... \n";
		bestname = logname.substr(0,logname.size()-4) + ".init_pop";
		sort(pop.begin(),pop.end(),SortFit());
	}
	else {
		s.out << "saving last pop... \n";
		bestname = logname.substr(0,logname.size()-4)+".last_pop";
		sort(pop.begin(),pop.end(),SortComplexity());
		stable_sort(pop.begin(),pop.end(),SortFit());
	}
	
	std::ofstream fout;
	fout.open(bestname);
	//boost::progress_timer timer;
	fout << "--- Population ---------------------------------------------------------------" << "\n";
	fout << "Corresponding Logfile: " + logname + "\n";
	fout << "Total Evaluations: " << s.totalevals() << "\n";

	for (int h = 0; h<pop.size(); ++h){
		fout << "--- Individual "<< h << " ------------------------------------------------------------" << "\n";
		fout << "f = " + pop.at(h).eqn + "\n";
		fout << "gline: ";
		for(unsigned int i =0;i<pop.at(h).line.size();++i)
		{
			if (pop.at(h).line[i].type=='n')
				fout <<pop.at(h).line[i].value << "\t";
			else if (pop.at(h).line[i].type=='v')
				fout << pop.at(h).line[i].varname << "\t";
			else
				fout << pop.at(h).line[i].type << "\t";
		}
		fout << endl;
		fout << "eline: ";
		for(unsigned int i =0;i<pop.at(h).line.size();++i){
			if (pop.at(h).line.at(i).on)
				fout <<"1\t";
			else
				fout <<"0\t";
		}
		fout << endl;
		fout << "size: " << pop.at(h).line.size() << "\n";
		fout << "eff size: " << pop.at(h).eff_size << "\n";
		fout << "complexity: " << pop.at(h).complexity << "\n";
		fout << "fitness: " << pop.at(h).fitness<< "\n";
		fout << "fitness_v: " << pop.at(h).fitness_v<< "\n";
		fout << "MAE: " << pop.at(h).abserror<< "\n";;
		fout << "MAE_v: " << pop.at(h).abserror_v<< "\n";
		fout << "correlation: " << pop.at(h).corr<< "\n";
		fout << "correlation_v: " << pop.at(h).corr_v<< "\n";
		fout << "VAF: " << pop.at(h).VAF<< "\n";
		fout << "VAF_v: " << pop.at(h).VAF_v<< "\n";
		fout <<  "rank: " << pop.at(h).rank << "\n";
		fout << "parent fitness: " << pop.at(h).parentfitness << "\n";
		fout << "origin: " << pop.at(h).origin << "\n";
		fout << "age: " << pop.at(h).age << "\n";
		fout << "eqn form: " << pop.at(h).eqn_form << "\n";
		/*fout << "output: ";
		for(unsigned int i =0;i<pop.at(h).output.size();++i)
		{
			fout << T.pop.at(h).output.at(i);
		}
		fout<<"\n";*/
		fout << "------------------------------------------------------------------" << "\n";
		}
}
void load_params(params &p, std::ifstream& fs)
{
	if (!fs.good()) 
	{
			cerr << "BAD PARAMETER FILE LOCATION" << "\n";
			exit(1);
	}

	string s;
	string varname;
	float tmpf;
	//int tmpi;
	string tmps;
	//bool tmpb;

	//string trash;
	//s.erase();
    //s.reserve(is.rdbuf()->in_avail());

    while(!fs.eof())
    {		
		getline(fs,s,'\n');
		istringstream ss(s);
		
		ss >> varname;

		//getline(is,s,'\t');

		if(varname.compare("g") == 0)
		{
			ss >> p.g;
			//p.g = tmp;
		}
		else if(varname.compare("popsize") == 0)
			ss >> p.popsize;
		else if(varname.compare("sel") == 0)
			ss>>p.sel;
		else if(varname.compare("tourn_size") == 0)
			ss>>p.tourn_size;
		else if(varname.compare("rt_rep") == 0)
			ss>>p.rt_rep;
		else if(varname.compare("rt_cross") == 0)
			ss>>p.rt_cross;
		else if(varname.compare("rt_mut") == 0)
			ss>>p.rt_mut;
		else if(varname.compare("cross") == 0)
			ss>>p.cross;
		else if(varname.compare("cross_ar") == 0)
			ss>>p.cross_ar;
		else if(varname.compare("mut_ar") == 0)
			ss>>p.mut_ar;
		//~ else if(varname.compare("stoperror") == 0)
			//~ ss>>p.stoperror;
		else if(varname.compare("init_validate_on") == 0)
			ss>>p.init_validate_on;
		else if(varname.compare(0,11,"resultspath") == 0)
		{
			int q=0;
			while (ss>>tmps)
			{
				if ( q > 0)
					p.resultspath = p.resultspath + ' ';
				p.resultspath.insert(p.resultspath.end(),tmps.begin(),tmps.end());
				++q;
			}
		}
		//~ else if(varname.compare(0,4,"loud") == 0)
			//~ ss>>p.loud;
		/*else if(varname.compare(0,8,"parallel") == 0)
			ss>>p.parallel;
		else if(varname.compare(0,8,"numcores") == 0)
			ss>>p.numcores;*/
		//~ else if(varname.compare(0,11,"sim_nom_mod") == 0)
			//~ ss>>p.sim_nom_mod;
		//~ else if(varname.compare(0,7,"nstates") == 0)
			//~ ss>>p.nstates;
		else if(varname.compare(0,7,"intvars") == 0)
		{
			while (ss>>tmps)
				p.intvars.push_back(tmps);
		}
		else if(varname.compare(0,7,"extvars") == 0)
		{
			while (ss>>tmps)
				p.extvars.push_back(tmps);
		}
		/*else if(varname.compare("cons") == 0)
		{
			while (ss>>tmps)
				p.cons.push_back(tmps);
		}*/
		else if(varname.compare("cvals") == 0)
		{
			while (ss>>tmpf){
				p.cvals.push_back(tmpf);
				p.cons.push_back(std::to_string(static_cast<long double>(tmpf)));
			}

		}
		else if(varname.compare("seeds") == 0)
		{
			while (ss>>tmps)
				p.seeds.push_back(tmps);
		}
		else if(varname.compare("ERC") == 0)
			ss>>p.ERC;
		else if(varname.compare("ERCints") == 0)
			ss>>p.ERCints;
		else if(varname.compare("maxERC") == 0)
			ss>>p.maxERC;
		else if(varname.compare("minERC") == 0)
			ss>>p.minERC;
		else if(varname.compare("numERC") == 0)
			ss>>p.numERC;
		else if(varname.compare("fit_type") == 0)
			ss>>p.fit_type;
		else if(varname.compare("max_fit") == 0)
			ss>>p.max_fit;
		else if(varname.compare("min_fit") == 0)
			ss>>p.min_fit;
		else if(varname.compare("op_list") == 0)
		{
			while (ss>>tmps)
				p.op_list.push_back(tmps);
		}
		else if(varname.compare("op_weight") == 0)
		{
			while(ss>>tmpf)
				p.op_weight.push_back(tmpf);
		}
		else if(varname.compare("weight_ops_on") == 0)
			ss>>p.weight_ops_on;
		else if(varname.compare("min_len") == 0)
			ss>>p.min_len;
		else if(varname.compare("max_len") == 0)
			ss>>p.max_len;
		//~ else if(varname.compare("max_dev_len") == 0)
			//~ ss>>p.max_dev_len;
		else if(varname.compare("complex_measure") == 0)
			ss>>p.complex_measure;
		//~ else if(varname.compare("precision") == 0)
			//~ ss>>p.precision;
		else if(varname.compare("lineHC_on") == 0)
			ss>>p.lineHC_on;
		else if(varname.compare("lineHC_its") == 0)
			ss>>p.lineHC_its;
		else if(varname.compare("pHC_on") == 0)
			ss>>p.pHC_on;
		else if(varname.compare("pHC_delay_on") == 0)
			ss>>p.pHC_delay_on;
		else if(varname.compare("pHC_size") == 0)
			ss>>p.pHC_size;
		else if(varname.compare("pHC_its") == 0)
			ss>>p.pHC_its;
		else if(varname.compare("pHC_gauss") == 0)
			ss>>p.pHC_gauss;
		else if(varname.compare("eHC_on") == 0)
			ss>>p.eHC_on;
		else if(varname.compare("eHC_its") == 0)
			ss>>p.eHC_its;
		else if(varname.compare("eHC_prob") == 0)
			ss>>p.eHC_prob;
		//~ else if(varname.compare("eHC_size") == 0)
			//~ ss>>p.eHC_size;
		//~ else if(varname.compare("eHC_cluster") == 0)
			//~ ss>>p.eHC_cluster;
		//~ else if(varname.compare("eHC_dev") == 0)
			//~ ss>>p.eHC_dev;
		//~ else if(varname.compare("eHC_best") == 0)
			//~ ss>>p.eHC_best;
		else if(varname.compare("eHC_init")==0)
			ss>>p.eHC_init;
		//~ else if(varname.compare("eHC_prob_scale") == 0)
			//~ ss>>p.eHC_prob_scale;
		//~ else if(varname.compare("eHC_max_prob") == 0)
			//~ ss>>p.eHC_max_prob;
		//~ else if(varname.compare("eHC_min_prob") == 0)
			//~ ss>>p.eHC_min_prob;
		else if(varname.compare("eHC_mut") == 0)
			ss>>p.eHC_mut;
		else if(varname.compare("eHC_slim") == 0)
			ss>>p.eHC_slim;
		else if(varname.compare("lexpool") == 0)
			ss>>p.lexpool;
		else if(varname.compare("lexage") == 0)
			ss>>p.lexage;
		else if(varname.compare("prto_arch_on") == 0)
			ss>>p.prto_arch_on;
		else if(varname.compare("prto_arch_size") == 0)
			ss>>p.prto_arch_size;
		else if(varname.compare("prto_sel_on") == 0)
			ss>>p.prto_sel_on;
		else if(varname.compare("islands") ==0)
			ss>>p.islands;
		else if(varname.compare("island_gens") ==0)
			ss>>p.island_gens;
		else if(varname.compare("train") ==0)
			ss>>p.train;
		else if(varname.compare("train_pct") ==0)
			ss>>p.train_pct;
		else if(varname.compare("printeverypop") == 0)
			ss>>p.printeverypop;
		else if(varname.compare("estimate_fitness") == 0)
			ss>>p.EstimateFitness;
		else if(varname.compare("FE_pop_size") == 0)
			ss>>p.FE_pop_size;
		else if(varname.compare("FE_ind_size") == 0)
			ss>>p.FE_ind_size;
		else if(varname.compare("FE_train_size") == 0)
			ss>>p.FE_train_size;
		else if(varname.compare("FE_train_gens") == 0)
			ss>>p.FE_train_gens;
		else if(varname.compare("FE_rank") == 0)
			ss>>p.FE_rank;
		else if(varname.compare("estimate_generality") == 0)
			ss>>p.estimate_generality;
		else if(varname.compare("G_sel") == 0)
			ss>>p.G_sel;
		else if(varname.compare("G_shuffle") == 0)
			ss>>p.G_shuffle;
		else if(varname.compare("norm_error") == 0)
			ss>>p.norm_error;
		else if(varname.compare("shuffle_data") == 0)
			ss>>p.shuffle_data;
		else if(varname.compare("init_trees") == 0)
			ss>>p.init_trees;
		else if(varname.compare("limit_evals") == 0)
			ss>>p.limit_evals;
		else if(varname.compare("max_evals") == 0)
			ss>>p.max_evals;
		else if(varname.compare("print_homology") == 0)
			ss>>p.print_homology;
		else if(varname.compare("print_log") == 0)
			ss>>p.print_log;
		else if(varname.compare("print_init_pop") == 0)
			ss>>p.print_init_pop;
		else if(varname.compare("print_genome") == 0)
			ss>>p.print_genome;
		else if(varname.compare("print_epigenome") == 0)
			ss>>p.print_epigenome;
		else if(varname.compare("num_log_pts") == 0)
			ss>>p.num_log_pts;
		else if(varname.compare("PS_sel") == 0)
			ss>>p.PS_sel;
		else if(varname.compare("pop_restart") == 0)
			ss>>p.pop_restart;
		else if(varname.compare("pop_restart_path") == 0)
			ss>>p.pop_restart_path;
		else{}
    }
	p.allvars = p.intvars;
	p.allvars.insert(p.allvars.end(), p.extvars.begin(), p.extvars.end());
	p.allblocks = p.allvars;
	p.allblocks.insert(p.allblocks.end(),p.cons.begin(),p.cons.end());
	p.allblocks.insert(p.allblocks.end(),p.seeds.begin(),p.seeds.end());

	p.seed = time(0);
	
	for (unsigned int i=0; i<p.op_list.size(); ++i)
	{
		if (p.op_list.at(i).compare("n")==0 )//&& ( p.ERC || !p.cvals.empty() ) )
		{
			p.op_choice.push_back(0);
			p.op_arity.push_back(0);
		}
		else if (p.op_list.at(i).compare("v")==0)
		{
			p.op_choice.push_back(1);
			p.op_arity.push_back(0);
		}
		else if (p.op_list.at(i).compare("+")==0)
		{
			p.op_choice.push_back(2);
			p.op_arity.push_back(2);
		}
		else if (p.op_list.at(i).compare("-")==0)
		{
			p.op_choice.push_back(3);
			p.op_arity.push_back(2);
		}
		else if (p.op_list.at(i).compare("*")==0)
		{
			p.op_choice.push_back(4);
			p.op_arity.push_back(2);
		}
		else if (p.op_list.at(i).compare("/")==0)
		{
			p.op_choice.push_back(5);
			p.op_arity.push_back(2);
		}
		else if (p.op_list.at(i).compare("sin")==0)
		{
			p.op_choice.push_back(6);
			p.op_arity.push_back(1);
		}
		else if (p.op_list.at(i).compare("cos")==0)
		{
			p.op_choice.push_back(7);
			p.op_arity.push_back(1);
		}
		else if (p.op_list.at(i).compare("exp")==0)
		{
			p.op_choice.push_back(8);
			p.op_arity.push_back(1);
		}
		else if (p.op_list.at(i).compare("log")==0)
		{
			p.op_choice.push_back(9);
			p.op_arity.push_back(1);
		}
		else 
			cout << "bad command (load params op_choice)" << "\n";
	}
	
	p.rep_wheel.push_back(p.rt_rep);
	p.rep_wheel.push_back(p.rt_cross);
	p.rep_wheel.push_back(p.rt_mut);

	partial_sum(p.rep_wheel.begin(), p.rep_wheel.end(), p.rep_wheel.begin());
	
	if(!p.seeds.empty()) // get seed stacks
	{
		p.op_choice.push_back(10); // include seeds in operation choices
		p.op_weight.push_back(1); // include opweight if used

		for (int i=0; i<p.seeds.size();++i)
		{
			p.seedstacks.push_back(vector<node>());

			Eqn2Line(p.seeds.at(i),p.seedstacks.at(i));
		}
	}
	//normalize fn weights
	if (p.weight_ops_on) 
	{
		float sumweight = accumulate(p.op_weight.begin(),p.op_weight.end(),0);
		for(unsigned int i=0;i<p.op_weight.size();++i)
                        p.op_weight.at(i) = p.op_weight.at(i)/sumweight;
	}
}
void load_data(data &d, std::ifstream& fs,params& p)
{
	string s;
	string varname;
	float tmpf;
	float tarf;
	//int tmpi;
	string tmps;
	//bool tmpb;
	//int varnum =0;
	int rownum=0;
	unsigned int index=0;
	bool pass = 0;
	// get variable names from first line / number of variables
	getline(fs,s,'\n');
	istringstream ss(s);
	ss>>tmps;

	while (ss>>tmps)
	{
		while (!pass && index<p.allvars.size())
		{
			if(tmps.compare(p.allvars.at(index))==0)
			{
				d.label.push_back(p.allvars.at(index)); // populate data labels from all vars based on column location
				pass=1;
				
			}
			else
			{
				++index;
			}
		}

		pass=0;
		index=0;
	}
	vector<int> shuffler;
    while(!fs.eof())
    {		
		getline(fs,s,'\n');
		istringstream ss2(s);
		// get target data
		ss2 >> tarf;
		/*if(boost::math::isnan(tarf)){
			cerr << "NaN value in data file. exiting..\n";
			exit(1);
		}*/
		d.target.push_back(tarf);
		d.vals.push_back(vector<float>());
		// get variable data
		while(ss2 >> tmpf)
		{
			d.vals[rownum].push_back(tmpf);
			//varnum++;
		}
		++rownum;
    }
	// pop end in case of extra blank lines in data file
	while(d.vals.back().empty())
	{
		d.vals.pop_back();
	}

	
	

	//d.dattovar.resize(p.allvars.size());
	//d.mapdata();
}
void load_lexdata(data &d, std::ifstream& fs,params& p)
{
	/*if (!fs.good())
	{
		cout << "BAD DATA FILE LOCATION" << "\n";
	}*/

	string s;
	string varname;
	float tmpf;
	float tarf;
	//int tmpi;
	string tmps;
	//bool tmpb;
	//int varnum =0;
	int rownum=0;
	int lexnum=0;
	unsigned int index=0;
	bool pass = 0;
	// get variable names first line / number of variables
	getline(fs,s,'\n');
	istringstream ss(s);
	ss>>tmps;

	while (ss>>tmps )
	{
		while (!pass && index<p.allvars.size())
		{
			if(tmps.compare(p.allvars.at(index))==0)
			{
				d.label.push_back(p.allvars.at(index)); // populate data labels from all vars based on column location
				pass=1;
				
			}
			else
			{
				++index;
			}
		}

		pass=0;
		index=0;
	}
	vector<int> shuffler;
	int c=0;
	d.lexvals.push_back(vector<vector<float>>());
	d.targetlex.push_back(vector<float>());
    while(!fs.eof())
    {		
		getline(fs,s,'\n');
		istringstream ss2(s);

		if (s.compare("case")==0 ){
			++c; //YAY!
			lexnum=0;
			d.lexvals.push_back(vector<vector<float>>());
			d.targetlex.push_back(vector<float>());
		}
		else{
			// get target data
			ss2 >> tarf;
			d.target.push_back(tarf);
			d.targetlex[c].push_back(tarf);
			d.vals.push_back(vector<float>());
			d.lexvals[c].push_back(vector<float>());
			// get variable data
			while(ss2 >> tmpf)
			{
				d.vals[rownum].push_back(tmpf);
				d.lexvals[c][lexnum].push_back(tmpf);
				//varnum++;
			}
			++rownum;
			++lexnum;
		}
    }
	// set number of cases
	p.numcases = c+1;
	// pop end in case of extra blank lines in data file
	while(d.vals.back().empty())
	{
		d.vals.pop_back();
	}

	//d.dattovar.resize(p.allvars.size());
	//d.mapdata();
}
int load_pop(vector<ind>& pop,params& p,state& s)
{
	// load population from file filename. 
	ifstream fs(p.pop_restart_path);
	if(!fs.is_open())
	{
		std::cerr << "Error: couldn't open population file " + p.pop_restart_path + ".\n";
		exit(1);
	}
	
	string s0;
	string varname;
	float tmpf;
	//int tmpi;
	string tmps;
	bool tmpb;

	//string trash;
	//s.erase();
    //s.reserve(is.rdbuf()->in_avail());
	int i=0;
	int j=0;
	while(!fs.eof() && i<p.popsize)
    {		
		getline(fs,s0,'\n');
		istringstream ss(s0);
		ss >> varname;

		if(varname.compare(0,6,"gline:") == 0)
		{
			//pop.push_back(ind());
			
			while (ss>>tmps){
				if (tmps.compare("+")==0 || tmps.compare("-")==0 || tmps.compare("*")==0 || tmps.compare("/")==0 || tmps.compare("s")==0 || tmps.compare("c")==0 || tmps.compare("e")==0 || tmps.compare("l")==0) // operator node
					pop[i].line.push_back(node(char(tmps[0])));
				else if (isdigit(tmps[0])) // constant node
					pop[i].line.push_back(node(std::stof(tmps)));
				else //variable node
					pop[i].line.push_back(node(tmps));

			}
			++i;
		}
		else if(p.eHC_on && varname.compare(0,6,"eline:") == 0)
		{
			while (ss>>tmpb){
				pop[i-1].line[j].on = tmpb;
				++j;
			}
			j=0;
		}
		/*else if(varname.compare(0,4,"age:") == 0)
			ss>>pop[i-1].age;*/
	}
	if (!fs.eof()) s.out << "WARNING: truncating population from file to first " + to_string(static_cast<long long>(p.popsize)) + " individuals...\n";
	return i; // size of loaded population
}

void shuffle_data(data& d, params& p, vector<Randclass>& r,state& s)
{
	vector<int> shuffler;
	vector<float> newtarget;
	vector<vector<float>> newvals;
	if (p.sel!=3) 
	{
		for(int i=0;i<d.vals.size();++i)
			shuffler.push_back(i);

		std::random_shuffle(shuffler.begin(),shuffler.end(),r[omp_get_thread_num()]);
		
		
		if (p.EstimateFitness){
			bool tmp = s.out.trials;
			s.out.trials=1; // keep output from going to console
			s.out << "data shuffle index: ";
			for (int i=0; i<shuffler.size(); ++i)
				s.out << shuffler.at(i) << " ";
			s.out << "\n";
			s.out.trials=tmp;
		}
		for(int i=0;i<d.vals.size();++i)
		{
			newtarget.push_back(d.target.at(shuffler.at(i)));
			newvals.push_back(d.vals.at(shuffler.at(i)));

		}
		using std::swap;
		std::swap(d.target,newtarget);
		std::swap(d.vals,newvals);

	}
	else // lexicase selection
	{
		
		vector<vector<float>> newtarlex(d.lexvals.size());
		vector<vector<vector<float>>> newlexvals(d.lexvals.size());
		for (int h=0; h<d.lexvals.size();++h)
		{
			shuffler.clear();
			for(int i=0;i<d.lexvals[h].size();++i)
				shuffler.push_back(i);

			std::random_shuffle(shuffler.begin(),shuffler.end(),r[omp_get_thread_num()]);
			
			for(int i=0;i<d.lexvals[h].size();++i)
			{
				//for (int j =0; j<d.lexvals[h][i].size();++j){
					newtarlex[h].push_back(d.targetlex[h][shuffler.at(i)]);
					newlexvals[h].push_back(d.lexvals[h][shuffler.at(i)]);

					newtarget.push_back(newtarlex[h].back());
					newvals.push_back(newlexvals[h].back());
				//}
			}
		}
		using std::swap;
		swap(d.target,newtarget);
		swap(d.vals,newvals);
		swap(d.targetlex,newtarlex);
		swap(d.lexvals,newlexvals);
	}
	
}
bool stopcondition(tribe& T,params p,data& d,state& s,FitnessEstimator& FE)
{
	if (!p.EstimateFitness){
		if (T.bestFit() <= 0.000001){
			s.out << "best fitness criterion achieved: " << T.bestFit() <<"\n";
			return true;
		}
		else
			return false;
	}
	else{
		vector<ind> best(1);
		T.getbestind(best[0]);
		
		p.EstimateFitness=0;
		Fitness(best,p,d,s,FE);
		p.EstimateFitness=1;
		if (best[0].fitness <= 0.000001){
			s.out << "best fitness criterion achieved: " << best[0].fitness <<"\n";
			return true;
		}
		else
			return false;
	}
}
int get_next_task(int& index,vector<int>& task_assignments)
{
	// 0: evolve solution 
	// 1: evolve fitness estimation
	// 2: print output (after all 0 tasks finish)
	// 3: update pareto archive (after all 0 tasks finish)
	if (index == task_assignments.size()-1)
		return -1;
	else{
		return task_assignments.at(++index);
	}
}
void runEllenGP(string paramfile, string datafile,bool trials,int trialnum)
{
	//string paramfile(param_in);
	//string datafile(data_in);
	/* steps:
	Initialize population
		make genotypes
		genotype to phenotype
		calculate fitness
		hill climb
	Next Generation
		Select Parents
		Create children genotypes
		genotype to phenotype
		calculate fitness
		hill climb
	Store Statistics
	Print Update
	
	INPUTS
	paramfile: parameter file
	datafile: data set: target in first column, dependent variables in second column
	*/
	struct params p; 
	struct data d;
	struct state s;

	vector <Randclass> r;
	
	// load parameter file
	ifstream fs(paramfile);
	if (!fs.is_open()){
		cerr << "Error: couldn't open parameter file " + paramfile << "\n";
		exit(1);
	}
	load_params(p, fs);
	// load data file
	ifstream ds(datafile);
	if (!ds.is_open()){
			cerr << "Error: couldn't open data file " + datafile << "\n";
			exit(1);
		}
	if (p.sel == 3) load_lexdata(d,ds,p);
	else load_data(d,ds,p);
	
	std::time_t t =  std::time(NULL);
    

#if defined(_WIN32)
	std::tm tm;
	localtime_s(&tm,&t);
	char tmplog[100];
	strftime(tmplog,100,"%Y-%m-%d_%H-%M-%S",&tm);
#else
	std::tm * tm = localtime(&t);
	char tmplog[100];
	strftime(tmplog,100,"%F_%H-%M-%S",tm);
#endif

	
   // string tmplog = "777";
	const char * c = p.resultspath.c_str();
	#if defined(_WIN32)
		_mkdir(c);
	#else
		mkdir(c, 0777); // notice that 777 is different than 0777
	#endif
	 int thrd = omp_get_thread_num();
	 string thread;
	 if (trials) thread = std::to_string(static_cast<long long>(trialnum));
	 else thread = std::to_string(static_cast<long long>(thrd));
#if defined(_WIN32)
	 string pname = paramfile.substr(paramfile.rfind('\\')+1,paramfile.size());
	 string dname = datafile.substr(datafile.rfind('\\')+1,datafile.size());
	 pname = pname.substr(0,pname.size()-4);
	 dname = dname.substr(0,dname.size()-4);
     string logname = p.resultspath + '\\' + "DevelEP_" + tmplog + "_" + pname + "_" + dname + "_" + thread + ".log";
#else
	 string pname = paramfile.substr(paramfile.rfind('/')+1,paramfile.size());
	 string dname = datafile.substr(datafile.rfind('/')+1,datafile.size());
	 pname = pname.substr(0,pname.size()-4);
	 dname = dname.substr(0,dname.size()-4);
	 string logname = p.resultspath + '/' + "DevelEP_" + tmplog + "_" + pname + "_" + dname + "_" + thread + ".log";
#endif



	 s.out.set(trials);
	 s.out.open(logname);
	 if (!s.out.is_open()){
		 cerr << "Write-to File " << logname << " did not open correctly.\n";
		 exit(1);
	 }
	 std::ofstream dfout;
	 initdatafile(dfout,logname);
	 s.out << "_______________________________________________________________________________ \n";
	 s.out << "                                    ellenGP                                     \n";
	 s.out << "_______________________________________________________________________________ \n";
	 //s.out << "Time right now is " << std::put_time(&tm, "%c %Z") << '\n';
	 s.out<< "Results Path: " << p.resultspath  << "\n";
	 s.out << "parameter file: " << paramfile << "\n";
	 s.out << "data file: " << datafile << "\n";
	 if(trials) {s.out << "Running in Trial Mode\n";
	 s.out << "Trial ID: " << trialnum << "\n";
	 }
	 s.out << "Settings: \n";
	 // get evolutionary method
	 s.out << "Evolutionary Method: ";
	 switch(p.sel){
	 case 1:
		 s.out << "Standard Tournament\n";
		 break;
	 case 2:
		 s.out << "Deterministic Crowding\n";
		 break;
	 case 3:
		 s.out << "Lexicase Selection\n";
		 break;
	 case 4:
		 s.out << "Age-Fitness Pareto\n";
		 break;
	 }
	 s.out << "ERC: " << p.ERC << "\n";
	 s.out << "Parameter Hill Climber: " << p.pHC_on <<"\n";
	 s.out << "Epigenetic Hill Climber: " << p.eHC_on <<"\n";
	 if(p.train) s.out << "Data split " << p.train_pct << "/" << 1-p.train_pct << " for training and validation.\n";
	 s.out << "Total Population Size: " << p.popsize << "\n";
	 if (p.limit_evals) s.out << "Maximum Point Evals: " << p.max_evals << "\n";
	 else s.out << "Maximum Generations: " << p.g << "\n";
	 s.out << "Number of log points: " << p.num_log_pts << " (0 means log all points)\n";
	 if (trials && p.islands){
		s.out << "WARNING: cannot run island populations in trial mode. This trial will run on one core.\n";
		p.islands = false;
	 }
	
	int nt=0;
	//int ntt=0;	

	#pragma omp parallel 
	{
		nt = omp_get_num_threads();
	}
	s.out << "Number of threads: " << nt << "\n";
	//s.out << "OMP Number of threads: " << omp_get_num_threads() << "\n";
	p.nt = nt;

	//initialize random number generator
	unsigned int seed1 = int(time(NULL));
	s.out << "seeds: \n";
	r.resize(omp_get_max_threads());
	#pragma omp parallel 
	{
			//cout << "seeder: " << seeder <<endl;
			//cout << "seed1: " << seed1*seeder <<endl;
			if(!trials){
				s.out << to_string(static_cast<long long>(seed1*(omp_get_thread_num()+1))) + "\n";
				//r.at(seeder).SetSeed(seed1*(seeder+1));
				r.at(omp_get_thread_num()).SetSeed(seed1*(omp_get_thread_num()+1));
			}
			else
				r.at(omp_get_thread_num()).SetSeed(seed1*(omp_get_thread_num()+1)*trialnum);
	}
	if(trials)
		s.out << (omp_get_thread_num()+1)*seed1*trialnum << "\n";

	//shuffle data for training
	if (p.shuffle_data) 
		shuffle_data(d,p,r,s);
	
	boost::timer time;
	paretoarchive A;
	if (p.prto_arch_on){
		A.resize(p.prto_arch_size);
		tribe FinalArchive(p.prto_arch_size,p.max_fit,p.min_fit);
		FinalArchive.pop = A.pop;
	}

	vector<FitnessEstimator> FE(1);
	vector<ind> trainers;

	if (p.islands)
	{
		//p.parallel=false;
		// determine number of threads
		
		//int num_islands=omp_get_max_threads(); //
		int num_islands=nt;
		int subpops = p.popsize/num_islands;
		vector<tribe> T;
		tribe World(subpops*num_islands,p.max_fit,p.min_fit); //total population of tribes
		s.out << num_islands << " islands of " << subpops << " individuals, total pop " << subpops*num_islands <<"\n";
		 

		for(int i=0;i<num_islands;++i)
			T.push_back(tribe(subpops,p.max_fit,p.min_fit));
		// run separate islands 
		if (p.init_validate_on)
		{
			s.out << "Initial validation..."; 
			
			if (p.EstimateFitness)
			{
				#pragma omp parallel for
				for(int i=0;i<num_islands;++i)
					InitPop(T.at(i).pop,p,r);
				
				// construct world population
				for(int j=0;j<T.size();++j){
					for(int k=0;k<T[0].pop.size();++k){
						World.pop.at(j*T[0].pop.size()+k)=T.at(j).pop.at(k);
						//makenew(World.pop[j*T[0].pop.size()+k]);	
					}
				}
				// initialize fitness estimation pop
				InitPopFE(FE,World.pop,trainers,p,r,d,s);

				//discard invalid individuals
				#pragma omp parallel for
				for(int i=0;i<num_islands;++i){
					float worstfit;
					float bestfit;
					vector<ind> tmppop;

					Fitness(T.at(i).pop,p,d,s,FE[0]);
					worstfit = T.at(i).worstFit();
					bestfit = T.at(i).bestFit();

					
					int counter=0;
					while(worstfit == p.max_fit && counter<100)
					{
						for (vector<ind>::iterator j=T.at(i).pop.begin();j!=T.at(i).pop.end();)
						{
							if ( (*j).fitness == p.max_fit)
							{
								j=T.at(i).pop.erase(j);
								tmppop.push_back(ind());
							}
							else
								++j;
						}

						InitPop(tmppop,p,r);
						Fitness(tmppop,p,d,s,FE[0]);
						T.at(i).pop.insert(T.at(i).pop.end(),tmppop.begin(),tmppop.end());
						tmppop.clear();
						worstfit = T.at(i).worstFit();
						counter++;
						if(counter==100)
							s.out << "initial population count exceeded. Starting evolution...\n";
					}
				}
				s.setgenevals();
				s.out << " number of evals: " << s.getgenevals() << "\n";
			}
			else{
				
				#pragma omp parallel for 
				for(int i=0;i<num_islands;++i)
				{
					float worstfit;
					float bestfit;
					vector<ind> tmppop;
					// s.out << "Initialize Population..." << "\n";
					InitPop(T.at(i).pop,p,r);
					// s.out << "Fitness..." << "\n";
					Fitness(T.at(i).pop,p,d,s,FE[0]);
					worstfit = T.at(i).worstFit();
					bestfit = T.at(i).bestFit();
					int counter=0;
					while(worstfit == p.max_fit && counter<100)
					{
						for (vector<ind>::iterator j=T.at(i).pop.begin();j!=T.at(i).pop.end();)
						{
							if ( (*j).fitness == p.max_fit)
							{
								j=T.at(i).pop.erase(j);
								tmppop.push_back(ind());
							}
							else
								++j;
						}

						InitPop(tmppop,p,r);
						Fitness(tmppop,p,d,s,FE[0]);
						T.at(i).pop.insert(T.at(i).pop.end(),tmppop.begin(),tmppop.end());
						tmppop.clear();
						worstfit = T.at(i).worstFit();
						counter++;
						if(counter==100)
							s.out << "initial population count exceeded. Starting evolution...\n";
					}
				
				}
			
				s.setgenevals();
				s.out << " number of evals: " << s.getgenevals() << "\n";
			}
		}
		else // normal population initialization
		{
			/*bool tmp = p.EstimateFitness;
			p.EstimateFitness=0;*/
			#pragma omp parallel for 
			for(int i=0;i<num_islands;++i)
			{
				InitPop(T.at(i).pop,p,r);
				
				// s.out << "Gen 2 Phen..." << "\n";
				// s.out << "Fitness..." << "\n";
				if(!p.EstimateFitness)
					Fitness(T.at(i).pop,p,d,s,FE[0]);
			}
			// construct world population
			for(int j=0;j<T.size();++j){
				for(int k=0;k<T[0].pop.size();++k){
					World.pop.at(j*T[0].pop.size()+k)=T.at(j).pop.at(k);
					//makenew(World.pop[j*T[0].pop.size()+k]);	
				}
			}
			if (p.EstimateFitness){
				InitPopFE(FE,World.pop,trainers,p,r,d,s);
				#pragma omp parallel for 
				for(int i=0;i<num_islands;++i)
					Fitness(T.at(i).pop,p,d,s,FE[0]);
			}
		}
		
		// use tmpFE for updating in parallel
		vector<FitnessEstimator> tmpFE = FE; 

		int gen=0;
		long long termits,term, print_trigger;
		if (p.limit_evals){
			termits = s.totalptevals();
			term = p.max_evals;
			if (p.num_log_pts ==0) print_trigger=0;
			else print_trigger = p.max_evals/p.num_log_pts;
		
		}
		else {
			print_trigger=0;
			termits=1;
			term = p.g;
		}
		bool pass=1;
		int mixtrigger=p.island_gens*(p.popsize+p.popsize*p.eHC_on+p.popsize*p.pHC_on);
		//int trainer_trigger=0;
		int trainer_trigger=p.FE_train_gens;//p.FE_train_gens*(p.popsize+p.popsize*p.eHC_on+p.popsize*p.pHC_on);
		
		bool migrate=false;
		// while(gen<=p.g && !stopcondition(World.best))
		// {
			int q;
//			int task_num;
			//int index=-1;
//			int cntr;
			/*
			vector<int> task_status;
			vector<int> task_assignments(num_islands,0);
			task_assignments.push_back(1);	*/

		//#pragma omp parallel private(q, task_num, cntr) shared(pass, index)
		//{
		//	q = omp_get_thread_num();

		//	while(gen<=p.g && pass)
		//	{
		//		#pragma omp critical
		//		{
		//			task_num = get_next_task(index,task_assignments);
		//		}
		//		while(task_num!=-1){
		//			switch(task_num){
		//			case 0: //evolve solution population
		//				if(pass){			
		//					Generation(T[q].pop,p,r,d,s,FE[0]);	

		//					if (stopcondition(T[q],p,d,s,FE[0]))
		//						pass=0;
		//				}

		//				if (pass) {
		//					if (p.pHC_on && p.ERC){
		//						for(int k=0; k<T[q].pop.size(); ++k)
		//							HillClimb(T[q].pop.at(k),p,r,d,s,FE[0]);
		//					}
		//					if (p.eHC_on){
		//						for(int m=0; m<T[q].pop.size(); m++)
		//							EpiHC(T[q].pop.at(m),p,r,d,s,FE[0]);
		//					}
		//	
		//					if (stopcondition(T[q],p,d,s,FE[0]))
		//						pass=0;
		//				}

		//				// construct world population
		//				cntr=0;
		//				for(int k=q*subpops;k<(q+1)*subpops;++k){
		//					World.pop.at(k)=T[q].pop.at(cntr);
		//					makenew(World.pop.at(k));
		//					cntr++;
		//				}
		//				break;
		//			case 1: //evolve fitness estimator
		//				s.out << "Evolving fitness estimators...\n";
		//				EvolveFE(World.pop,tmpFE,trainers,p,d,s,r);
		//								
		//				break;
			//		case 2: // print out 
			//			printstatsP(World,gen,s,p,A);
			//			if (p.printeverypop) printpop(World.pop,p,s,logname,2);
			//			s.out << "Total Time: " << (int)floor(time.elapsed()/3600) << " hr " << ((int)time.elapsed() % 3600)/60 << " min " << (int)time.elapsed() % 60 << " s\n";
			//			s.out << "Total Evals: " << s.totalevals() << "\n";
			//			s.out << "Point Evals: " << s.totalptevals() << "\n";
			//			s.out << "Average evals per second: " << (float)s.totalevals()/time.elapsed() << "\n";
			//			s.out << "Average point evals per second: " << (float)s.totalptevals()/time.elapsed() << "\n";
			//			break;
			//			}// switch(task_num)
			//			#pragma omp critical
			//			{
			//				task_num = get_next_task(index,task_assignments);
			//			}
		//		} // while(task_num!=-1)
		//		#pragma omp barrier
		//		#pragma omp single  // mix island populations at regular intervals
		//		{
		//			s.setgenevals();
		//			if(s.totalevals()>mixtrigger) 
		//			{
		//				//shuffle population	
		//				std::random_shuffle(World.pop.begin(),World.pop.end(),r[q]);
		//				//redistribute populations to islands
		//				s.out << "Shuffling island populations...\n";
		//				migrate = true;
		//				mixtrigger+=p.island_gens*(p.popsize+p.popsize*p.eHC_on+p.popsize*p.pHC_on);
		//			}
		//			else
		//				migrate=false;
		//		}
		//				
		//		#pragma omp single // pick new trainers for fitness estimation at regular intervals
		//		{
		//			if (p.EstimateFitness){
		//				FE.assign(tmpFE.begin(),tmpFE.end());
		//				if(s.totalevals()>trainer_trigger) {
		//					s.out << "Picking trainers...\n";
		//					PickTrainers(World.pop,FE,trainers,p,d,s);
		//					trainer_trigger+=p.island_gens*(p.popsize+p.popsize*p.eHC_on+p.popsize*p.pHC_on);
		//					//trainer_trigger=0;
		//				}
		//			}
		//		}
		//		#pragma omp single // print status each generation
		//		{
		//			//update pareto archive
		//			A.update(World.pop); 
		//			printpop(A.pop,p,s,logname,1);
		//			
		//			ge++n;
		//			if (p.EstimateFitness)
		//			{
		//				s.out << "Best FE fit: " << tmpFE[0].fitness <<"\n";
		//				s.out << "Fitness Estimators:\n";
		//				for (int a=0;a<tmpFE.size();a++){
		//					for (int b=0;b<tmpFE[a].FEpts.size();b++)
		//						s.out << tmpFE[a].FEpts[b] << " ";
		//					s.out << "\n";		
		//				}
		//			}	
		//		}

		//		if (migrate)					
		//			T[q].pop.assign(World.pop.begin()+q*subpops,World.pop.begin()+(q+1)*subpops);
				

		//		if (gen>p.g) pass=0;
		//	} // while gen<=p.g
		//} // pragma omp parallel
		//print initial population
		// construct world population
		int cntr=0;
		for (int q0 = 0; q0<nt;++q0){ 
			for(int k=q0*subpops;k<(q0+1)*subpops;++k){
				World.pop.at(k)=T.at(q0).pop.at(cntr);
				//makenew(World.pop.at(k));
				++cntr;
			}
			cntr=0;
		}
		if (p.print_init_pop) printpop(World.pop,p,s,logname,3);
		if (p.print_genome) printGenome(World,0,logname,d,p);
		#pragma omp parallel private(q) shared(pass)
		{
		
			q = omp_get_thread_num();

			while(termits<=term && pass)
			{
			

				if(pass){			
					Generation(T[q].pop,p,r,d,s,FE[0]);	

					if (stopcondition(T[q],p,d,s,FE[0]))
						pass=0;
				}

				if (pass) {
					if (p.pHC_on && p.ERC)
					{
							for(int k=0; k<T[q].pop.size(); ++k)
								HillClimb(T[q].pop.at(k),p,r,d,s,FE[0]);
					}
					if (p.eHC_on && !p.eHC_mut) 
					{
							for(int m=0; m<T[q].pop.size(); m++)
								EpiHC(T[q].pop.at(m),p,r,d,s,FE[0]);
					}
			
					if (stopcondition(T[q],p,d,s,FE[0]))
						pass=0;
				}

			
					// construct world population
					int cntr=0;
					//std::vector<ind>::iterator it = T[q].begin();
					//World.pop.assign(
					for(int k=q*subpops;k<(q+1)*subpops;++k){
						World.pop.at(k)=T[q].pop.at(cntr);
						//makenew(World.pop.at(k));
						cntr++;
					}
					
					#pragma omp barrier
					
					#pragma omp single  nowait //coevolve fitness estimators
					{
					
						if (p.EstimateFitness){

						
							float aveFEfit=0;
							for (int u=0;u<FE.size();u++)
								aveFEfit+=FE[u].fitness;
							aveFEfit /= FE.size();
							if (aveFEfit==0)
								std::random_shuffle(tmpFE.begin(),tmpFE.end(),r[omp_get_thread_num()]);
							else
								EvolveFE(World.pop,tmpFE,trainers,p,d,s,r);
							if (!p.limit_evals || s.totalptevals() >= print_trigger){
								s.out << "Evolving fitness estimators...\n";
								s.out << "Best FE fit: " << FE[0].fitness <<"\n";
								if (p.estimate_generality) s.out << "Best FE genty: " << FE[0].genty <<"\n";
								s.out << "Ave FE fit: " << aveFEfit << "\n";
								s.out << "Current Fitness Estimator:\n";
						
								for (int b=0;b<FE[0].FEpts.size();b++)
									s.out << FE[0].FEpts[b] << " ";
								s.out << "\n";
							}
						
						}
					}
					#pragma omp single 
					{
						s.setgenevals();
						if(s.totalevals()>mixtrigger) 
						{
							//shuffle population	
							std::random_shuffle(World.pop.begin(),World.pop.end(),r[q]);
							//redistribute populations to islands
							s.out << "Shuffling island populations...\n";
							migrate = true;
							mixtrigger+=p.island_gens*(p.popsize+p.popsize*p.eHC_on+p.popsize*p.pHC_on);
						}
						else
							migrate=false;
					}
					#pragma omp single 
					{
						if(p.prto_arch_on){
							A.update(World.pop);
							printpop(A.pop,p,s,logname,1);
						}
						if (!p.limit_evals || s.totalptevals() >= print_trigger){
							printdatafile(World,s,p,r,dfout,gen);
							if (p.printeverypop) printpop(World.pop,p,s,logname,2);
							if (p.print_log) {
								printstats(World,gen,s,p,A);
								s.out << "Total Time: " << (int)floor(time.elapsed()/3600) << " hr " << ((int)time.elapsed() % 3600)/60 << " min " << (int)time.elapsed() % 60 << " s\n";
								s.out << "Total Evals: " << s.totalevals() << "\n";
								s.out << "Point Evals: " << s.totalptevals() << "\n";
								s.out << "Average evals per second: " << (float)s.totalevals()/time.elapsed() << "\n";
								s.out << "Average point evals per second: " << (float)s.totalptevals()/time.elapsed() << "\n";
							}
							if (print_trigger!=0) print_trigger += p.max_evals/p.num_log_pts;
							if (p.print_genome) printGenome(World,gen,logname,d,p);
						}
						++gen;
						if (p.limit_evals) termits = s.totalptevals();
						else ++termits;
					
						if (!p.EstimateFitness && p.estimate_generality && p.G_shuffle)
							shuffle_data(d,p,r,s);	
									
					}
					
					#pragma omp single
					{
						if (p.EstimateFitness){
							// assign tmpFE to FE
							FE.assign(tmpFE.begin(),tmpFE.end());
							/*float avefit=0;
							for (int u=0;u<FE.size();u++)
								avefit+=FE[u].fitness;
							avefit /= FE.size();
							if (avefit>1)
								cout <<"avefit error\n";*/
							if(gen>trainer_trigger) { //pick new trainers when FE pop has converged or when it has been enough generations
							//if(gen>trainer_trigger) {
								s.out << "Picking trainers...\n";
								PickTrainers(World.pop,FE,trainers,p,d,s);
								trainer_trigger = gen + p.FE_train_gens; //p.island_gens*(p.popsize+p.popsize*p.eHC_on+p.popsize*p.pHC_on);
								//trainer_trigger=0;
							}
						}
					}

					if (migrate)					
						T[q].pop.assign(World.pop.begin()+q*subpops,World.pop.begin()+(q+1)*subpops);

					//if (gen>p.g) pass=0;
								
			
			}  s.out << "exited while loop...\n";
		} s.out << "exited parallel region ...\n";
		

		if (p.EstimateFitness){// assign real fitness values to final population and archive
			p.EstimateFitness=0;
			Fitness(World.pop,p,d,s,FE[0]);
			if (p.prto_arch_on) Fitness(A.pop,p,d,s,FE[0]);
			p.EstimateFitness=1;
		}
		printdatafile(World,s,p,r,dfout,gen);
		printbestind(World,p,s,logname);
		printpop(World.pop,p,s,logname,0);
		if (p.prto_arch_on)
			printpop(A.pop,p,s,logname,1);
	}
	else //no islands
	{
		tribe T(p.popsize,p.max_fit,p.min_fit);
		if (p.pop_restart) // initialize population from file
		{
			s.out << "loading pop from " + p.pop_restart_path + "...\n";
			int tmp = load_pop(T.pop,p,s);
			if (tmp < p.popsize){
				s.out << "WARNING: population loaded from file is smaller than set pop size. Randomly initiated individuals will be added...\n";
				vector<ind> tmppop(p.popsize-tmp);
				InitPop(tmppop,p,r);
				swap_ranges(T.pop.end()-(p.popsize-tmp),T.pop.end(),tmppop.begin());
			}
			Fitness(T.pop,p,d,s,FE[0]);
		}
		else{
			if (p.init_validate_on)
			{
				s.out << "Initial validation..."; 
				bool tmp = p.EstimateFitness;
				p.EstimateFitness=0;
				float worstfit;
				int cnt=0;
				//float bestfit;
				vector<ind> tmppop;
				// s.out << "Initialize Population..." << "\n";
				InitPop(T.pop,p,r);
				assert (T.pop.size()== p.popsize) ;
				// s.out << "Gen 2 Phen..." << "\n";
				// s.out << "Fitness..." << "\n";
				Fitness(T.pop,p,d,s,FE[0]);
				assert (T.pop.size()== p.popsize) ;
			
				worstfit = T.worstFit();
				while(worstfit == p.max_fit && cnt<100)
				{
					for (vector<ind>::iterator j=T.pop.begin();j!=T.pop.end();)
					{
						if ( (*j).fitness == p.max_fit)
						{
							j=T.pop.erase(j);
							tmppop.push_back(ind());
						}
						else
							++j;
					}
				
					s.out << "\ntmppop size: " << tmppop.size();
					InitPop(tmppop,p,r);
					Fitness(tmppop,p,d,s,FE[0]);
				
					T.pop.insert(T.pop.end(),tmppop.begin(),tmppop.end());
					tmppop.clear();
					worstfit = T.worstFit();
					cnt++;
					if(cnt==100)
						s.out << "initial population count exceeded. Starting evolution...\n";
				}
				p.EstimateFitness=tmp;
			}
			else // normal population initialization
			{
				InitPop(T.pop,p,r);	

				bool tmp = p.EstimateFitness;
				p.EstimateFitness=0;
				Fitness(T.pop,p,d,s,FE[0]);
				p.EstimateFitness=tmp;
			}
		}
		s.setgenevals();
		s.out << " number of evals: " << s.getgenevals() << "\n";
		int its = 1;
		long long termits;
		int gits=1;
		
		if (p.limit_evals) termits = s.totalptevals();
		else termits=1;
		int trigger=0;
		int trainer_trigger=p.FE_train_gens;//*(p.popsize+p.popsize*p.eHC_on+p.popsize*p.pHC_on);
		

		int gen=0;
		int counter=0;
		//if(p.sel==2) // if using deterministic crowding, increase gen size
		//{
		//	gen = p.g*(p.popsize*(p.rt_mut+p.rt_rep) + p.popsize*p.rt_cross/2);
		//	trigger = p.popsize*(p.rt_mut+p.rt_rep)+p.popsize*p.rt_cross/2;
		//}
		//else
			gen=p.g;
		long long term, print_trigger;
		bool printed=false;
		if (p.limit_evals) {
			term = p.max_evals;
			if (p.num_log_pts ==0) print_trigger=0;
			else print_trigger = p.max_evals/p.num_log_pts;
		}
		else{
			print_trigger=0;
			term = gen;
		}
		if (p.EstimateFitness)
			InitPopFE(FE,T.pop,trainers,p,r,d,s);
		float etmp;
		//print initial population
		if (p.print_init_pop) printpop(T.pop,p,s,logname,3);
		if (p.print_genome) printGenome(T,0,logname,d,p);

		while (termits<=term && !stopcondition(T,p,d,s,FE[0]))
		{
			
			 etmp = s.numevals[omp_get_thread_num()];

			 if (!p.EstimateFitness && p.estimate_generality && p.G_shuffle)
				shuffle_data(d,p,r,s);
			 assert (T.pop.size() == p.popsize);
			 Generation(T.pop,p,r,d,s,FE[0]);
			 assert (T.pop.size()== p.popsize);
			 //s.out << "Generation evals = " + to_string(static_cast<long long>(s.numevals[omp_get_thread_num()]-etmp)) + "\n";
			 

			 if (its>trigger)
			 {
				 if (p.pHC_on && p.ERC)
				 {
					 etmp = s.numevals[omp_get_thread_num()];
					//#pragma omp parallel for
		 			for(int k=0; k<T.pop.size(); ++k)
		 				HillClimb(T.pop.at(k),p,r,d,s,FE[0]);
		 			 //s.out << "Hill climb evals = " + to_string(static_cast<long long>(s.numevals[omp_get_thread_num()]-etmp)) + "\n";	
					

		 		 }
				
				 if (p.eHC_on&& !p.eHC_mut) 
				 {
					 etmp = s.numevals[omp_get_thread_num()];
					// boost::progress_timer tm1;
					//#pragma omp parallel for
					for(int m=0; m<T.pop.size(); m++)
						EpiHC(T.pop.at(m),p,r,d,s,FE[0]);
					 //s.out << "EHC evals = " + to_string(static_cast<long long>(s.numevals[omp_get_thread_num()]-etmp)) + "\n";
				 } 

				s.setgenevals();
				//s.out << "Elapsed time: \n";
				if (p.prto_arch_on){
					A.update(T.pop);
					printpop(A.pop,p,s,logname,1);
				}
				if (!p.limit_evals || s.totalptevals() >= print_trigger){ 
					printdatafile(T,s,p,r,dfout,counter);
					if (p.printeverypop) printpop(T.pop,p,s,logname,2);
					if (p.print_log){
						printstats(T,counter,s,p,A);
						s.out << "Total Time: " << (int)floor(time.elapsed()/3600) << " hr " << ((int)time.elapsed() % 3600)/60 << " min " << (int)time.elapsed() % 60 << " s\n";
						s.out << "Total Evals: " << s.totalevals() << "\n";
						s.out << "Point Evals: " << s.totalptevals() << "\n";
						s.out << "Average evals per second: " << (float)s.totalevals()/time.elapsed() << "\n";
						s.out << "Average point evals per second: " << (float)s.totalptevals()/time.elapsed() << "\n";
					}
					if (print_trigger!=0) print_trigger += p.max_evals/p.num_log_pts;
					if (p.print_genome) printGenome(T,counter,logname,d,p);
					printed=true;
				}

				//if (p.sel==2)
					//trigger+=p.popsize*(p.rt_mut+p.rt_rep)+p.popsize*p.rt_cross/2;
				++counter;
		
			 }
			
			if (p.EstimateFitness){
				EvolveFE(T.pop,FE,trainers,p,d,s,r);
				if (!p.limit_evals || printed){ 
					s.out << "Evolving fitness estimators...\n";
					s.out << "Best FE fit: " << FE[0].fitness <<"\n";
					if (p.estimate_generality) s.out << "Best FE genty: " << FE[0].genty <<"\n";
					s.out << "Current Fitness Estimator:\n";
					
					for (int b=0;b<FE[0].FEpts.size();b++)
						s.out << FE[0].FEpts[b] << " ";
					s.out << "\n";
				}
			}
			
				
			if (p.EstimateFitness){
				if(counter>trainer_trigger) {
					s.out << "Picking trainers...\n";
					PickTrainers(T.pop,FE,trainers,p,d,s);
					trainer_trigger= counter + p.FE_train_gens; //*(p.popsize+p.popsize*p.eHC_on+p.popsize*p.pHC_on);
					//trainer_trigger=0;
				}
			}
				

			if (p.limit_evals) termits = s.totalptevals(); 
			else termits++;
			its++;
			printed=false;
		}
		

		if (p.EstimateFitness){// assign real fitness values to final population and archive
			p.EstimateFitness=0;
			Fitness(T.pop,p,d,s,FE[0]);
			if (p.prto_arch_on) Fitness(A.pop,p,d,s,FE[0]);
			p.EstimateFitness=1;
		}

		printdatafile(T,s,p,r,dfout,counter);
		printbestind(T,p,s,logname);
		printpop(T.pop,p,s,logname,0);
		if (p.prto_arch_on)
			printpop(A.pop,p,s,logname,1);
	}

	
	s.out << "\n Program finished sucessfully.\n";
	



}
