// Devilep.cpp : Defines the entry point for the console application.
// This is the main program for running Developmental Linear Epigenetic Programming, or Devilep. 

#include "stdafx.h"
// mine
#include "runDevelep.h"
#include <exception>
//#include "logger.h"
//#include "pop.h"
//#include "data.h"
//#include "InitPop.h"
//#include "Fitness.h"
//#include "stdafx.h"
//#include "state.h"
//#include "Generation.h"
//#include "HillClimb.h"
//#include "instructionset.h"
//#include "evaluator.h"
//#include "strdist.h"
//#include <ctime>
//#include <cstring>

using namespace std;

//// global parameters structure
//struct params p; 
//struct data d;
//struct state s;
//class logger fcout;
////struct evaluator e;
//vector <Randclass> r;


//void load_params(params &p, std::ifstream& is);
//void load_data(data &d, std::ifstream& is);
//bool stopcondition(float&);
//void printstats(tribe&,int&,int&);

//template<typename T>
//int _tmain(int argc, _TCHAR *argv[])
int main(int argc, char** argv)
{
	try 
	{
		string paramfile(argv[1]);
		string datafile(argv[2]);
		runDevelep(paramfile,datafile,0);
	}
	catch(exception& er) 
	{
		cout << "Error: " << er.what() << endl;

	}
	catch(...)
	{
		cout << "Exception Occurred."<<endl;
	}
	return 0;
}
//void runDevelep(string& paramfile, string& datafile)
//{
//	/* steps:
//	Initialize population
//		make genotypes
//		genotype to phenotype
//		calculate fitness
//		hill climb
//	Next Generation
//		Select Parents
//		Create children genotypes
//		genotype to phenotype
//		calculate fitness
//		hill climb
//	Store Statistics
//	Print Update
//	
//	INPUTS
//	paramfile: parameter file
//	datafile: data set: target in first column, dependent variables in second column
//	*/
//	// load parameter file
//	ifstream fs(paramfile);
//	load_params(p, fs);
//	// load data file
//	ifstream ds(datafile);
//	load_data(d,ds);
//
//	std::time_t t =  std::time(NULL);
//    std::tm tm    = *std::localtime(&t);
//	char tmplog[100];
//	strftime(tmplog,100,"%Y-%m-%d_%H-%M-%S",&tm);
//	const char * c = p.resultspath.c_str();
//	_mkdir(c);
//    string logname = p.resultspath + '\\' + "Develep_" + tmplog + ".txt";
//	fcout.open(logname);
//	fcout << "_________________________________________________ \n";
//	fcout << "              Develep v 1.0                      \n";
//	fcout << "_________________________________________________ \n";
//	fcout << "Time right now is " << std::put_time(&tm, "%c %Z") << '\n';
//	//fcout<< "Results Path: " << logname  << "\n";
//	fcout << "parameter file: " << paramfile << "\n";
//	fcout << "data file: " << datafile << "\n";
//	fcout << "Settings: \n";
//	fcout << "ERC: " << p.ERC << "\n";
//	fcout << "Parameter Hill Climber: " << p.pHC_on <<"\n";
//	fcout << "Epigenetic Hill Climber: " << p.eHC_on <<"\n";
//	fcout << "Total Population Size: " << p.popsize << "\n";
//	fcout << "Max Generations: " << p.g << "\n";
//	//fcout.setfile(logname);
//	//vector <ind> pop(p.popsize);
//	//ofstream jtbs(logname);
//	//fcout << "Hello Saved World" << "\n";
//	//fcout << "Why won't you write? " << "\n";
//	//jtbs << "Am i nuts or what?\n";
//	// setup equation parser 
//	/*if(p.parallel)
//	{
//		e.resize(omp_get_max_threads());
//		e.set_vars();
//	}
//	else
//	{
//		e.resize(1);
//		e.set_vars();
//	}*/
//	//struct evaluator e;
//
//
//	//initialize random number generator
//	unsigned seed1 = time(0);
//	if(p.parallel)
//	{
//		r.resize(omp_get_max_threads());
//		for(int seeder=0;seeder<omp_get_max_threads();seeder++)
//			{
//				seed1=seed1*1.987/7.891;
//				r.at(seeder).SetSeed(seed1);
//		}
//	}
//	else
//	{
//		r.resize(1);
//		r[0].SetSeed(seed1);
//	}
//
//	
//	boost::timer time;
//	
//	if (p.islands)
//	{
//		p.parallel=false;
//		
//		int num_islands=omp_get_max_threads();
//		int subpops = p.popsize/num_islands;
//		vector<tribe> T;
//		tribe World(subpops*num_islands); //total population of tribes
//		if ( p.loud ) fcout << num_islands << " islands of " << subpops << " individuals, total pop " << subpops*num_islands <<"\n";
//
//		for(int i=0;i<num_islands;i++)
//			 T.push_back(tribe(subpops));
//		// run separate islands 
//		if (p.init_validate_on)
//		{
//			fcout << "Initial validation..."; 
//			#pragma omp parallel for 
//			for(int i=0;i<num_islands;i++)
//			{
//				float worstfit;
//				vector<ind> tmppop;
//				//fcout << "Initialize Population..." << "\n";
//				InitPop(T.at(i).pop);
//				//fcout << "Gen 2 Phen..." << "\n";
//				Gen2Phen(T.at(i).pop);
//				//fcout << "Fitness..." << "\n";
//				Fitness(T.at(i).pop);
//				worstfit = T.at(i).worstFit();
//				while(worstfit == p.max_fit)
//				{
//					for (int j=0;j<T.at(i).pop.size(); j++)
//					{
//						if ( T.at(i).pop.at(j).fitness == p.max_fit)
//						{
//							T.at(i).pop.erase(T.at(i).pop.begin()+j);
//							tmppop.push_back(ind());
//						}
//					}
//					InitPop(tmppop);
//					Gen2Phen(tmppop);
//					Fitness(tmppop);
//					T.at(i).pop.insert(T.at(i).pop.end(),tmppop.begin(),tmppop.end());
//					tmppop.clear();
//					worstfit = T.at(i).worstFit();
//				}
//				
//			}
//			s.setgenevals();
//			fcout << " number of evals: " << s.getgenevals() << "\n";
//			//totalevals=0;
//			
//		}
//		else
//		{
//			#pragma omp parallel for
//			for(int i=0;i<num_islands;i++)
//			{
//				InitPop(T.at(i).pop);
//				//fcout << "Gen 2 Phen..." << "\n";
//				Gen2Phen(T.at(i).pop);
//				//fcout << "Fitness..." << "\n";
//				Fitness(T.at(i).pop);
//			}
//		}
//		int gen=0;
//		bool pass=1;
//		int run=0;
//		while(gen<=p.g && !stopcondition(World.best))
//		{
//			
//			#pragma omp parallel for
//			for(int i=0;i<num_islands;i++)
//			{
//				if(pass)
//				{
//					for(int j=0;j<p.island_gens;j++)
//					{
//						if(pass)
//						{
//							Generation(T.at(i).pop);
//							run=1;
//							if (stopcondition(T.at(i).best))
//							{
//								pass=0;
//							}
//						}
//
//						//cout << "thread " << omp_get_thread_num() << " numevals: " << s[i].numevals << "\n";
//					}
//				}
//
//			// construct world population
//				swap_ranges(World.pop.begin()+i*subpops,World.pop.begin()+(i+1)*subpops,T.at(i).pop.begin());
//			}
//
//				
//
//
//		//shuffle population	
//		std::random_shuffle(World.pop.begin(),World.pop.end(),r[omp_get_thread_num()]);
//		//redistribute populations to islands
//		#pragma omp parallel for
//		for(int i=0;i<num_islands;i++)
//		{
//			T.at(i).pop.assign(World.pop.begin()+i*subpops,World.pop.begin()+(i+1)*subpops);
//		}
//		s.setgenevals();
//		printstats(World,gen,s.genevals.back());
//		fcout << "Total Time: " << (int)floor(time.elapsed()/3600) << " hr " << ((int)time.elapsed() % 3600)/60 << " min " << (int)time.elapsed() % 60 << " s\n";
//		fcout << "Total Evals: " << s.totalevals() << "\n";
//		fcout << "Evals per second: " << (float)s.totalevals()/time.elapsed() << "\n";
//		//World.pop.clear();
//		gen++;
//		}
//	}
//	else
//	{
//		tribe T(p.popsize);
//		if (p.init_validate_on)
//		{//while loop 
//			fcout << "Initialize Population..." << "\n";
//			InitPop(T.pop);
//			fcout << "Gen 2 Phen..." << "\n";
//			Gen2Phen(T.pop);
//			fcout << "Fitness..." << "\n";
//			Fitness(T.pop);
//		}
//		int i=1;
//		while (i<=p.g && !stopcondition(T.best))
//		{
//			
//			fcout << "Elapsed time: \n";
//			Generation(T.pop);
//			fcout << "     print stats...\n";
//			printstats(T,i,s.genevals.back());
//
//			i++;
//			
//			fcout << "Total Time: " << (int)floor(time.elapsed()/3600) << " hr " << ((int)time.elapsed() % 3600)/60 << " min " << (int)time.elapsed() % 60 << " s\n";
//			fcout << "Total Evals: " << s.totalevals() << "\n";
//			fcout << "Evals per second: " << (float)s.totalevals()/time.elapsed() << "\n";
//		}
//	}
//	char key;
//	fcout << "Program finished. Press return to exit.\n";
//	key = getchar();
//}
//
//void printstats(tribe& T,int &i,int& totalevals)
//{
//	//boost::progress_timer timer;
//	fcout << "--- Generation " << i << "-------------------------------------------------------------" << "\n";
//fcout << "Number of evals: " << totalevals << "\n";
//fcout << "Best Fitness: " << T.bestFit() <<"\n";
//fcout << "Median Fitness: " << T.medFit()<<"\n";
//fcout << "Mean Size: " << T.meanSize() << "\n";
//fcout << "Median Size: " << T.medSize() << "\n";
//
//fcout << "Equation" << "\t \t \t \t \t" << "Abs Error" << "\n";
//vector <ind> besteqns;
//T.topTen(besteqns);
//for(unsigned int i=0;i<besteqns.size();i++)
//{
//	fcout << besteqns.at(i).eqn <<"\t \t \t \t \t" <<besteqns.at(i).abserror << "\n";
//}
//fcout << "-------------------------------------------------------------------------------" << "\n";
//}
//void load_params(params &p, std::ifstream& fs)
//{
//	if (!fs.good()) 
//	{
//			fcout << "BAD PARAMETER FILE LOCATION" << "\n";
//	}
//
//	string s;
//	string varname;
//	float tmpf;
//	//int tmpi;
//	string tmps;
//	//bool tmpb;
//
//	//string trash;
//	//s.erase();
//    //s.reserve(is.rdbuf()->in_avail());
//
//    while(!fs.eof())
//    {		
//		getline(fs,s,'\n');
//		istringstream ss(s);
//		
//		ss >> varname;
//
//		//getline(is,s,'\t');
//
//		if(varname.compare("g") == 0)
//		{
//			ss >> p.g;
//			//p.g = tmp;
//		}
//		else if(varname.compare("popsize") == 0)
//			ss >> p.popsize;
//		else if(varname.compare("numits") == 0)
//			ss>>p.numits;
//		else if(varname.compare("sel") == 0)
//			ss>>p.sel;
//		else if(varname.compare("tourn_size") == 0)
//			ss>>p.tourn_size;
//		else if(varname.compare("rt_rep") == 0)
//			ss>>p.rt_rep;
//		else if(varname.compare("rt_cross") == 0)
//			ss>>p.rt_cross;
//		else if(varname.compare("rt_mut") == 0)
//			ss>>p.rt_mut;
//		else if(varname.compare("cross") == 0)
//			ss>>p.cross;
//		else if(varname.compare("cross_ar") == 0)
//			ss>>p.cross_ar;
//		else if(varname.compare("mut_ar") == 0)
//			ss>>p.mut_ar;
//		else if(varname.compare("stoperror") == 0)
//			ss>>p.stoperror;
//		else if(varname.compare("init_validate_on") == 0)
//			ss>>p.init_validate_on;
//		else if(varname.compare(0,11,"resultspath") == 0)
//		{
//			int q=0;
//			while (ss>>tmps)
//			{
//				if ( q > 0)
//					p.resultspath = p.resultspath + ' ';
//				p.resultspath.insert(p.resultspath.end(),tmps.begin(),tmps.end());
//				q++;
//			}
//		}
//		else if(varname.compare(0,4,"loud") == 0)
//			ss>>p.loud;
//		else if(varname.compare(0,8,"parallel") == 0)
//			ss>>p.parallel;
//		else if(varname.compare(0,8,"numcores") == 0)
//			ss>>p.numcores;
//		else if(varname.compare(0,11,"sim_nom_mod") == 0)
//			ss>>p.sim_nom_mod;
//		else if(varname.compare(0,7,"nstates") == 0)
//			ss>>p.nstates;
//		else if(varname.compare(0,7,"intvars") == 0)
//		{
//			while (ss>>tmps)
//				p.intvars.push_back(tmps);
//		}
//		else if(varname.compare(0,7,"extvars") == 0)
//		{
//			while (ss>>tmps)
//				p.extvars.push_back(tmps);
//		}
//		/*else if(varname.compare("cons") == 0)
//		{
//			while (ss>>tmps)
//				p.cons.push_back(tmps);
//		}*/
//		else if(varname.compare("cvals") == 0)
//		{
//			while (ss>>tmpf){
//				p.cvals.push_back(tmpf);
//				p.cons.push_back(std::to_string(static_cast<long double>(tmpf)));
//			}
//
//		}
//		else if(varname.compare("seeds") == 0)
//		{
//			while (ss>>tmps)
//				p.seeds.push_back(tmps);
//		}
//		else if(varname.compare("ERC") == 0)
//			ss>>p.ERC;
//		else if(varname.compare("ERCints") == 0)
//			ss>>p.ERCints;
//		else if(varname.compare("maxERC") == 0)
//			ss>>p.maxERC;
//		else if(varname.compare("minERC") == 0)
//			ss>>p.minERC;
//		else if(varname.compare("numERC") == 0)
//			ss>>p.numERC;
//		else if(varname.compare("fit_type") == 0)
//			ss>>p.fit_type;
//		else if(varname.compare("max_fit") == 0)
//			ss>>p.max_fit;
//		else if(varname.compare("min_fit") == 0)
//			ss>>p.min_fit;
//		else if(varname.compare("op_list") == 0)
//		{
//			while (ss>>tmps)
//				p.op_list.push_back(tmps);
//		}
//		else if(varname.compare("op_weight") == 0)
//		{
//			while(ss>>tmpf)
//				p.op_weight.push_back(tmpf);
//		}
//		else if(varname.compare("weight_ops_on") == 0)
//			ss>>p.weight_ops_on;
//		else if(varname.compare("min_len") == 0)
//			ss>>p.min_len;
//		else if(varname.compare("max_len") == 0)
//			ss>>p.max_len;
//		else if(varname.compare("max_dev_len") == 0)
//			ss>>p.max_dev_len;
//		else if(varname.compare("complex_measure") == 0)
//			ss>>p.complex_measure;
//		else if(varname.compare("precision") == 0)
//			ss>>p.precision;
//		else if(varname.compare("lineHC_on") == 0)
//			ss>>p.lineHC_on;
//		else if(varname.compare("lineHC_its") == 0)
//			ss>>p.lineHC_its;
//		else if(varname.compare("pHC_on") == 0)
//			ss>>p.pHC_on;
//		else if(varname.compare("pHC_delay_on") == 0)
//			ss>>p.pHC_delay_on;
//		else if(varname.compare("pHC_size") == 0)
//			ss>>p.pHC_size;
//		else if(varname.compare("pHC_its") == 0)
//			ss>>p.pHC_its;
//		else if(varname.compare("pHC_gauss") == 0)
//			ss>>p.pHC_gauss;
//		else if(varname.compare("eHC_on") == 0)
//			ss>>p.eHC_on;
//		else if(varname.compare("eHC_its") == 0)
//			ss>>p.eHC_its;
//		else if(varname.compare("eHC_prob") == 0)
//			ss>>p.eHC_prob;
//		else if(varname.compare("eHC_size") == 0)
//			ss>>p.eHC_size;
//		else if(varname.compare("eHC_cluster") == 0)
//			ss>>p.eHC_cluster;
//		else if(varname.compare("eHC_dev") == 0)
//			ss>>p.eHC_dev;
//		else if(varname.compare("eHC_best") == 0)
//			ss>>p.eHC_best;
//		else if(varname.compare("eHC_prob_scale") == 0)
//			ss>>p.eHC_prob_scale;
//		else if(varname.compare("eHC_max_prob") == 0)
//			ss>>p.eHC_max_prob;
//		else if(varname.compare("eHC_min_prob") == 0)
//			ss>>p.eHC_min_prob;
//		else if(varname.compare("prto_arch_on") == 0)
//			ss>>p.prto_arch_on;
//		else if(varname.compare("prto_arch_size") == 0)
//			ss>>p.prto_arch_size;
//		else if(varname.compare("prto_sel_on") == 0)
//			ss>>p.prto_sel_on;
//		else if(varname.compare("islands") ==0)
//			ss>>p.islands;
//		else if(varname.compare("island_gens") ==0)
//			ss>>p.island_gens;
//		else{}
//    }
//	p.allvars = p.intvars;
//	p.allvars.insert(p.allvars.end(), p.extvars.begin(), p.extvars.end());
//	p.allblocks = p.allvars;
//	p.allblocks.insert(p.allblocks.end(),p.cons.begin(),p.cons.end());
//	p.allblocks.insert(p.allblocks.end(),p.seeds.begin(),p.seeds.end());
//
//	p.seed = time(0);
//	//normalize fn weights
//	if (p.weight_ops_on) 
//	{
//		float sumweight = accumulate(p.op_weight.begin(),p.op_weight.end(),0);
//		for(unsigned int i=0;i<p.op_weight.size();i++)
//                        p.op_weight.at(i) = p.op_weight.at(i)/sumweight;
//	}
//
//	for (unsigned int i=0; i<p.op_list.size(); i++)
//	{
//		if (p.op_list.at(i).compare("insert")==0)
//			p.op_choice.push_back(0);
//		else if (p.op_list.at(i).compare("absf")==0)
//			p.op_choice.push_back(1);
//		else if (p.op_list.at(i).compare("add")==0)
//			p.op_choice.push_back(2);
//		else if (p.op_list.at(i).compare("cosf")==0)
//			p.op_choice.push_back(3);
//		else if (p.op_list.at(i).compare("DEL0")==0)
//			p.op_choice.push_back(4);
//		else if (p.op_list.at(i).compare("DEL1")==0)
//			p.op_choice.push_back(5);
//		else if (p.op_list.at(i).compare("divL")==0)
//			p.op_choice.push_back(6);
//		else if (p.op_list.at(i).compare("divR")==0)
//			p.op_choice.push_back(7);
//		else if (p.op_list.at(i).compare("DNL")==0)
//			p.op_choice.push_back(8);
//		else if (p.op_list.at(i).compare("DNR")==0)
//			p.op_choice.push_back(9);
//		else if (p.op_list.at(i).compare("FLIP")==0)
//			p.op_choice.push_back(10);
//		else if (p.op_list.at(i).compare("mul")==0)
//			p.op_choice.push_back(11);
//		else if (p.op_list.at(i).compare("NOOP")==0)
//			p.op_choice.push_back(12);
//		else if (p.op_list.at(i).compare("sinf")==0)
//			p.op_choice.push_back(13);
//		else if (p.op_list.at(i).compare("subL")==0)
//			p.op_choice.push_back(14);
//		else if (p.op_list.at(i).compare("subR")==0)
//			p.op_choice.push_back(15);
//		else if (p.op_list.at(i).compare("totheL")==0)
//			p.op_choice.push_back(16);
//		else if (p.op_list.at(i).compare("totheR")==0)
//			p.op_choice.push_back(17);
//		else if (p.op_list.at(i).compare("UP")==0)
//			p.op_choice.push_back(18);
//		else if (p.op_list.at(i).compare("UP2")==0)
//			p.op_choice.push_back(19);
//		else if (p.op_list.at(i).compare("UP3")==0)
//			p.op_choice.push_back(20);
//		else 
//			fcout << "bad command (develep line 306)" << "\n";
//	}
//	
//	p.rep_wheel[0]=(p.rt_rep);
//	p.rep_wheel[1]=(p.rt_cross);
//	p.rep_wheel[2]=(p.rt_mut);
//
//	partial_sum(p.rep_wheel.begin(), p.rep_wheel.end(), p.rep_wheel.begin());
//	
//}
//
//void load_data(data &d, std::ifstream& fs)
//{
//	if (!fs.good()) 
//	{
//			fcout << "BAD DATA FILE LOCATION" << "\n";
//	}
//
//	string s;
//	string varname;
//	float tmpf;
//	float tarf;
//	//int tmpi;
//	string tmps;
//	//bool tmpb;
//	//int varnum =0;
//	int rownum=0;
//	unsigned int index=0;
//	bool pass = 0;
//	// get variable names from first line / number of variables
//	getline(fs,s,'\n');
//	istringstream ss(s);
//	ss>>tmps;
//
//	while (ss>>tmps)
//	{
//		while (!pass && index<p.allvars.size())
//		{
//			if(tmps.compare(p.allvars.at(index))==0)
//			{
//				d.label.push_back(p.allvars.at(index)); // populate data labels from all vars based on column location
//				pass=1;
//				
//			}
//			else
//			{
//				index++;
//			}
//		}
//
//		pass=0;
//		index=0;
//	}
//    while(!fs.eof())
//    {		
//		getline(fs,s,'\n');
//		istringstream ss(s);
//		// get target data
//		ss >> tarf;
//		d.target.push_back(tarf);
//		d.vals.push_back(vector<float>());
//		// get variable data
//		while(ss >> tmpf)
//		{
//			d.vals[rownum].push_back(tmpf);
//			//varnum++;
//		}
//		rownum++;
//    }
//	d.dattovar.resize(p.allvars.size());
//	//symbol_table.add_constants();
//	//for (unsigned int i=0; i<p.cons.size();i++)
//	//	d.symbol_table.add_constant(p.cons.at(i),p.cvals.at(i));
//
//	
//}
//bool stopcondition(float &bestfit)
//{
//	if (bestfit <= 0.0001)
//		return true;
//	else
//		return false;
//}
