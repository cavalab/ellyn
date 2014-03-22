#include "stdafx.h"
// mine

#include "pop.h"
#include "data.h"
#include "logger.h"
#include "state.h"
#include "InitPop.h"
#include "Fitness.h"
#include "Generation.h"
#include "instructionset.h"
#include "Gen2Phen.h"
#include "Generationfns.h"
#include "strdist.h"
#include <time.h>
#include <cstring>

using namespace std;

// global parameters structure


void load_params(params &p, std::ifstream& is);
void load_data(data &d, std::ifstream& is,params&);
bool stopcondition(float&);
void printstats(tribe& T,int &i,state& s,params& p);
void printbestind(tribe& T,params& p,string& logname);

void runDevelep(string& paramfile, string& datafile,bool trials)
{
	
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
	//class logger fcout;
	vector <Randclass> r;
	
	// load parameter file
	ifstream fs(paramfile);
	load_params(p, fs);
	// load data file
	ifstream ds(datafile);
	load_data(d,ds,p);
	s.out.set(trials);

	std::time_t t =  std::time(NULL);
    std::tm tm    = *std::localtime(&t);
	char tmplog[100];
	strftime(tmplog,100,"%Y-%m-%d_%H-%M-%S",&tm);
	const char * c = p.resultspath.c_str();
	_mkdir(c);
	 int thrd = omp_get_thread_num();
	 string thread = std::to_string(static_cast<long long>(thrd));
	 string pname = paramfile.substr(paramfile.rfind('\\')+1,paramfile.size());
	 pname = pname.substr(0,pname.size()-4);
	 string dname = datafile.substr(datafile.rfind('\\')+1,datafile.size());
	 dname = dname.substr(0,dname.size()-4);
     string logname = p.resultspath + '\\' + "Develep_" + tmplog + "_" + pname + "_" + dname + "_" + thread + ".log";
	 s.out.open(logname);
	 s.out << "_______________________________________________________________________________ \n";
	 s.out << "                                    Develep                                     \n";
	 s.out << "_______________________________________________________________________________ \n";
	 s.out << "Time right now is " << std::put_time(&tm, "%c %Z") << '\n';
	// s.out<< "Results Path: " << logname  << "\n";
	 s.out << "parameter file: " << paramfile << "\n";
	 s.out << "data file: " << datafile << "\n";
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
	 s.out << "Total Population Size: " << p.popsize << "\n";
	 s.out << "Max Generations: " << p.g << "\n";

	//initialize random number generator
	unsigned int seed1 = unsigned int(time(NULL));
	s.out << "seeds: \n";
	r.resize(omp_get_max_threads());
	#pragma omp parallel for 
	for(int seeder=0;seeder<omp_get_max_threads();seeder++)
	{
			//cout << "seeder: " << seeder <<endl;
			//cout << "seed1: " << seed1*seeder <<endl;
			if(!trials) s.out << to_string(static_cast<long long>(seed1*(seeder+1))) + "\n";
			r.at(seeder).SetSeed(seed1*(seeder+1));
	}
	if(trials)
		s.out << (omp_get_thread_num()+1)*seed1 << "\n";


	
	boost::timer time;
	
	if (p.islands)
	{
		//p.parallel=false;
		
		int num_islands=omp_get_max_threads();
		int subpops = p.popsize/num_islands;
		vector<tribe> T;
		tribe World(subpops*num_islands,p.max_fit,p.min_fit); //total population of tribes
		 s.out << num_islands << " islands of " << subpops << " individuals, total pop " << subpops*num_islands <<"\n";

		for(int i=0;i<num_islands;i++)
			T.push_back(tribe(subpops,p.max_fit,p.min_fit));
		// run separate islands 
		if (p.init_validate_on)
		{
			 s.out << "Initial validation..."; 
			#pragma omp parallel for 
			for(int i=0;i<num_islands;i++)
			{
				float worstfit;
				float bestfit;
				vector<ind> tmppop;
				// s.out << "Initialize Population..." << "\n";
				InitPop(T.at(i).pop,p,r);
				// s.out << "Gen 2 Phen..." << "\n";
				Gen2Phen(T.at(i).pop,p);
				// s.out << "Fitness..." << "\n";
				Fitness(T.at(i).pop,p,d,s);
				worstfit = T.at(i).worstFit();
				bestfit = T.at(i).bestFit();
				int counter=0;
				while(worstfit == p.max_fit && counter<100)
				{
					for (int j=0;j<T.at(i).pop.size(); j++)
					{
						if ( T.at(i).pop.at(j).fitness == p.max_fit)
						{
							T.at(i).pop.erase(T.at(i).pop.begin()+j);
							tmppop.push_back(ind());
						}
					}

					InitPop(tmppop,p,r);
					Gen2Phen(tmppop,p);
					Fitness(tmppop,p,d,s);
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
			//totalevals=0;
			
		}
		else
		{
			#pragma omp parallel for
			for(int i=0;i<num_islands;i++)
			{
				InitPop(T.at(i).pop,p,r);
				// s.out << "Gen 2 Phen..." << "\n";
				Gen2Phen(T.at(i).pop,p);
				// s.out << "Fitness..." << "\n";
				Fitness(T.at(i).pop,p,d,s);
			}
		}
		int gen=0;
		bool pass=1;
		int mixtrigger=0;
		while(gen<=p.g && !stopcondition(World.best))
		{
			
			#pragma omp parallel for
			for(int i=0;i<num_islands;i++)
			{
				if(pass)
				{
					if (p.sel==2)
					{
						for(int j=0;j<p.island_gens;j++)
						{
							if(pass)
							{
								Generation(T.at(i).pop,p,r,d,s);
								if (stopcondition(T.at(i).best))
								{
									pass=0;
								}
							}

							//cout << "thread " << omp_get_thread_num() << " numevals: " << s[i].numevals << "\n";
						}
					}
					else
					{
						if(pass)
							{
								Generation(T.at(i).pop,p,r,d,s);								
								if (stopcondition(T.at(i).best))
								{
									pass=0;
								}
							}
					}
					if (p.pHC_on && p.ERC)
					{
							for(int k=0; k<T.at(i).pop.size(); k++)
								HillClimb(T.at(i).pop.at(k),p,r,d,s);
					}
					if (p.eHC_on) 
					{
							for(int m=0; m<T.at(i).pop.size(); m++)
								EpiHC(T.at(i).pop.at(m),p,r,d,s);
					}
				}

			// construct world population
				swap_ranges(World.pop.begin()+i*subpops,World.pop.begin()+(i+1)*subpops,T.at(i).pop.begin());
			}

				
		s.setgenevals();
		if(s.genevals.back()>mixtrigger) 
		{
			//shuffle population	
			std::random_shuffle(World.pop.begin(),World.pop.end(),r[omp_get_thread_num()]);
			//redistribute populations to islands
			#pragma omp parallel for
			for(int i=0;i<num_islands;i++)
			{
				T.at(i).pop.assign(World.pop.begin()+i*subpops,World.pop.begin()+(i+1)*subpops);
			}
			mixtrigger+=p.island_gens*subpops;
		}
		 
		 printstats(World,gen,s,p);
		 s.out << "Total Time: " << (int)floor(time.elapsed()/3600) << " hr " << ((int)time.elapsed() % 3600)/60 << " min " << (int)time.elapsed() % 60 << " s\n";
		 s.out << "Total Evals: " << s.totalevals() << "\n";
		 s.out << "Average evals per second: " << (float)s.totalevals()/time.elapsed() << "\n";
		//World.pop.clear();
		gen++;
		}
		printbestind(World,p,logname);
	}
	else //no islands
	{
		tribe T(p.popsize,p.max_fit,p.min_fit);
		if (p.init_validate_on)
		{		
			float worstfit;
			int counter=0;
			//float bestfit;
			vector<ind> tmppop;
			// s.out << "Initialize Population..." << "\n";
			InitPop(T.pop,p,r);
			// s.out << "Gen 2 Phen..." << "\n";
			Gen2Phen(T.pop,p);
			// s.out << "Fitness..." << "\n";
			Fitness(T.pop,p,d,s);
			worstfit = T.worstFit();
			while(worstfit == p.max_fit && counter<100)
			{
				for (int j=0;j<T.pop.size(); j++)
				{
					if ( T.pop.at(j).fitness == p.max_fit)
					{
						T.pop.erase(T.pop.begin()+j);
						tmppop.push_back(ind());
					}
				}
				InitPop(tmppop,p,r);
				Gen2Phen(tmppop,p);
				Fitness(tmppop,p,d,s);
				T.pop.insert(T.pop.end(),tmppop.begin(),tmppop.end());
				tmppop.clear();
				worstfit = T.worstFit();
				counter++;
				if(counter==100)
					s.out << "initial population count exceeded. Starting evolution...\n";
			}
		}
		s.setgenevals();
		s.out << " number of evals: " << s.getgenevals() << "\n";
		int i=1;
		int trigger=0;
		int gen=0;
		int counter=0;
		if(p.sel==2) // if using deterministic crowding, increase gen size
		{
			gen = p.g*(p.popsize*(p.rt_mut+p.rt_rep) + p.popsize*p.rt_cross/2);
			trigger = p.popsize*(p.rt_mut+p.rt_rep)+p.popsize*p.rt_cross/2;
		}
		else
			gen=p.g;
		while (i<=gen && !stopcondition(T.best))
		{
			
			 
			 Generation(T.pop,p,r,d,s);


			 if (i>trigger)
			 {
				 if (p.pHC_on && p.ERC)
				 {
		 			for(int k=0; k<T.pop.size(); k++)
		 				HillClimb(T.pop.at(k),p,r,d,s);
		 		 }
				 if (p.eHC_on) 
				 {
					for(int m=0; m<T.pop.size(); m++)
						EpiHC(T.pop.at(m),p,r,d,s);
				 } 

				 s.setgenevals();
				 s.out << "Elapsed time: \n";
				 printstats(T,counter,s,p);
				
				s.out << "Total Time: " << (int)floor(time.elapsed()/3600) << " hr " << ((int)time.elapsed() % 3600)/60 << " min " << (int)time.elapsed() % 60 << " s\n";
				s.out << "Total Evals: " << s.totalevals() << "\n";
				s.out << "Average evals per second: " << (float)s.totalevals()/time.elapsed() << "\n";

				if (p.sel==2)
					trigger+=p.popsize*(p.rt_mut+p.rt_rep)+p.popsize*p.rt_cross/2;
				counter++;
			 }
			i++;
		}
		printbestind(T,p,logname);
	}

	if(!trials)
	{
		char key;
		s.out << "Program finished. Press return to exit.\n";
		key = getchar();
	}
	else
	{
		s.out << "Program finished. \n";
	}

}

void printstats(tribe& T,int &i,state& s,params& p)
{
	//boost::progress_timer timer;
	s.out << "--- Generation " << i << "---------------------------------------------------------------" << "\n";
s.out << "Number of evals: " << s.genevals.back() << "\n";
s.out << "Best Fitness: " << T.bestFit() <<"\n";
s.out << "Median Fitness: " << T.medFit()<<"\n";
s.out << "Mean Size: " << T.meanSize() << "\n";
s.out << "Mean Eff Size: " << T.meanEffSize() << "\n";
if(p.pHC_on)
	s.out << "Parameter updates: " << s.getpHCupdates() << "\n";
if(p.eHC_on)
	s.out << "Epigenetic updates: " << s.geteHCupdates() << "\n";
s.out << "Beneficial Genetics: " << s.getGoodCrossPct() << "\%\n";
s.out << "Neutral Genetics: " << s.getNeutCrossPct() << "\%\n";
s.out << "Bad Genetics: " << s.getBadCrossPct() << "%\n";
s.clearCross();
s.out << "Equation" << "\t \t \t \t \t" << "Abs Error" << "\n";
vector <ind> besteqns;
T.topTen(besteqns);
for(unsigned int i=0;i<besteqns.size();i++)
{
	s.out << besteqns.at(i).eqn <<"\t \t \t \t \t" <<besteqns.at(i).abserror << "\n";
}
s.out << "-------------------------------------------------------------------------------" << "\n";
}
void printbestind(tribe& T,params& p,string& logname)
{
	ind best;
	T.getbestind(best);
	string bestname = logname.substr(0,logname.size()-4)+"_best_ind.log";
	std::ofstream fout;
	fout.open(bestname);
	//boost::progress_timer timer;
	fout << "--- Best Individual ---------------------------------------------------------------" << "\n";
	fout << "Corresponding Logfile: " + logname + "\n";
	fout << "f = " + best.eqn + "\n";
	fout << " line: ";
	for(unsigned int i =0;i<best.line.size();i++)
		fout << best.line[i] << "\t";
	fout << endl;
	fout << "eline: ";
	for(unsigned int i =0;i<best.epiline.size();i++)
		fout << best.epiline[i] << "\t";
	fout << endl;
	fout << "args: ";
	for(unsigned int i =0;i<best.args.size();i++)
		fout << best.args.at(i) << "\t";
	fout << endl;
	fout << "size: " << best.line.size() << "\n";
	fout << "eff size: " << best.eff_size << "\n";
	fout << "abs error: " << best.abserror<< "\n";;
	fout << "correlation: " << best.corr<< "\n";;
	fout << "fitness: " << best.fitness<< "\n";;

	fout << "parent fitness: " << best.parentfitness << "\n";;
	fout << "origin: " << best.origin << "\n";
	fout << "age: " << best.age << "\n";
	fout << "eqn form: " << best.eqn_form << "\n";
}
void load_params(params &p, std::ifstream& fs)
{
	if (!fs.good()) 
	{
			cout << "BAD PARAMETER FILE LOCATION" << "\n";
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
		else if(varname.compare("numits") == 0)
			ss>>p.numits;
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
		else if(varname.compare("stoperror") == 0)
			ss>>p.stoperror;
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
				q++;
			}
		}
		else if(varname.compare(0,4,"loud") == 0)
			ss>>p.loud;
		else if(varname.compare(0,8,"parallel") == 0)
			ss>>p.parallel;
		else if(varname.compare(0,8,"numcores") == 0)
			ss>>p.numcores;
		else if(varname.compare(0,11,"sim_nom_mod") == 0)
			ss>>p.sim_nom_mod;
		else if(varname.compare(0,7,"nstates") == 0)
			ss>>p.nstates;
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
		else if(varname.compare("max_dev_len") == 0)
			ss>>p.max_dev_len;
		else if(varname.compare("complex_measure") == 0)
			ss>>p.complex_measure;
		else if(varname.compare("precision") == 0)
			ss>>p.precision;
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
		else if(varname.compare("eHC_size") == 0)
			ss>>p.eHC_size;
		else if(varname.compare("eHC_cluster") == 0)
			ss>>p.eHC_cluster;
		else if(varname.compare("eHC_dev") == 0)
			ss>>p.eHC_dev;
		else if(varname.compare("eHC_best") == 0)
			ss>>p.eHC_best;
		else if(varname.compare("eHC_init")==0)
			ss>>p.eHC_init;
		else if(varname.compare("eHC_prob_scale") == 0)
			ss>>p.eHC_prob_scale;
		else if(varname.compare("eHC_max_prob") == 0)
			ss>>p.eHC_max_prob;
		else if(varname.compare("eHC_min_prob") == 0)
			ss>>p.eHC_min_prob;
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
		else{}
    }
	p.allvars = p.intvars;
	p.allvars.insert(p.allvars.end(), p.extvars.begin(), p.extvars.end());
	p.allblocks = p.allvars;
	p.allblocks.insert(p.allblocks.end(),p.cons.begin(),p.cons.end());
	p.allblocks.insert(p.allblocks.end(),p.seeds.begin(),p.seeds.end());

	p.seed = time(0);
	//normalize fn weights
	if (p.weight_ops_on) 
	{
		float sumweight = accumulate(p.op_weight.begin(),p.op_weight.end(),0);
		for(unsigned int i=0;i<p.op_weight.size();i++)
                        p.op_weight.at(i) = p.op_weight.at(i)/sumweight;
	}

	for (unsigned int i=0; i<p.op_list.size(); i++)
	{
		if (p.op_list.at(i).compare("insert")==0)
			p.op_choice.push_back(0);
		else if (p.op_list.at(i).compare("absf")==0)
			p.op_choice.push_back(1);
		else if (p.op_list.at(i).compare("add")==0)
			p.op_choice.push_back(2);
		else if (p.op_list.at(i).compare("cosf")==0)
			p.op_choice.push_back(3);
		else if (p.op_list.at(i).compare("DEL0")==0)
			p.op_choice.push_back(4);
		else if (p.op_list.at(i).compare("DEL1")==0)
			p.op_choice.push_back(5);
		else if (p.op_list.at(i).compare("divL")==0)
			p.op_choice.push_back(6);
		else if (p.op_list.at(i).compare("divR")==0)
			p.op_choice.push_back(7);
		else if (p.op_list.at(i).compare("DNL")==0)
			p.op_choice.push_back(8);
		else if (p.op_list.at(i).compare("DNR")==0)
			p.op_choice.push_back(9);
		else if (p.op_list.at(i).compare("FLIP")==0)
			p.op_choice.push_back(10);
		else if (p.op_list.at(i).compare("mul")==0)
			p.op_choice.push_back(11);
		else if (p.op_list.at(i).compare("NOOP")==0)
			p.op_choice.push_back(12);
		else if (p.op_list.at(i).compare("sinf")==0)
			p.op_choice.push_back(13);
		else if (p.op_list.at(i).compare("subL")==0)
			p.op_choice.push_back(14);
		else if (p.op_list.at(i).compare("subR")==0)
			p.op_choice.push_back(15);
		else if (p.op_list.at(i).compare("totheL")==0)
			p.op_choice.push_back(16);
		else if (p.op_list.at(i).compare("totheR")==0)
			p.op_choice.push_back(17);
		else if (p.op_list.at(i).compare("UP")==0)
			p.op_choice.push_back(18);
		else if (p.op_list.at(i).compare("UP2")==0)
			p.op_choice.push_back(19);
		else if (p.op_list.at(i).compare("UP3")==0)
			p.op_choice.push_back(20);
		else if (p.op_list.at(i).compare("expf")==0)
			p.op_choice.push_back(21);
		else if (p.op_list.at(i).compare("logf")==0)
			p.op_choice.push_back(22);
		else 
			cout << "bad command (develep line 306)" << "\n";
	}
	
	p.rep_wheel.push_back(p.rt_rep);
	p.rep_wheel.push_back(p.rt_cross);
	p.rep_wheel.push_back(p.rt_mut);

	partial_sum(p.rep_wheel.begin(), p.rep_wheel.end(), p.rep_wheel.begin());
	
}

void load_data(data &d, std::ifstream& fs,params& p)
{
	if (!fs.good()) 
	{
			cout << "BAD DATA FILE LOCATION" << "\n";
	}

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
				index++;
			}
		}

		pass=0;
		index=0;
	}
    while(!fs.eof())
    {		
		getline(fs,s,'\n');
		istringstream ss(s);
		// get target data
		ss >> tarf;
		d.target.push_back(tarf);
		d.vals.push_back(vector<float>());
		// get variable data
		while(ss >> tmpf)
		{
			d.vals[rownum].push_back(tmpf);
			//varnum++;
		}
		rownum++;
    }
	d.dattovar.resize(p.allvars.size());
	//symbol_table.add_constants();
	//for (unsigned int i=0; i<p.cons.size();i++)
	//	d.symbol_table.add_constant(p.cons.at(i),p.cvals.at(i));

	
}
bool stopcondition(float &bestfit)
{
	if (bestfit <= 0.0001)
		return true;
	else
		return false;
}
