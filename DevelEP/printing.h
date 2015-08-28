#pragma once
#ifndef LOAD_PRINTING_H
#define LOAD_PRINTING_H

using namespace std;


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
	if (p.classification && p.class_m3gp){
		fout << "M:\n";
		fout << best.M << "\n";
		for (unsigned i = 0; i<p.number_of_classes; ++i){
			fout << "C[" << i << "]:\n";
			fout << best.C[i] << "\n";
		}	
	}
	fout << "output:\n";
	for(unsigned int i =0;i<best.output.size();++i)
	{
		fout << best.output.at(i) << "\n";
	}
	//fout<<"\n";
}
void initdatafile(std::ofstream& dfout,string & logname,params& p)
{
	string dataname = logname.substr(0,logname.size()-4)+".data";
	dfout.open(dataname,std::ofstream::out | std::ofstream::app);
	//dfout.open(dataname,std::ofstream::app);
	//dfout << "pt_evals \t best_eqn \t best_fit \t best_fit_v \t med_fit \t med_fit_v \t best_MAE \t best_MAE_v \t best_R2 \t best_R2_v \t best_VAF \t best_VAF_v \t size \t eff_size \t pHC_pct \t eHC_pct \t good_g_pct \t neut_g_pct \t bad_g_pct \t tot_hom \t on_hom \t off_hom\n";
	dfout << "gen \t pt_evals \t best_eqn \t best_fit \t best_fit_v \t med_fit \t med_fit_v \t best_MAE \t best_MAE_v \t best_R2 \t best_R2_v \t best_VAF \t best_VAF_v \t size \t eff_size \t pHC_pct \t eHC_pct \t good_g_pct \t neut_g_pct \t bad_g_pct ";
	if (p.print_homology)
		dfout << "\t tot_hom \t on_hom \t off_hom";
	if (p.classification && p.class_m3gp) 
		dfout << "\t dimension" ;
	if (p.print_novelty)
		dfout << "\t novelty" ;
	if (p.print_protected_operators)
		dfout << "\t best_eqn_matlab";
	dfout << "\n";
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
	dfout << gen << "\t" << s.totalptevals() << "\t" << best_ind.eqn << "\t" << T.bestFit() << "\t" << T.bestFit_v() << "\t" << T.medFit() << "\t" << T.medFit_v() << "\t" << best_ind.abserror << "\t" << best_ind.abserror_v << "\t" << best_ind.corr << "\t" << best_ind.corr_v << "\t" << best_ind.VAF << "\t" << best_ind.VAF_v << "\t" << T.meanSize() << "\t" << T.meanEffSize() << "\t" << s.current_pHC_updates/float(p.popsize)*100.0 << "\t" << s.current_eHC_updates/float(p.popsize)*100.0 << "\t" <<  s.good_cross_pct << "\t" << s.neut_cross_pct << "\t" << s.bad_cross_pct;
	if (p.print_homology){
		float tot_hom, on_hom, off_hom;
		T.hom(r,tot_hom,on_hom,off_hom);
		dfout << "\t" << tot_hom << "\t" << on_hom << "\t" << off_hom;
	}
	if (p.classification && p.class_m3gp) 
		dfout << "\t" << best_ind.dim ;
	if (p.print_novelty){
		float novelty;
		T.novelty(novelty);
		dfout << "\t" << novelty;
	}
	if (p.print_protected_operators)
		dfout << "\t " + best_ind.eqn_matlab;

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
		//bestname = "lastpop.last_pop";
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
		if (p.classification && p.class_m3gp) fout << "dimensions: " << pop[h].dim << "\n";
		/*fout << "output: ";
		for(unsigned int i =0;i<pop.at(h).output.size();++i)
		{
			fout << T.pop.at(h).output.at(i);
		}
		fout<<"\n";*/
		fout << "------------------------------------------------------------------" << "\n";
		}
}

#endif