#pragma once
#ifndef LOAD_DATA_H
#define LOAD_DATA_H

using namespace std;

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
	ss>>d.target_var;
	vector<int> d_col; // keep indices; columns of data that are to be used (they are variables in p.allvars)
	int d_index=0;
	bool useall=false;
	// if intvars is not specified, set it to all variables in data file
	if (p.allvars.empty())
		useall = true;
	
	while (ss>>tmps)
	{
		while (!pass && (useall || index<p.allvars.size()))
		{
			if(useall || tmps.compare(p.allvars.at(index))==0)
			{
				d.label.push_back(tmps); // populate data labels from all vars based on column location
				pass=1;
				d_col.push_back(d_index);
				++d_index;
			}
			else
			{
				++index;
			}
			
		}
		if (!pass) ++d_index;
		pass=0;
		index=0;
	}
	if (useall) // assign data labels to p.allvars
		p.allvars = d.label;
	int varcount=0;
	int cur_col=0;
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
		
		while(ss2 >> tmpf && cur_col<d_col.size())
		{
			if (varcount == d_col[cur_col]){ 
				d.vals[rownum].push_back(tmpf);
				++cur_col;
			}
			++varcount;
		}
		++rownum;
		cur_col=0;
		varcount=0;
    }
	// pop end in case of extra blank lines in data file
	while(d.vals.back().empty())
	{
		d.vals.pop_back();
	}

	if (p.classification){
		// set p.num_cases based on unique elements in target
		vector<float> tmp = d.target;
		sort(tmp.begin(),tmp.end());
		tmp.erase(unique(tmp.begin(),tmp.end()),tmp.end());
		p.number_of_classes = tmp.size();
		// set lowest class equal to zero
		int offset  = *std::min_element(d.target.begin(),d.target.end());
		if (offset>0){
			for (unsigned i=0;i<d.target.size();++i) 
				d.target[i] -= offset;
		}
	}
	
	if (p.AR){ // make auto-regressive variables
		// add data columns to d.vals 
		int numvars = d.vals[0].size();
		for (unsigned i = 0; i<p.AR_n; ++i){
			for (unsigned j = 0; j<d.vals.size(); ++j){
				for (unsigned k = 0; k<numvars;++k){
					if(j>i) 
						d.vals[j].push_back(d.vals[j-i-1][k]);
					else 
						d.vals[j].push_back(0.0);
				}
			}
		}
		// add data labels to d.label and p.allvars
			int tmp = d.label.size();
			for (unsigned i = 0; i<p.AR_n; ++i){ 
				d.vals.erase(d.vals.begin()); // erase zero padding from data
				d.target.erase(d.target.begin()); //erase zero padding from target
				for (unsigned k=0; k<tmp;++k){ // add delay variable names
					d.label.push_back(d.label[k] + "_" + to_string(static_cast<long long>(i+1)));
					p.allvars.push_back(d.label[k] + "_" + to_string(static_cast<long long>(i+1)));
				}
			}
			// load AR variables if necessary
		
			for (int i=0;i<p.AR_n;++i){
				p.allvars.push_back(d.target_var + "_" + to_string(static_cast<long long>(i+1)));
				d.label.push_back(d.target_var + "_" + to_string(static_cast<long long>(i+1)));
			}
		
	}
	//d.dattovar.resize(p.allvars.size());
	//d.mapdata();
}

#endif