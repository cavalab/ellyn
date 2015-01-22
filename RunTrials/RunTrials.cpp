// RunTrials.cpp : Run trials of EllenGP on openMP threads.
// input is a file containing the number of trials, the parameter file for those trials, and the data for those trials. 
// the settings are passed on to run Develep.  parameter file names and data are then passed on to Develep to run. 

#include "stdafxRT.h"
#include "runEllenGP.h"
#include "pop.h"
using namespace std;

void getTrialSetup(ifstream& fs,int&,vector<int>&,vector<string>&,vector<string>&);

int main(int argc, char** argv)
{
	//_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	//_CrtMemState s1;
	//_CrtMemCheckpoint( &s1 );

	/* run multiple trials of Develep 
	   input: trial text file
	   first column: number of trials to run
	   second column: parameter file for trials
	   third column: data for trials
	*/
	try
	{
		string trialsetup(argv[1]);

		int totaltrials = 0;
		vector<int> trialset; 
		vector<string> paramfile;
		vector<string> datafile;

		ifstream fs(trialsetup);
		getTrialSetup(fs,totaltrials,trialset,paramfile,datafile);
		

		cout << "Running Trials: \n";

		#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<totaltrials;++i)
		{
#if defined(_WIN32)
			cout << to_string(static_cast<long long>(i)) + ": " + paramfile.at(i).substr(paramfile.at(i).rfind('\\')+1,paramfile.at(i).size()) + ", " + datafile.at(i).substr(datafile.at(i).rfind('\\')+1,datafile.at(i).size())  + "\n";
#else
			cout << to_string(static_cast<long long>(i)) + ": " + paramfile.at(i).substr(paramfile.at(i).rfind('/')+1,paramfile.at(i).size()) + ", " + datafile.at(i).substr(datafile.at(i).rfind('/')+1,datafile.at(i).size())  + "\n";
#endif
			runEllenGP(paramfile.at(i),datafile.at(i),1,i); //run develep
			//_CrtMemDumpStatistics( &s1 );
		}
		cout << "All trials completed. Exiting." << endl;
		//key = getchar();
	}
	catch(const std::bad_alloc&)
	{
		cout << "bad allocation error. \n";
		exit(1);
	}
	catch(exception& er) 
	{
		cout << "Error: " << er.what() << endl;
		exit(1);

	}
	catch(...)
	{
		cout << "Exception Occurred. Ending..."<<endl;
		exit(1);
	}
	return 0;
}

void getTrialSetup(ifstream& fs,int& totaltrials,vector<int>& trialset,vector<string>& paramfile,vector<string>& datafile)
{
	if (!fs.good()) 
	{
			cerr << "error: can't find input file." << "\n";
			exit(1);
	}
	
	vector<string> paramset;
	vector<string> dataset;
	string s;
	string tmps;
	int tmpi;
	//string trash;
	//s.erase();
    //s.reserve(is.rdbuf()->in_avail());
	//#if defined(_WIN32)
	//	using namespace std;
	//#else
	//	using namespace boost;
	//#endif
	//regex re("~");
	int n=0;
    while(!fs.eof())
    {		
		getline(fs,s,'\n');
		// in case of extra blank lines in data file
		if (!s.empty()){
		
			istringstream ss(s);
		
			ss >> tmpi;
			trialset.push_back(tmpi);
			ss >> tmps;
			//paramset.push_back(regex_replace(tmps,re," "));
			paramset.push_back(tmps);
			ss >> tmps;
			//dataset.push_back(regex_replace(tmps,re," "));
			dataset.push_back(tmps);
			++n;
		}
	}
	for(unsigned int i=0;i<trialset.size();++i)
	{
		for (unsigned int j = 0; j<trialset[i];++j)
		{
			paramfile.push_back(paramset[i]);
			datafile.push_back(dataset[i]);
		}
		totaltrials+=trialset[i];
	}
}
