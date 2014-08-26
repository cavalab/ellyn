// RunTrials.cpp : Run trials of Develep.
// input is a file containing the number of trials, the parameter file for those trials, and the data for those trials. 
// the settings are passed on to run Develep.  parameter file names and data are then passed on to Develep to run. 

#include "stdafxRT.h"
#include "runDevelep.h"
#include "pop.h"
#include "mpi.h"
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
		//MPI stuff
		MPI::Init();
		int size = MPI::COMM_WORLD.Get_size();
		int rank = MPI::COMM_WORLD.Get_rank();

#if defined(_WIN32)
			cout << "running process " + to_string(static_cast<long long>(rank)) + "of" + to_string(static_cast<long long>(size)) + ": " + paramfile.at(i).substr(paramfile.at(i).rfind('\\')+1,paramfile.at(i).size()) + ", " + datafile.at(i).substr(datafile.at(i).rfind('\\')+1,datafile.at(i).size())  + "\n";
#else
			cout << "running process " + to_string(static_cast<long long>(rank)) + "of" + to_string(static_cast<long long>(size)) + ": " + paramfile.at(i).substr(paramfile.at(i).rfind('/')+1,paramfile.at(i).size()) + ", " + datafile.at(i).substr(datafile.at(i).rfind('/')+1,datafile.at(i).size())  + "\n";
#endif

		runDevelep(paramfile.at(i).c_str(),datafile.at(i).c_str(),1); //run develep

		MPI::Finalize();
		char key;
		cout << "All trials completed. Press return to exit." << endl;
		key = getchar();
	}
	catch(const std::bad_alloc&)
	{
		cout << "bad allocation error. \n";
		abort();
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
	boost::regex re("~");
	int n=0;
    while(!fs.eof())
    {		
		getline(fs,s,'\n');
		istringstream ss(s);
		
		ss >> tmpi;
		trialset.push_back(tmpi);
		ss >> tmps;
		paramset.push_back(boost::regex_replace(tmps,re," "));
		ss >> tmps;
		dataset.push_back(boost::regex_replace(tmps,re," "));
		n++;
	}
	for(unsigned int i=0;i<trialset.size();i++)
	{
		for (unsigned int j = 0; j<trialset[i];j++)
		{
			paramfile.push_back(paramset[i]);
			datafile.push_back(dataset[i]);
		}
		totaltrials+=trialset[i];
	}
}
