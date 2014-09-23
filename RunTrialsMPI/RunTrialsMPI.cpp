// RunTrials.cpp : Run trials of Develep.
// input is a file containing the number of trials, the parameter file for those trials, and the data for those trials. 
// the settings are passed on to run Develep.  parameter file names and data are then passed on to Develep to run. 

#include "stdafxRTMPI.h"
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
       // cout << "trialsetup: " + trialsetup + "\n";
		int totaltrials = 0;
		vector<int> trialset; 
		vector<string> paramfile;
		vector<string> datafile;

		ifstream fs(trialsetup);
		getTrialSetup(fs,totaltrials,trialset,paramfile,datafile);
		
		int numsent=0;

		//MPI stuff
		int master=0;
		int ierr;
		MPI::Init();
		int numprocs = MPI::COMM_WORLD.Get_size();
		int myid = MPI::COMM_WORLD.Get_rank();
		//cout << "I am process " + to_string(static_cast<long long>(myid)) + " of " + to_string(static_cast<long long>(numprocs)) + "\n";
		MPI::Status status;
		//const char * pbuff,dbuff;
        if (myid==master){
        	//cout << "total trials: " + to_string(static_cast<long long>(totaltrials)) + "\n";
        	//cout << "In master loop\n";
        	cout << "Running Trials: \n";
        	// schedule tasks from master node
        	for (int i=0;i<min(numprocs-1,totaltrials);i++){
        		//cout << "sending " + paramfile.at(i) + " to process " + to_string(static_cast<long long>(i)) + "\n";
        		MPI::COMM_WORLD.Send(paramfile.at(i).c_str(),paramfile.at(i).length(),MPI::CHAR,i+1,i+1);
        		//cout << "sending " + datafile.at(i) + " to process " + to_string(static_cast<long long>(i)) + "\n";
				MPI::COMM_WORLD.Send(datafile.at(i).c_str(),datafile.at(i).length(),MPI::CHAR,i+1,i+1);
				numsent++;
				//cout << "numsent: " + to_string(static_cast<long long>(numsent)) + "\n";
        	}
        	//int curnumsent=numsent;
        	int stops =0;
        	while(numsent<=totaltrials && stops<numprocs-1){
        		int ans;
        		MPI::COMM_WORLD.Recv(&ans,1,MPI::INT,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
        		const int sender = status.Get_source();
        		//int anstype = status.Get_tag();
        		if (numsent < totaltrials){
        			MPI::COMM_WORLD.Send(paramfile.at(numsent).c_str(),paramfile.at(numsent).length(),MPI::CHAR,sender,numsent+1);
					MPI::COMM_WORLD.Send(datafile.at(numsent).c_str(),datafile.at(numsent).length(),MPI::CHAR,sender,numsent+1);
					++numsent;
        		}
        		else{
        			cout << "sending stop command to process " + to_string(static_cast<long long>(sender)) + "\n";
        			MPI::COMM_WORLD.Send(MPI::BOTTOM,0,MPI::CHAR,sender,0);
        			++stops;
        		}

        	} cout << "out of master while loop\n";
        }
        else{
        	//cout << "in slave task \n";
        	// receive tasks and send completion messages to master
        	//cout << "in slave task. myid is " + to_string(static_cast<long long>(myid)) + " and totaltrials is " + to_string(static_cast<long long>(totaltrials)) + "\n";
        	bool cont = true;
        	while (cont){
				if (myid < totaltrials){
					//char * pbuff,dbuff;
					//cout << "probe master status\n";
					MPI::COMM_WORLD.Probe(master, MPI::ANY_TAG, status);
					int l1 = status.Get_count(MPI::CHAR);
					char * pbuff = new char[l1];
					//cout << "Receive packet\n";
					MPI::COMM_WORLD.Recv(pbuff,l1,MPI::CHAR,master, MPI::ANY_TAG,status);
					//cout << "received pbuff value: " + string(pbuff) + "\n";
					if(status.Get_tag() !=0 ){

						MPI::COMM_WORLD.Probe(master, MPI::ANY_TAG, status);
						int l2 = status.Get_count(MPI::CHAR);
						char * dbuff = new char[l2];
						MPI::COMM_WORLD.Recv(dbuff,l2,MPI::CHAR,master, MPI::ANY_TAG,status);
						//cout << "received dbuff value: " + string(dbuff) + "\n";
						if(status.Get_tag() !=0 ){
							int tag = status.Get_tag();
							string pfile(pbuff,l1);
							string dfile(dbuff,l2);
							cout << "running process " + to_string(static_cast<long long>(tag)) + " of " + to_string(static_cast<long long>(totaltrials)) + " on processor " + to_string(static_cast<long long>(myid)) + " : " + pfile.substr(pfile.rfind('/')+1,pfile.size()) + ", " + dfile.substr(dfile.rfind('/')+1,dfile.size())  + "\n";
							//run develep
							runDevelep(pbuff,dbuff,1);

							// send message when finished
							int tmp = 1;
							MPI::COMM_WORLD.Send(&tmp,1,MPI::INT,master,myid);
						}
						else{
							cout << "status tag is zero on process " + to_string(static_cast<long long>(myid)) + "\n";
							cont=false;

						}


						delete [] dbuff;
					}
					else{
						cout << "status tag is zero on process " + to_string(static_cast<long long>(myid)) + "\n";
						cont=false;
					}

					delete [] pbuff;

				}
        	}


        }





		MPI::Finalize();
		char key;
		if(myid==master)
			cout << "All trials completed. Exiting..." << endl;
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
		cout << "Exception Occurred."<<endl;
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
	boost::regex re("~");
	int n=0;
	while(!fs.eof())
	{
		getline(fs,s,'\n');
		istringstream ss(s);
		if (s.compare("")!=0){
			//cout << "s: " << s << endl;
			ss >> tmpi;
			trialset.push_back(tmpi);
			//cout << "trialset " << n << ": " << trialset.back() << endl;
			ss >> tmps;
			paramset.push_back(tmps);
			//cout << "paramset " << n << ": " << paramset.back() << endl;
			ss >> tmps;
			dataset.push_back(tmps);
			//cout << "dataset " << n << ": " << dataset.back() << endl;
			n++;
		}
	}
    //cout << "trialset size: " + to_string(static_cast<long long>(trialset.size())) + "\n";

	for(unsigned int i=0;i<trialset.size();i++)
	{
		//cout << "trialset[" + to_string(static_cast<long long>(i)) + "] = " + to_string(static_cast<long long>(trialset[i])) + "\n";
		for (unsigned int j = 0; j<trialset[i];j++)
		{
			paramfile.push_back(paramset[i]);
			datafile.push_back(dataset[i]);
		}
		totaltrials+=trialset[i];
	}
}
