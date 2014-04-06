// RPN_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "pop.h"
//#include "InitPop.h"
#include "rnd.h"
#include <unordered_map>
#include "op_node.h"
using namespace std;

string Line2Eqn(vector<shared_ptr<node>>& line)
{
	vector<string> eqnstack;

	for(unsigned int i=0;i<line.size();i++)
	{
		if(line.at(i)->type=='n')
			eqnstack.push_back(to_string(static_cast<long double>(static_pointer_cast<n_num>(line.at(i))->value)));
		else if(line.at(i)->type=='v')
			eqnstack.push_back(static_pointer_cast<n_sym>(line.at(i))->varname);
		else{
			string sop;
			char tmp = line.at(i)->type;
			if(line.at(i)->arity==1 && eqnstack.size()>=1){
				switch (line.at(i)->type){
				case 'l':
					sop="log";
					break;
				case 's':
					sop="sin";
					break;
				case 'c':
					sop="cos";
					break;
				case 'e':
					sop="exp";
					break;
				default:
					cout << "error";
				}

				string s1 = eqnstack.back(); eqnstack.pop_back();
				s1.insert(0,"(");
				s1+=")";
				eqnstack.push_back(sop + s1);
			}
			else if (line.at(i)->arity==2 && eqnstack.size()>=2){
				string s1 = eqnstack.back(); eqnstack.pop_back();
				string s2 = eqnstack.back(); eqnstack.pop_back();
				eqnstack.push_back("("+s2+line.at(i)->type+s1+")");
			}
			//else
				//cout <<"arity screwed up.\n";
		}			
	}
	if (eqnstack.empty())
	{
		cout << "equation stack empty.\n";
		return "unwriteable";
	}
	
	return eqnstack.front();
}
int _tmain(int argc, _TCHAR* argv[])
{

	// Initialize population and test fitness
	int size=1000;
	float maxf=2000000;
	float minf=0.00000001;
	//tribe T(size,maxf,minf);
	unordered_map <string,float*> datatable;
	float x = 15.198;
	float y = -1.3450;
	bool ERCints = 0;

	datatable.insert(pair<string,float*>("x",&x));
	datatable.insert(pair<string,float*>("y",&y));
	vector <string> label;
	label.push_back("x");
	label.push_back("y");

	Randclass r;
	vector<int> op_choice;

	for (int i=0;i<10;i++)
		op_choice.push_back(i);
	vector<ind> pop(100);
	//vector<node*> line;
	int choice;
	string varchoice;
	for(int i=0;i<pop.size();i++)
	{
		
		for(int j=0;j<20;j++)
		{
			choice = r.rnd_int(0,op_choice.size()-1);
		
			switch (op_choice.at(choice))
			{
			case 0: //load number
				if(ERCints)
					pop.at(i).line.push_back(shared_ptr<node>(new n_num((float)r.rnd_int(-10,10))));
				else
					pop.at(i).line.push_back(shared_ptr<node>(new n_num(r.rnd_flt(-10,10))));
				break;
			case 1: //load variable
				varchoice = label.at(r.rnd_int(0,label.size()-1));
				pop.at(i).line.push_back(shared_ptr<node>(new n_sym(datatable.at(varchoice),varchoice)));
				break;
			case 2: // +
				pop.at(i).line.push_back(shared_ptr<node>(new n_add()));
				break;
			case 3: // -
				pop.at(i).line.push_back(shared_ptr<node>(new n_sub()));
				break;
			case 4: // *
				pop.at(i).line.push_back(shared_ptr<node>(new n_mul()));
				break;
			case 5: // /
				pop.at(i).line.push_back(shared_ptr<node>(new n_div()));
				break;
			case 6: // sin
				pop.at(i).line.push_back(shared_ptr<node>(new n_sin()));
				break;
			case 7: // cos
				pop.at(i).line.push_back(shared_ptr<node>(new n_cos()));
				break;
			case 8: // exp
				pop.at(i).line.push_back(shared_ptr<node>(new n_exp()));
				break;
			case 9: // log
				pop.at(i).line.push_back(shared_ptr<node>(new n_log()));
				break;
			}
		}
		//if (choice==0) // number reference to argument from all blocks in next position in line
		//{
		//	if(ERCints)
		//		line.push_back(ops(choice,r.rnd_int((int)minf,(int)maxf)));
		//	else
		//		line.push_back(ops(choice,r.rnd_flt(minf,maxf)));
		//}
		//else if (choice==1) // pick pointer values from mapped variables in data struct
		//{
		//	unordered_map<string,float*>::iterator it = datatable.begin();
		//	it+=r.rnd_int(0,datatable.size());
		//	line.push_back(ops(choice,it.second));

		//}
		//else 
		//	line.push_back(choice);
		

		pop.at(i).eqn = Line2Eqn(pop.at(i).line);
		cout << "Equation" << i << ": f=" << pop.at(i).eqn << "\n";
		vector<float> outstack;
		for (int j=0;j<1000;j++)
		{
			x = float(j)/1000;
			y = float(j);

			for(int k=0;k<pop.at(i).line.size();k++)
				pop.at(i).line.at(k)->eval(outstack);

			if(!outstack.empty())
				pop.at(i).output.push_back(outstack.front());
		/*if(boost::math::isinf(outstack.front())){
			char blah;
			cout << "inifinite result. waiting on you...\n";
			cin >> blah;
		}*/
			//cout << "value: " << outstack.front() << endl;
			outstack.clear();
		}
	}
	return 0;
}


