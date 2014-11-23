#include "stdafx.h"
#include "op_node.h"

using namespace std;

string Line2Eqn(vector<node>& line)
{
	vector<string> eqnstack;

	for(unsigned int i=0;i<line.size();i++)
	{
		if(line.at(i).on){
			if(line.at(i).type=='n')
				eqnstack.push_back(to_string(static_cast<long double>(line.at(i).value)));
			else if(line.at(i).type=='v')
				eqnstack.push_back(line.at(i).varname);
			else{
				string sop;
				char tmp = line.at(i).type;
				if(line.at(i).arity==1 && eqnstack.size()>=1){
					switch (line.at(i).type){
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
				else if (line.at(i).arity==2 && eqnstack.size()>=2){
					string s1 = eqnstack.back(); eqnstack.pop_back();
					string s2 = eqnstack.back(); eqnstack.pop_back();
					eqnstack.push_back("("+s2+line.at(i).type+s1+")");
				}
				//else
					//cout <<"arity screwed up.\n";
			}	
		}
	}
	if (eqnstack.empty())
	{
		//cout << "equation stack empty.\n";
		return "unwriteable";
	}
	
	return eqnstack.back();
}