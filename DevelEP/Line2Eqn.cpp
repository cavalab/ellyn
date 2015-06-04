#include "stdafx.h"
#include "op_node.h"
#include "params.h"
using namespace std;

string Line2Eqn(vector<node>& line,string& eqnForm,params& p)
{
	vector<string> eqn_floatstack;
	vector<string> eqn_binstack;
	vector<string> form_floatstack;
	vector<string> form_binstack;

	for(unsigned int i=0;i<line.size();++i)
	{
		if(line.at(i).on){
			if(line.at(i).type=='n'){
				eqn_floatstack.push_back(to_string(static_cast<long double>(line.at(i).value)));
				form_floatstack.push_back("n");
			}
			else if(line.at(i).type=='v'){
				eqn_floatstack.push_back(line.at(i).varname);
				form_floatstack.push_back(line.at(i).varname);
			}
			else{
				string sop;
				string sop2; // for protected functions
				char tmp = line.at(i).type;
				if(line.at(i).arity==1 && eqn_floatstack.size()>=1){
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
					case 'q':
						sop="sqrt";
						break;
					default:
						cout << "error";
					}

					string s1 = eqn_floatstack.back(); eqn_floatstack.pop_back();
					string s1f = form_floatstack.back(); form_floatstack.pop_back();
					s1.insert(0,"(");
					s1f.insert(0,"(");
					s1+=")";
					s1f+=")";
					//if (sop.compare("log")==0 ){ //|| sop.compare("sqrt")==0
					//	eqn_floatstack.push_back(sop + "(abs" + s1 + ")");
					//	form_floatstack.push_back(sop + "(abs" + s1f + ")");
					//}
					//else{
					eqn_floatstack.push_back(sop + s1);
					form_floatstack.push_back(sop + s1f);
					//}
				}
				else if (line.at(i).arity==2 && ((!line.at(i).op_takes_bool && eqn_floatstack.size()>=2) || (line.at(i).op_takes_bool && eqn_floatstack.size()>=1 && eqn_binstack.size()>=1))){
					string s1, s2, s1f, s2f;
					if (line[i].op_takes_bool){
						s1 = eqn_binstack.back(); eqn_binstack.pop_back();
						s1f = form_binstack.back(); form_binstack.pop_back();
					}
					else{
						s1 = eqn_floatstack.back(); eqn_floatstack.pop_back();
						s1f = form_floatstack.back(); form_floatstack.pop_back();
					}

					s2 = eqn_floatstack.back(); eqn_floatstack.pop_back();
					s2f = form_floatstack.back(); form_floatstack.pop_back();
										
					string sop(&line.at(i).type); sop = sop[0];
					switch (line.at(i).type){
					case '{':
						sop = "<=";
						break;
					case '}':
						sop = ">=";
						break;
					case 'i':
						sop = "IF-THEN";
						break;
					default:
						break;
					}
					

					if (line.at(i).type=='i'){
						eqn_floatstack.push_back("(IF-THEN("+s1+","+s2+"))");
						form_floatstack.push_back("(IF-THEN("+s1f+","+s2f+"))");
					}
					else{
						if (line[i].op_returns_bool){
							eqn_binstack.push_back("("+s2+sop+s1+")");
							form_binstack.push_back("("+s2f+sop+s1f+")");
						}
						else{
							eqn_floatstack.push_back("("+s2+sop+s1+")");
							form_floatstack.push_back("("+s2f+sop+s1f+")");
						}
					}
				}
				else if (line.at(i).arity==3 && eqn_binstack.size()>=1 && eqn_floatstack.size()>=2){
					string s1 = eqn_binstack.back(); eqn_binstack.pop_back();
					string s2 = eqn_floatstack.back(); eqn_floatstack.pop_back();
					string s3 = eqn_floatstack.back(); eqn_floatstack.pop_back();
					string s1f = form_binstack.back(); form_binstack.pop_back();
					string s2f = form_floatstack.back(); form_floatstack.pop_back();
					string s3f = form_floatstack.back(); form_floatstack.pop_back();
					
					switch (line.at(i).type){
					case 't':
						eqn_floatstack.push_back("(IF-THEN-ELSE("+s1+","+s3+","+s2+"))");
						form_floatstack.push_back("(IF-THEN-ELSE("+s1f+","+s3f+","+s2f+"))");
						break;
					default:
						std:cerr << "unspecified arity 3 operator in Line2Eqn";
						throw;
						break;
					}
				}
					
					//cout <<"arity screwed up.\n";
			}	
		}
	}
	if ((eqn_floatstack.empty() && !p.classification) || (eqn_binstack.empty() && p.classification) )
	{
		eqnForm = "unwriteable";
		//cout << "equation stack empty.\n";
		return "unwriteable";
	}
	if (p.classification){
		eqnForm = form_binstack.back();
		return eqn_binstack.back();
	}
	else{
		eqnForm = form_floatstack.back();
		return eqn_floatstack.back();
	}
}