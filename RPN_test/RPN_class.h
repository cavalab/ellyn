#include "stdafx.h"
#include <functional>

using namespace std::placeholders;

void e_num(vector<float> & outstack,float& val)
{
	outstack.push_back(val);
}
void e_var(vector<float>& outstack,float* valpt)
{
	outstack.push_back(*valpt);
}
float rec_exponent(float a1, float a2)
{
	return a2<0 ? 0 : (a2==0?1:a1*rec_exponent(a1, a2-1));
}

void e_exponent(vector<float>& outstack)
{
	if(outstack.size()>=2){
		float a2 = outstack.back(); outstack.pop_back();
		float a1 = outstack.back(); outstack.pop_back();
	
		outstack.push_back(rec_exponent(a1,a2));
	}	
}

void e_mul(vector<float>& outstack) 
{
	if(outstack.size()>=2){
		float a2 = outstack.back(); outstack.pop_back();
		float a1 = outstack.back(); outstack.pop_back();
		outstack.push_back(a1*a2);
	}
}
void e_div(vector<float>& outstack) 
{
	if(outstack.size()>=2){
		float a2 = outstack.back(); outstack.pop_back();
		float a1 = outstack.back(); outstack.pop_back();
		if(!a2) {
			outstack.push_back(0);
		}
		outstack.push_back(a1/a2);
	}
}
void e_mod(vector<float>& outstack) 
{
	if(outstack.size()>=2){
		float a2 = outstack.back(); outstack.pop_back();
		float a1 = outstack.back(); outstack.pop_back();
		if(!a2) {
			//fprintf(stderr, "ERROR: Division by zero\n");
			//exit(EXIT_FAILURE);
			outstack.push_back(0.0);
		}
		outstack.push_back(float(int(a1)%int(a2)));
	}
}
void e_add(vector<float>& outstack) 
{
	if(outstack.size()>=2){
		float a2 = outstack.back(); outstack.pop_back();
		float a1 = outstack.back(); outstack.pop_back();
		outstack.push_back(a1+a2);
	}
}
void e_sub(vector<float>& outstack) 
{
	if(outstack.size()>=2){
		float a2 = outstack.back(); outstack.pop_back();
		float a1 = outstack.back(); outstack.pop_back();
		outstack.push_back(a1-a2);
	}
}
void e_sin(vector<float>& outstack)
{
	if(outstack.size()>=1){
		float a1 = outstack.back(); outstack.pop_back();
		outstack.push_back(sin(a1));
	}
}
void e_cos(vector<float>& outstack)
{
	if(outstack.size()>=1){
		float a1 = outstack.back(); outstack.pop_back();
		outstack.push_back(cos(a1));
	}
}
void e_exp(vector<float>& outstack)
{
	if(outstack.size()>=1){
		float a1 = outstack.back(); outstack.pop_back();
		outstack.push_back(exp(a1));
	}
}
void e_log(vector<float>& outstack)
{
	if(outstack.size()>=1){
		float a1 = outstack.back(); outstack.pop_back();
		outstack.push_back(log(a1));
	}
}

struct ops {
	char type;
	bool on;
	float val;
	float* valpt;
	//typedef std::function<void> evalfn;
	//evalfn eval;
	void (*eval)(vector<float>&);
	void (*loadn)(vector<float>&,float&);
	void (*loadv)(vector<float>&,float*);
	ops(int setfn)
	{
		on=true;

		switch (setfn) { //choose instruction to point pointer function to.
		case 0:
			//eval = std::bind(e_num,_1,val);
			loadn = e_num;
			type = 'n';
			break;
		case 1:
			//eval = std::bind(e_var,_1,valpt);
			loadv = e_var;
			type = 'v';
			break;
		case 2:
			//eval = std::bind(e_exponent,_1);
			eval = e_exponent;
			type = '^';
			break;
		case 3:
			//eval = std::bind(e_mul,_1);
			eval=e_mul;
			type = '*';
			break;
		case 4:
			//eval = std::bind(e_div,_1);
			eval = e_div;
			type = '/';
			break;
		case 5:
			//eval = std::bind(e_mod,_1);
			eval=e_mod;
			type = '%';
			break;
		case 6:
			//eval = std::bind(e_add,_1);
			eval=e_add;
			type = '+';
			break;
		case 7:
			//eval = std::bind(e_sub,_1);
			eval=e_sub;
			type = '-';
			break;
		case 8: 
			//eval = std::bind(e_sin,_1);
			eval=e_sin;
			type = 's';
			break;
		case 9: 
			//eval = std::bind(e_cos,_1);
			eval=e_cos;
			type = 'c';
			break;
		case 10: 
			//eval = std::bind(e_exp,_1);
			eval = e_exp;
			type = 'e';
			break;
		case 11: 
			//eval = std::bind(e_log,_1);
			eval=e_log;
			type = 'l';
			break;
		}
	}
	ops(int setfn,int v)
	{
		on=true;
		if (setfn==0)
			val = (float)v;
		else
			cout << "Wrong operator type for float argument. \n";
	}
	ops(int setfn,float& v)
	{
		on=true;
		if (setfn==0)
			val = v;
		else
			cout << "Wrong operator type for float argument. \n";
	}
	ops(int setfn,float* v)
	{
		on=true;
		if (setfn==1)
			valpt = v;
		else
			cout << "Wrong operator type for float pointer argument. \n";
	}
	~ops() {}

	
};