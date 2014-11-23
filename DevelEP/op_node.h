#pragma once
#ifndef OP_NODE_H
#define OP_NODE_H
//#include "stdafx.h"
//#include <vector>
using namespace std;
//#include <boost/ptr_container/ptr_vector.hpp>
//#include <boost/utility.hpp>
//enum {NUM,OP,SYM};

class node{
public:
	char type;
	bool on;
	int arity;
	float value;
	float* valpt;
	string varname;
	node() {type=0; on=1; arity=0;}
	//node(int set) {type=set;on=1;}
	//operator with specified arity
	node(char stype,int sarity){type=stype;arity=sarity;on=1;}
	//operator with arity lookup
	node(char stype)
	{
		type=stype;
		if (type=='s' || type=='c' || type=='e' || type=='l')
			arity = 1;
		else
			arity = 2;

		on=1;
	}
	//number
	node(float svalue){type='n'; value=svalue; on=1; arity=0;}
	//variable
	node(string& vname){type='v';varname=vname;on=1;}
	//set pointer for variables
	void setpt(float* set)
	{
		if (set==NULL)
			cout<<"bad set";
		valpt=set;
	}
	~node() {}
	//void eval(vector<float> & outstack)=0;
	/*node* clone() const
	{
		return do_clone();
	}*/
//private:
	/*virtual node* do_clone() const
	{
		return new node(*this); 
	}*/
};

//inline node* new_clone(const node& a)
//{
//	return a.clone();
//}

class n_num : public node{
public:
	float value;
	n_num(){on=1;arity=0;}
	n_num(float set) {type='n'; value=set; on=1; arity=0;}
	~n_num(){}
	virtual void eval(vector<float> & outstack)	
	{
		//if (on)
			outstack.push_back(value);
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_num(*this); 
//	}
};

class n_sym : public node{
public: 
	float* valpt;
	string varname;
	n_sym(){on=1; arity=0;};
	n_sym(string& vname)
	{
		type='v';
		varname=vname;
		on=1;
	}
	~n_sym(){}
	void setpt(float* set)
	{
		if (set==NULL)
			cout<<"bad set";

		valpt=set;
	}
	virtual void eval(vector<float> & outstack)
	{
		if (valpt==NULL)
			cout<<"problem";
		else
			outstack.push_back(*valpt);
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_sym(*this); 
//	}
};

class n_add: public node{
public:
	n_add(){type='+';arity=2;on=1;}
	~n_add(){}
	virtual void eval(vector<float> & outstack)
	{
		//if (on){
			if(outstack.size()>=2){
				float n1 = outstack.back(); outstack.pop_back();
				float n2 = outstack.back(); outstack.pop_back();

				outstack.push_back(n2+n1);
			}
		//}
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_add(*this); 
//	}
};

class n_sub: public node{
public:
	n_sub(){type='-';arity=2;on=1;}
	~n_sub(){}
	virtual void eval(vector<float> & outstack)
	{
		//if (on){
			if(outstack.size()>=2){
				float n1 = outstack.back(); outstack.pop_back();
				float n2 = outstack.back(); outstack.pop_back();

				outstack.push_back(n2-n1);
			}
		//}
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_sub(*this); 
//	}
};
class n_mul: public node{
public:
	n_mul(){type='*';arity=2;on=1;}
	~n_mul(){}
	virtual void eval(vector<float> & outstack)
	{
		//if (on){
			if(outstack.size()>=2){
				float n1 = outstack.back(); outstack.pop_back();
				float n2 = outstack.back(); outstack.pop_back();

				outstack.push_back(n2*n1);
			}
		//}
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_mul(*this); 
//	}
};
class n_div: public node{
public:
	n_div(){type='/';on=1;arity=2;}
	~n_div(){}
	virtual void eval(vector<float> & outstack)
	{
		//if (on){
			if(outstack.size()>=2){
				float n1 = outstack.back(); outstack.pop_back();
				float n2 = outstack.back(); outstack.pop_back();
				if(abs(n1)<0.0001)
					outstack.push_back(0);
				else
					outstack.push_back(n2/n1);
			}
		//}
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_div(*this); 
//	}
};

class n_sin: public node{
public:
	n_sin(){type='s';arity=1;on=1;}
	~n_sin(){}
	virtual void eval(vector<float> & outstack)
	{
		//if (on){
			if(outstack.size()>=1){
				float n1 = outstack.back(); outstack.pop_back();
				outstack.push_back(sin(n1));
			}
		//}
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_sin(*this); 
//	}

};

class n_cos: public node{
public:
	n_cos(){type='c';arity=1;on=1;}
	~n_cos(){}
	virtual void eval(vector<float> & outstack)
	{
		//if (on){
			if(outstack.size()>=1){
				float n1 = outstack.back(); outstack.pop_back();
				outstack.push_back(cos(n1));
			}
		//}
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_cos(*this); 
//	}
};

class n_exp: public node{
public:
	n_exp(){type='e';arity=1;on=1;}
	~n_exp(){}
	virtual void eval(vector<float> & outstack)
	{
		//if (on){
			if(outstack.size()>=1){
				float n1 = outstack.back(); outstack.pop_back();
				outstack.push_back(exp(n1));
			}
		//}
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_exp(*this); 
//	}
};

class n_log: public node{
public:
	n_log(){type='l';arity=1;on=1;}
	~n_log(){}
	virtual void eval(vector<float> & outstack)
	{
		//if (on){
			if(outstack.size()>=1){
				float n1 = outstack.back(); outstack.pop_back();
				if (abs(n1)<0.0001)
					outstack.push_back(0);
				else
					outstack.push_back(log(n1));
			}
		//}
	}
//private:
//	virtual node* do_clone() const
//	{
//		return new n_log(*this); 
//	}
};

#endif