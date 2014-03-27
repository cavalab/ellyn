#include "stdafx.h"

enum {NUM,OP,SYM};

class node{
public:
	int type;
	node() {type=0;}
	node(int set) {type=set;}
	~node() {}
	virtual void eval(vector<float> & outstack);
};

class n_sym : public node{
	float* valpt;
public: 
	n_sym(float& set)
	{
		valpt = &set;
	}
	~n_sym(){}
	virtual void eval(vector<float> & outstack)
	{
		outstack.push_back(*valpt);
	}
};

class n_num : public node{
public:
	float value;
	n_num(float set) {value=set;}
	virtual void eval(vector<float> & outstack)	{outstack.push_back(value);}
};

class n_add: public node{
public:
	n_add(){}
	~n_add(){}
	virtual void eval(vector<float> & outstack)
	{
		float n1 = outstack.back(); outstack.pop_back();
		float n2 = outstack.back(); outstack.pop_back();

		outstack.push_back(n1+n2);
	}
};

class n_sub: public node{
public:
	n_sub(){}
	~n_sub(){}
	virtual void eval(vector<float> & outstack)
	{
		float n1 = outstack.back(); outstack.pop_back();
		float n2 = outstack.back(); outstack.pop_back();

		outstack.push_back(n1-n2);
	}
};
class n_mul: public node{
public:
	n_mul(){}
	~n_mul(){}
	virtual void eval(vector<float> & outstack)
	{
		float n1 = outstack.back(); outstack.pop_back();
		float n2 = outstack.back(); outstack.pop_back();

		outstack.push_back(n1*n2);
	}
};
class n_div: public node{
public:
	n_div(){}
	~n_div(){}
	virtual void eval(vector<float> & outstack)
	{
		float n1 = outstack.back(); outstack.pop_back();
		float n2 = outstack.back(); outstack.pop_back();

		outstack.push_back(n1+n2);
	}
};

class n_sin: public node{
public:
	n_sin(){}
	~n_sin(){}
	virtual void eval(vector<float> & outstack)
	{
		float n1 = outstack.back(); outstack.pop_back();
		outstack.push_back(sin(n1));
	}
};

class n_cos: public node{
public:
	n_cos(){}
	~n_cos(){}
	virtual void eval(vector<float> & outstack)
	{
		float n1 = outstack.back(); outstack.pop_back();
		outstack.push_back(cos(n1));
	}
};

class n_exp: public node{
public:
	n_exp(){}
	~n_exp(){}
	virtual void eval(vector<float> & outstack)
	{
		float n1 = outstack.back(); outstack.pop_back();
		outstack.push_back(exp(n1));
	}
};

class n_log: public node{
public:
	n_log(){}
	~n_log(){}
	virtual void eval(vector<float> & outstack)
	{
		float n1 = outstack.back(); outstack.pop_back();
		outstack.push_back(log(n1));
	}
};