#include "stdafx.h"

enum {NUM,OP,SYM};

class node{
public:
	int type;
	node(int set) {type=set;}
	~node() {}
	virtual void eval(vector<float> & outstack);
};

class symbol : public node{
	float* valpt;
public: 
	symbol(float& set)
	{
		valpt = &set;
	}
	virtual void eval(vector<float> & outstack)
	{
		outstack.push_back(*valpt);
	}
};

class number : public node{
public:
	float value;
	number(float set) {value=set;}
	virtual void eval(vector<float> & outstack)	{outstack.push_back(value);}
};

class addition: public node{
public:
	virtual void eval(vector<float> & outstack)
	{
		float n1 = outstack.back(); outstack.pop_back();
		float n2 = outstack.back(); outstack.pop_back();

		outstack.push_back(n1+n2);
	}
};


