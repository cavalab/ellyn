// header file for ind struct
#pragma once
#ifndef DATA_H
#define DATA_H
#include <iostream>
#include <string>
#include <vector>
#include <random>
//#include "exprtk.hpp"
//#include <unordered_map>

using namespace std;


//typedef exprtk::parser_error::type error_t;

struct data {

	vector<string> label; // variables corresponding to those defined in parameters 
	vector<vector<float>> vals; //2D vector of all data, can be accessed data[row][column]
	//vector<float> dattovar;
	vector<float> target;
	/*unordered_map <string,float*> datatable;*/
	//mymap.
	/*mymap.insert(pair<string,int*>("alpha",10));
	mymap.insert(pair<string,int>("beta",20));
	mymap.insert(pair<string,int>("G1",30));*/

	void clear()
	{
		label.clear();
		vals.clear();
	//	dattovar.clear();
		target.clear();
	//	datatable.clear();
	}
	/*void mapdata()
	{
		for (unsigned int i=0;i<label.size(); i++)
			datatable.insert(pair<string,float*>(label[i],&dattovar[i]));
	}*/
	~data() {}
};

#endif