// header file for ind struct
#pragma once
#ifndef DATA_H
#define DATA_H
#include <iostream>
#include <string>
#include <vector>
#include <random>
using namespace std;


//typedef exprtk::parser_error::type error_t;

struct data {

	vector<string> label; // variables corresponding to those defined in parameters 
	vector<vector<float>> vals; //2D vector of all data, can be accessed data[row][column]
	vector<float> dattovar;
	vector<float> target;

	void clear()
	{
		label.clear();
		vals.clear();
		dattovar.clear();
		target.clear();
	}

};

#endif