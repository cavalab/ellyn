#pragma once
#ifndef GENERAL_FNS_H
#define GENERAL_FNS_H
#include "stdafx.h"
#include "pop.h"
#include "params.h"

bool is_number(const std::string& s);
void NewInstruction(ind& newind,int loc,params&,vector<Randclass>& r,data& d);
void makenewcopy(ind& newind);
void copystack(vector<shared_ptr<node>>& line, vector<shared_ptr<node>>& newline);


#endif