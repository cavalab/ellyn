#pragma once
#ifndef GENERAL_FNS_H
#define GENERAL_FNS_H
#include "pop.h"
#include "params.h"

bool is_number(const std::string& s);
void NewInstruction(ind& newind,int loc,params&,vector<Randclass>& r);
#endif