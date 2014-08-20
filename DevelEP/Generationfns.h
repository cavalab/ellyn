#pragma once
#ifndef GENERATIONFNS_H
#define GENERATIONFNS_H
//#include "FitnessEstimator.h"
void Tournament(vector<ind>&,vector<unsigned int>&,params&,vector<Randclass>& r);
void DC(vector<ind>&,params&,vector<Randclass>& r,data& d,state&,FitnessEstimator& FE);
void ApplyGenetics(vector<ind>&,vector<unsigned int>&,params& p,vector<Randclass>& r,data& d,state& s,FitnessEstimator& FE);
void Crossover(ind&,ind&,vector<ind>&,params&,vector<Randclass>& r);
void CrossoverP(ind&,ind&,ind&,ind&,params&,vector<Randclass>& r);
void Mutate(ind&,vector<ind>&,params&,vector<Randclass>& r,data& d);
void MutateP(ind&,ind& tmpind,params&,vector<Randclass>& r);
void HillClimb(ind&,params&,vector<Randclass>& r,data& d,state& s,FitnessEstimator& FE);
void EpiHC(ind&,params&,vector<Randclass>& r,data& d,state& s,FitnessEstimator& FE);
void AgeBreed(vector<ind>& pop,params& p,vector<Randclass>& r,data& d,state&,FitnessEstimator& FE);
void AgeFitSelect(vector<ind>& pop,params& p,vector<Randclass>& r);
void LexicaseSelect(vector<ind>& pop,vector<unsigned int>& parloc,params& p,vector<Randclass>& r);
#endif