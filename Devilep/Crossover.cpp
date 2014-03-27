#include "stdafx.h"
#include "pop.h"

void Crossover(ind& p1,ind& p2,vector<ind>& tmppop,params& p,vector<Randclass>& r)
{
	vector<ind> parents;
	parents.push_back(p1);
	parents.push_back(p2);

	vector<ind> kids(2);

	int r2,off1,off2,offset,head,it;
	//std::vector<int>::iterator it;
	int tmpinssize = 0;
	
	if(p.cross==1)
	{
		for (int r1=0; r1 < 2; r1++)
		{
			if(r1==0)
				r2=1;
			else
				r2=0;
		
			if(parents[r1].line.size()>parents[r2].line.size())
			{
				off1 = r[omp_get_thread_num()].rnd_int(0,parents[r1].line.size()-parents[r2].line.size()-1);
				off2 = 0;
			}
			else if (parents[r2].line.size()>parents[r1].line.size())
			{
				off1 = 0;
				off2 = r[omp_get_thread_num()].rnd_int(0,parents[r2].line.size()-parents[r1].line.size()-1);
			}
			else
			{
				off1=0; 
				off2=0;
			}
			head = r1;
			offset=off1;
			it = 0;
			// assign beginning of parent to kid if it is longer
			while(kids.at(r1).line.size() < offset)
			{
				if (parents[r1].line.at(it)>=100)
				{
					kids.at(r1).args.push_back(parents[r1].args.at(parents[r1].line.at(it)-100));
					kids.at(r1).line.push_back(kids.at(r1).args.size()+99);
					
				}
				else
				{
					kids.at(r1).line.push_back(parents[r1].line.at(it));
				}
				if (p.eHC_on)
						kids.at(r1).epiline.push_back(parents[r1].epiline.at(it));
				++it;
			}


			for (unsigned int i=0;i<std::min(parents[r1].line.size(),parents[r2].line.size());i++)
			{
				if (r[omp_get_thread_num()].rnd_flt(0,1)<p.cross_ar)
				{
					if(head==r1)
					{
						head=r2;
						offset=off2;
					}
					else
					{
						head=r1;
						offset=off1;
					}
				}
					
				if (p.eHC_on)
						kids.at(r1).epiline.push_back(parents[r1].epiline.at(it));
			
				if(parents[head].line.at(i+offset)>99)
				{// grab arguments
					kids.at(r1).args.push_back(parents[head].args.at(parents[head].line.at(i+offset)-100));
					kids.at(r1).line.push_back(kids.at(r1).args.size()+99);
				}
				else
					kids.at(r1).line.push_back(parents[head].line.at(i+offset));
			
			/*tmpinssize=0;
			for(unsigned int t=0; t<kids[r1].line.size();t++)
				if(kids[r1].line.at(t)>99) tmpinssize++;
			if(tmpinssize!=kids[r1].args.size())
				cout << "size mismatch" << endl;*/
			}
		
			
			while(kids.at(r1).line.size() < parents[r1].line.size())
			{
				if (parents[r1].line.at(kids.at(r1).line.size())>=100)
				{
					kids.at(r1).args.push_back(parents[r1].args.at(parents[r1].line.at(kids.at(r1).line.size())-100));
					kids.at(r1).line.push_back(kids.at(r1).args.size()+99);
				}
				else
					kids.at(r1).line.push_back(parents[r1].line.at(kids.at(r1).line.size()));
				
				if (p.eHC_on)
						kids.at(r1).epiline.push_back(parents[r1].epiline.at(it));
			}
			if (p.eHC_on && kids[r1].epiline.size() != kids[r1].line.size())
				cout << "size mismatch between epiline and line\n";
		}
	}
	else if (p.cross==2) // two-point crossover
	{

	}
	//tmpinssize=0;
	//for(unsigned int t=0; t<kids[0].line.size();t++)
	//	if(kids[0].line.at(t)>99) tmpinssize++;
	//if(tmpinssize!=kids[0].args.size())
	//	cout << "size mismatch" << endl;

	tmpinssize=0;

	for(unsigned int t=0; t<kids[1].line.size();t++)
		if(kids[1].line.at(t)>99) tmpinssize++;
	if(p.eHC_on && tmpinssize!=kids[1].args.size())
		cout << "size mismatch" << endl;

	kids[0].origin = 'c';
	kids[0].parentfitness = parents[0].fitness;
	kids[1].origin = 'c';
	kids[1].parentfitness = parents[1].fitness;

	if(!p.eHC_on) // set epiline to all ones 
	{
		kids[0].epiline.assign(kids[0].line.size(),1);
		kids[1].epiline.assign(kids[1].line.size(),1);
	}
	tmppop.push_back(kids[0]);
	tmppop.push_back(kids[1]);

	if (kids[0].epiline.size() != kids[0].line.size())
				cout << "size mismatch between epiline and line\n";
	if (kids[1].epiline.size() != kids[1].line.size())
				cout << "size mismatch between epiline and line\n";
}