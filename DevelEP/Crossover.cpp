#include "stdafx.h"
#include "pop.h"
#include "params.h"
#include "rnd.h"
#include "data.h"
#include "general_fns.h"

void Crossover(ind& p1,ind& p2,vector<ind>& tmppop,params& p,vector<Randclass>& r)
{
	vector<ind> parents;
	parents.push_back(p1);
	parents.push_back(p2);

	//makenew(parents[0]);
	//makenew(parents[1]);

	vector<ind> kids(2);

	int r2,off1,off2,offset,head,it;
	//std::vector<int>::iterator it;
	int tmpinssize = 0;
	vector<int> psize; 
	psize.push_back(p1.line.size());
	psize.push_back(p2.line.size());

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
			// assign beginning of parent to kid if it is longer
			kids.at(r1).line.insert(kids.at(r1).line.end(),parents[r1].line.begin(),parents[r1].line.begin()+offset);

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
					
				kids.at(r1).line.push_back(parents[head].line.at(i+offset));
			}
			if(kids.at(r1).line.size() < parents[r1].line.size())
			{
				int gap = kids.at(r1).line.size() < parents[r1].line.size()+1;
				kids.at(r1).line.insert(kids.at(r1).line.end(),parents[r1].line.end()-gap,parents[r1].line.end());
			}
				//kids.at(r1).line.push_back(parents[r1].line.at(kids.at(r1).line.size()));
		}
	}
	else if (p.cross==2) // one-point crossover
	{
		int point1 = r[omp_get_thread_num()].rnd_int(0,min(p1.line.size(),p2.line.size()));
		
		kids[0].line.assign(parents[0].line.begin(),parents[0].line.begin()+point1);
		kids[0].line.insert(kids[0].line.end(),parents[1].line.begin()+point1,parents[1].line.end());

		kids[1].line.assign(parents[1].line.begin(),parents[1].line.begin()+point1);
		kids[1].line.insert(kids[1].line.end(),parents[0].line.begin()+point1,parents[0].line.end());

	}
	//tmpinssize=0;
	//for(unsigned int t=0; t<kids[0].line.size();t++)
	//	if(kids[0].line.at(t)>99) tmpinssize++;
	//if(tmpinssize!=kids[0].args.size())
	//	cout << "size mismatch" << endl;

	tmpinssize=0;

	kids[0].origin = 'c';
	kids[0].parentfitness = parents[0].fitness;
	kids[1].origin = 'c';
	kids[1].parentfitness = parents[1].fitness;

	//makenew(kids[0]);
	//makenew(kids[1]);

	tmppop.push_back(kids[0]);
	tmppop.push_back(kids[1]);

	kids.clear();
	parents.clear();
}