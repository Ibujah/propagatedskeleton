/*
Copyright (c) 2016 Bastien Durix

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


/**
 *  \file ThetaMedialAxis.cpp
 *  \brief Defines scale axis transform, a pruning algorithm
 *  \author Bastien Durix
 */

#include "ThetaMedialAxis.h"

unsigned int connectedTo(const skeleton::GraphSkel2d::Ptr grsimp, unsigned int ind, std::map<unsigned int,bool> &connected)
{
	std::list<unsigned int> list_ind;
	grsimp->getAllNodes(list_ind);
	for(std::list<unsigned int>::iterator it = list_ind.begin(); it != list_ind.end(); it++)
	{
		connected.insert(std::make_pair(*it,false));
	}

	std::list<unsigned int> front;
	front.push_back(ind);
	while(!front.empty())
	{
		unsigned int cur = *(front.begin());
		front.pop_front();
		if(!connected[cur])
		{
			connected[cur] = true;
			std::list<unsigned int> l_nod;
			// get the neighbors of the node
			grsimp->getNeighbors(cur,l_nod);
			for(std::list<unsigned int>::iterator it = l_nod.begin(); it != l_nod.end(); it++)
			{
				front.push_back(*it);
			}
		}
	}
	unsigned int cpt = 0;
	for(std::map<unsigned int,bool>::iterator it = connected.begin(); it != connected.end(); it++)
		if(it->second)
			cpt++;
	return cpt;
}

double computeTheta(const Eigen::Vector2d &C1, double r1, const Eigen::Vector2d &C2, double r2)
{
	double dist = (C1 - C2).norm();
	double cosphi = (dist*dist + r1*r1 - r2*r2)/(2.0*r1*dist);
	double phi = acos(cosphi);
	double theta = 2.0*phi;
	return theta;
}

skeleton::GraphSkel2d::Ptr algorithm::pruning::ThetaMedialAxis(const skeleton::GraphSkel2d::Ptr grskel, const double &theta)
{
	typename skeleton::GraphSkel2d::Ptr grsimp(new skeleton::GraphSkel2d(*grskel));

	std::list<unsigned int> list_ind;
	grsimp->getAllNodes(list_ind);

	std::list<std::tuple<unsigned int,Eigen::Vector2d,double> > list_ind_size;
	for(std::list<unsigned int>::iterator it = list_ind.begin(); it != list_ind.end(); it++)
	{
		Eigen::Vector3d vec = grskel->getNode(*it);
		list_ind_size.push_back(std::make_tuple(*it,vec.block<2,1>(0,0),vec(2)));
	}

	//sort nodes in decreasing order
	list_ind_size.sort([](const std::tuple<unsigned int,Eigen::Vector2d,double> &p1, const std::tuple<unsigned int,Eigen::Vector2d,double> &p2)
			{
				return std::get<2>(p1) > std::get<2>(p2);
			});
	
	bool fini = true;
	do{
		fini = true;
		std::list<std::tuple<unsigned int,Eigen::Vector2d,double> >::iterator it = list_ind_size.begin();
		while(it != list_ind_size.end())
		{
			try
			{
				unsigned int refkey = std::get<0>(*it);
				std::list<unsigned int> l_nod;

				// get the neighbors of the node
				grsimp->getNeighbors(refkey,l_nod);

				bool removed = false;
				for(std::list<unsigned int>::iterator itno = l_nod.begin(); itno != l_nod.end() && !removed; itno++)
				{
					Eigen::Vector3d vec = grskel->getNode(*itno);
					Eigen::Vector2d C2 = vec.block<2,1>(0,0);
					double r2 = vec(2);
					double ang = computeTheta(std::get<1>(*it),std::get<2>(*it),C2,r2);
					if(ang < theta)
					{
						grsimp->remEdge(refkey,*itno);
						std::map<unsigned int,bool> connected;
						unsigned int nb = connectedTo(grsimp,*itno,connected);
						
						if(!connected[refkey] && nb == 1)
						{
							for(std::map<unsigned int,bool>::iterator itm = connected.begin(); itm != connected.end(); itm++)
							{
								if(itm->second)
									grsimp->remNode(itm->first);
							}
							removed = true;
						}
						else
						{
							grsimp->addEdge(refkey,*itno);
						}
					}
				}
				if(!removed)
					it++;
				else
					fini = false;
			}
			catch(...)
			{
				it++;
			}
		}
	}while(!fini);

	return grsimp;
}
