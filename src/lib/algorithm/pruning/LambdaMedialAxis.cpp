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
 *  \file LambdaMedialAxis.cpp
 *  \brief Defines scale axis transform, a pruning algorithm
 *  \author Bastien Durix
 */

#include "LambdaMedialAxis.h"

std::list<unsigned int> closInd(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C)
{
	std::list<unsigned int> inds;
	float distmin = -1.0;
	for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
	{
		Eigen::Vector2d pt = disbnd->getCoordinates(i);
		float dist = (pt - C).norm();
		
		if(dist < distmin || distmin == -1.0)
		{
			distmin = dist;
		}
	}
	
	for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
	{
		Eigen::Vector2d pt = disbnd->getCoordinates(i);
		float dist = (pt - C).norm();

		if(dist == distmin)// || std::abs(dist-distmin) <= 10.0*std::abs(std::min(dist,distmin))*std::numeric_limits<float>::epsilon())
		{
			inds.push_back(i);
		}
	}
	
	return inds;
}

double computeLambdaCir(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C, double rad)
{
	std::list<unsigned int> inds = closInd(disbnd,C);
	std::vector<Eigen::Vector2d> pts(0);
	pts.reserve(inds.size());
	
	for(std::list<unsigned int>::iterator it = inds.begin(); it != inds.end(); it++)
		pts.push_back(disbnd->getCoordinates(*it));
	
	// 1: find the most distant points
	double distmax = 0.0;
	std::pair<unsigned int,unsigned int> pind;
	for(unsigned int i = 0; i < pts.size(); i++)
		for(unsigned int j = i+1; j < pts.size(); j++)
		{
			double dist = (pts[i] - pts[j]).norm();
			if(dist > distmax)
			{
				distmax = dist;
				pind.first = i;
				pind.second = j;
			}
		}

	Eigen::Vector2d ctr = (pts[pind.first] + pts[pind.second])*0.5;
	double lambdarad = distmax/2.0;
	
	bool valid = true;
	for(unsigned int i = 0; i < pts.size() && valid; i++)
	{
		if(i != pind.first && i != pind.second)
		{
			double dist = (ctr - pts[i]).norm();
			valid = (dist < lambdarad);
		}
	}
	
	if(!valid)
		lambdarad = rad;
	return lambdarad;
}

skeleton::GraphSkel2d::Ptr algorithm::pruning::LambdaMedialAxis(const skeleton::GraphSkel2d::Ptr grskel, const boundary::DiscreteBoundary<2>::Ptr disbnd, const double &lambda)
{
	skeleton::GraphSkel2d::Ptr grsimp(new skeleton::GraphSkel2d(*grskel));

	std::list<unsigned int> list_ind;
	grsimp->getAllNodes(list_ind);

	std::list<std::pair<unsigned int,double> > list_ind_size;
	for(std::list<unsigned int>::iterator it = list_ind.begin(); it != list_ind.end(); it++)
	{
		Eigen::Vector3d vec = grskel->getNode(*it);
		double lambdarad = computeLambdaCir(disbnd,vec.block<2,1>(0,0),vec(2));
		list_ind_size.push_back(std::pair<unsigned int,double>(*it,lambdarad));
	}

	//sort nodes in decreasing order
	list_ind_size.sort([](const std::pair<unsigned int,double> &p1, const std::pair<unsigned int,double> &p2)
			{
				return p1.second < p2.second;
			});
	
	bool fini = true;
	do{
		fini = true;
		std::list<std::pair<unsigned int,double> >::iterator it = list_ind_size.begin();
		while(it != list_ind_size.end())
		{
			try
			{
				unsigned int refkey = it->first;
				double lambdarad = it->second;
				bool removed = false;
				if(lambdarad < lambda)
				{
					unsigned int deg = grsimp->getNodeDegree(refkey);
					if(deg == 1)
					{
						grsimp->remNode(refkey);
						removed = true;
					}
				}
				
				if(!removed)
					it++;
				else
				{
					it = list_ind_size.erase(it);
					fini = false;
				}
			}
			catch(...)
			{
				it++;
			}
		}
	}while(!fini);

	return grsimp;
}
