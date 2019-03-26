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
 *  \file DistanceField.cpp
 *  \brief Defines functions related to distance field
 *  \author Bastien Durix
 */

#include "DistanceField.h"
#include <boost/math/distributions/normal.hpp>
#include <set>
#include <nlopt.hpp>
#include <iostream>

unsigned int algorithm::skeletonization::closestInd(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C)
{
	double distmin = 0.0;
	unsigned int ind = 0;
	for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
	{
		Eigen::Vector2d vec = disbnd->getCoordinates(i) - C;
		double dist = vec.norm();
		
		if(dist < distmin || i == 0)
		{
			distmin = dist;
			ind = i;
		}
	}

	return ind;
}

double algorithm::skeletonization::fieldValue(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C, double noise)
{
	boost::math::normal norm(0.0,noise);
	
	std::vector<double> dists(disbnd->getNbVertices());
	
	double distmin = 0.0;
	for(unsigned int i = 0; i < dists.size(); i++)
	{
		Eigen::Vector2d vec = disbnd->getCoordinates(i) - C;
		double dist = vec.norm();
		
		if(dist < distmin || i == 0) distmin = dist;
		
		dists[i] = dist;
	}
	
	double denom = 0.0;
	double numer = 0.0;
	
	for(unsigned int i = 0; i < dists.size(); i++)
	{
		double diff = dists[i] - distmin;
		double proba = 2.0*(1.0-boost::math::cdf(norm,diff));
		
		numer += proba * dists[i];
		
		denom += proba;
	}
	
	double val = numer/denom;
	
	if(val - distmin > 2.0)
	std::cout << val - distmin << std::endl;

	return val;
}

double algorithm::skeletonization::fieldValue(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C, Eigen::Vector2d &grad, double noise)
{
	boost::math::normal norm(0.0,noise);
	
	std::vector<double> dists(disbnd->getNbVertices());
	
	double distmin = 0.0;
	for(unsigned int i = 0; i < dists.size(); i++)
	{
		Eigen::Vector2d vec = disbnd->getCoordinates(i) - C;
		double dist = vec.norm();
		
		if(dist < distmin || i == 0) distmin = dist;
		
		dists[i] = dist;
	}

	double denom = 0.0;
	double numer = 0.0;
	Eigen::Vector2d numvec(0.0,0.0);
	
	for(unsigned int i = 0; i < dists.size(); i++)
	{
		Eigen::Vector2d vec = disbnd->getCoordinates(i) - C;
		double diff = dists[i] - distmin;
		double proba = 2.0*(1.0-boost::math::cdf(norm,diff));
		
		numer += proba * dists[i];
		
		numvec += proba*vec.normalized();
		
		denom += proba;
	}
	
	double val = numer/denom;

	if(val - distmin > 2.0)
	std::cout << val - distmin << std::endl;
	
	grad = numvec * (1.0/denom);
	
	return val;
}

void algorithm::skeletonization::tangencyBoundary(const boundary::DiscreteBoundary<2>::Ptr disbnd,
												  const Eigen::Vector2d &C,
												  const double &rad,
												  double noise,
												  std::vector<std::list<unsigned int> > &v_ind,
												  double fac)
{
	std::vector<bool> used(disbnd->getNbVertices(),false);
	std::vector<double> val(disbnd->getNbVertices());
	std::vector<bool> isin(disbnd->getNbVertices(),false);
	for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
	{
		val[i] = (disbnd->getVertex(i) - C).norm() - rad;
		if(val[i] < noise) isin[i] = true;
	}
	
	bool nochange = false;
	
	do
	{
		nochange = true;
		for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
		{
			if(!isin[i])
			{
				unsigned int prev = disbnd->getPrev(i);
				unsigned int next = disbnd->getNext(i);
				if((isin[prev] || isin[next]) && val[i] < fac)
				{
					isin[i] = true;
					nochange = false;
				}
			}
		}
	}while(!nochange);
	
	std::list<std::list<unsigned int> > l_ind;
	for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
	{
		if(isin[i] && !used[i])
		{
			std::list<unsigned int> cind;
			
			used[i] = true;
			unsigned int indn = i;
			do
			{
				cind.push_back(indn);
				used[indn] = true;
				indn = disbnd->getNext(indn);
			}while(isin[indn] && indn != *(cind.begin()));

			unsigned int indp = i;
			cind.pop_front();
			do
			{
				cind.push_front(indp);
				used[indp] = true;
				indp = disbnd->getPrev(indp);
			}while(isin[indp] && indp != *(cind.rbegin()));
			if(indp == *(cind.rbegin()))
				cind.push_front(indp);
			
			while(val[*(cind.begin())] > noise)
				cind.pop_front();
			while(val[*(cind.rbegin())] > noise)
				cind.pop_back();
			if(cind.size() > 1)
				l_ind.push_back(cind);
		}
	}
	
	l_ind.sort(
		[disbnd,C](const std::list<unsigned int> &c1, const std::list<unsigned int> &c2)
		{
			Eigen::Vector2d p1 = disbnd->getVertex(*(c1.begin())) - C;
			Eigen::Vector2d p2 = disbnd->getVertex(*(c2.begin())) - C;
			double ang1 = atan2(p1.y(),p1.x());
			double ang2 = atan2(p2.y(),p2.x());
			return ang1 < ang2;
		});
	v_ind = std::vector<std::list<unsigned int> >(l_ind.begin(),l_ind.end());
}

struct DataClosestCenter2d
{
	const boundary::DiscreteBoundary<2>::Ptr &disbnd;
	double noise;
	
	DataClosestCenter2d(
		const boundary::DiscreteBoundary<2>::Ptr &disbnd_,
		double noise_) :
		disbnd(disbnd_), noise(noise_) {};
};

double minClosestCenter2d(const std::vector<double> &vt, std::vector<double> &, void *dataFun)
{
	DataClosestCenter2d *data = (DataClosestCenter2d*) dataFun;
	
	Eigen::Vector2d C(vt[0],vt[1]);
	
	Eigen::Vector2d grad;
	algorithm::skeletonization::fieldValue(data->disbnd, C, grad, data->noise);
	//return -algorithm::skeletonization::recFieldValue(data->vdisbnd, data->vdisshp, data->vcam, C, data->noise);
	return grad.squaredNorm();
}

bool algorithm::skeletonization::closestCenter(const boundary::DiscreteBoundary<2>::Ptr disbnd, Eigen::Vector2d &C, double noise)
{
	unsigned int nbiter = 10000;
	
	bool fini = false;
	do
	{
		Eigen::Vector2d grad;
		fieldValue(disbnd,C,grad,noise);
		fini = grad.norm() < 0.01*noise;
		if(!fini) C -= grad * noise * 0.1;
		nbiter--;
	}while(nbiter != 0 && !fini);
	
	std::vector<double> lb(2);
	std::vector<double> ub(2);
	std::vector<double> C_init(2);

	Eigen::Vector2d grad;
	double rad = fieldValue(disbnd,C,grad,noise);
	
	C_init[0] = C.x();
	C_init[1] = C.y();
	
	lb[0] = C.x() - rad*0.2;
	lb[1] = C.y() - rad*0.2;
	
	ub[0] = C.x() + rad*0.2;
	ub[1] = C.y() + rad*0.2;
	
	DataClosestCenter2d data(disbnd,noise);

	nlopt::opt opt(nlopt::LN_COBYLA, 2);
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
	
	opt.set_min_objective(minClosestCenter2d, &data);
	
	opt.set_xtol_rel(1e-2);
	
	double res;
	opt.optimize(C_init, res);
	C.x() = C_init[0];
	C.y() = C_init[1];

	return true;
}

void algorithm::skeletonization::closestCenterOnArc(const boundary::DiscreteBoundary<2>::Ptr disbnd,
													const Eigen::Vector2d &C,
													const std::pair<double,double> &pang,
													double dist,
													Eigen::Vector2d &Cmov,
													double noise)
{
	double step = noise / (dist*dist);
	
	double ang1 = pang.first;
	double ang2 = pang.second;

	Eigen::Vector2d P1 = C + dist*Eigen::Vector2d(cos(ang1),sin(ang1));
	Eigen::Vector2d P2 = C + dist*Eigen::Vector2d(cos(ang2),sin(ang2));
	
	while(ang2 - ang1 > step)
	{
		double angmid = (ang1 + ang2) / 2.0;
		
		Eigen::Vector2d grad;
		Eigen::Vector2d vec(cos(angmid),sin(angmid));
		Cmov = C + dist * vec;
		fieldValue(disbnd,Cmov,grad,noise);
		
		//double der = dist*(grad.y() * vec.x() - grad.x() * vec.y());
		//if(der < 0.0) ang1 = angmid;
		//if(der > 0.0) ang2 = angmid;
		
		Eigen::Vector2d vec1 = (P1 - Cmov).normalized();
		Eigen::Vector2d vec2 = (P2 - Cmov).normalized();
		double sim1 = grad.dot(vec1);
		double sim2 = grad.dot(vec2);
		if(sim1 > sim2) ang1 = angmid;
		if(sim2 > sim1) ang2 = angmid;
	}

	double angmid = (ang1 + ang2) / 2.0;

	Eigen::Vector2d vec(cos(angmid),sin(angmid));
	Cmov = C + dist * vec;
}

void algorithm::skeletonization::discontinuitiesOnArc(const boundary::DiscreteBoundary<2>::Ptr disbnd,
													  const Eigen::Vector2d &C,
													  const std::pair<double,double> &pang,
													  double dist,
													  std::vector<Eigen::Vector2d> &vecC,
													  std::vector<std::pair<unsigned int,unsigned int> > &vecInd,
													  double noise)
{
	std::list<Eigen::Vector2d> lisC;
	std::list<std::pair<unsigned int,unsigned int> > lisInd;
	double step = noise / (dist*dist);
	
	double ang2 = pang.second;
	Eigen::Vector2d P2 = C + dist*Eigen::Vector2d(cos(ang2),sin(ang2));
	unsigned int ind2 = closestInd(disbnd,P2);
	
	double ang1 = pang.first;
	Eigen::Vector2d P1 = C + dist*Eigen::Vector2d(cos(ang1),sin(ang1));
	unsigned int ind1 = closestInd(disbnd,P1);
	
	//std::cout << "disc " << ind1 << " ";

	while(ind1 != ind2)
	{
		while(ang2 - ang1 > step)
		{
			double angmid = (ang1 + ang2) / 2.0;
			Eigen::Vector2d vec(cos(angmid),sin(angmid));
			Eigen::Vector2d Cmov = C + dist * vec;
			
			unsigned int indmid = closestInd(disbnd,Cmov);

			if(indmid == ind1)
			{
				ang1 = angmid;
			}
			else
			{
				ang2 = angmid;
			}
		}
		
		Eigen::Vector2d P1 = C + dist*Eigen::Vector2d(cos(ang1),sin(ang1));
		unsigned int ind1b = closestInd(disbnd,P1);
		Eigen::Vector2d P2 = C + dist*Eigen::Vector2d(cos(ang2),sin(ang2));
		unsigned int ind2b = closestInd(disbnd,P2);

		double angmid = (ang1 + ang2) / 2.0;
		Eigen::Vector2d vec(cos(angmid),sin(angmid));
		Eigen::Vector2d Cmov = C + dist * vec;
		
		lisC.push_back(Cmov);
		lisInd.push_back(std::make_pair(ind1b,ind2b));
		
		ind1 = ind2b;
		ang1 = ang2;
		ang2 = pang.second;
		//std::cout << ind1 << " ";
	}
	//std::cout << std::endl;
	
	vecC = std::vector<Eigen::Vector2d>(lisC.begin(),lisC.end());
	vecInd = std::vector<std::pair<unsigned int,unsigned int> >(lisInd.begin(),lisInd.end());
}
