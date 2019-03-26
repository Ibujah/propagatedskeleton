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
 *  \file MovingCenter.cpp
 *  \brief Defines object containing informations about moving centers
 *  \author Bastien Durix
 */

#include "MovingCenter.h"
#include "DistanceField.h"
#include <set>
#include <iostream>

algorithm::skeletonization::MovingCenter::MovingCenter() {}

algorithm::skeletonization::MovingCenter::MovingCenter(const Eigen::Vector2d &center) : m_center(center) {}

void algorithm::skeletonization::MovingCenter::correction(const boundary::DiscreteBoundary<2>::Ptr disbnd)
{
	std::list<Eigen::Vector2d> lpt;
	for(unsigned int i = 0; i < m_tgtinds.size(); i++)
	{
		double distmin = -1.0;
		unsigned int indmin;
		Eigen::Vector2d pt;
		for(std::list<unsigned int>::const_iterator it = m_tgtinds[i].begin(); it != m_tgtinds[i].end(); it++)
		{
			Eigen::Vector2d p1 = disbnd->getCoordinates(*it);
			double dist = (p1-m_center).norm();
			//std::cout << "(" << *it << "," << dist << ") ";
			if(distmin == -1.0 || dist < distmin)
			{
				distmin = dist;
				pt = p1;
				indmin = *it;
			}
		}
		lpt.push_back(pt);
	}

	if(lpt.size() == 2)
	{
		Eigen::Vector2d P1 = *(lpt.begin());
		Eigen::Vector2d P2 = *(lpt.rbegin());
		double lambda = 0.5*((P2-P1).norm()) - (m_center - P1).dot((P2 - P1).normalized());
		m_center = m_center + lambda*((P2-P1).normalized());
	}
	else if(lpt.size() >= 3)
	{
		for(unsigned int i = 0; i < 10; i++)
		{
			for(std::list<Eigen::Vector2d>::iterator it = lpt.begin(); it != lpt.end(); it++)
			{
				for(std::list<Eigen::Vector2d>::iterator it2 = std::next(it); it2 != lpt.end(); it2++)
				{
					Eigen::Vector2d P1 = *(it);
					Eigen::Vector2d P2 = *(it2);
					double lambda = 0.5*((P2-P1).norm()) - (m_center - P1).dot((P2 - P1).normalized());
					m_center = m_center + lambda*((P2-P1).normalized());
				}
			}
		}
	}
}

void algorithm::skeletonization::MovingCenter::computeTangencyData(const boundary::DiscreteBoundary<2>::Ptr disbnd, double noise, double alpha)
{
	m_rad = fieldValue(disbnd,m_center,noise);
	tangencyBoundary(disbnd,m_center,m_rad,noise,m_tgtinds,alpha);
	m_ang.resize(m_tgtinds.size());
	for(unsigned int i = 0; i < m_tgtinds.size(); i++)
	{
		/*Eigen::Vector2d pt;
		double distmin = -1.;
		for(std::list<unsigned int>::const_iterator it = m_tgtinds[i].begin(); it != m_tgtinds[i].end(); it++)
		{
			Eigen::Vector2d p1 = disbnd->getCoordinates(*it);
			double dist = (p1-m_center).norm();
			if(distmin == -1.0 || dist < distmin)
			{
				distmin = dist;
				pt = p1;
			}
		}
		double ang = atan2(pt.y() - m_center.y(), pt.x() - m_center.x());
		m_ang[i].first  = ang;
		m_ang[i].second = ang;*/
		
		Eigen::Vector2d p1 = disbnd->getCoordinates(*(m_tgtinds[i].begin()));
		Eigen::Vector2d p2 = disbnd->getCoordinates(*(m_tgtinds[i].rbegin()));
		double ang1 = atan2(p1.y() - m_center.y(), p1.x() - m_center.x());
		double ang2 = atan2(p2.y() - m_center.y(), p2.x() - m_center.x());
		if(ang2 < ang1) ang2 += 2.0*M_PI;
		m_ang[i].first  = ang1;
		m_ang[i].second = ang2;
	}
}

void algorithm::skeletonization::MovingCenter::computeTangencyData(const MovingCenter &mov, unsigned int dir, const boundary::DiscreteBoundary<2>::Ptr disbnd, double noise, double alpha)
{
	m_rad = fieldValue(disbnd,m_center,noise);
	tangencyBoundary(disbnd,m_center,m_rad,noise,m_tgtinds,alpha);
	
	double ang1 = mov.getAng()[(dir+1)%mov.getAng().size()].second;
	double ang2 = mov.getAng()[dir].first;
	while(ang2 < ang1) ang2 += 2.0 * M_PI;
	
	for(std::vector<std::list<unsigned int> >::iterator it = m_tgtinds.begin(); it != m_tgtinds.end();)
	{
		bool isin = true;
		for(std::list<unsigned int>::iterator itp = it->begin(); itp != it->end() && isin; itp++)
		{
			Eigen::Vector2d vec = disbnd->getCoordinates(*itp) - mov.getCenter();
			double angpt = atan2(vec.y(),vec.x());
			while(angpt < ang1) angpt += 2.0 * M_PI;
			if(angpt > ang2)
				isin = false;
		}
		
		if(isin)
		{
			it = m_tgtinds.erase(it);
		}
		else
		{
			it++;
		}
	}

	m_ang.resize(m_tgtinds.size());
	for(unsigned int i = 0; i < m_tgtinds.size(); i++)
	{
		/*Eigen::Vector2d pt;
		double distmin = -1.;
		for(std::list<unsigned int>::const_iterator it = m_tgtinds[i].begin(); it != m_tgtinds[i].end(); it++)
		{
			Eigen::Vector2d p1 = disbnd->getCoordinates(*it);
			double dist = (p1-m_center).norm();
			if(distmin == -1.0 || dist < distmin)
			{
				distmin = dist;
				pt = p1;
			}
		}
		double ang = atan2(pt.y() - m_center.y(), pt.x() - m_center.x());
		m_ang[i].first  = ang;
		m_ang[i].second = ang;*/
		
		Eigen::Vector2d p1 = disbnd->getCoordinates(*(m_tgtinds[i].begin()));
		Eigen::Vector2d p2 = disbnd->getCoordinates(*(m_tgtinds[i].rbegin()));
		double ang1 = atan2(p1.y() - m_center.y(), p1.x() - m_center.x());
		double ang2 = atan2(p2.y() - m_center.y(), p2.x() - m_center.x());
		if(ang2 < ang1) ang2 += 2.0*M_PI;
		m_ang[i].first  = ang1;
		m_ang[i].second = ang2;
	}
}

void algorithm::skeletonization::MovingCenter::computeTangencyDataaccurate(const MovingCenter &mov, unsigned int dir, const boundary::DiscreteBoundary<2>::Ptr disbnd, double noise, double alpha)
{
	m_rad = fieldValue(disbnd,m_center,noise);
	tangencyBoundary(disbnd,m_center,m_rad,noise,m_tgtinds,alpha);
	
	bool fini = false;
	unsigned int dirmov;
	unsigned int offset;
	for(unsigned int j = 0; j < m_tgtinds.size() - 1 && !fini; j++)
	{
		for(unsigned int i = 0; i < m_tgtinds.size() && !fini; i++)
		{
			fini = intersect(*this,i,mov,dir,j);
			if(fini)
			{
				dirmov = i;
				offset = j;
			}
		}
	}

	if(fini) // supprimer les tangences strictement entre dirmov et dirmov+offset
	{
		std::list<std::list<unsigned int> > ltgtinds;
		unsigned int dirend = (dirmov + offset)%m_tgtinds.size();
		for(unsigned int i = 0; i < m_tgtinds.size(); i++)
		{
			if(dirmov < dirend)
			{
				if(i <= dirmov || i >= dirend)
					ltgtinds.push_back(m_tgtinds[i]);
			}
			else
			{
				if(i >= dirmov || i <= dirend)
					ltgtinds.push_back(m_tgtinds[i]);
			}
		}
		m_tgtinds = std::vector<std::list<unsigned int> >(ltgtinds.begin(), ltgtinds.end());
	}

	m_ang.resize(m_tgtinds.size());
	for(unsigned int i = 0; i < m_tgtinds.size(); i++)
	{
		/*Eigen::Vector2d pt;
		double distmin = -1.;
		for(std::list<unsigned int>::const_iterator it = m_tgtinds[i].begin(); it != m_tgtinds[i].end(); it++)
		{
			Eigen::Vector2d p1 = disbnd->getCoordinates(*it);
			double dist = (p1-m_center).norm();
			if(distmin == -1.0 || dist < distmin)
			{
				distmin = dist;
				pt = p1;
			}
		}
		double ang = atan2(pt.y() - m_center.y(), pt.x() - m_center.x());
		m_ang[i].first  = ang;
		m_ang[i].second = ang;*/
		
		Eigen::Vector2d p1 = disbnd->getCoordinates(*(m_tgtinds[i].begin()));
		Eigen::Vector2d p2 = disbnd->getCoordinates(*(m_tgtinds[i].rbegin()));
		double ang1 = atan2(p1.y() - m_center.y(), p1.x() - m_center.x());
		double ang2 = atan2(p2.y() - m_center.y(), p2.x() - m_center.x());
		if(ang2 < ang1) ang2 += 2.0*M_PI;
		m_ang[i].first  = ang1;
		m_ang[i].second = ang2;
	}
}

bool algorithm::skeletonization::MovingCenter::propagate(const boundary::DiscreteBoundary<2>::Ptr disbnd,
														 unsigned int dir,
														 double noise,
														 algorithm::skeletonization::MovingCenter &mov,
														 std::list<unsigned int> &lnext,
														 double alpha) const
{
	std::pair<double,double> pang;
	pang.first  = m_ang[dir].second;
	pang.second = m_ang[(dir+1)%m_ang.size()].first;
	while(pang.second < pang.first) pang.second += 2.0*M_PI;
	
	double dist0 = 0.0, dist1 = m_rad;

	Eigen::Vector2d Cmov0;
	unsigned int nbdir = 0;
	while(dist1 - dist0 > noise)
	{
		Eigen::Vector2d Cmov;
		double dist = (dist1 + dist0)/2.0;
		bool fini = false;
		closestCenterOnArc(disbnd,m_center,pang,dist,Cmov,noise);
		MovingCenter mov(Cmov);
		mov.computeTangencyData(*this,dir,disbnd,noise,alpha);
		for(unsigned int i = 0; i < mov.m_tgtinds.size() && !fini; i++)
		{
			fini = intersect(*this,dir,mov,i);
		}
		if(fini)
		{
			Cmov0 = Cmov;
			dist0 = dist;
			nbdir = mov.getNbdir();
		}
		else
			dist1 = dist;
	}
	
	double dist = dist0;
	if(nbdir != 2)
	{
		double ddist0 = 0.0, ddist1 = dist0;
		Eigen::Vector2d grad0;
		fieldValue(disbnd,m_center,grad0,noise);
		Eigen::Vector2d grad1;
		fieldValue(disbnd,Cmov0,grad1,noise);
		while(ddist1 - ddist0 > noise)
		{
			Eigen::Vector2d Cmov;
			double ddist = (ddist1 + ddist0)/2.0;
			closestCenterOnArc(disbnd,m_center,pang,ddist,Cmov,noise);
			MovingCenter mov(Cmov);
			mov.computeTangencyData(*this,dir,disbnd,noise,alpha);
			Eigen::Vector2d grad;
			fieldValue(disbnd,Cmov,grad,noise);
			if(mov.getNbdir() != nbdir)
			{
				//grad0 = grad;
				ddist0 = ddist;
				Cmov0 = Cmov;
			}
			else
			{
				//if(grad.dot(grad0) > grad.dot(grad1))
				if((grad - grad0).squaredNorm() < (grad - grad1).squaredNorm())
				{
					//grad0 = grad;
					ddist0 = ddist;
					Cmov0 = Cmov;
				}
				else
				{
					//grad1 = grad;
					ddist1 = ddist;
				}
			}
		}
		dist = ddist1;
	}
	bool fini = false;
	unsigned int dirmov;
	mov = MovingCenter(Cmov0);
	mov.computeTangencyData(*this,dir,disbnd,noise,alpha);
	for(unsigned int i = 0; i < mov.m_tgtinds.size() && !fini; i++)
	{
		fini = intersect(*this,dir,mov,i);
		if(fini) dirmov = i;
	}
	for(unsigned int i = 0; i < mov.m_tgtinds.size(); i++)
	{
		if(i != dirmov) lnext.push_back(i);
	}
	
	/*if(!fini || dist == 0.0)
	{
		fini = propagateaccurate(disbnd,dir,noise,mov,lnext,alpha);
	}*/

	//mov.correction(disbnd);

	return fini;
}


bool algorithm::skeletonization::MovingCenter::propagateaccurate(const boundary::DiscreteBoundary<2>::Ptr disbnd,
														 unsigned int dir,
														 double noise,
														 algorithm::skeletonization::MovingCenter &mov,
														 std::list<unsigned int> &lnext,
														 double alpha) const
{
	//std::cout << "accurate" << std::endl;
	//std::cout << m_center.transpose() << std::endl;
	//std::cout << "512 " << (m_center - disbnd->getCoordinates(512)).norm() << std::endl;
	std::pair<double,double> pang;
	Eigen::Vector2d pt;
	double distmin = -1.0;
	unsigned int indmin1, indmin2;
	for(std::list<unsigned int>::const_iterator it = m_tgtinds[dir].begin(); it != m_tgtinds[dir].end(); it++)
	{
		Eigen::Vector2d p1 = disbnd->getCoordinates(*it);
		double dist = (p1-m_center).norm();
		//std::cout << "(" << *it << "," << dist << ") ";
		if(distmin == -1.0 || dist < distmin)
		{
			distmin = dist;
			pt = p1;
			indmin1 = *it;
		}
	}
	//std::cout << std::endl;
	//std::cout << indmin1 << std::endl;
	double ang = atan2(pt.y() - m_center.y(), pt.x() - m_center.x());
	pang.first  = ang;
	
	distmin = -1.0;
	for(std::list<unsigned int>::const_iterator it = m_tgtinds[(dir+1)%m_tgtinds.size()].begin(); it != m_tgtinds[(dir+1)%m_tgtinds.size()].end(); it++)
	{
		Eigen::Vector2d p1 = disbnd->getCoordinates(*it);
		double dist = (p1-m_center).norm();
		//std::cout << "(" << *it << "," << dist << ") ";
		if(distmin == -1.0 || dist < distmin)
		{
			distmin = dist;
			pt = p1;
			indmin2 = *it;
		}
	}
	//std::cout << std::endl;
	//std::cout << indmin2 << std::endl;
	ang = atan2(pt.y() - m_center.y(), pt.x() - m_center.x());
	pang.second = ang;
	while(pang.second < pang.first) pang.second += 2.0*M_PI;
	
	double dist0 = 0.0, dist1 = m_rad;

	Eigen::Vector2d Cmov0;
	while(dist1 - dist0 > noise)
	{
		//std::cout << dist0 << " " << dist1 << std::endl;
		Eigen::Vector2d Cmov;
		double dist = (dist1 + dist0)/2.0;
		bool fini = false;
		std::vector<Eigen::Vector2d> vecC;
		std::vector<std::pair<unsigned int,unsigned int> > vecInd;
		discontinuitiesOnArc(disbnd,m_center,pang,dist,vecC,vecInd,noise);
		
		std::list<Eigen::Vector2d> lisfC;
		std::list<std::pair<unsigned int,unsigned int> > lisfInd;
		
		//std::cout << "vecInd " << vecInd.size() << std::endl;
		if(vecInd.size() != 0 && vecInd[0].first == indmin1 && vecInd[vecInd.size()-1].second == indmin2)
		{
			// filter
			for(unsigned int i = 0; i < vecC.size(); i++)
			{
				MovingCenter mov(vecC[i]);
				mov.computeTangencyDataaccurate(*this,dir,disbnd,noise,alpha);
				
				bool filter = false;
				for(unsigned int j = 0; j < mov.m_tgtinds.size() && !filter; j++)
				{
					if(std::find(mov.m_tgtinds[j].begin(),mov.m_tgtinds[j].end(),vecInd[i].first) != mov.m_tgtinds[j].end() &&
					   std::find(mov.m_tgtinds[j].begin(),mov.m_tgtinds[j].end(),vecInd[i].second) != mov.m_tgtinds[j].end())
					{
						filter = true;
					}
				}
				if(!filter)
				{
					lisfC.push_back(vecC[i]);
					lisfInd.push_back(vecInd[i]);
				}
			}
			//std::cout << "lisfC " << lisfC.size() << std::endl;
			
			if(lisfC.size() > 1)
			{
				dist1 = dist;
			}
			else if(lisfC.size() == 0)
			{
				//new filter
				lisfC.clear();
				lisfInd.clear();
				for(unsigned int i = 0; i < vecC.size(); i++)
				{
					MovingCenter mov(vecC[i]);
					mov.computeTangencyDataaccurate(*this,dir,disbnd,noise,alpha);
					for(unsigned int i = 0; i < mov.m_tgtinds.size() && !fini; i++)
					{
						fini = intersect(*this,dir,mov,i);
					}
					if(fini)
					{
						lisfC.push_back(vecC[i]);
						lisfInd.push_back(vecInd[i]);
					}
				}
				
				if(lisfC.size() != 0)
				{
					Cmov0 = *(lisfC.begin());
					dist0 = dist;
				}
				else
					dist1 = dist;
			}
			else
			{
				Cmov = *(lisfC.begin());
				MovingCenter mov(Cmov);
				mov.computeTangencyDataaccurate(*this,dir,disbnd,noise,alpha);
	

				
				//for(unsigned int i = 0; i < mov.m_tgtinds.size(); i++)
				//{
				//	std::cout << "inds ";
				//	for(std::list<unsigned int>::const_iterator it = mov.m_tgtinds[i].begin(); it != mov.m_tgtinds[i].end(); it++)
				//	{
				//		std::cout << *it << " ";
				//	}
				//	std::cout << std::endl;
				//}



				for(unsigned int i = 0; i < mov.m_tgtinds.size() && !fini; i++)
				{
					fini = intersect(*this,dir,mov,i);
				}
				if(fini)
				{
					Cmov0 = Cmov;
					dist0 = dist;
				}
				else
					dist1 = dist;
			}
		}
		else
		{
			double ang = (pang.first + pang.second)/2.0;
			Cmov0 = m_center + dist*Eigen::Vector2d(cos(ang),sin(ang));
			dist0 = dist;
		}
	}

	bool fini = false;
	unsigned int dirmov;
	mov = MovingCenter(Cmov0);
	mov.computeTangencyDataaccurate(*this,dir,disbnd,noise,alpha);
	for(unsigned int i = 0; i < mov.m_tgtinds.size() && !fini; i++)
	{
		fini = intersect(*this,dir,mov,i);
		if(fini) dirmov = i;
	}
	for(unsigned int i = 0; i < mov.m_tgtinds.size(); i++)
	{
		if(i != dirmov) lnext.push_back(i);
	}

	//std::cout << mov.m_center.transpose() << std::endl;
	if(dist0 == 0.0 || !fini)
	{
		//std::cout << "arf" << std::endl;
		throw std::logic_error("arf...");
	}
	mov.correction(disbnd);
	
	return fini && dist0 != 0.0;
}

unsigned int algorithm::skeletonization::MovingCenter::getNbdir() const
{
	if(m_tgtinds.size() == 1) if(*(m_tgtinds[0].begin()) == *(m_tgtinds[0].rbegin())) return 0;
	return m_tgtinds.size();
}

const Eigen::Vector2d& algorithm::skeletonization::MovingCenter::getCenter() const
{
	return m_center;
}

double algorithm::skeletonization::MovingCenter::getRadius() const
{
	return m_rad;
}

const std::vector<std::list<unsigned int> >& algorithm::skeletonization::MovingCenter::getTgt() const
{
	return m_tgtinds;
}

const std::vector<std::pair<double,double> >& algorithm::skeletonization::MovingCenter::getAng() const
{
	return m_ang;
}

bool algorithm::skeletonization::MovingCenter::isInNext(unsigned int dir, const Eigen::Vector2d &C) const
{
	double ang1 = m_ang[dir].second;
	double ang2 = m_ang[(dir+1)%m_ang.size()].first;
	while(ang2 < ang1) ang2 += 2.0 * M_PI;
	
	Eigen::Vector2d vec = C - m_center;
	
	double angpt = atan2(vec.y(),vec.x());
	while(angpt < ang1) angpt += 2.0 * M_PI;
	return angpt < ang2;
}

bool algorithm::skeletonization::MovingCenter::intersect(const algorithm::skeletonization::MovingCenter &mov1, unsigned int dir1, const algorithm::skeletonization::MovingCenter &mov2, unsigned int dir2, unsigned int offset1)
{
	std::set<unsigned int> set11(mov1.m_tgtinds[dir1].begin(),
								 mov1.m_tgtinds[dir1].end()),
						   set12(mov1.m_tgtinds[(dir1+offset1)%mov1.m_tgtinds.size()].begin(),
								 mov1.m_tgtinds[(dir1+offset1)%mov1.m_tgtinds.size()].end()),
						   set21(mov2.m_tgtinds[dir2].begin(),
								 mov2.m_tgtinds[dir2].end()),
						   set22(mov2.m_tgtinds[(dir2+1)%mov2.m_tgtinds.size()].begin(),
								 mov2.m_tgtinds[(dir2+1)%mov2.m_tgtinds.size()].end());

	if(mov1.m_tgtinds.size() == 1)
	{
		std::list<unsigned int>::const_iterator it = mov1.m_tgtinds[0].begin();
		for(unsigned int i = 0; i < mov1.m_tgtinds[0].size()/2; i++)
			it++;
		set11.clear();
		set11.insert(it,mov1.m_tgtinds[0].end());
		set12.clear();
		set12.insert(mov1.m_tgtinds[0].begin(),it);
	}

	if(mov2.m_tgtinds.size() == 1)
	{
		std::list<unsigned int>::const_iterator it = mov2.m_tgtinds[0].begin();
		for(unsigned int i = 0; i < mov2.m_tgtinds[0].size()/2; i++)
			it++;
		set21.clear();
		set21.insert(it,mov2.m_tgtinds[0].end());
		set22.clear();
		set22.insert(mov2.m_tgtinds[0].begin(),it);
	}
	
	std::set<unsigned int> set1, set2;
	std::set_intersection(set11.begin(),set11.end(),set22.begin(),set22.end(),std::inserter(set1,set1.begin()));
	std::set_intersection(set12.begin(),set12.end(),set21.begin(),set21.end(),std::inserter(set2,set2.begin()));
	
	return set1.size() && set2.size();
}
