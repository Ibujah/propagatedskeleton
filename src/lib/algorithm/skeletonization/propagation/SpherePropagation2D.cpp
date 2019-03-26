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
 *  \file SpherePropagation2D.cpp
 *  \brief Defines functions to compute 2d skeleton with sphere propagation algorithm
 *  \author Bastien Durix
 */

#include "SpherePropagation2D.h"
#include "DistanceField.h"
#include "MovingCenter.h"
#include <time.h>

#include <iostream>
#include <fstream>

/*void infosLocales(const boundary::DiscreteBoundary<2>::Ptr disbnd,
				  const algorithm::skeletonization::MovingCenter &mov,
				  double noise,
				  std::ofstream &ofs,
				  bool showgrad = false)
{
	for(unsigned int i= 0; i < mov.getTgt().size(); i++)
	{
		Eigen::Vector2d p1 = disbnd->getCoordinates(*(mov.getTgt()[i].begin()));
		Eigen::Vector2d p2 = disbnd->getCoordinates(*(mov.getTgt()[i].rbegin()));
		Eigen::Vector2d C = mov.getCenter();
		ofs << "plot([" << p1.x() + 0.5 << ", " << p2.x() + 0.5 << "], [" << p1.y() + 0.5 << ", " << p2.y() + 0.5 << "],'color','red');" << std::endl;
		ofs << "plot([" << C.x() + 0.5 << ", " << p1.x() + 0.5 << "], [" << C.y() + 0.5 << ", " << p1.y() + 0.5 << "],'color','red');" << std::endl;
		ofs << "plot([" << C.x() + 0.5 << ", " << p2.x() + 0.5 << "], [" << C.y() + 0.5 << ", " << p2.y() + 0.5 << "],'color',[1. 0. 1.]);" << std::endl;
	}
	
	Eigen::Vector2d grad;
	algorithm::skeletonization::fieldValue(disbnd,mov.getCenter(),grad,noise);
	ofs << "plot([" << mov.getCenter().x() + 0.5 << ", " << mov.getCenter().x() + 0.5 + 5.0*grad.x() << "], ["
					<< mov.getCenter().y() + 0.5 << ", " << mov.getCenter().y() + 0.5 + 5.0*grad.y() << "],'color','green');" << std::endl;
	
	if(showgrad)
	{
		unsigned int nbdiv = (unsigned int)(2.0*M_PI*mov.getRadius()/noise);
		for(unsigned int i = 0; i < nbdiv; i++)
		{
			double ang = (double)i*2.0*M_PI/(double)nbdiv;
			Eigen::Vector2d vec(cos(ang),sin(ang));
			Eigen::Vector2d C = mov.getCenter();
			Eigen::Vector2d Cmov = C + noise*vec;
			algorithm::skeletonization::fieldValue(disbnd,Cmov,noise);

			ofs << "plot([" << C.x() + 0.5 + 10.0*vec.x() << ", " << C.x() + 0.5 + 10.0*vec.x() + 10.0*grad.x() << "] , ["
				<< C.y() + 0.5 + 10.0*vec.y() << ", " << C.y() + 0.5 + 10.0*vec.y() + 10.0*grad.y() << "], 'blue');" << std::endl;
		}
	}
}


void makeMovie(const skeleton::GraphSkel2d::Ptr grskl,
			   const boundary::DiscreteBoundary<2>::Ptr disbnd,
			   std::ofstream &ofs,
			   const algorithm::skeletonization::OptionsSphProp &options)
{
	ofs << "fig = figure('units','normalized','outerposition',[0 0 1 1]);" << std::endl;
	ofs << "dist = dist - min(dist(:));" << std::endl;
	ofs << "dist(:) = uint8(63*dist(:)/max(dist(:)))+1;" << std::endl;
	ofs << "dist = ind2rgb(dist,jet);" << std::endl;
	ofs << "dist(img(:) < 125) = 1.0;" << std::endl;
	//ofs << "imagesc(img);" << std::endl;
	ofs << "imagesc(dist);" << std::endl;
	ofs << "hold on" << std::endl;
	ofs << "axis equal" << std::endl;
	ofs << "axis off" << std::endl;
	ofs << "hold off" << std::endl;
	//ofs << "imagesc(img);" << std::endl;
	ofs << "imagesc(dist);" << std::endl;
	ofs << "hold on" << std::endl;
	ofs << "axis equal" << std::endl;
	ofs << "axis off" << std::endl;
	ofs << "axis([ 0 size(img,2) 0 size(img,1) ]);" << std::endl;
	ofs << "ax = gca;" << std::endl;
	ofs << "mov = moviein(" << grskl->getNbNodes() + 10 << ",ax);" << std::endl;
	
	std::list<unsigned int> nods;
	grskl->getAllNodes(nods);
	std::vector<bool> valid(disbnd->getNbVertices(),false);

	std::list<std::pair<unsigned int,unsigned int> > ledg;
	grskl->getAllEdges(ledg);

	unsigned int cpt = 1;
	for(std::list<unsigned int>::iterator it = nods.begin(); it != nods.end(); it++)
	{
		ofs << "figure(fig);" << std::endl;
		ofs << "hold off" << std::endl;
		//ofs << "imagesc(img);" << std::endl;
		ofs << "imagesc(dist);" << std::endl;
		ofs << "hold on" << std::endl;
		ofs << "axis equal" << std::endl;
		ofs << "axis off" << std::endl;
		ofs << "axis([ 0 size(img,2) 0 size(img,1) ]);" << std::endl;
		ofs << "ax = gca;" << std::endl;
		
		Eigen::Vector3d cir = grskl->getNode(*it);
		algorithm::skeletonization::MovingCenter mov(cir.block<2,1>(0,0));
		mov.computeTangencyData(disbnd,options.noise,options.alpha);
		for(unsigned int i = 0; i < mov.getTgt().size(); i++)
		{
			for(auto itp = mov.getTgt()[i].begin(); itp != mov.getTgt()[i].end(); itp++)
			{
				valid[*itp] = true;
			}
		}

		for(unsigned int i = 0; i < valid.size(); i++)
		{
			if(valid[i])
				ofs << "plot(" << disbnd->getVertex(i).getCoords().x() + 0.5 << ", " << disbnd->getVertex(i).getCoords().y() + 0.5 << ", '.', 'color', 'red');" << std::endl;
		}
		
		for(auto ite = ledg.begin(); ite != ledg.end(); ite++)
		{
			if(ite->first <= *it && ite ->second <= *it)
			{
				Eigen::Vector3d pt1 = grskl->getNode(ite->first);
				Eigen::Vector3d pt2 = grskl->getNode(ite->second);
				
				ofs << "plot([" << pt1.x() + 0.5 << " " << pt2.x() + 0.5 << "], [" << pt1.y() + 0.5 << " " << pt2.y() + 0.5 << "], 'color', 'black', 'linewidth', 2);" << std::endl;
				
				if(grskl->getNodeDegree(ite->first) == 1)
					ofs << "draw_circle(" << pt1.y() + 0.5 << ", " << pt1.x() + 0.5 << ", " << pt1.z() << ",'black');" << std::endl;
				
				if(grskl->getNodeDegree(ite->second) == 1)
					ofs << "draw_circle(" << pt2.y() + 0.5 << ", " << pt2.x() + 0.5 << ", " << pt2.z() << ",'black');" << std::endl;
			}
			else if(ite->first <= *it)
			{
				Eigen::Vector3d cirit = grskl->getNode(ite->first);
				ofs << "draw_circle(" << cirit.y() + 0.5 << ", " << cirit.x() + 0.5 << ", " << cirit.z() << ",'black');" << std::endl;
			}
			else if(ite->second <= *it)
			{
				Eigen::Vector3d cirit = grskl->getNode(ite->second);
				ofs << "draw_circle(" << cirit.y() + 0.5 << ", " << cirit.x() + 0.5 << ", " << cirit.z() << ",'black');" << std::endl;
			}

		}
		
		ofs << "draw_circle(" << mov.getCenter().y() + 0.5 << ", " << mov.getCenter().x() + 0.5 << ", " << mov.getRadius() << ",'black');" << std::endl;
		
		ofs << "mov(:," << cpt << ") = getframe(ax);" << std::endl;
		cpt++;
	}

	ofs << "figure(fig);" << std::endl;
	ofs << "hold off" << std::endl;
	//ofs << "imagesc(img);" << std::endl;
	ofs << "imagesc(dist);" << std::endl;
	ofs << "hold on" << std::endl;
	ofs << "axis equal" << std::endl;
	ofs << "axis off" << std::endl;
	ofs << "axis([ 0 size(img,2) 0 size(img,1) ]);" << std::endl;
	ofs << "ax = gca;" << std::endl;

	for(unsigned int i = 0; i < valid.size(); i++)
	{
		if(valid[i])
			ofs << "plot(" << disbnd->getVertex(i).getCoords().x() + 0.5 << ", " << disbnd->getVertex(i).getCoords().y() + 0.5 << ", '.', 'color', 'red');" << std::endl;
	}

	for(auto ite = ledg.begin(); ite != ledg.end(); ite++)
	{
		Eigen::Vector3d pt1 = grskl->getNode(ite->first);
		Eigen::Vector3d pt2 = grskl->getNode(ite->second);

		ofs << "plot([" << pt1.x() + 0.5 << " " << pt2.x() + 0.5 << "], [" << pt1.y() + 0.5 << " " << pt2.y() + 0.5 << "], 'color', 'black', 'linewidth', 2);" << std::endl;
	}

	for(unsigned int i = 0; i < 10; i++)
	{
		ofs << "mov(:," << cpt << ") = getframe(ax);" << std::endl;
		cpt++;
	}

	ofs << "myVideo = VideoWriter('myfile.avi','Uncompressed AVI');" << std::endl;
	ofs << "myVideo.FrameRate = 5;" << std::endl;
	ofs << "myVideo.Quality = 100;" << std::endl;
	ofs << "open(myVideo);" << std::endl;
	ofs << "writeVideo(myVideo,mov);" << std::endl;
	ofs << "close(myVideo);" << std::endl;
}*/

skeleton::GraphSkel2d::Ptr algorithm::skeletonization::SpherePropagation2D_old(const boundary::DiscreteBoundary<2>::Ptr disbnd, const shape::DiscreteShape<2>::Ptr disshp, const OptionsSphProp &options)
{
	srand(time(NULL));
	skeleton::GraphSkel2d::Ptr skel(new skeleton::GraphSkel2d(skeleton::model::Classic<2>()));

	/*std::ofstream ofs("draw.m",std::ofstream::out);
	
	ofs << "clear all" << std::endl;
	ofs << "close all" << std::endl;
	
	ofs << "img = imread('cur.png');" << std::endl;

	ofs << "dist = zeros(size(img,1),size(img,2));" << std::endl;
	for(unsigned int l = 0; l < disshp->getHeight(); l++)
		for(unsigned int c = 0; c < disshp->getWidth(); c++)
		{
			Eigen::Vector2d P(c,l);
			if(disshp->isIn(mathtools::affine::Point<2>(P)))
			{
				double val = fieldValue(disbnd,P,options.noise);
				ofs << "dist(" << l+1 << ", " << c+1 << ") = " << val << ";" << std::endl;
			}
			else 
				ofs << "dist(" << l+1 << ", " << c+1 << ") = 0.0;" << std::endl;
		}
	
	ofs << "figure(2)" << std::endl;
	ofs << "imagesc(img)" << std::endl;
	ofs << "hold on" << std::endl;*/

	/*ofs << "figure(2)" << std::endl;
	ofs << "hold off" << std::endl;
	ofs << "imagesc(Dist);" << std::endl;
	ofs << "axis equal" << std::endl;
	ofs << "axis off" << std::endl;
	ofs << "hold on" << std::endl;
	ofs << "colorbar" << std::endl;
	
	ofs << "figure(3)" << std::endl;
	ofs << "hold off" << std::endl;
	ofs << "imagesc(Dist2);" << std::endl;
	ofs << "axis equal" << std::endl;
	ofs << "axis off" << std::endl;
	ofs << "hold on" << std::endl;
	ofs << "colorbar" << std::endl;*/
	
	/*ofs << "figure(2)" << std::endl;
	ofs << "[dx,dy] = gradient(dist);" << std::endl;
	ofs << "dist = dist - min(dist(:));" << std::endl;
	ofs << "dist(:) = uint8(63*dist(:)/max(dist(:)))+1;" << std::endl;
	ofs << "dist = ind2rgb(dist,jet);" << std::endl;
	ofs << "dist(img(:) < 125) = 1.0;" << std::endl;
	ofs << "imagesc(dist);" << std::endl;
	ofs << "hold on" << std::endl;
	ofs << "axis equal" << std::endl;
	ofs << "axis off" << std::endl;
	ofs << "[x,y] = meshgrid(1:size(dist,2),1:size(dist,1));" << std::endl;
	ofs << "mask = repmat([false],size(dist,1),size(dist,2));" << std::endl;
	ofs << "mask(10:10:end,10:10:end) = true;" << std::endl;
	ofs << "mask = (sum(img,3)>125) & mask;" << std::endl;
	ofs << "quiver(x(mask), y(mask), dx(mask), dy(mask),'black')" << std::endl;

	ofs << "disteuc = bwdist(sum(img,3) < 125);" << std::endl;
	ofs << "figure(3)" << std::endl;
	ofs << "[dx,dy] = gradient(disteuc);" << std::endl;
	ofs << "disteuc = disteuc - min(disteuc(:));" << std::endl;
	ofs << "disteuc(:) = uint8(63*disteuc(:)/max(disteuc(:)))+1;" << std::endl;
	ofs << "disteuc = ind2rgb(disteuc,jet);" << std::endl;
	ofs << "disteuc(img(:) < 125) = 1.0;" << std::endl;
	ofs << "imagesc(disteuc);" << std::endl;
	ofs << "hold on" << std::endl;
	ofs << "axis equal" << std::endl;
	ofs << "axis off" << std::endl;
	ofs << "[x,y] = meshgrid(1:size(disteuc,2),1:size(disteuc,1));" << std::endl;
	ofs << "mask = repmat([false],size(disteuc,1),size(disteuc,2));" << std::endl;
	ofs << "mask(10:10:end,10:10:end) = true;" << std::endl;
	ofs << "mask = (sum(img,3)>125) & mask;" << std::endl;
	ofs << "quiver(x(mask), y(mask), dx(mask), dy(mask),'black')" << std::endl;

	ofs << "figure;" << std::endl;
	ofs << "imagesc(img);" << std::endl;
	ofs << "hold on" << std::endl;
	ofs << "axis equal" << std::endl;
	ofs << "axis off" << std::endl;*/
	

	/**
	 *  First center computation
	 */
	Eigen::Vector2d C;
	bool fini = false;
	do
	{
		C.x() = rand()%disshp->getWidth();
		C.y() = rand()%disshp->getHeight();
		if(disshp->isIn(mathtools::affine::Point<2>(C)) && fieldValue(disbnd,C,options.noise) > 3.0*options.noise)
			fini = true;
	}while(!fini);
	/*C.x() = 301;
	C.y() = 231;*/

	/*Eigen::Vector2d coords = disbnd->getCoordinates(330);
	ofs << "plot(" << coords.x() + 0.5 << ", " << coords.y() + 0.5 << ",'+');" << std::endl;
	coords = disbnd->getCoordinates(332);
	ofs << "plot(" << coords.x() + 0.5 << ", " << coords.y() + 0.5 << ",'*');" << std::endl;
	coords = disbnd->getCoordinates(353);
	ofs << "plot(" << coords.x() + 0.5 << ", " << coords.y() + 0.5 << ",'*');" << std::endl;
	ofs << "plot(" << 360.032 + 0.5 << ", " << 143.163 + 0.5 << ", '*', 'color', 'red');" << std::endl;*/
	 
	//ofs << "plot(" << C.x() + 0.5 << ", " << C.y() + 0.5 << ", '*', 'color', 'red');" << std::endl;
	//C.x() = 425.395;
	//C.y() = 157.712;
	//C.x() = 425.155;
	//C.y() = 157.512;
	if(closestCenter(disbnd,C,options.noise))
	{
		std::cout << C.x() << " " << C.y() << std::endl;
		//ofs << "plot(" << C.x() + 0.5 << ", " << C.y() + 0.5 << ", '*', 'color', 'red');" << std::endl;
		MovingCenter mov(C);
		mov.computeTangencyData(disbnd,options.noise,options.alpha);

		double rad = mov.getRadius();
		//ofs << "draw_circle(" << C.y() + 0.5 << ", " << C.x() + 0.5 << ", " << rad << ",'green');" << std::endl;
		unsigned int nbdir = mov.getNbdir();
		
		std::list<std::tuple<unsigned int,MovingCenter,unsigned int> > lctr;
		
		unsigned int ind = skel->addNode(Eigen::Vector3d(C.x(),C.y(),rad));
		
		for(unsigned int i = 0; i < nbdir; i++)
		{
			std::tuple<unsigned int,MovingCenter,unsigned int> cdir;
			cdir = std::make_tuple(ind,mov,i);
			lctr.push_back(cdir);
		}
		
		//if(nbdir != 2)
		//	infosLocales(vbnd,C,rad,options.noise,disshp,v_ind,ofs);
		
		unsigned int cpt = 100000;
		if(lctr.size() != 0)
		do
		{
			std::tuple<unsigned int,MovingCenter,unsigned int> cdir = *(lctr.begin());
			lctr.pop_front();
			
			bool skip = false;
			for(std::list<std::tuple<unsigned int,MovingCenter,unsigned int> >::iterator it = lctr.begin(); it != lctr.end();)
			{
				if(std::get<1>(*it).isInNext(std::get<2>(*it),std::get<1>(cdir).getCenter()))
				{
					if(MovingCenter::intersect(std::get<1>(*it),std::get<2>(*it),std::get<1>(cdir),std::get<2>(cdir)))
					{
						//ofs << "plot([" << std::get<1>(cdir).getCenter().x() + 0.5 << ", " << std::get<1>(*it).getCenter().x() + 0.5 << "] , [" 
						//				<< std::get<1>(cdir).getCenter().y() + 0.5 << ", " << std::get<1>(*it).getCenter().y() + 0.5 << "] ,'blue');" << std::endl;
						skel->addEdge(std::get<0>(cdir),std::get<0>(*it));
						skip = true;
						it = lctr.erase(it);
					}
					else
					{
						it++;
					}
				}
				else
				{
					it++;
				}
			}
			
			if(!skip)
			{
				MovingCenter movn;
				std::list<unsigned int> lnext;

				//if(std::get<1>(cdir).propagateaccurate(disbnd,std::get<2>(cdir),options.noise,movn,lnext,options.alpha))
				if(std::get<1>(cdir).propagate(disbnd,std::get<2>(cdir),options.noise,movn,lnext,options.alpha))
				{
					unsigned int indcur = skel->addNode(Eigen::Vector3d(movn.getCenter().x(),movn.getCenter().y(),movn.getRadius()));
					skel->addEdge(indcur,std::get<0>(cdir));

					/*unsigned int nb = 0;
					double var = 0.0;
					for(unsigned int i = 0; i < movn.getTgt().size(); i++)
					{
						for(std::list<unsigned int>::const_iterator it = movn.getTgt()[i].begin(); it != movn.getTgt()[i].end(); it++)
						{
							double dist = (movn.getCenter() - disbnd->getVertex(*it).getCoords()).norm();
							var += (dist - movn.getRadius())*(dist - movn.getRadius());
							nb++;
						}
					}
					var /= (double)nb;
					std::cout << var << " " << movn.getTgt().size() << std::endl;
					
					if(sqrt(var) > options.noise)
						ofs << "draw_circle(" << movn.getCenter().y() + 0.5 << ", " << movn.getCenter().x() + 0.5 << ", " << movn.getRadius() << ",'red');" << std::endl;
					else*/
						//ofs << "draw_circle(" << movn.getCenter().y() + 0.5 << ", " << movn.getCenter().x() + 0.5 << ", " << movn.getRadius() << ",'green');" << std::endl;
					//ofs << "plot(" << movn.getCenter().x() + 0.5 << ", " << movn.getCenter().y() + 0.5 << ",'+','color','blue');" << std::endl;
					/*ofs << "plot([" << std::get<1>(cdir).getCenter().x() + 0.5 << ", " << movn.getCenter().x() + 0.5 << "] , [" 
						<< std::get<1>(cdir).getCenter().y() + 0.5 << ", " << movn.getCenter().y() + 0.5 << "] ,'blue');" << std::endl;*/
					for(std::list<unsigned int>::iterator it = lnext.begin(); it != lnext.end(); it++)
					{
						std::tuple<unsigned int,MovingCenter,unsigned int> cdirn;
						cdirn = std::make_tuple(indcur,movn,*it);
						lctr.push_back(cdirn);
					}

					/*if(v_indcur.size() != 2)
					{
						std::cout << v_indcur.size() << " sides" << std::endl;
						infosLocales(vbnd,Cmov,radmov,options.noise,disshp,v_indcur,ofs);
						ofs << "draw_circle(" << Cmov.y() + 0.5 << ", " << Cmov.x() + 0.5 << ", " << radmov << ",'blue');" << std::endl;
					}*/
					/*if(v_indcur.size() == 2)
						ofs << "plot(" << Cmov.x() + 0.5 << " , " << Cmov.y() + 0.5 << ",'+','color','blue');" << std::endl;
					else if(v_indcur.size() > 2)
						ofs << "plot(" << Cmov.x() + 0.5 << " , " << Cmov.y() + 0.5 << ",'*','color','green');" << std::endl;*/
				}
				else
				{
					//infosLocales(disbnd,std::get<1>(cdir),options.noise,ofs);
					//ofs << "draw_circle(" << std::get<1>(cdir).getCenter().y() + 0.5 << ", " << std::get<1>(cdir).getCenter().x() + 0.5 << ", " << std::get<1>(cdir).getRadius() << ",'blue');" << std::endl;

					////movn.computeTangencyData(disbnd,options.noise,options.alpha);
					//ofs << "draw_circle(" << movn.getCenter().y() + 0.5 << ", " << movn.getCenter().x() + 0.5 << ", " << movn.getRadius() << ",'red');" << std::endl;
					//ofs << "plot([" << std::get<1>(cdir).getCenter().x() + 0.5 << ", " << movn.getCenter().x() + 0.5 << "] , [" 
					//	<< std::get<1>(cdir).getCenter().y() + 0.5 << ", " << movn.getCenter().y() + 0.5 << "] ,'red');" << std::endl;
					//infosLocales(disbnd,movn,options.noise,ofs);
				}
			}
			cpt--;
		}while(!lctr.empty() && cpt != 0);
	}
	
	std::list<std::pair<unsigned int,unsigned int> > ledg;
	skel->getAllEdges(ledg);
	
	/*for(std::list<std::pair<unsigned int, unsigned int> >::iterator it = ledg.begin(); it != ledg.end(); it++)
	{
		Eigen::Vector3d v1 = skel->getNode(it->first);
		Eigen::Vector3d v2 = skel->getNode(it->second);
		ofs << "plot([" << v1.x() + 0.5 << ", " << v2.x() + 0.5 << "] , [" 
						<< v1.y() + 0.5 << ", " << v2.y() + 0.5 << "] ,'blue');" << std::endl;
	}

	ofs << "axis([ 0 size(img,2) 0 size(img,1) ]);" << std::endl;

	makeMovie(skel,disbnd,ofs,options);
	ofs.close();*/

	return Subdiv(skel,disbnd,disshp,options);
}

skeleton::GraphSkel2d::Ptr algorithm::skeletonization::SpherePropagation2D(const boundary::DiscreteBoundary<2>::Ptr disbnd, const shape::DiscreteShape<2>::Ptr &disshp, const OptionsSphProp &options)
{
	srand(time(NULL));
	skeleton::GraphSkel2d::Ptr skel(new skeleton::GraphSkel2d(skeleton::model::Classic<2>()));

	/**
	 *  First center computation
	 */
	//Eigen::Vector2d C(6.0,4.0);
	Eigen::Vector2d C;
	bool fini = false;
	do
	{
		C.x() = rand()%disshp->getWidth();
		C.y() = rand()%disshp->getHeight();
		if(disshp->isIn(mathtools::affine::Point<2>(C)) && fieldValue(disbnd,C,options.noise) > 3.0*options.noise)
			fini = true;
	}while(!fini);
	/*C.x() = 301;
	C.y() = 231;*/
	
	if(closestCenter(disbnd,C,options.noise))
	{
		MovingCenter mov(C);

		mov.computeTangencyData(disbnd,options.noise,options.alpha);
		unsigned int indc = closestInd(disbnd,C);

		unsigned int nbdir = mov.getNbdir();
		
		std::list<std::tuple<unsigned int,MovingCenter,unsigned int> > lctr;
		
		unsigned int ind = skel->addNode(Eigen::Vector3d(mov.getCenter().x(),mov.getCenter().y(),mov.getRadius()));
		
		for(unsigned int i = 0; i < nbdir; i++)
		{
			std::tuple<unsigned int,MovingCenter,unsigned int> cdir;
			cdir = std::make_tuple(ind,mov,i);
			lctr.push_back(cdir);
		}
		
		//if(nbdir != 2)
		//	infosLocales(vbnd,C,rad,options.noise,disshp,v_ind,ofs);
		
		unsigned int cpt = 100000;
		if(lctr.size() != 0)
		do
		{
			std::tuple<unsigned int,MovingCenter,unsigned int> cdir = *(lctr.begin());
			lctr.pop_front();
			
			bool skip = false;
			for(std::list<std::tuple<unsigned int,MovingCenter,unsigned int> >::iterator it = lctr.begin(); it != lctr.end();)
			{
				if(std::get<1>(*it).isInNext(std::get<2>(*it),std::get<1>(cdir).getCenter()))
				{
					if(MovingCenter::intersect(std::get<1>(*it),std::get<2>(*it),std::get<1>(cdir),std::get<2>(cdir)))
					{
						skel->addEdge(std::get<0>(cdir),std::get<0>(*it));
						skip = true;
						it = lctr.erase(it);
					}
					else
					{
						it++;
					}
				}
				else
				{
					it++;
				}
			}
			
			if(!skip)
			{
				MovingCenter movn;
				std::list<unsigned int> lnext;

				//if(std::get<1>(cdir).propagateaccurate(disbnd,std::get<2>(cdir),options.noise,movn,lnext,options.alpha))
				if(std::get<1>(cdir).propagate(disbnd,std::get<2>(cdir),options.noise,movn,lnext,options.alpha))
				{
					unsigned int indcur = skel->addNode(Eigen::Vector3d(movn.getCenter().x(),movn.getCenter().y(),movn.getRadius()));
					unsigned int indc = closestInd(disbnd,movn.getCenter());
					skel->addEdge(indcur,std::get<0>(cdir));

					/*unsigned int nb = 0;
					double var = 0.0;
					for(unsigned int i = 0; i < movn.getTgt().size(); i++)
					{
						for(std::list<unsigned int>::const_iterator it = movn.getTgt()[i].begin(); it != movn.getTgt()[i].end(); it++)
						{
							double dist = (movn.getCenter() - disbnd->getVertex(*it).getCoords()).norm();
							var += (dist - movn.getRadius())*(dist - movn.getRadius());
							nb++;
						}
					}
					var /= (double)nb;
					std::cout << var << " " << movn.getTgt().size() << std::endl;
					
					if(sqrt(var) > options.noise)
						ofs << "draw_circle(" << movn.getCenter().y() + 0.5 << ", " << movn.getCenter().x() + 0.5 << ", " << movn.getRadius() << ",'red');" << std::endl;
					else*/
					for(std::list<unsigned int>::iterator it = lnext.begin(); it != lnext.end(); it++)
					{
						std::tuple<unsigned int,MovingCenter,unsigned int> cdirn;
						cdirn = std::make_tuple(indcur,movn,*it);
						lctr.push_back(cdirn);
					}

					/*if(v_indcur.size() != 2)
					{
						std::cout << v_indcur.size() << " sides" << std::endl;
						infosLocales(vbnd,Cmov,radmov,options.noise,disshp,v_indcur,ofs);
						ofs << "draw_circle(" << Cmov.y() + 0.5 << ", " << Cmov.x() + 0.5 << ", " << radmov << ",'blue');" << std::endl;
					}*/
					/*if(v_indcur.size() == 2)
						ofs << "plot(" << Cmov.x() + 0.5 << " , " << Cmov.y() + 0.5 << ",'+','color','blue');" << std::endl;
					else if(v_indcur.size() > 2)
						ofs << "plot(" << Cmov.x() + 0.5 << " , " << Cmov.y() + 0.5 << ",'*','color','green');" << std::endl;*/
				}
			}
			cpt--;
		}while(!lctr.empty() && cpt != 0);
	}

	/*std::ofstream ofs("draw.m",std::ofstream::out);
	
	ofs << "clear all" << std::endl;
	ofs << "close all" << std::endl;

	std::list<Eigen::Vector2d> lpt;
	disbnd->getVerticesVector(lpt);
	
	ofs << "B=[";

	for(std::list<Eigen::Vector2d>::iterator it = lpt.begin(); it != lpt.end(); it++)
	{
		ofs << it->x() << " " << it->y() << std::endl;
	}
	ofs << lpt.begin()->x() << " " << lpt.begin()->y() << "];" << std::endl;
	
	ofs << "figure" << std::endl;
	ofs << "plot(B(:,1),B(:,2));" << std::endl;
	ofs << "axis equal" << std::endl;
	ofs << "hold on" << std::endl;
	ofs << "plot(B(:,1),B(:,2),'+');" << std::endl;
	
	std::list<std::pair<unsigned int,unsigned int> > ledg;
	skel->getAllEdges(ledg);
	
	for(std::list<std::pair<unsigned int,unsigned int> >::iterator it = ledg.begin(); it != ledg.end(); it++)
	{
		Eigen::Vector3d v1 = skel->getNode(it->first);
		Eigen::Vector3d v2 = skel->getNode(it->second);
		ofs << "plot([" << v1.x() + 0.5 << ", " << v2.x() + 0.5 << "] , [" 
						<< v1.y() + 0.5 << ", " << v2.y() + 0.5 << "] ,'blue');" << std::endl;
	}*/

	
	return skel;//Subdiv(skel,disbnd,disshp,options);
}

skeleton::GraphSkel2d::Ptr algorithm::skeletonization::Subdiv(const skeleton::GraphSkel2d::Ptr grskel, const boundary::DiscreteBoundary<2>::Ptr disbnd, const shape::DiscreteShape<2>::Ptr disshp, const OptionsSphProp &options)
{
	skeleton::GraphSkel2d::Ptr grskelcpy(new skeleton::GraphSkel2d(*grskel));
	std::list<std::pair<unsigned int,unsigned int> > ledg;
	grskelcpy->getAllEdges(ledg);
	
	for(std::list<std::pair<unsigned int, unsigned int> >::iterator it = ledg.begin(); it != ledg.end(); it++)
	{
		Eigen::Vector2d v1 = grskel->getNode(it->first).block<2,1>(0,0);
		Eigen::Vector2d v2 = grskel->getNode(it->second).block<2,1>(0,0);
		
		Eigen::Vector2d ctr = (v1 + v2)*0.5;
		Eigen::Vector2d dir = Eigen::Vector2d((v2.y() - v1.y()),-(v2.x() - v1.x())).normalized();
		
		double rad = algorithm::skeletonization::fieldValue(disbnd,ctr,options.noise);
		
		double dist0 = -rad*0.5, dist1 = rad*0.5;
		if(dist0 < -2.0*options.noise) dist0 = -2.0*options.noise;
		if(dist1 > 2.0*options.noise) dist1 = 2.0*options.noise;
		
		Eigen::Vector2d grad0;
		fieldValue(disbnd,ctr+dist0*dir,grad0,options.noise);
		Eigen::Vector2d grad1;
		fieldValue(disbnd,ctr+dist1*dir,grad1,options.noise);

		while(dist1 - dist0 > 0.1*options.noise)
		{
			double dist = (dist1 + dist0)/2.0;
			Eigen::Vector2d Cmov = ctr + dist*dir;
			Eigen::Vector2d grad;
			fieldValue(disbnd,Cmov,grad,options.noise);
			if((grad0 - grad).squaredNorm() < (grad1 - grad).squaredNorm())
				dist0 = dist;
			else
				dist1 = dist;
		}
		double dist = (dist1 + dist0)/2.0;
		Eigen::Vector2d Ccur = ctr + dist*dir;
		double radcur = fieldValue(disbnd,Ccur,options.noise);
		
		unsigned int indcur = grskelcpy->addNode(Eigen::Vector3d(Ccur.x(),Ccur.y(),radcur));
		grskelcpy->remEdge(it->first,it->second);
		grskelcpy->addEdge(it->first,indcur);
		grskelcpy->addEdge(it->second,indcur);
		
		
	}
	
	return grskelcpy;
}
