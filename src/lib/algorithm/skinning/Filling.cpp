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
 *  \file Filling.h
 *  \brief Computes a discrete shape associated to a continuous skeleton
 *  \author Bastien Durix
 */

#include "Filling.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <mathtools/geometry/euclidian/HyperSphere.h>
#include <mathtools/affine/Point.h>

void algorithm::skinning::Filling(shape::DiscreteShape<2>::Ptr shape, const skeleton::GraphSkel2d::Ptr grskl)
{
	cv::Mat im_shape(shape->getHeight(),shape->getWidth(),CV_8U,&shape->getContainer()[0]);
	
	std::list<unsigned int> lind;
	grskl->getAllNodes(lind);

	for(std::list<unsigned int>::iterator it = lind.begin(); it != lind.end(); it++)
	{
		mathtools::geometry::euclidian::HyperSphere<2> sph = grskl->getNode<mathtools::geometry::euclidian::HyperSphere<2> >(*it);
		
    	cv::circle(im_shape,cv::Point2i(sph.getCenter().getCoords(shape->getFrame()).x()+0.5,sph.getCenter().getCoords(shape->getFrame()).y()+0.5),sph.getRadius(),255,-1);
	}
}
