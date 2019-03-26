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
 *  \file DistanceField.h
 *  \brief Defines functions related to distance field
 *  \author Bastien Durix
 */

#ifndef _DISTANCEFIELD_H_
#define _DISTANCEFIELD_H_

#include <boundary/DiscreteBoundary2.h>
#include <shape/DiscreteShape.h>

/**
 *  \brief Lots of algorithms
 */
namespace algorithm
{
	/**
	 *  \brief skeletonization algorithms
	 */
	namespace skeletonization
	{
		unsigned int closestInd(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C);
		/**
		 *  \brief Estimates value of a distance field, up to a noise
		 *
		 *  \param disbnd  boundary of the shape
		 *  \param C       actual point inside the shape
		 *  \param noise   amount of noise
		 *
		 *  \return distance field value, corresponding to the distance to the boundary
		 */
		double fieldValue(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C, double noise);

		/**
		 *  \brief Estimates value and gradient of a distance field, up to a noise
		 *
		 *  \param disbnd  boundary of the shape
		 *  \param C       actual point inside the shape
		 *  \param grad    output estimated gradient at C
		 *  \param noise   amount of noise
		 *
		 *  \return distance field value, corresponding to the distance to the boundary
		 */
		double fieldValue(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C, Eigen::Vector2d &grad, double noise);
		
		/**
		 *  \brief Computes tangency angles between circle and boundaries
		 *
		 *  \param disbnd  discrete boundary of the shape
		 *  \param C       input initialisation and output found point
		 *  \param rad     output found distance to the boundary
		 *  \param noise   amount of noise
		 *  \param v_ind   output vector of indices
		 */
		void tangencyBoundary(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C, const double &rad, double noise, std::vector<std::list<unsigned int> > &v_ind, double fac);

		/**
		 *  \brief Computes closest point on medial axis to an initialisation
		 *
		 *  \param disbnd  discrete boundary of the shape
		 *  \param C       input initialisation and output found point
		 *  \param noise   amount of noise
		 *
		 *  \return false if the minimization has diverged
		 */
		bool closestCenter(const boundary::DiscreteBoundary<2>::Ptr disbnd, Eigen::Vector2d &C, double noise);

		/**
		 *  \brief Computes intersection of medial axis and circle arc
		 *
		 *  \param disbnd  vector of coordinates of the boundary points
		 *  \param C       center of circle
		 *  \param pang    delimitation of circle arc
		 *  \param dist    distance to C
		 *  \param Cmov    output center
		 *  \param noise   amount of noise
		 */
		void closestCenterOnArc(const boundary::DiscreteBoundary<2>::Ptr disbnd, const Eigen::Vector2d &C, const std::pair<double,double> &pang, double dist, Eigen::Vector2d &Cmov, double noise);

		/**
		 *  \brief Computes discontinuities on arc
		 *
		 *  \param disbnd  vector of coordinates of the boundary points
		 *  \param C       center of circle
		 *  \param pang    delimitation of circle arc
		 *  \param dist    distance to C
		 *  \param vecC    output centers
		 *  \param noise   amount of noise
		 */
		void discontinuitiesOnArc(const boundary::DiscreteBoundary<2>::Ptr disbnd,
								  const Eigen::Vector2d &C,
								  const std::pair<double,double> &pang,
								  double dist,
								  std::vector<Eigen::Vector2d> &vecC,
								  std::vector<std::pair<unsigned int,unsigned int> > &vecInd,
								  double noise);
	}
}

#endif //_DISTANCEFIELD_H_
