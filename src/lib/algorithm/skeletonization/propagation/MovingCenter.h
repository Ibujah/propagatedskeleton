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
 *  \file MovingCenter.h
 *  \brief Defines object containing informations about moving centers
 *  \author Bastien Durix
 */

#ifndef _MOVINGCENTER_H_
#define _MOVINGCENTER_H_

#include <Eigen/Dense>
#include <boundary/DiscreteBoundary2.h>

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
		/**
		 *  \brief Contains properties about moving centers
		 */
		class MovingCenter
		{
			protected:
				/**
				 *  \brief Position of the center
				 */
				Eigen::Vector2d m_center;

				/**
				 *  \brief Radius of the circle at the center
				 */
				double m_rad;
				
				/**
				 *  \brief Vector of tangency points to the boundary
				 */
				std::vector<std::list<unsigned int> > m_tgtinds;

				/**
				 *  \brief Vector of angles
				 */
				std::vector<std::pair<double, double> > m_ang;
				
				void correction(const boundary::DiscreteBoundary<2>::Ptr disbnd);

			public:
				/**
				 *  \brief Default constructor
				 */
				MovingCenter();

				/**
				 *  \brief Constructor
				 *
				 *  \param center  position of the center
				 */
				MovingCenter(const Eigen::Vector2d &center);

				/**
				 *  \brief Infers data from boundary
				 *
				 *  \param disbnd  boundary to study
				 *  \param noise   amout of noise
				 */
				void computeTangencyData(const boundary::DiscreteBoundary<2>::Ptr disbnd, double noise, double alpha);

				/**
				 *  \brief Infers data from boundary and previous center
				 *  
				 *  \param mov     previous moving center
				 *  \param dir     direction in which mov was propagated
				 *  \param disbnd  boundary to study
				 *  \param noise   amout of noise
				 */
				void computeTangencyData(const MovingCenter &mov, unsigned int dir, const boundary::DiscreteBoundary<2>::Ptr disbnd, double noise, double alpha);
				void computeTangencyDataaccurate(const MovingCenter &mov, unsigned int dir, const boundary::DiscreteBoundary<2>::Ptr disbnd, double noise, double alpha);

				/**
				 *  \brief Propagates circle center
				 *
				 *  \param disbnd  boundary to study
				 *  \param dir     direction in which propagate the center
				 *  \param noise   amout of noise
				 *  \param mov     output propagated center
				 *  \param lnext   output possible directions from thise moving center
				 *
				 *  \return true if center has been computed correctly
				 */
				bool propagate(const boundary::DiscreteBoundary<2>::Ptr disbnd, unsigned int dir, double noise, MovingCenter &mov, std::list<unsigned int> &lnext, double alpha) const;

				bool propagateaccurate(const boundary::DiscreteBoundary<2>::Ptr disbnd, unsigned int dir, double noise, MovingCenter &mov, std::list<unsigned int> &lnext, double alpha) const;

				/**
				 *  \brief Number of directions getter
				 *
				 *  \return number of directions
				 */
				unsigned int getNbdir() const;

				/**
				 *  \brief Center getter
				 *
				 *  \return center
				 */
				const Eigen::Vector2d& getCenter() const;

				/**
				 *  \brief Radius getter
				 *
				 *  \return radius
				 */
				double getRadius() const;
				
				/**
				 *  \brief Tangents getter
				 *
				 *  \return vector of tangents
				 */
				const std::vector<std::list<unsigned int> >& getTgt() const;
				
				/**
				 *  \brief Tangency angles getter
				 *
				 *  \return vector of tangency angles
				 */
				const std::vector<std::pair<double,double> >& getAng() const;

				/**
				 *  \brief Verifies if a point is in the cône of next possible
				 *
				 *  \param dir  direction of the cone
				 *  \param C    coordinates of the point
				 */
				bool isInNext(unsigned int dir, const Eigen::Vector2d &C) const;
				
				/**
				 *  \brief Verifies intersection between two moving centers
				 *
				 *  \param mov1  first moving center
				 *  \param dir1  first direction
				 *  \param mov2  second moving center
				 *  \param dir2  second direction
				 *
				 *  \return true if mov1 and mov2 intersect
				 */
				static bool intersect(const MovingCenter &mov1, unsigned int dir1, const MovingCenter &mov2, unsigned int dir2, unsigned int offset1 = 1);
		};
	}
}

#endif //_MOVINGCENTER_H_
