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
 *  \file SpherePropagation2D.h
 *  \brief Defines functions to compute 2d skeleton with sphere propagation algorithm
 *  \author Bastien Durix
 */

#ifndef _SPHEREPROPAGATION_H_
#define _SPHEREPROPAGATION_H_

#include <skeleton/Skeletons.h>
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
		/**
		 *  \brief Sphere propagation options structure
		 */
		struct OptionsSphProp
		{
			/**
			 *  \brief Boundary noise
			 */
			double noise;

			/**
			 *  \brief Hausdorff distance
			 */
			double alpha;

			/**
			 *  \brief Default constructor
			 */
			OptionsSphProp(double noise_ = 1.0, double alpha_=2.5) : noise(noise_), alpha(alpha_)
			{}
		};
		
		/**
		 *  \brief 2D skeletonization, by sphere propagation
		 *
		 *  \param disbnd   discrete boundary of the shape
		 *  \param disshp   discrete shape
		 *  \param options  options of the algorithm
		 *
		 *  \return pointer to the computed 2d graph skeleton
		 */
		skeleton::GraphSkel2d::Ptr SpherePropagation2D_old(const boundary::DiscreteBoundary<2>::Ptr disbnd, const shape::DiscreteShape<2>::Ptr disshp, const OptionsSphProp &options = OptionsSphProp());

		skeleton::GraphSkel2d::Ptr SpherePropagation2D(const boundary::DiscreteBoundary<2>::Ptr disbnd, const shape::DiscreteShape<2>::Ptr &dissh, const OptionsSphProp &options = OptionsSphProp());

		skeleton::GraphSkel2d::Ptr Subdiv(const skeleton::GraphSkel2d::Ptr grskel, const boundary::DiscreteBoundary<2>::Ptr disbnd, const shape::DiscreteShape<2>::Ptr disshp, const OptionsSphProp &options = OptionsSphProp());
	}
}

#endif //_SPHEREPROPAGATION_H_
