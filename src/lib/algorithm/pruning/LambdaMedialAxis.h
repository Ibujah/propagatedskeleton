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
 *  \file LambdaMedialAxis.h
 *  \brief Defines scale axis transform, a pruning algorithm
 *  \author Bastien Durix
 */

#ifndef _LAMBDAMEDIALAXIS_H_
#define _LAMBDAMEDIALAXIS_H_

#include <skeleton/Skeletons.h>
#include <boundary/DiscreteBoundary2.h>

/**
 *  \brief Lots of algorithms
 */
namespace algorithm
{
	/**
	 *  \brief Defines pruning algorithms
	 */
	namespace pruning
	{
		/**
		 *  \brief Scale Axis Transform pruning of classic skeleton
		 *
		 *  \param grskel : Skeleton to prune
		 *  \param scale : scale parameter
		 *
		 *  \return Pruned skeleton
		 */
		skeleton::GraphSkel2d::Ptr LambdaMedialAxis(const skeleton::GraphSkel2d::Ptr grskel, const boundary::DiscreteBoundary<2>::Ptr disbnd, const double &lambda);
	}
}

#endif //_LAMBDAMEDIALAXIS_H_
