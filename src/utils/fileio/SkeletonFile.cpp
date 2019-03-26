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
 *  \file SkeletonFile.cpp
 *  \brief Defines skeleton file writer
 *  \author Bastien Durix
 */

#include "SkeletonFile.h"

#include <string>
#include <fstream>
#include <sstream>

void fileio::WriteSkeleton2D(const skeleton::GraphSkel2d::Ptr skel, const std::string &filename)
{
	std::ofstream file(filename);

	if(file)
	{
		std::list<unsigned int> lpts;
		skel->getAllNodes(lpts);
		std::vector<unsigned int> vpts(lpts.begin(),lpts.end());
		std::map<unsigned int,unsigned int> mind;
		
		std::list<std::pair<unsigned int,unsigned int> > ledg;
		skel->getAllEdges(ledg);
		
		file << vpts.size() << " " << ledg.size() << std::endl;
		file << "\%points" << std::endl;
		for(unsigned int i = 0; i < vpts.size(); i++)
		{
			Eigen::Vector3d pt = skel->getNode(vpts[i]);
			file << pt.x() << " " << pt.y() << " " << pt.z() << std::endl;
			mind[vpts[i]] = i;

		}
		file << std::endl << std::endl;
		
		file << "\%edges" << std::endl;
		for(std::list<std::pair<unsigned int,unsigned int> >::iterator it = ledg.begin(); it != ledg.end(); it++)
		{
			file << mind[it->first] << " " << mind[it->second] << std::endl;;
		}
		file << std::endl;

		file.close();
	}
}
