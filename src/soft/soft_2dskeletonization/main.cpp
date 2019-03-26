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
 *  \brief 2D skeletonization
 *  \author Bastien Durix
 */

#include <boost/program_options.hpp>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>

#include <shape/DiscreteShape.h>
#include <boundary/DiscreteBoundary.h>
#include <skeleton/Skeletons.h>

#include <algorithm/extractboundary/NaiveBoundary.h>
#include <algorithm/skeletonization/propagation/SpherePropagation2D.h>
#include <algorithm/skeletonization/voronoi/VoronoiSkeleton2D.h>
#include <algorithm/skinning/Filling.h>
#include <algorithm/evaluation/ReprojError.h>
#include <algorithm/pruning/ScaleAxisTransform.h>
#include <algorithm/pruning/LambdaMedialAxis.h>
#include <algorithm/pruning/ThetaMedialAxis.h>

#include <displayopencv/DisplayShapeOCV.h>
#include <displayopencv/DisplayBoundaryOCV.h>
#include <displayopencv/DisplaySkeletonOCV.h>

std::tuple<double,double,int,int> EvalSkel(const shape::DiscreteShape<2>::Ptr dissh,
									   const boundary::DiscreteBoundary<2>::Ptr disbnd,
									   const skeleton::GraphSkel2d::Ptr skel)
{
	shape::DiscreteShape<2>::Ptr shp(new shape::DiscreteShape<2>(dissh->getWidth(),dissh->getHeight()));
	algorithm::skinning::Filling(shp,skel);
	
	double res = algorithm::evaluation::SymDiffArea(dissh,shp);
	double res2 = algorithm::evaluation::HausDist(skel,disbnd,dissh->getFrame());
	
	std::list<unsigned int> lnod;
	skel->getAllNodes(lnod);
	unsigned int nbbr = 0;
	for(std::list<unsigned int>::iterator it = lnod.begin(); it != lnod.end(); it++)
	{
		unsigned int deg = skel->getNodeDegree(*it);
		if(deg != 2)
			nbbr += deg;
	}
	nbbr /= 2;
	std::tuple<double,double,int,int> result = std::make_tuple(res*100.0,res2,skel->getNbNodes(),nbbr);
	
	return result;
}

int main(int argc, char** argv)
{
	std::string imgfile, fileskl, fileimg, filebnd;
	bool output = false;
	bool compare = false;
	bool eval = false;
	double alpha, noise;

	boost::program_options::options_description desc("OPTIONS");
	
	desc.add_options()
		("help", "Help message")
		("imgfile", boost::program_options::value<std::string>(&imgfile)->default_value("mask"), "Binary image file (*.png)")
		("output", boost::program_options::value<bool>(&output)->implicit_value(true), "Returns output images")
		("sigma", boost::program_options::value<double>(&noise)->default_value(1.0), "Shape noise (sigma parameter)")
		("alpha", boost::program_options::value<double>(&alpha)->default_value(2.1), "Skeleton precision (alpha parameter)")
		("fileimg", boost::program_options::value<std::string>(&fileimg)->default_value("skelpropagortho"), "Skeleton img file")
		("compare", boost::program_options::value<bool>(&compare)->implicit_value(true), "Compare result with pruning methods")
		;
	
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
	boost::program_options::notify(vm);
	
	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 0;
	}

	time_t start,end;
	double diff;

	cv::Mat shpimggray = cv::imread(imgfile,cv::ImreadModes::IMREAD_GRAYSCALE);
	cv::Mat shpimg;
	cv::threshold(shpimggray,shpimg,125,255,cv::THRESH_BINARY);
	
	// topological closure
	cv::Mat shpdil;
	cv::Mat element = cv::getStructuringElement(cv::MORPH_RECT,cv::Size(3,3),cv::Point(1,1));
	cv::dilate(shpimg,shpdil,element);
	cv::erode(shpdil,shpimg,element);
	shape::DiscreteShape<2>::Ptr dissh = shape::DiscreteShape<2>::Ptr(new shape::DiscreteShape<2>(shpimg.cols,shpimg.rows));
	cv::Mat cpymat(shpimg.rows,shpimg.cols,CV_8U,&dissh->getContainer()[0]);
	shpimg.copyTo(cpymat);
	cv::Mat image(shpimg.rows,shpimg.cols,CV_8UC3,cv::Scalar(255,255,255));
	
	displayopencv::DisplayDiscreteShape(dissh,image,dissh->getFrame(),cv::Scalar(0,0,255));
	
	boundary::DiscreteBoundary<2>::Ptr disbnd = algorithm::extractboundary::NaiveBoundary(dissh);
	
	auto start0 = std::chrono::steady_clock::now();
	algorithm::skeletonization::OptionsSphProp options(noise,alpha);
	skeleton::GraphSkel2d::Ptr grskelpropag = algorithm::skeletonization::SpherePropagation2D(disbnd,dissh,options);
	auto duration0 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start0);
	std::tuple<double,double,int,int> respropag = EvalSkel(dissh,disbnd,grskelpropag);
	int t0 = duration0.count();
	double A0 = std::get<0>(respropag); // sym area diff
	double H0 = std::get<1>(respropag); // Hausdorff dist
	int N0 = std::get<2>(respropag); // nb nodes
	int B0 = std::get<3>(respropag); // nb branches

	if(H0 > alpha)
		std::cerr << "Problem while computing propagation" << std::endl;
	
	std::cout << "Propagation skeleton computation: " << duration0.count() << "ms." << std::endl;

	shape::DiscreteShape<2>::Ptr shppropag(new shape::DiscreteShape<2>(dissh->getWidth(),dissh->getHeight()));
	algorithm::skinning::Filling(shppropag,grskelpropag);

	cv::Mat imagepropag;
	image.copyTo(imagepropag);
	displayopencv::DisplayDiscreteShape(shppropag,imagepropag,shppropag->getFrame(),cv::Scalar(0,255,0));
	displayopencv::DisplayDiscreteBoundary(disbnd,imagepropag,dissh->getFrame(),cv::Scalar(0,0,0));
	displayopencv::DisplayGraphSkeleton(grskelpropag,imagepropag,dissh->getFrame(),cv::Scalar(255,0,0));
	
	if(output)
	{
		std::ostringstream oss;
		oss.setf( std::ios::fixed, std:: ios::floatfield );
		oss.precision(1);
		oss << fileimg << ".png";
		cv::imwrite(oss.str(), imagepropag);
	}

	if(compare)
	{
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
		std::cout.precision(2);

		auto startv = std::chrono::steady_clock::now();
		skeleton::GraphSkel2d::Ptr grskelvoro = algorithm::skeletonization::VoronoiSkeleton2d(disbnd);
		auto durationv = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startv);
		std::tuple<double,double,int,int> resv = EvalSkel(dissh,disbnd,grskelvoro);
		int tv = durationv.count();
		double Av = std::get<0>(resv); // sym area diff
		double Hv = std::get<1>(resv); // Hausdorff dist
		int Nv = std::get<2>(resv); // nb nodes
		int Bv = std::get<3>(resv); // nb branches
	
		int t1;
		double A1;
		double H1;
		int N1;
		int B1;
		double sat = 1.2;
		auto start1 = std::chrono::steady_clock::now();
		skeleton::GraphSkel2d::Ptr grskelsat = algorithm::pruning::ScaleAxisTransform(grskelvoro,sat);
		auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start1);
		std::tuple<double,double,int,int> res1 = EvalSkel(dissh,disbnd,grskelsat);
		t1 = tv + duration1.count();
		A1 = std::get<0>(res1); // sym area diff
		H1 = std::get<1>(res1); // Hausdorff dist
		N1 = std::get<2>(res1); // nb nodes
		B1 = std::get<3>(res1); // nb branches

		if(output)
		{
			shape::DiscreteShape<2>::Ptr shpskel(new shape::DiscreteShape<2>(dissh->getWidth(),dissh->getHeight()));
			algorithm::skinning::Filling(shpskel,grskelsat);

			cv::Mat imagerec;
			image.copyTo(imagerec);
			displayopencv::DisplayDiscreteShape(shpskel,imagerec,shpskel->getFrame(),cv::Scalar(0,255,0));
			displayopencv::DisplayDiscreteBoundary(disbnd,imagerec,dissh->getFrame(),cv::Scalar(0,0,0));
			displayopencv::DisplayGraphSkeleton(grskelsat,imagerec,dissh->getFrame(),cv::Scalar(255,0,0));

			std::ostringstream oss;
			oss.setf( std::ios::fixed, std:: ios::floatfield );
			oss.precision(1);
			oss << "sat.png";
			cv::imwrite(oss.str(), imagerec);
		}

		int t2;
		double A2;
		double H2;
		int N2;
		int B2;
		double lambda = 2.0;
		auto start2 = std::chrono::steady_clock::now();
		skeleton::GraphSkel2d::Ptr grskellambda = algorithm::pruning::LambdaMedialAxis(grskelvoro,disbnd,lambda);
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start2);
		std::tuple<double,double,int,int> res2 = EvalSkel(dissh,disbnd,grskellambda);
		t2 = tv + duration2.count();
		A2 = std::get<0>(res2); // sym area diff
		H2 = std::get<1>(res2); // Hausdorff dist
		N2 = std::get<2>(res2); // nb nodes
		B2 = std::get<3>(res2); // nb branches

		if(output)
		{
			shape::DiscreteShape<2>::Ptr shpskel(new shape::DiscreteShape<2>(dissh->getWidth(),dissh->getHeight()));
			algorithm::skinning::Filling(shpskel,grskellambda);

			cv::Mat imagerec;
			image.copyTo(imagerec);
			displayopencv::DisplayDiscreteShape(shpskel,imagerec,shpskel->getFrame(),cv::Scalar(0,255,0));
			displayopencv::DisplayDiscreteBoundary(disbnd,imagerec,dissh->getFrame(),cv::Scalar(0,0,0));
			displayopencv::DisplayGraphSkeleton(grskellambda,imagerec,dissh->getFrame(),cv::Scalar(255,0,0));

			std::ostringstream oss;
			oss.setf( std::ios::fixed, std:: ios::floatfield );
			oss.precision(0);
			oss << "lambda.png";
			cv::imwrite(oss.str(), imagerec);
		}

		int t3;
		double A3;
		double H3;
		int N3;
		int B3;
		double theta = 100*M_PI/180.0;
		auto start3 = std::chrono::steady_clock::now();
		skeleton::GraphSkel2d::Ptr grskeltheta = algorithm::pruning::ThetaMedialAxis(grskelvoro,theta);
		auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start3);
		std::tuple<double,double,int,int> res3 = EvalSkel(dissh,disbnd,grskeltheta);
		t3 = tv + duration3.count();
		A3 = std::get<0>(res3); // sym area diff
		H3 = std::get<1>(res3); // Hausdorff dist
		N3 = std::get<2>(res3); // nb nodes
		B3 = std::get<3>(res3); // nb branches

		if(output)
		{
			shape::DiscreteShape<2>::Ptr shpskel(new shape::DiscreteShape<2>(dissh->getWidth(),dissh->getHeight()));
			algorithm::skinning::Filling(shpskel,grskeltheta);

			cv::Mat imagerec;
			image.copyTo(imagerec);
			displayopencv::DisplayDiscreteShape(shpskel,imagerec,shpskel->getFrame(),cv::Scalar(0,255,0));
			displayopencv::DisplayDiscreteBoundary(disbnd,imagerec,dissh->getFrame(),cv::Scalar(0,0,0));
			displayopencv::DisplayGraphSkeleton(grskeltheta,imagerec,dissh->getFrame(),cv::Scalar(255,0,0));

			std::ostringstream oss;
			oss.setf( std::ios::fixed, std:: ios::floatfield );
			oss.precision(2);
			oss << "theta.png";
			cv::imwrite(oss.str(), imagerec);
		}
		
		std::cout << "\t Time \t SymArea  Hausdorff  Branches \t Nodes" << std::endl;
		std::cout << "Propag : " << t0 << " \t " << A0 << " \t  " << H0 << " \t     " << B0 << " \t " << N0 << std::endl;
		std::cout << "Voro   : " << tv << " \t " << Av << " \t  " << Hv << " \t     " << Bv << " \t " << Nv << std::endl;
		std::cout << "SAT    : " << t1 << " \t " << A1 << " \t  " << H1 << " \t     " << B1 << " \t " << N1 << std::endl;
		std::cout << "Lambda : " << t2 << " \t " << A2 << " \t  " << H2 << " \t     " << B2 << " \t " << N2 << std::endl;
		std::cout << "Theta  : " << t3 << " \t " << A3 << " \t  " << H3 << " \t     " << B3 << " \t " << N3 << std::endl;
	}

	return 0;
}
