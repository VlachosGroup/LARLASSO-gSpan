/*
    $Id: main.cpp,v 1.4 2004/05/21 05:50:13 taku-ku Exp $;
 
   Copyright (C) 2017 Geun Ho Gu, All rights reserved.
     This is free software with ABSOLUTELY NO WARRANTY.
  
   This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
    
   This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
   You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
     02111-1307, USA
*/
#include "gspan.h"
#include "cv.h"
#include "getopt.h"
#include <time.h>
#include <iomanip>
// TODO: parallelization

void usage (void)
{
	std::cout << "gspan+LarsLasso implementation by Geun Ho Gu" << std::endl;
	std::cout << std::endl;
	std::cout << "usage: gspan [-s minsup] [-D] [-e] [-w] [-m minnodes] [-M maxnodes] [-p maxparam] [-l minlam] [-c cvmethod] [-C cvparam] [-n cvrepeat] [-S]" << std::endl;
	std::cout << std::endl;
	std::cout << "options" << std::endl;
	std::cout << "  -h, show this usage help" << std::endl;
	std::cout << "  -s minsup, set the minimum support (absolute count)" << std::endl;
	std::cout << "  -D, use directed edges, default: undirected" << std::endl;
	std::cout << "  -e, output selected subgraphs" << std::endl;
	std::cout << "  -m minnodes, the minimum number of nodes in substructes (default: 0)" << std::endl; // Nope maximu
	std::cout << "  -M maxnodes, the maximum number of nodes in substructes (default: infinite)" << std::endl;
	std::cout << "  -p maxparam, maximum number of parameter mined (default: number of data - 1)" << std::endl;
	std::cout << "  -l minlam, minimum lambda considered. (default: 0)" << std::endl;
	std::cout << "  -c cvmethod, method of cross-validation. 0: no CV 1: hold-out; 2: k-fold; 3: leave-ont-out. (default: 0)" << std::endl;
	std::cout << "  -C cvparam, for hold-out, ratio of number of test set to number of entire data set, and, for k-fold, k value (default: 0)" << std::endl;
	std::cout << "  -n cvrepeat, repeatation for cv method. shuffle performed. (default: 1)" << std::endl;
	std::cout << "  -b search through DFS code tree in breadth first manner" << std::endl;
	
	/* List of constructor:
	* cvpartition(unsigned int n, unsigned int 0, float param)
	*		hold out method, where param is the ratio between test and training set
	*		n is equal to the number of data
	* cvpartition(unsigned int n, unsigned int 1, unsigned int param)
	*		KFold method, where param is the number of fold (k).
	* cvpartition(unsigned int n, unsigned int 2)
	*		Leave one out method
	*/
	std::cout << std::endl;

	std::cout << "The graphs are read from stdin, and have to be in this format:" << std::endl;
	std::cout << "t" << std::endl;
	std::cout << "y <continuous dependent variable>" << std::endl;
	std::cout << "v <vertex-index> <vertex-label>" << std::endl;
	std::cout << "..." << std::endl;
	std::cout << "e <edge-from> <edge-to> <edge-label>" << std::endl;
	std::cout << "..." << std::endl;
	std::cout << "<empty line>" << std::endl;
	std::cout << std::endl;

	std::cout << "Indices start at zero, labels are arbitrary unsigned integers." << std::endl;
	std::cout << std::endl;
}

int main (int argc, char **argv)
{
	unsigned int minsup = 1;
	unsigned int maxpat_min = 0;
	unsigned int maxpat_max = INT_MAX;
	unsigned int nmaxparam = 0xffffffff;
	unsigned int cvmethod = 0;
	float cvparam = 0;
	double minlambda = 0;
	unsigned int cvrepeat = 1;
	bool enc = false;
	bool directed = false;
	bool breadth_first = false;

	int opt;
	while ((opt = getopt(argc, argv, "edwbs:m:L:Dhn:M:p:c:C:n:l:")) != -1) { // : means it require argument. :: means that the argument is optional
		switch(opt) {
		case 's':
			minsup = atoi (optarg);
			break;
		case 'M':
			maxpat_max = atoi (optarg);
			break;
		case 'm':
			maxpat_min = atoi (optarg);
			break;
		case 'd': // same as original gSpan
		case 'e':
			enc = true;
			break;
		case 'D':
			directed = true;
			break;
		case 'h':
			usage();
			return -1;
		case 'p':
			nmaxparam = atoi (optarg);
			break;
		case 'c':
			cvmethod = atoi (optarg);
			break;
		case 'C':
			cvparam = atof (optarg);
			break;
		case 'n':
			cvrepeat = atoi(optarg);
			break;
		case 'l':
			minlambda = atof(optarg);
			break;
		case 'b':
			breadth_first = true;
			break;
		default:
			usage ();
			return -1;
		}
	}
	// Error check
	if (cvmethod > 3) { 
		std::cerr << " Invalid CV method."; 
		return -1;
	}
	// read data
	GSPAN::InputData input_data;
	GSPAN::read(std::cin, directed, input_data);
	/*
	   Initial run to get regression info with full regression set. 
	   This initial run achieves two key factor:
	       1. Build up cache, so in case of parallelized CV run, each CV run will have pre-built cache to be performed at fast speed.
		   2. It records lambdas to be surveyed in CV run. Lambda at each event of this initial run is used as surveying point. 
	*/
	std::vector<bool> initial_run_logical(input_data.size(), false); // No test set
	std::vector<double> * lambdas = new std::vector<double>; // lambdas at each event is recorded here
	
	GSPAN::gSpan gspan;
	gspan.breadth_first = breadth_first;
	gspan.run(input_data, std::cout, minsup, maxpat_min, maxpat_max, enc, nmaxparam, minlambda, directed, NULL, lambdas, NULL, NULL); // lambdas are recorded at each event.
	
	if (cvmethod != 0) {
		// perform CV
		gspan.report = false;
		gspan.cv_run = true; // when true, MSE is surveyed at inputted lambdas. 
		nmaxparam = UINT_MAX;
		CV::cvpartition * CV = NULL;
		std::vector< std::vector<double> > CVMSEs;
		std::vector < std::vector< std::vector<double> > > CVEs; // repeat, lambda, data point
		// initialize cross-validation setting
		if (cvmethod == 1) CV = new CV::cvpartition(input_data.size(), (unsigned)0, cvparam);                                   // hold-out
		else if (cvmethod == 2) CV = new CV::cvpartition(input_data.size(), (unsigned)1, static_cast<unsigned int>(cvparam));   // k-fold
		else if (cvmethod == 3) {                                                                                               // leave-ont-out
			CV = new CV::cvpartition(input_data.size(), (unsigned)2); 
			cvparam = input_data.size();
		}                                       

		//perform CV
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "Cross Validation" << std::endl;

		for (unsigned int i = 0; i < cvrepeat; i++) { // over repeatation
			CV->repartition();
			std::vector < std::vector<double> >  * CVE = new std::vector< std::vector<double> >(lambdas->size(), std::vector<double>(input_data.size(), 0)); // lambda, data point
			std::cout << "Repeat " << i + 1 << ": Test Set";
			for (unsigned int j = 0; j < cvparam; j++) { // over each test set
				std::cout << " " << j + 1;
				std::vector<double>  * CVMSE = new std::vector<double>;

				gspan.run(input_data, std::cout, minsup, maxpat_min, maxpat_max, enc, nmaxparam, minlambda, directed, &CV->at(j).test_logical_index, lambdas, CVMSE, CVE); // supplied lambdas are used to compute CV
				CVMSEs.push_back(*CVMSE);

			}
			CVEs.push_back(*CVE);
			std::cout << std::endl;
		}
		
		std::cout << "Lambda       TestMSE" << std::endl;
		unsigned int maxnlambda=UINT_MAX;
		for (unsigned int j = 0; j < CVMSEs.size(); j++) {
			if (CVMSEs.at(j).size() < maxnlambda) {
				maxnlambda = CVMSEs.at(j).size();
			}
		}
		for (unsigned int i = 0; i < maxnlambda; i++) {
			double AvgMSE = 0;
			for (unsigned int j = 0; j < CVMSEs.size(); j++) {
				AvgMSE += CVMSEs.at(j).at(i);
			}
			AvgMSE = AvgMSE / CVMSEs.size();
			std::cout << lambdas->at(i) << " " << AvgMSE << std::endl;
		}
		std::cout << "Number of descriptor searched: " << gspan.ndfs << std::endl;
		std::cout << "Lambda";
		for (unsigned int j = 0; j < maxnlambda; j++) {
			std::cout << " " << lambdas->at(j);
		}
		std::cout << std::endl;

		for (unsigned int i = 0; i <input_data.size(); i++) {
			std::cout << std::setw(5) << i;
			for (unsigned int j = 0; j < maxnlambda; j++) {
				double AvgE = 0;
				for (unsigned int k = 0; k < cvrepeat; k++) {
					AvgE += CVEs.at(k).at(j).at(i);
				}
				AvgE = AvgE / cvrepeat;
				std::cout << " " << AvgE;
			}
			std::cout << std::endl;
		}
	}
}
