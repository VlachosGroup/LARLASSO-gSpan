/*
$Id: cv.cpp,v 1.0 2017/05/17 Geun Ho Gu, Univeristy of Delaware $;

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
#include "cv.h"
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>

namespace CV {
	//round float to unsigned int
	unsigned int round(float num) {
		return (unsigned int)num + 0.5;
	}

	// pair used to keep track of index.
	typedef std::pair<unsigned int, unsigned int> RIpair;
	bool RIpairComparator(const RIpair& l, const RIpair& r)
	{
		return l.first < r.first;
	}

	// constructor for Hold out
	cvpartition::cvpartition(unsigned int n, unsigned int method, float param) {
		if (method == 2) std::invalid_argument("Leave one out method does not require third parameter");
		else if (method == 1) std::invalid_argument("K Fold method requires unsigned int third parameter");
		if (param<0 || param>1) std::invalid_argument("p must be between 0 and 1");
		// initialize random number generator
		std::srand((unsigned)std::time(0));
		// saving in case of repartitioning
		this->n = n;
		this->method = method;
		this->param = round(n * param); // test set size is saved.
										// initializing vector
		this->reserve(1);
		partition part;
		part.test_set_size = this->param;
		part.training_set_size = n - this->param;
		part.test_set.reserve(part.test_set_size);
		part.training_set.reserve(part.training_set_size);
		part.test_logical_index.resize(n);
		this->push_back(part);
		repartition();
	}

	// constructor for KFold
	cvpartition::cvpartition(unsigned int n, unsigned int method, unsigned int param) {
		if (method == 2) std::invalid_argument("Leave one out method does not require third parameter");
		else if (method == 0) std::invalid_argument("Hold out method requires float third parameter");
		if (param>n) std::invalid_argument("K cannot be bigger than number of data");
		// initialize random number generator
		std::srand((unsigned)std::time(0));
		// saving in case of repartitioning
		this->n = n;
		this->method = method;
		this->param = param; // k is saved
							 // initializing vector
		this->reserve(param);
		unsigned int quotient = n / param;
		unsigned int remainder = n % param;
		for (unsigned int i = 0; i < this->param; i++) {
			partition part;
			if (i < remainder) 	part.test_set_size = quotient + (unsigned int)1;
			else part.test_set_size = quotient;
			part.training_set_size = n - part.test_set_size;
			part.test_set.reserve(part.test_set_size);
			part.training_set.reserve(part.training_set_size);
			part.test_logical_index.resize(n);
			this->push_back(part);
		}
		repartition();
	}



	// constructor for leave out
	cvpartition::cvpartition(unsigned int n, unsigned int method) {
		if (method != 0) std::invalid_argument("KFold and Hold out method requires one more parameter");
		// saving in case of repartitioning
		this->n = n;
		this->method = method;

		// initializing vector
		this->reserve(n);
		for (unsigned int i = 0; i < n; i++) {
			partition part;
			// initialize
			part.test_set_size = 1;
			part.training_set_size = n - 1;
			part.test_set.reserve(1);
			part.training_set.reserve(n - 1);
			part.test_logical_index.resize(n);
			// append test_set and training_set
			for (unsigned int j = 0; j < n; j++) {
				if (i == j) {
					part.test_set.push_back(i);
					part.test_logical_index.at(j) = true;
				}
				else if (i != j) {
					part.training_set.push_back(j);
				}
			}
			this->push_back(part);
		}
	}

	void cvpartition::repartition(void) {
		// leave one out method does not require repartitioning.
		if (method == 2) return;
		// clear vectors
		for (std::vector<CV::partition>::iterator it = this->begin(); it != this->end(); it++) {
			it->test_set.clear();
			it->training_set.clear();
			std::fill(it->test_logical_index.begin(), it->test_logical_index.end(), false);
		}

		// get random permutation
		std::vector<RIpair> randvec;
		randvec.reserve(this->n);

		for (unsigned int i = 0; i < this->n; i++) {
			RIpair pair((unsigned)std::rand(), i);
			randvec.push_back(pair);
		}
		std::sort(randvec.begin(), randvec.end(), RIpairComparator);
		//hold out
		if (method == 0) {

			// first p*100 of the data set is split
			for (unsigned int i = 0; i < this->n; i++) {
				if (i < this->at(0).test_set_size)
				{
					this->at(0).test_set.push_back(randvec[i].second);
					this->at(0).test_logical_index.at(randvec[i].second) = true;
				}
				else {
					this->at(0).training_set.push_back(randvec[i].second);
				}
			}
		}
		//k fold
		else if (method == 1) {
			// vector partitoning randvec into test set. 
			std::vector<unsigned int> PartIndex;
			PartIndex.reserve(this->n);
			for (unsigned int i = 0; i < this->n; i++) {
				PartIndex.push_back(i % this->param);
			}
			// partitioning
			for (unsigned int i = 0; i < this->size(); i++) {
				for (unsigned int j = 0; j < this->n; j++) {
					if (PartIndex[j] == i)
					{
						this->at(i).test_set.push_back(randvec[j].second);
						this->at(i).test_logical_index.at(randvec[j].second) = true;
					}
					else if (PartIndex[j] != i) this->at(i).training_set.push_back(randvec[j].second);
				}
			}
		}

	}


}

