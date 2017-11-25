#pragma once
/*
$Id: cv.h,v 1.0 2017/05/17 Geun Ho Gu, Univeristy of Delaware $;

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
#include <vector>

namespace CV {

	class partition {
	public:
		std::vector<unsigned int> test_set;
		std::vector<unsigned int> training_set;
		std::vector<bool> test_logical_index;
		unsigned int test_set_size;
		unsigned int training_set_size;
		partition() {};
	};

	// each element in vector represents a set of test and training set.
	class cvpartition : public std::vector<partition> {

	public:
		/* List of constructor:
		* cvpartition(unsigned int n, unsigned int 0, float param)
		*		hold out method, where param is the ratio between test and training set
		*		n is equal to the number of data
		* cvpartition(unsigned int n, unsigned int 1, unsigned int param)
		*		KFold method, where param is the number of fold (k).
		* cvpartition(unsigned int n, unsigned int 2)
		*		Leave one out method
		*/
		cvpartition(unsigned int n, unsigned int method);
		cvpartition(unsigned int n, unsigned int method, unsigned int param);
		cvpartition(unsigned int n, unsigned int method, float param);
		unsigned int n;
		unsigned int method;
		unsigned int param;

		// randomization function
		void repartition(void);
	};
};


