/*
	$Id: gspan.cpp,v 1.8 2004/05/21 09:27:17 taku-ku Exp $;

   Copyright (C) 2004 Taku Kudo, All rights reserved.
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
#include <iterator>
#include <stdlib.h>
#include <ostream>
#include <fstream>
#include <iostream>
#include <strstream>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <float.h> // MAX
#include <iomanip> //setw
#include <deque> //iterative method

namespace GSPAN {
gSpan::gSpan(void) {
	initial_graph_search = true;
	report = true;
	cv_run = false; // when true, MSE is surveyed at gSpan::lambda
	single_node_enumerated = false;
	breadth_first = false;
	/*Things done to save memory:
	  1. cache no-longer saves graph. graph is generated on the fly when the graph is selected.*/
	
}

void gSpan::reset(void) {
	// clear up vectors
	beta.clear();
	gamma.clear();
	test_logical_index.clear();
	sign.clear();
	
	w.clear();
	v.clear();
	// initial graph search enabled.
	initial_graph_search = true;
	// reset caches.
	for (std::vector<DFSCodeCache *>::iterator it = linearly_dependent_caches.begin(); it != linearly_dependent_caches.end(); ++it) {
		(*it)->linearlydependent = false;
	}
	for (std::vector<DFSCodeCache *>::iterator it = X_cache_pointers.begin(); it != X_cache_pointers.end(); ++it) {
		(*it)->ActiveSet = false;
	}
	linearly_dependent_caches.clear();
	X_cache_pointers.clear();
}

std::vector<double> gSpan::compute_residual(std::vector<double> Beta, bool TestOrTrain) {
	// initialize residual vector
	std::vector<double> residual;
	if (TestOrTrain) residual.reserve(test_set_size);
	else residual.reserve(train_set_size);
	// compute residual
	for (unsigned int i = 0; i < Data_Set.size(); i++) {
		// check whether the index is test or train set
		if ((TestOrTrain && test_logical_index.at(i)) || (!TestOrTrain && !test_logical_index.at(i))) {
			// compute predicted y
			double ypi = 0;
			for (unsigned int j = 0; j < Beta.size(); j++) {
				ypi += Beta.at(j)*X_cache_pointers[j]->x.at(i);
			}
			// subtract from actual y and append
			residual.push_back(y.at(i) - ypi);
		}
	}
	return residual;
}

double gSpan::compute_MSE_from_residual(std::vector<double> residual) {
	// sum up squared error
	double MSE = 0;
	for (unsigned int i = 0; i < residual.size(); i++) {
		MSE += residual.at(i)*residual.at(i);
	}
	// take average
	MSE = MSE / residual.size();
	return MSE;
}


std::vector<double> gSpan::complete_regression (void) {
	std::vector<double> TempBeta;
	double d = rho0 / eta0;
	for (unsigned int i = 0; i < beta.size(); i++) {
		TempBeta.push_back(beta[i] + d*gamma[i]);
	}
	return TempBeta;
}

double gSpan::compute_covariance_r_xtrain(FrequencyVector x) {
	std::vector<double>::iterator wIT = w.begin();
	double covariance = 0;
	for (unsigned int i = 0; i < Data_Set.size(); i++) {
		// check whether the index is test or train set
		if (!test_logical_index.at(i)) {
			covariance += (*wIT)*x.at(i);
			wIT++;
		}
	}
	return covariance;
}

void gSpan::compute_v(void) {
	v.clear();
	//g^t*x
	for (unsigned int i = 0; i < Data_Set.size(); i++) {
		// check whether the index is test or train set
		if (!test_logical_index.at(i)){
			// compute predicted y
			double gtxi = 0;
			for (unsigned int j = 0; j < gamma.size(); j++) {
				gtxi += gamma.at(j)*X_cache_pointers[j]->x.at(i);
			}
			v.push_back(gtxi);
		}
	}
}

double gSpan::compute_covariance_g_xtrain(FrequencyVector x) {
	std::vector<double>::iterator vIT = v.begin();
	double covariance = 0;
	for (unsigned int i = 0; i < Data_Set.size(); i++) {
		// check whether the index is test or train set
		if (!test_logical_index.at(i)) {
			covariance += (*vIT)*x.at(i);
			vIT++;
		}
	}
	return covariance;
}

double gSpan::compute_dt(FrequencyVector x) {
	// can this be optimized?
	// compute dt
	double wx = compute_covariance_r_xtrain(x);
	double vx = compute_covariance_g_xtrain(x);
	double dt1 = (rho0 - wx) / (eta0 - vx);
	double dt2 = (rho0 + wx) / (eta0 + vx);
	// find positive minimum
	double dt;
	if (dt1 > 0 && dt1 <= dt2) // a is positive and smaller than or equal to b
		dt = dt1;
	else if (dt2 > 0) // b is positive and either smaller than a or a is negative
		dt = dt2;
	else //both values are not positive values (could be negative or zero)
		dt = DBL_MAX;
	return dt;
}

unsigned int gSpan::compute_d2(void) {
	double dt;
	unsigned int d2index = UINT_MAX;
	for (unsigned int i = 0; i < beta.size(); i++) {
		dt = -beta.at(i) / gamma.at(i);
		if (dt > 0 && d2 > dt) {
			d2 = dt;
			d2index = i; // saving the index in case we need to remove it. 
		}
	}
	return d2index;
}

// sign function
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}



unsigned int gSpan::support (std::vector<PDFS * > PDFSs)
{
	unsigned int oid = 0xffffffff;
	unsigned int size = 0;

	for (std::vector<PDFS * >::iterator cur = PDFSs.begin(); cur != PDFSs.end(); ++cur) { // the total number of PDFS is returned. but why is this related to support?
		if (oid != (*cur)->graph->id) {
			++size;
		}
		oid = (*cur)->graph->id;
	}

	return size;
}

double gSpan::compute_gain(FrequencyVector x) {
	double gain_t = 0;
	for (unsigned int i = 0; i < Data_Set.size(); i++) {
		// check whether the index is test or train set
		if (!test_logical_index.at(i)) gain_t += y.at(i)*x.at(i);
	}
	return fabs(gain_t);
}

bool gSpan::gain_prune_condition(FrequencyVector x) {
	double gain_p1 = 0;
	double gain_p2 = 0;
	for (unsigned int i = 0; i < Data_Set.size(); i++) {
		if (!test_logical_index.at(i)) {
			if (y.at(i) < 0) {
				gain_p1 += fabs(y.at(i))*x.at(i);
			}
			else if (y.at(i) > 0) {
				gain_p2 += (y.at(i))*x.at(i);
			}
		}
	}
	return largest_gain > std::max(gain_p1, gain_p2);

}

bool gSpan::d1_prune_condition(FrequencyVector x) {
	// solving pruning condition equation.
	double bw1 = 0;
	double bw2 = 0;
	std::vector<double>::iterator wIT = w.begin();
	for (unsigned int i = 0; i < Data_Set.size(); i++) {
		if (!test_logical_index.at(i)) {
			if (*wIT < 0) {
				bw1 += fabs(*wIT)*x.at(i);
			}
			else if (*wIT > 0) {
				bw2 += *wIT*x.at(i);
			}
			wIT++;
		}
	}
	double bv1 = 0;
	double bv2 = 0;
	std::vector<double>::iterator vIT = v.begin();
	for (unsigned int i = 0; i < Data_Set.size(); i++) {
		if (!test_logical_index.at(i)) {
			if (*vIT < 0) {
				bv1 += fabs(*vIT)*x.at(i);
			}
			else if (*vIT > 0) {
				bv2 += *vIT*x.at(i);
			}
			vIT++;
		}
	}
	return (std::max(bw1, bw2) + d1*std::max(bv1, bv2)) < (fabs(rho0) - d1*fabs(eta0)); //RHS can be calculated when new d1 is calculated
}


Eigen::MatrixXd gSpan::AppendXtoGram(Eigen::MatrixXd TempGram, FrequencyVector x) {
	// must be executed before appending selected_cache to X_cache_pointers. 
	// resizing
	TempGram.conservativeResize(X_cache_pointers.size() + 1, X_cache_pointers.size() + 1);
	// off-diagonal 
	for (unsigned int i = 0; i < X_cache_pointers.size(); i++) {
		double num = 0;
		for (unsigned int k = 0; k < Data_Set.size(); k++) {
			if (!test_logical_index.at(k)) num += X_cache_pointers[i]->x.at(k) * x.at(k);
		}
		TempGram(X_cache_pointers.size(), i) = num;
		TempGram(i, X_cache_pointers.size()) = num;
	}
	//diagonal
	double num = 0;
	for (unsigned int k = 0; k < Data_Set.size(); k++) {
		if (!test_logical_index.at(k)) num += x.at(k) * x.at(k);
	}
	TempGram(X_cache_pointers.size(), X_cache_pointers.size()) = num;
	return TempGram;
}

Eigen::MatrixXd gSpan::RemoveXfromGram(Eigen::MatrixXd Gram, unsigned int dimensionToRemove) {
	unsigned int Dim = Gram.rows() - 1;

	if (dimensionToRemove < Dim) {
		Gram.block(dimensionToRemove, 0, Dim - dimensionToRemove, Dim + 1) = Gram.block(dimensionToRemove + 1, 0, Dim - dimensionToRemove, Dim + 1);
		Gram.block(0, dimensionToRemove, Dim , Dim - dimensionToRemove) = Gram.block(0, dimensionToRemove + 1, Dim, Dim - dimensionToRemove);
	}

	Gram.conservativeResize(Dim, Dim);
	return Gram;
}

////////////////////////////////////// new vector < unsigned int > -> sort -> find method
// This is the slowest part of the algorithm for the large graph data.
// Tried several thing to optimize:
//   1. Tried edge pointer for isomorphism. didn't work well.
//   2. Tried unordered_set for unsigned int container which is the fastest container for find to avoid using sort. Not the fastest
//   3. Sort -> find is the fastest.
// here are some related discussion
// https://stackoverflow.com/questions/26431593/compare-two-integer-arrays-and-check-if-they-have-equal-values
// https://stackoverflow.com/questions/17394149/how-to-efficiently-compare-vectors-with-c
// https://stackoverflow.com/questions/6985572/which-is-the-fastest-stl-container-for-find
FrequencyVector gSpan::GetFrequencyVector(std::vector<PDFS * > PDFSs) {
	// initialize zero vector
	FrequencyVector x(Data_Set.size(), 0);

	unsigned int oid = 0xffffffff;
	// get size for reserving 
	unsigned int n_edge = 0;
	for (PDFS *p = (*PDFSs.begin()); p; p = p->prev) n_edge++;
	std::vector <std::vector<unsigned int> > UniquePDFS; // uniquePDFS given a data set.
													//std::unique(SortedPDFS.begin(), SortedPDFS.begin(),)
													//std::cout << projected.size() << std::endl;
													//int i = 0;
	// edge id pointer
	std::vector<unsigned int> * edge_ids = new std::vector<unsigned int>;
	edge_ids->reserve(n_edge);
	for (std::vector<PDFS *>::iterator cur = PDFSs.begin(); cur != PDFSs.end(); ++cur) {
		//i += 1;
		//std::cout << i << std::endl;
		if (oid != (*cur)->graph->id && cur != PDFSs.begin()) { // or the next iteration involves different data point (or graph)

			x.at(oid) = UniquePDFS.size();
			UniquePDFS.clear();
		}
		
		// push ids
		for (PDFS *p = (*cur); p; p = p->prev) edge_ids->push_back(p->edge->id);
		// sort.
		std::sort(edge_ids->begin(), edge_ids->end());
		// compare with existing unique vector of edge_pts.
		if (std::find(UniquePDFS.begin(), UniquePDFS.end(), *edge_ids) != UniquePDFS.end()) 
			//isomorphic graph.
			edge_ids->clear(); 
		else {
			//non-isomorphic graph. add to unique list.
			UniquePDFS.push_back(*edge_ids); 
			delete edge_ids;
			edge_ids = new std::vector<unsigned int>;
		}

		// If iterator reaches the end or find next graph, it assigns x
		oid = (*cur)->graph->id;
	}
	delete edge_ids;
	// for the last graph repeat
	x.at(oid) = UniquePDFS.size();
	return x;
}

bool gSpan::DFS_check_condition(std::vector<PDFS*> PDFSs) {
	// Check if the pattern is frequent enough. Support prunning is currently done for test set and training set combined. 
	unsigned int sup = support(PDFSs);
	if (sup < minsup) 
		return false;

	// The minimal DFS code check is more expensive than the support check,
	// hence it is done now, after checking the support.
	if (is_min() == false) 
		return false;

	// In case we have a valid upper bound and our graph already exceeds it,
	// return.  Note: we do not check for equality as the DFS exploration may
	// still add edges within an existing subgraph, without increasing the
	// number of nodes.
	if (DFS_CODE.nodeCount() > maxpat_max) 
		return false;

	return true;
}

// This code below is an example of recursive function version.
void gSpan::project_depth_first(Projected * const &projected)
{	
	DFSCodeCache *DFS_Code_Cache = projected->DFS_Code_Cache;
	/* GGDEBUG
	if (DFS_Code_Cache->g == debug_graph) {
		std::cout << "GG";
	}
	GGDEBUG*/
	// compute gain or dt and prune
	if (!DFS_Code_Cache->maxpatmincheck) {
		if (initial_graph_search) { // initial descriptor search
									// Line 3 of algorithm 3 
			double gain_t = compute_gain(DFS_Code_Cache->x);
			// Line 10,11,12.
			gain_t = fabs(gain_t);
			if (gain_t > largest_gain) {
				largest_gain = gain_t;
				selected_cache = DFS_Code_Cache;
			}
			// Line 13,14,15
			if (gain_prune_condition(DFS_Code_Cache->x_upper)) return;

		}
		else { // regular run
			double dt = compute_dt(DFS_Code_Cache->x); // Line 9 of algorithm 2
			if (!DFS_Code_Cache->ActiveSet && !DFS_Code_Cache->linearlydependent && dt < d1 && dt != 0) {
				d1 = dt; // Line 11
				selected_cache = DFS_Code_Cache;
			}
			if (d1_prune_condition(DFS_Code_Cache->x_upper)) return;
		}
	}
	// enumerate next extension if it has not been done.
	if (!DFS_Code_Cache->next_extension_enumerated) {
		// Here we enumerate supergraphs as it has passed all the test.
		const RMPath &rmpath = DFS_CODE.buildRMPath();
		int minlabel = DFS_CODE[0].fromlabel;
		DFS_Code_Cache->maxtoc = DFS_CODE[rmpath[0]].to;
		EdgeList edges;

		// Enumerate all possible one edge extensions of the current substructure.
		for (unsigned int n = 0; n < projected->PDFSs.size(); ++n) {

			Graph *graph = projected->PDFSs.at(n)->graph;
			History history(*graph, projected->PDFSs.at(n));

			// backward
			for (int i = (int)rmpath.size() - 1; i >= 1; --i) {
				Edge *e = get_backward(*graph, history[rmpath[i]], history[rmpath[0]], history);
				if (e)
					DFS_Code_Cache->new_bck_DFS[DFS_CODE[rmpath[i]].from][e->elabel].push(graph, e, projected->PDFSs.at(n));
			}
			// pure forward
			if (get_forward_pure(*graph, history[rmpath[0]], minlabel, history, edges))
				for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
					DFS_Code_Cache->new_fwd_DFS[DFS_Code_Cache->maxtoc][(*it)->elabel][(*graph)[(*it)->to].label].push(graph, *it, projected->PDFSs.at(n));

			// backtracked forward
			for (int i = 0; i < (int)rmpath.size(); ++i)
				if (get_forward_rmpath(*graph, history[rmpath[i]], minlabel, history, edges))
					for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
						DFS_Code_Cache->new_fwd_DFS[DFS_CODE[rmpath[i]].from][(*it)->elabel][(*graph)[(*it)->to].label].push(graph, *it, projected->PDFSs.at(n));
		}

		// create caches for next extensions. Backward
		for (Projected_iterator2 to = DFS_Code_Cache->new_bck_DFS.begin();
			to != DFS_Code_Cache->new_bck_DFS.end(); )// from each root
		{
			for (Projected_iterator1 elabel = to->second.begin();
				elabel != to->second.end(); )
			{
				DFS_CODE.push(DFS_Code_Cache->maxtoc, to->first, -1, elabel->first, -1);
				if (DFS_check_condition(elabel->second.PDFSs)) {
					elabel->second.DFS_Code_Cache = new DFSCodeCache; // create cache
					elabel->second.DFS_Code_Cache->Parent_DFS_Code_Cache = DFS_Code_Cache;
					// record graph
					elabel->second.DFS_Code_Cache->g = GRAPH_IS_MIN; // graph is built when is_min() is called.
					// minimum check does not stop algorithm as DFS code needs to be continued to be expanded before its size can exceed minimum.
					elabel->second.DFS_Code_Cache->maxpatmincheck = DFS_CODE.nodeCount() < maxpat_min;
					// Get frequency vector
					elabel->second.DFS_Code_Cache->x = GetFrequencyVector(elabel->second.PDFSs);
					elabel->second.DFS_Code_Cache->x_upper = elabel->second.DFS_Code_Cache->x;
					// Update upper bound of the parent DFS_Code
					for (DFSCodeCache* p = DFS_Code_Cache; p; p = p->Parent_DFS_Code_Cache) {
						for (unsigned int i = 0; i < elabel->second.DFS_Code_Cache->x_upper.size(); ++i) {
							if (p->x_upper.at(i) < elabel->second.DFS_Code_Cache->x.at(i)) {
								p->x_upper.at(i) = elabel->second.DFS_Code_Cache->x.at(i);
							}
						}
					}
					++elabel;
					ndfs += 1;
				}
				else {
					// erase PDFSs
					for (std::vector< PDFS* >::iterator it = elabel->second.PDFSs.begin(); it != elabel->second.PDFSs.end(); ++it) delete (*it);
					elabel->second.PDFSs.clear();
					std::vector<PDFS*>().swap(elabel->second.PDFSs);
					// delete map key
					to->second.erase(elabel++);
				}
				DFS_CODE.pop();
			}
			if (DFS_Code_Cache->new_bck_DFS.size() != 0) ++to;
			else DFS_Code_Cache->new_bck_DFS.erase(to++);
		}

		// create caches for next extensions. Forward
		for (Projected_iterator3 fromlabel = DFS_Code_Cache->new_fwd_DFS.begin();
			fromlabel != DFS_Code_Cache->new_fwd_DFS.end();)// from each root
		{
			for (Projected_iterator2 elabel = fromlabel->second.begin(); // given std::map<X, Y> where X is key and Y is value, std::map<X, Y>->second will give Y
				elabel != fromlabel->second.end();)
			{
				for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end();)
				{

					if (DFS_CODE.size() == 0) DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);
					else DFS_CODE.push(fromlabel->first, DFS_Code_Cache->maxtoc + 1, -1, elabel->first, tolabel->first);

					if (DFS_check_condition(tolabel->second.PDFSs)) {
						tolabel->second.DFS_Code_Cache = new DFSCodeCache; // create cache
						tolabel->second.DFS_Code_Cache->Parent_DFS_Code_Cache = DFS_Code_Cache;
						// record graph
						tolabel->second.DFS_Code_Cache->g = GRAPH_IS_MIN; // graph is built when is_min() is called.
						// minimum check does not stop algorithm as DFS code needs to be continued to be expanded before its size can exceed minimum.
						tolabel->second.DFS_Code_Cache->maxpatmincheck = DFS_CODE.nodeCount() < maxpat_min;
						// Get frequency vector
						tolabel->second.DFS_Code_Cache->x = GetFrequencyVector(tolabel->second.PDFSs);
						tolabel->second.DFS_Code_Cache->x_upper = tolabel->second.DFS_Code_Cache->x;
						// Update upper bound of the parent DFS_Code
						for (DFSCodeCache* p = DFS_Code_Cache; p; p = p->Parent_DFS_Code_Cache) {
							for (unsigned int i = 0; i < tolabel->second.DFS_Code_Cache->x_upper.size(); ++i) {
								if (p->x_upper.at(i) < tolabel->second.DFS_Code_Cache->x.at(i)) {
									p->x_upper.at(i) = tolabel->second.DFS_Code_Cache->x.at(i);
								}
							}
						}
						++tolabel;
						ndfs += 1;
					}
					else {
						// erase PDFSs
						for (std::vector< PDFS* >::iterator it = tolabel->second.PDFSs.begin(); it != tolabel->second.PDFSs.end(); ++it) delete (*it);
						tolabel->second.PDFSs.clear();
						std::vector<PDFS*>().swap(tolabel->second.PDFSs);
						// delete map key
						elabel->second.erase(tolabel++);
					}
					DFS_CODE.pop();
				}
				if (fromlabel->second.size() != 0) 	++elabel; 
				else fromlabel->second.erase(elabel++);
			}
			if (DFS_Code_Cache->new_fwd_DFS.size() != 0) ++fromlabel;
			else DFS_Code_Cache->new_fwd_DFS.erase(fromlabel++);
		}
		// clear to save memory.
		projected->PDFSs.clear();
		std::vector<PDFS*>().swap(projected->PDFSs);
		DFS_Code_Cache->next_extension_enumerated = true;
	}

	// Test all extended substructures.
	// backward
	for (Projected_iterator2 to = DFS_Code_Cache->new_bck_DFS.begin(); to != DFS_Code_Cache->new_bck_DFS.end(); ++to) {
		for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
			DFS_CODE.push(DFS_Code_Cache->maxtoc, to->first, -1, elabel->first, -1);
			project_depth_first(&(elabel->second));
			DFS_CODE.pop();
		}
	}
	// forward
	for (Projected_riterator3 from = DFS_Code_Cache->new_fwd_DFS.rbegin();
		from != DFS_Code_Cache->new_fwd_DFS.rend(); ++from)
	{
		for (Projected_iterator2 elabel = from->second.begin();
			elabel != from->second.end(); ++elabel)
		{
			for (Projected_iterator1 tolabel = elabel->second.begin();
				tolabel != elabel->second.end(); ++tolabel)
			{
				if (DFS_CODE.size() == 0) DFS_CODE.push(0, 1, from->first, elabel->first, tolabel->first);
				else DFS_CODE.push(from->first, DFS_Code_Cache->maxtoc + 1, -1, elabel->first, tolabel->first);
				project_depth_first(&(tolabel->second));
				DFS_CODE.pop();
			}
		}
	}
	
	return;
}

void gSpan::project_breadth_first(void)
{
	// stack for the iterative for-loop
	std::deque<Projected *> Stack;
	// initial stack
	for (std::map<int, Projected>::iterator SingleNodeIt = Single_Node_DFS_Code_Cache_Map.begin();
		SingleNodeIt != Single_Node_DFS_Code_Cache_Map.end(); ++SingleNodeIt) Stack.push_back(&(SingleNodeIt->second));
	// enter iterative loop
	while (!Stack.empty()) {
		// Pick the very first element from the stack.
		Projected *projected = Stack.front();
		Stack.pop_front();
		DFSCodeCache *DFS_Code_Cache = projected->DFS_Code_Cache;
		// Get DFS Code as well
		DFS_CODE = DFS_Code_Cache->DFS_Code;

		// compute gain or dt and prune
		if (!DFS_Code_Cache->maxpatmincheck) {
			if (initial_graph_search) { // initial descriptor search
										// Line 3 of algorithm 3 
				double gain_t = compute_gain(DFS_Code_Cache->x);
				// Line 10,11,12.
				gain_t = fabs(gain_t);
				if (gain_t > largest_gain) {
					largest_gain = gain_t;
					selected_cache = DFS_Code_Cache;
				}
				// Line 13,14,15
				if (gain_prune_condition(DFS_Code_Cache->x)) continue;

			}
			else { // regular run
				double dt = compute_dt(DFS_Code_Cache->x); // Line 9 of algorithm 2
				if (!DFS_Code_Cache->ActiveSet && !DFS_Code_Cache->linearlydependent && dt < d1 && dt != 0) {
					d1 = dt; // Line 11
					selected_cache = DFS_Code_Cache;
				}
				if (d1_prune_condition(DFS_Code_Cache->x)) 	continue;
			}
		}
		// enumerate next extension if it has not been done.
		if (!DFS_Code_Cache->next_extension_enumerated) {
			// Here we enumerate supergraphs as it has passed all the test.
			const RMPath &rmpath = DFS_CODE.buildRMPath();
			int minlabel = DFS_CODE[0].fromlabel;
			DFS_Code_Cache->maxtoc = DFS_CODE[rmpath[0]].to;
			EdgeList edges;

			// Enumerate all possible one edge extensions of the current substructure.
			for (unsigned int n = 0; n < projected->PDFSs.size(); ++n) {

				Graph *graph = projected->PDFSs.at(n)->graph;
				History history(*graph, projected->PDFSs.at(n));

				// XXX: do we have to change something here for directed edges?

				// backward
				for (int i = (int)rmpath.size() - 1; i >= 1; --i) {
					Edge *e = get_backward(*graph, history[rmpath[i]], history[rmpath[0]], history);
					if (e)
						DFS_Code_Cache->new_bck_DFS[DFS_CODE[rmpath[i]].from][e->elabel].push(graph, e, projected->PDFSs.at(n));
				}

				// pure forward
				// FIXME: here we pass a too large e->to (== history[rmpath[0]]->to
				// into get_forward_pure, such that the assertion fails.
				//
				// The problem is:
				// history[rmpath[0]]->to > Data_Set[id].size()
				if (get_forward_pure(*graph, history[rmpath[0]], minlabel, history, edges))
					for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
						DFS_Code_Cache->new_fwd_DFS[DFS_Code_Cache->maxtoc][(*it)->elabel][(*graph)[(*it)->to].label].push(graph, *it, projected->PDFSs.at(n));

				// backtracked forward
				for (int i = 0; i < (int)rmpath.size(); ++i)
					if (get_forward_rmpath(*graph, history[rmpath[i]], minlabel, history, edges))
						for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
							DFS_Code_Cache->new_fwd_DFS[DFS_CODE[rmpath[i]].from][(*it)->elabel][(*graph)[(*it)->to].label].push(graph, *it, projected->PDFSs.at(n));
			}

			// create caches for next extensions. Backward
			for (Projected_iterator2 to = DFS_Code_Cache->new_bck_DFS.begin();
				to != DFS_Code_Cache->new_bck_DFS.end(); )// from each root
			{
				for (Projected_iterator1 elabel = to->second.begin();
					elabel != to->second.end(); )
				{
					DFS_CODE.push(DFS_Code_Cache->maxtoc, to->first, -1, elabel->first, -1);
					if (DFS_check_condition(elabel->second.PDFSs)) {
						elabel->second.DFS_Code_Cache = new DFSCodeCache; // create cache
						elabel->second.DFS_Code_Cache->Parent_DFS_Code_Cache = DFS_Code_Cache;
						// record graph
						elabel->second.DFS_Code_Cache->g = GRAPH_IS_MIN; // graph is built when is_min() is called.
																		 // minimum check does not stop algorithm as DFS code needs to be continued to be expanded before its size can exceed minimum.
						elabel->second.DFS_Code_Cache->DFS_Code = DFS_CODE;
						elabel->second.DFS_Code_Cache->maxpatmincheck = DFS_CODE.nodeCount() < maxpat_min;
						// Get frequency vector
						elabel->second.DFS_Code_Cache->x = GetFrequencyVector(elabel->second.PDFSs);
						elabel->second.DFS_Code_Cache->x_upper = elabel->second.DFS_Code_Cache->x;
						// Update upper bound of the parent DFS_Code
						for (DFSCodeCache* p = DFS_Code_Cache; p; p = p->Parent_DFS_Code_Cache) {
							for (unsigned int i = 0; i < elabel->second.DFS_Code_Cache->x_upper.size(); ++i) {
								if (p->x_upper.at(i) < elabel->second.DFS_Code_Cache->x.at(i)) {
									p->x_upper.at(i) = elabel->second.DFS_Code_Cache->x.at(i);
								}
							}
						}
						++elabel;
						ndfs += 1;
					}
					else {
						// erase PDFSs
						for (std::vector< PDFS* >::iterator it = elabel->second.PDFSs.begin(); it != elabel->second.PDFSs.end(); ++it) delete (*it);
						elabel->second.PDFSs.clear();
						std::vector<PDFS*>().swap(elabel->second.PDFSs);
						// delete map key
						to->second.erase(elabel++);
					}
					DFS_CODE.pop();
				}
				if (DFS_Code_Cache->new_bck_DFS.size() != 0) ++to;
				else DFS_Code_Cache->new_bck_DFS.erase(to++);
			}

			// create caches for next extensions. Forward
			for (Projected_iterator3 fromlabel = DFS_Code_Cache->new_fwd_DFS.begin();
				fromlabel != DFS_Code_Cache->new_fwd_DFS.end();)// from each root
			{
				for (Projected_iterator2 elabel = fromlabel->second.begin(); // given std::map<X, Y> where X is key and Y is value, std::map<X, Y>->second will give Y
					elabel != fromlabel->second.end();)
				{
					for (Projected_iterator1 tolabel = elabel->second.begin();
						tolabel != elabel->second.end();)
					{

						if (DFS_CODE.size() == 0) DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);
						else DFS_CODE.push(fromlabel->first, DFS_Code_Cache->maxtoc + 1, -1, elabel->first, tolabel->first);

						if (DFS_check_condition(tolabel->second.PDFSs)) {
							tolabel->second.DFS_Code_Cache = new DFSCodeCache; // create cache
							tolabel->second.DFS_Code_Cache->Parent_DFS_Code_Cache = DFS_Code_Cache;
							// record graph
							tolabel->second.DFS_Code_Cache->g = GRAPH_IS_MIN; // graph is built when is_min() is called.
																			  // minimum check does not stop algorithm as DFS code needs to be continued to be expanded before its size can exceed minimum.
							tolabel->second.DFS_Code_Cache->DFS_Code = DFS_CODE;
							tolabel->second.DFS_Code_Cache->maxpatmincheck = DFS_CODE.nodeCount() < maxpat_min;
							// Get frequency vector
							tolabel->second.DFS_Code_Cache->x = GetFrequencyVector(tolabel->second.PDFSs);
							tolabel->second.DFS_Code_Cache->x_upper = tolabel->second.DFS_Code_Cache->x;
							// Update upper bound of the parent DFS_Code
							for (DFSCodeCache* p = DFS_Code_Cache; p; p = p->Parent_DFS_Code_Cache) {
								for (unsigned int i = 0; i < tolabel->second.DFS_Code_Cache->x_upper.size(); ++i) {
									if (p->x_upper.at(i) < tolabel->second.DFS_Code_Cache->x.at(i)) {
										p->x_upper.at(i) = tolabel->second.DFS_Code_Cache->x.at(i);
									}
								}
							}
							++tolabel;
							ndfs += 1;
						}
						else {
							// erase PDFSs
							for (std::vector< PDFS* >::iterator it = tolabel->second.PDFSs.begin(); it != tolabel->second.PDFSs.end(); ++it) delete (*it);
							tolabel->second.PDFSs.clear();
							std::vector<PDFS*>().swap(tolabel->second.PDFSs);
							// delete map key
							elabel->second.erase(tolabel++);
						}
						DFS_CODE.pop();
					}
					if (fromlabel->second.size() != 0) 	++elabel;
					else fromlabel->second.erase(elabel++);
				}
				if (DFS_Code_Cache->new_fwd_DFS.size() != 0) ++fromlabel;
				else DFS_Code_Cache->new_fwd_DFS.erase(fromlabel++);
			}
			// clear to save memory.
			projected->PDFSs.clear();
			std::vector<PDFS*>().swap(projected->PDFSs);
			DFS_Code_Cache->next_extension_enumerated = true;
		}

		// Test all extended substructures.
		// backward
		for (Projected_iterator2 to = DFS_Code_Cache->new_bck_DFS.begin(); to != DFS_Code_Cache->new_bck_DFS.end(); ++to) {
			for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
				Stack.push_back(&(elabel->second));
			}
		}
		// forward
		for (Projected_riterator3 from = DFS_Code_Cache->new_fwd_DFS.rbegin();
			from != DFS_Code_Cache->new_fwd_DFS.rend(); ++from)
		{
			for (Projected_iterator2 elabel = from->second.begin();
				elabel != from->second.end(); ++elabel)
			{
				for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end(); ++tolabel)
				{
					Stack.push_back(&(tolabel->second));
				}
			}
		}
	}
}


void gSpan::run(InputData _input_data, std::ostream &_os,
	unsigned int _minsup, unsigned int _maxpat_min, unsigned int _maxpat_max,
	bool _enc, unsigned int _nmaxparam, double _minlambda, bool _directed,
	std::vector<bool> * _test_logical_index, std::vector<double> * lambdas, std::vector<double> * CVMSE, std::vector < std::vector<double> >  * CVE)
{	/************************* Begin Initialization *************************/
	// Error check
	if (_maxpat_max < _maxpat_min) {
		std::cerr << " Maximum number of nodes cannot be lower than minimum number of nodes.";
		return;
	}
	reset();
	// Reading Input
	os = &_os;					// output stream
	ID = 0;						// mined subgraph index for report function
	minsup = _minsup;			// minimum support
	maxpat_min = _maxpat_min;	// minimum number of nodes in mined pattern
	maxpat_max = _maxpat_max;	// maximum number of nodes in mined pattern
	enc = _enc;					// print out selected descriptor
	directed = _directed;		// edge direction specification
	unsigned int nmaxparam;		// maximum number of pattern mined
	if (_nmaxparam == 0xffffffff) nmaxparam = _input_data.size()-1;
	else nmaxparam = _nmaxparam;
	double minlambda = _minlambda;
	// initialize
	Eigen::MatrixXd Gram;
	std::cout << std::scientific;
	std::vector<double>::iterator lambdasIT = lambdas->begin();   // this is only used for CV run. 
	unsigned int nevent = 1;
	ilambda = 0;
	std::ofstream XmatFIO;
	XmatFIO.open("Xmat.txt");
	// Reading Data
	Data_Set = _input_data;
	y.reserve(Data_Set.size());
	for (std::vector < Graph >::iterator cur = Data_Set.begin(); cur != Data_Set.end(); ++cur) {
		y.push_back(cur->y);
	}

	// Set up response vector y.
	if (cv_run) {
		test_logical_index = *_test_logical_index;
		test_set_size = std::count(test_logical_index.begin(), test_logical_index.end(), true);
		train_set_size = std::count(test_logical_index.begin(), test_logical_index.end(), false);
	}
	else {
		test_logical_index = std::vector<bool>(Data_Set.size(), false);
		test_set_size = 0;
		train_set_size = Data_Set.size();
	}

	/************************** End Initialization **************************/
	/********************** Begin Pick First Descriptor **********************/
	if (report) std::cout << "   #: Event         (N graph): Lambda       TrainMSE     TestMSE" << std::endl;
	// line 1  in Algorithm 1
	// set up beta
	beta.push_back(0.0);

	// initial descriptor search
	largest_gain = 0;

	// enumerate single node DFS code
	if (!single_node_enumerated) Enumerate_single_node_DFScode();
	single_node_enumerated = true;

	// DFS search through each vertex, and pick a optimal graph

	if (breadth_first) {
		project_breadth_first();
	}
	else {
		for (std::map<int, Projected>::iterator SingleNodeIt = Single_Node_DFS_Code_Cache_Map.begin();
			SingleNodeIt != Single_Node_DFS_Code_Cache_Map.end(); ++SingleNodeIt) project_depth_first(&(SingleNodeIt->second));
	}
	
	// Append first graph
	selected_cache->ActiveSet = true;
	X_cache_pointers.push_back(selected_cache);
	// compute train set MSE
	std::vector < double > residual = compute_residual(beta, false);
	double trainMSE = compute_MSE_from_residual(residual);
	// compute test set MSE
	residual = compute_residual(beta, true);
	double testMSE = compute_MSE_from_residual(residual);
	Lambda = largest_gain / train_set_size;
	
	// record lambdas (non-CV run), or record MSE (CV run)
	if (!cv_run) {
		lambdas->push_back(Lambda);
		XmatFIO << "Event " << nevent << "; Lambda = " << Lambda << std::endl;
		for (unsigned int i = 0; i < X_cache_pointers[0]->x.size(); i++) {
			for (unsigned int j = 0; j < X_cache_pointers.size(); j++) {
				XmatFIO << X_cache_pointers[j]->x[i] << " ";
			}
			XmatFIO << std::endl;
		}
		XmatFIO << std::endl;
	}
	else if (cv_run){
		residual = compute_residual(std::vector<double>(1, 0.0), false);
		double InitialtrainMSE = compute_MSE_from_residual(residual);
		residual = compute_residual(std::vector<double>(1, 0.0), true);
		double InitialtestMSE = compute_MSE_from_residual(residual);
		

		while (lambdasIT != lambdas->end() && *lambdasIT > Lambda ) {
			if (report) {

				std::cout << std::setw(4) << nevent << ": Step to Lambda   (" << std::setw(4) << X_cache_pointers.size() << "): ";
				std::cout << *lambdasIT << " " << InitialtrainMSE << " ";
				if (test_set_size != 0) std::cout << InitialtestMSE << std::endl;
				else std::cout << "n/a" << std::endl;
			}
			//If lambda of interests are below computed lambda, 
			CVMSE->push_back(InitialtestMSE);
			//record error
			std::vector<double>::iterator rit = residual.begin();
			for (unsigned int i = 0; i < Data_Set.size(); i++) {
				if (test_logical_index.at(i)) {
					CVE->at(ilambda).at(i) = (*rit);
					rit++;
				}
			}
			ilambda++;
			lambdasIT++;
			nevent++;
		}
	}
	
	// Report
	if (report) {
		std::cout << std::setw(4) << nevent << ": Adding a graph   (   1): " << Lambda << " " << trainMSE << " ";
		if (test_set_size != 0) std::cout << testMSE << std::endl;
		else std::cout << "n/a" << std::endl;
	}
	nevent++;

	// initialize Gram Matrix
	Gram.resize(1, 1);
	double num = 0;
	for (unsigned int i = 0; i < Data_Set.size(); i++) {
		// pick out training set molecules
		if (!test_logical_index.at(i)) {
			num += selected_cache->x.at(i) * selected_cache->x.at(i);
		}
	}
	Gram(0, 0) = num;
	w = compute_residual(beta, false);
	// set up sign vector and gamma vector
	sign.push_back(sgn(compute_covariance_r_xtrain(selected_cache->x)));
	gamma.push_back(sign[0]/Gram(0,0));

	// turn off initial search
	initial_graph_search = false;
	/*********************** End Pick First Descriptor ***********************/
	/************************ Begin Perform LarsLasso ************************/
	// computation of w, v, rho0, and eta0. These are function of active set.
	w = compute_residual(beta, false);
	rho0 = compute_covariance_r_xtrain(X_cache_pointers[0]->x);
	compute_v();
	eta0 = compute_covariance_g_xtrain(X_cache_pointers[0]->x);

	/* GG DEBUG
	debug_graph.resize(2);
	debug_graph[0].label = 0;
	debug_graph[1].label = 0;
	debug_graph[0].push(0, 1, 1);
	debug_graph[1].push(1, 0, 1);

	debug_graph.resize(3);
	debug_graph[2].label = 0;
	debug_graph[1].push(1, 2, 1);
	debug_graph[2].push(2, 1, 1);

	debug_graph[2].push(2, 0, 3);
	debug_graph[0].push(0, 2, 3);

	//debug_graph.resize(4);
	//debug_graph[3].label = 0;
	//debug_graph[2].push(2, 3, 1);
	//debug_graph[3].push(3, 2, 1);
	GGDEBUG	*/

	unsigned int nlindep = 0;
	while (X_cache_pointers.size() < nmaxparam && X_cache_pointers.size() < train_set_size && Lambda > minlambda) {
		// Line 4 algorithm 1
		d1 = DBL_MAX;
		d2 = DBL_MAX;
		ID = 0;
		selected_cache = NULL;
		/* GGDEBUG
		if (nevent == 31) {
			std::cout << "wow";
		}
		GGDEBUG*/

		// DFS search through each vertex, and pick a optimal graph
		if (breadth_first) {
			project_breadth_first();
		}
		else {
			for (std::map<int, Projected>::iterator SingleNodeIt = Single_Node_DFS_Code_Cache_Map.begin();
				SingleNodeIt != Single_Node_DFS_Code_Cache_Map.end(); ++SingleNodeIt) project_depth_first(&(SingleNodeIt->second));
		}
		/* GGDEBUG
		// This means that no graph has selected. perhaps every graph has been added as active set. terminate
		if (nevent == 31) {
			std::cout << "wow";
		}
		*/
		if (selected_cache == NULL) {
			if (report) std::cout << std::setw(4) << nevent <<": Could not find a new descriptor." << std::endl;
			break;
		}
		
		// Checking the rank when the new descriptor is added.
		Eigen::MatrixXd TempGram = AppendXtoGram(Gram, selected_cache->x);
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QRGram(TempGram);
		unsigned int rank = QRGram.rank(); // is there a faster method to compute rank?
		if (rank != X_cache_pointers.size() + 1) {
			// This descriptor is dropped for good.
			nlindep++;
			selected_cache->linearlydependent = true; 
			linearly_dependent_caches.push_back(selected_cache);
			continue;
		}
		unsigned d2index = compute_d2(); // line 6 of algorithm 1.

		// Lambda of interest requires smaller step than d1 and d2, then the step is temporarily taken and MSE are computed.
		if (cv_run) {
			double d1lambda = d1 / train_set_size;
			double d2lambda = d2 / train_set_size;
			while (lambdasIT != lambdas->end() && (Lambda - *lambdasIT) < d1lambda && (Lambda - *lambdasIT) < d2lambda) {
				// compute a new beta
				std::vector<double> tempbeta = beta;
				for (unsigned int i = 0; i < tempbeta.size(); i++) tempbeta[i] += (Lambda - *lambdasIT)*train_set_size*gamma[i];
				// compute MSEs
				residual = compute_residual(tempbeta, false);
				double trainMSE = compute_MSE_from_residual(residual);
				residual = compute_residual(tempbeta, true);
				double testMSE = compute_MSE_from_residual(residual);
				// report
				if (report) {
					std::cout << std::setw(4) << nevent << ": Step to Lambda   (" << std::setw(4) << X_cache_pointers.size() << "): ";
					std::cout << *lambdasIT << " " << trainMSE << " ";
					if (test_set_size != 0) std::cout << testMSE << std::endl;
					else std::cout << "n/a" << std::endl;
				}
				// save MSE, look at next lambda
				CVMSE->push_back(testMSE);
				//record error
				std::vector<double>::iterator rit = residual.begin();
				for (unsigned int i = 0; i < Data_Set.size(); i++) {
					if (test_logical_index.at(i)) {
						CVE->at(ilambda).at(i) = (*rit);
						rit++;
					}
				}
				ilambda++;
				lambdasIT++;
				nevent++;
			}
			if (lambdasIT == lambdas->end()) {
				break;
			}
		}
		/*
		if (nevent == 30) {
			std::cout << "wow";
		}
		*/
		// line 7-8
		double d = std::min(d1, d2);
		for (unsigned int i = 0; i < beta.size(); i++) {
			beta[i] += d*gamma[i];
		}
		// line 9-10
		// Update w, v, rho0, and eta0 for next iteration and for covariance calculation.
		w = compute_residual(beta, false);
		rho0 = compute_covariance_r_xtrain(X_cache_pointers[0]->x);
		if (d1 < d2) { // append new graph.
			Gram = AppendXtoGram(Gram, selected_cache->x);
			sign.push_back(sgn(compute_covariance_r_xtrain(selected_cache->x)));
			selected_cache->ActiveSet = true;
			X_cache_pointers.push_back(selected_cache);
			beta.push_back(0.0);
			gamma.push_back(0.0);
			if (report) std::cout << std::setw(4) << nevent << ": Adding a graph   (" << std::setw(4) << X_cache_pointers.size() << "): ";
		}
		else { // remove graph
			Gram = RemoveXfromGram(Gram, d2index);
			sign.erase(sign.begin() + d2index);
			X_cache_pointers[d2index]->ActiveSet = false;
			X_cache_pointers.erase(X_cache_pointers.begin() + d2index);
			beta.erase(beta.begin() + d2index);
			gamma.erase(gamma.begin() + d2index);
			if (report) std::cout << std::setw(4) << nevent << ": Removing a graph (" << std::setw(4) << X_cache_pointers.size() << "): ";
			// See if any descriptor is no-longer linearly dependent.
			for (std::vector<DFSCodeCache *>::iterator it = linearly_dependent_caches.begin();
				it != linearly_dependent_caches.end();) {
				Eigen::MatrixXd TempGram = AppendXtoGram(Gram, (*it)->x);
				Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QRGram(TempGram);
				unsigned int rank = QRGram.rank(); 
				if (rank == X_cache_pointers.size() + 1) {
					// This descriptor is dropped for good.
					nlindep--;
					(*it)->linearlydependent = false;
					it = linearly_dependent_caches.erase(it);
					continue;
				}
				else {
					++it;
				}
			}
			
		}

		// line 11 inverse
		Eigen::VectorXd TempSign(sign.size());
		for (unsigned int i = 0; i < sign.size(); i++) {
			TempSign(i) = sign.at(i);
		}
		// Cholesky inverse solve.
		Eigen::VectorXd Tempgamma(gamma.size());
		Tempgamma = Gram.llt().solve(TempSign);
		for (unsigned int i = 0; i < sign.size(); i++) {
			gamma[i] = Tempgamma(i);
		}
		// Update w, v, rho0, and eta0 for next iteration and for covariance calculation.
		compute_v();
		eta0 = compute_covariance_g_xtrain(X_cache_pointers[0]->x);

		// Report. Regression is completed with current set. 
		residual = compute_residual(beta, false);
		double trainMSE = compute_MSE_from_residual(residual);
		residual = compute_residual(beta, true);
		double testMSE = compute_MSE_from_residual(residual);
		Lambda = fabs(rho0)/train_set_size;
		if (report) {
			std::cout << Lambda << " " << trainMSE << " ";
			if (test_set_size != 0) std::cout << testMSE << std::endl;
			else std::cout << "n/a" << std::endl;
		}
		if (!cv_run) {
			lambdas->push_back(Lambda);
			XmatFIO << "Event " << nevent << "; Lambda = " << Lambda << std::endl;
			for (unsigned int i = 0; i < X_cache_pointers[0]->x.size(); i++) {
				for (unsigned int j = 0; j < X_cache_pointers.size(); j++) {
					XmatFIO << X_cache_pointers[j]->x[i] << " ";
				}
				XmatFIO << std::endl;
			}
			XmatFIO << std::endl;
		}
		nevent++;

	}
	/************************* End Perform LarsLasso *************************/
	/***************************** Begin Report *****************************/
	
	beta = complete_regression();
	// compute MSEs
	residual = compute_residual(beta, false);
	trainMSE = compute_MSE_from_residual(residual);
	residual = compute_residual(beta, true);
	testMSE = compute_MSE_from_residual(residual);
	// compute w, v, rho0, and eta0
	w = compute_residual(beta, false);
	rho0 = compute_covariance_r_xtrain(X_cache_pointers[0]->x);
	Lambda = fabs(rho0) / train_set_size;
	if (report) {
		std::cout << std::setw(4) << nevent << ": Final step       (" << std::setw(4) << X_cache_pointers.size() <<
			"): " << Lambda << " " << trainMSE << " ";
		if (test_set_size != 0) std::cout << testMSE << std::endl;
		else std::cout << "n/a" << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "Number of descriptor searched: " << ndfs << std::endl;
		std::cout << "Number of dropped Linearly dependent graph: " << nlindep << std::endl;
		std::cout << "Regression Coefficients" << std::endl;
		for (unsigned int i = 0; i < beta.size(); i++) {
			std::cout << std::setw(13) << beta[i] << std::endl;
		}
		if (cv_run) {
			std::cout << "Test Set Residual" << std::endl;
			std::vector<double> residual = compute_residual(beta, true);
		}
		for (unsigned int i = 0; i < residual.size(); i++) {
			std::cout << residual[i] << std::endl;
		}
		if (enc) {
			std::cout << "Selected Graphs" << std::endl;
			for (unsigned int i = 0; i < X_cache_pointers.size(); i++) {
				std::cout << "{";
				for (unsigned int j = 0; j < X_cache_pointers[i]->x.size(); j++) {
					std::cout << " " << X_cache_pointers[i]->x[j];
				}
				std::cout << " }" << std::endl;
				report_graph(X_cache_pointers[i]->g, X_cache_pointers[i]->x);
				
			}
		}
	}
	/****************************** End Report ******************************/
	
	
	return;
}

void gSpan::Enumerate_single_node_DFScode()
{
	if (1 > maxpat_max) return;
	/* Do single node handling, as the normal gspan DFS code based processing
	* cannot find subgraphs of size |subg|==1.  Hence, we find frequent node
	* labels explicitly.
	*/
	std::map<unsigned int, std::map<unsigned int, unsigned int> > singleVertex;
	std::map<unsigned int, unsigned int> singleVertexLabel;

	// All the nodes are organized
	for (unsigned int id = 0; id < Data_Set.size(); ++id) {// for each graph
		for (unsigned int nid = 0; nid < Data_Set[id].size(); ++nid) { // for each vertex
			if (singleVertex[id][Data_Set[id][nid].label] == 0) {
				// number of graphs it appears in
				singleVertexLabel[Data_Set[id][nid].label] += 1; // record support
			}

			singleVertex[id][Data_Set[id][nid].label] += 1; // record how much each graph has.
		}
	}
	/* All minimum support node labels are frequent 'subgraphs'.
		* singleVertexLabel[nodelabel] gives the number of graphs it appears
		* in.
		*/
	// Organized nodes are processed further more and put into cache.
	for (std::map<unsigned int, unsigned int>::iterator it =
		singleVertexLabel.begin(); it != singleVertexLabel.end(); ++it)
	{
		ndfs += 1;
		// new cache 
		DFSCodeCache * Single_Node_DFS_Code_Cache = new DFSCodeCache;

		Single_Node_DFS_Code_Cache_Map[(*it).first].DFS_Code_Cache = Single_Node_DFS_Code_Cache;
		// set flags
		Single_Node_DFS_Code_Cache->next_extension_enumerated = true;
		Single_Node_DFS_Code_Cache->maxtoc = 0;
		// record Max pattern check
		Single_Node_DFS_Code_Cache->maxpatmincheck = 1 < maxpat_min;

		// append graph. For book keeping.
		Graph g(directed);
		g.resize(1);
		g[0].label = (*it).first;
		Single_Node_DFS_Code_Cache->g = g;
		// append frequency vector
		FrequencyVector x(Data_Set.size(), 0);
		std::vector<unsigned int> counts(Data_Set.size());
		for (std::map<unsigned int, std::map<unsigned int, unsigned int> >::iterator it2 =
			singleVertex.begin(); it2 != singleVertex.end(); ++it2)
		{
			x.at((*it2).first) = (*it2).second[(*it).first];
		}
		Single_Node_DFS_Code_Cache->x = x;	
		Single_Node_DFS_Code_Cache->x_upper = x;
	}

	// Enumerate first edge projected DFS.

	EdgeList edges; //vector of edges
					// enumerate first edges 
	for (unsigned int id = 0; id < Data_Set.size(); ++id) {
		Graph &g = Data_Set[id];
		for (unsigned int from = 0; from < g.size(); ++from) {
			if (get_forward_root(g, g[from], edges)) { // forward edges are appended here, and return true if edge list is not empty. 
				for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
					// append PDFS to the projected. Here Null is passed as an memory address for PDFS
					Single_Node_DFS_Code_Cache_Map[g[from].label].DFS_Code_Cache->new_fwd_DFS[0][(*it)->elabel][g[(*it)->to].label].push(&g, *it, 0);
			}
		}
	}

	// create a list of pointers to caches

	for (std::map<int, Projected>::iterator SingleNodeIt = Single_Node_DFS_Code_Cache_Map.begin();
		SingleNodeIt != Single_Node_DFS_Code_Cache_Map.end(); ++SingleNodeIt) 
	{
		for (Projected_iterator2 elabel = SingleNodeIt->second.DFS_Code_Cache->new_fwd_DFS[0].begin(); // given std::map<X, Y> where X is key and Y is value, std::map<X, Y>->second will give Y
			elabel != SingleNodeIt->second.DFS_Code_Cache->new_fwd_DFS[0].end();)
		{
			for (Projected_iterator1 tolabel = elabel->second.begin();
				tolabel != elabel->second.end();)
			{
				DFS_CODE.push(0, 1, SingleNodeIt->first, elabel->first, tolabel->first);
				

				if (DFS_check_condition(tolabel->second.PDFSs)) {
					tolabel->second.DFS_Code_Cache = new DFSCodeCache; // create cache
					tolabel->second.DFS_Code_Cache->Parent_DFS_Code_Cache = SingleNodeIt->second.DFS_Code_Cache;
					// record graph
					tolabel->second.DFS_Code_Cache->g = GRAPH_IS_MIN; // graph is built when is_min() is called.
																	  // minimum check does not stop algorithm as DFS code needs to be continued to be expanded before its size can exceed minimum.
					tolabel->second.DFS_Code_Cache->DFS_Code = DFS_CODE;
					tolabel->second.DFS_Code_Cache->maxpatmincheck = DFS_CODE.nodeCount() < maxpat_min;
					// Get frequency vector
					tolabel->second.DFS_Code_Cache->x = GetFrequencyVector(tolabel->second.PDFSs);
					tolabel->second.DFS_Code_Cache->x_upper = tolabel->second.DFS_Code_Cache->x;
					// Update upper bound of the parent DFS_Code
					for (DFSCodeCache* p = SingleNodeIt->second.DFS_Code_Cache; p; p = p->Parent_DFS_Code_Cache) {
						for (unsigned int i = 0; i < tolabel->second.DFS_Code_Cache->x_upper.size(); ++i) {
							if (p->x_upper.at(i) < tolabel->second.DFS_Code_Cache->x.at(i)) {
								p->x_upper.at(i) = tolabel->second.DFS_Code_Cache->x.at(i);
							}
						}
					}
					++tolabel;
					ndfs += 1;
				}
				else {
					// erase PDFSs
					for (std::vector< PDFS* >::iterator it = tolabel->second.PDFSs.begin(); it != tolabel->second.PDFSs.end(); ++it) delete (*it);
					tolabel->second.PDFSs.clear();
					std::vector<PDFS*>().swap(tolabel->second.PDFSs);
					// delete map key
					elabel->second.erase(tolabel++);
				}
				DFS_CODE.pop();
			}
			if (SingleNodeIt->second.DFS_Code_Cache->new_fwd_DFS[0].size() != 0) ++elabel;
			else SingleNodeIt->second.DFS_Code_Cache->new_fwd_DFS[0].erase(elabel++);
		}
	}
}

		
}
