/*
    $Id: gspan.h,v 1.6 2004/05/21 05:50:13 taku-ku Exp $;
 
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
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <Eigen/Core>


// Namespaces provide a method for preventing name conflicts in large projects.
namespace GSPAN {

// template is a most generic way to define things. Below just swaps two values.. wow. C++ needs a separate function for this.. lol.
template <class T> inline void _swap (T &x, T &y) { T z = x; x = y; y = z; }


struct Edge {
	unsigned int from;	    // Index of from vertex
	unsigned int to;        // Index of to vertex
	unsigned int elabel;    // edge label
	unsigned int id;        // Index of edge
	Edge(): from(0), to(0), elabel(0), id(0) {};
};

class Vertex
{
public:
	
	typedef std::vector<Edge>::iterator edge_iterator;// edge_iterator = std::vector<Edge>::iterator

	int label;
	std::vector<Edge> edge; // vector of sturcture edges. 

	void push (int from, int to, int elabel) // this adds a new edge to the vertex
	{
		edge.resize (edge.size()+1);
		edge[edge.size()-1].from = from;
		edge[edge.size()-1].to = to;
		edge[edge.size()-1].elabel = elabel;
		return;
	}
};

class Graph: public std::vector<Vertex> { // Graph is a child class of vector that is a vector of vertex
private:
	unsigned int edge_size_;
public:
	typedef std::vector<Vertex>::iterator vertex_iterator;

	Graph (bool _directed) // constructor
	{
		directed = _directed;
	};
	bool directed; // attribute
	unsigned int id;
	double y; // response/dependent variable.
	//  int y; // class label
	unsigned int edge_size ()   { return edge_size_; } //definition
	unsigned int vertex_size () { return (unsigned int)size(); } //definition
	void buildEdge (); //definition
	void _setEdgesize(unsigned int size); //definition
	std::istream &read (std::istream &); // read
	std::ostream &write (std::ostream &); // write
	void check (void); //definition
	friend bool operator == (const Graph &g1, const Graph &g2)
	{
		// check vertex size
		if (g1.size() != g2.size()) return false;
		// check edge size and vertex label
		for (unsigned int i = 0; i != g1.size(); i++) {
			if (g1[i].edge.size() != g2[i].edge.size()) return false;
			if (g1[i].label != g2[i].label) return false;
		}
		// check edges
		for (unsigned int i = 0; i != g1.size(); i++) {
			for (unsigned int j = 0; j != g1[i].edge.size(); j++) {
				if (g1[i].edge[j].from != g2[i].edge[j].from || g1[i].edge[j].to != g2[i].edge[j].to || g1[i].edge[j].elabel != g2[i].edge[j].elabel) return false;
			}
		}
		return true;
	}

	Graph(): edge_size_(0), directed(false) {};
};

class DFS { // same definition as edge in the paper. 
public:
	int from;
	int to;
	int fromlabel;
	int elabel;
	int tolabel;
	friend bool operator == (const DFS &d1, const DFS &d2)
	{
		return (d1.from == d2.from && d1.to == d2.to && d1.fromlabel == d2.fromlabel
			&& d1.elabel == d2.elabel && d1.tolabel == d2.tolabel);
	}
	friend bool operator != (const DFS &d1, const DFS &d2) { return (! (d1 == d2)); }
	DFS(): from(0), to(0), fromlabel(0), elabel(0), tolabel(0) {};
};

typedef std::vector<int> RMPath;

struct DFSCode: public std::vector <DFS> {
private:
	RMPath rmpath;
public:
	const RMPath& buildRMPath ();

	/* Convert current DFS code into a graph.
	 */
	bool toGraph (Graph &);

	/* Clear current DFS code and build code from the given graph.
	 */
	void fromGraph (Graph &g);

	/* Return number of nodes in the graph.
	 */
	unsigned int nodeCount (void);

	DFSCode operator=(const DFSCode & rhs);

	void push (int from, int to, int fromlabel, int elabel, int tolabel)
	{
		resize (size() + 1);
		DFS &d = (*this)[size()-1];

		d.from = from;
		d.to = to;
		d.fromlabel = fromlabel;
		d.elabel = elabel;
		d.tolabel = tolabel;
	}
	void pop () { resize (size()-1); }

	std::ostream &write (std::ostream &); // write
};

struct PDFS { // This is a DFS code that is projected to individual graphs in the data set
	Graph       *graph;	// original input graph
	Edge        *edge;  // this is the edge supplied by user, not edges of DFS tree. This is basically the edge added for this particular DFS node compared to its parent DFS node
	PDFS        *prev;
	PDFS(): graph(0), edge(0), prev(0) {};
	~PDFS() { 
		graph = 0;
		edge = 0;
		prev = 0;
	};
};

class History: public std::vector<Edge*> {
private:
	std::vector<int> edge;
	std::vector<int> vertex;

public:
	bool hasEdge   (unsigned int id) { return (bool)edge[id]; }
	bool hasVertex (unsigned int id) { return (bool)vertex[id]; }
	void build     (Graph &, PDFS *);
	History() {};
	History (Graph& g, PDFS *p) { build (g, p); }

};

class DFSCodeCache;

class Projected{
public:
	std::vector<PDFS * > PDFSs;
	DFSCodeCache* DFS_Code_Cache;

	void push (Graph *graph, Edge *edge, PDFS *prev)
	{
		PDFS * d = new PDFS;
		d->graph = graph; d->edge = edge; d->prev = prev;
		PDFSs.push_back(d);
	}
};



typedef std::vector <Edge*> EdgeList;

bool  get_forward_pure   (Graph&, Edge *,  int,    History&, EdgeList &);
bool  get_forward_rmpath (Graph&, Edge *,  int,    History&, EdgeList &);
bool  get_forward_root   (Graph&, Vertex&, EdgeList &);
Edge *get_backward       (Graph&, Edge *,  Edge *, History&);


//GHG, 3/31/2017 - Newly added caching function for faster re-runs
//Below two were originally private in gspan class. Moved to public for caching
typedef std::map<int, std::map <int, std::map <int, Projected> > >			Projected_map3;
typedef std::map<int, std::map <int, Projected> >							Projected_map2;
typedef std::map<int, Projected>                                            Projected_map1; // as in DFS projected to data set
typedef std::vector<unsigned int>													FrequencyVector;


// caching function to run faster
class DFSCodeCache {
public:
	Graph g;
	DFSCode DFS_Code;
	DFSCodeCache* Parent_DFS_Code_Cache;
	bool maxpatmincheck; // min number of node pattern check
	int maxtoc;
	Projected_map3 new_fwd_DFS; // cached right-most extension result
	Projected_map2 new_bck_DFS;
	FrequencyVector x;
	FrequencyVector x_upper;
	bool next_extension_enumerated;
	bool ActiveSet;
	bool linearlydependent;
	DFSCodeCache() {
		next_extension_enumerated = false;
		Parent_DFS_Code_Cache = NULL;
		ActiveSet = false;
		linearlydependent = false;
	};
};



typedef std::vector<DFSCodeCache*> DFSCodeCaches;
typedef std::map<int, Projected> SingleNodeDFSCodeCacheMap;
typedef std::vector<DFSCodeCache*>::iterator DFS_Code_Cache_iterator;

typedef std::vector < Graph > InputData;
std::istream  &read(std::istream &is, bool directed, InputData &input_data);

// typedefs
typedef std::map<int, std::map <int, std::map <int, Projected> > >::iterator Projected_iterator3;
typedef std::map<int, std::map <int, Projected> >::iterator                  Projected_iterator2;
typedef std::map<int, Projected>::iterator                                   Projected_iterator1;
typedef std::map<int, std::map <int, std::map <int, Projected> > >::reverse_iterator Projected_riterator3;

class gSpan {

private:

	// Variables
	std::vector < Graph >       Data_Set;
	DFSCode                     DFS_CODE;
	DFSCode                     DFS_CODE_IS_MIN;
	SingleNodeDFSCodeCacheMap	Single_Node_DFS_Code_Cache_Map; // contains cached information about single node. 
	Graph                       GRAPH_IS_MIN;
	unsigned int                ID;			// Sub graph index
	unsigned int                minsup;		// Minimum support
	unsigned int                maxpat_min;	// lower bound on node count
	unsigned int                maxpat_max;	// upper bound on node count
	bool                        where;					// Switch for whether or not to print out which data each subgraphs are located
	bool                        enc;
	bool                        directed;
	bool                        single_node_enumerated;
	std::ostream*               os;
	std::vector<DFSCodeCache *> linearly_dependent_caches;
	unsigned int                ilambda;
	// Definitions
	bool DFS_check_condition(std::vector<PDFS*> PDFSs);
	void report_graph(Graph &g, FrequencyVector fv);
	bool is_min ();
	bool project_is_min (Projected &);
	unsigned int support (std::vector<PDFS*> PDFSs);
	void project_depth_first (Projected * const &);
	void project_breadth_first(void);
	void Enumerate_single_node_DFScode();
	FrequencyVector GetFrequencyVector(std::vector<PDFS * > PDFSs);
	Eigen::MatrixXd AppendXtoGram(Eigen::MatrixXd Gram, FrequencyVector x);//
	Eigen::MatrixXd RemoveXfromGram(Eigen::MatrixXd Gram, unsigned int dimensionToRemove);
	double compute_gain(FrequencyVector x);
	bool gain_prune_condition(FrequencyVector x);
	double compute_covariance_r_xtrain(FrequencyVector x);					//
	void compute_v(void);											//
	double compute_covariance_g_xtrain(FrequencyVector x);					//
	double compute_dt(FrequencyVector x);							//
	bool d1_prune_condition(FrequencyVector x);						//
	unsigned int compute_d2(void);											//
	std::vector<double> complete_regression(void);
	std::vector<double> compute_residual(std::vector<double> Beta, bool TestOrTrain);
	double compute_MSE_from_residual(std::vector<double> residual);
	
public:

	// Variables
	std::vector<double>	        y;
	std::vector<double>         beta;
	std::vector<double>         gamma;
	unsigned int                ndfs;
	std::vector<bool>           test_logical_index;
	unsigned int                test_set_size;
	unsigned int                train_set_size;
	std::vector<int>            sign;
	std::vector<DFSCodeCache *> X_cache_pointers;
	double                      largest_gain;
	double                      d1;
	double                      d2;
	std::vector<double>         w;
	std::vector<double>         v;
	double                      rho0;
	double                      eta0;
	bool                        initial_graph_search;
	bool                        cv_run;
	bool                        report;
	bool                        breadth_first;
	double                      Lambda;
	DFSCodeCache *              selected_cache;
	//Graph debug_graph;
	// Definitions
	gSpan(void);
	void run(InputData input_data, std::ostream &_os,
		unsigned int _minsup, unsigned int _maxpat_min, unsigned int _maxpat_max,
		bool _enc, unsigned int _nmaxparam, double _minlambda, bool _directed,
		std::vector<bool> * test_logical_index, std::vector<double> * lambdas,
		std::vector<double> * CVMSE, std::vector < std::vector<double> >  * CVE);
	void reset(void);
	//GHG, 3/31/2017 
};


};


