#include "gspan.h"
#include <iterator>
#include <ostream>

namespace GSPAN {

	std::istream  &read(std::istream &is, bool directed, InputData &input_data)
	{
		unsigned int id = 0;
		Graph g(directed);
		while (true) {
			g.read(is);
			g.id = id;
			if (g.empty()) break;
			input_data.push_back(g); // vector of graphs
			id++;
		}
		return is;
	}
	/* Special report function for single node graphs.
	*/
	void gSpan::report_graph(Graph &g, FrequencyVector x) {
		unsigned int sup = 0;
		

		if (where == false) {
			for (std::vector<unsigned int>::iterator it = x.begin();
				it != x.end(); ++it)
			{
				sup += *it;
			}
			*os << "t # " << ID << " * " << sup;
		}
		*os << '\n';

		g.write(*os);
		*os << '\n';
	}


}
