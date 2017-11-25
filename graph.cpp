/*
    $Id: graph.cpp,v 1.4 2004/05/21 05:50:13 taku-ku Exp $;
 
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
#include <cstring>
#include <string>
#include <iterator>
#include <strstream>
#include <set>
#include <assert.h>
// debug


namespace GSPAN {

template <class T, class Iterator>
void tokenize (const char *str, Iterator iterator)
{
	std::istrstream is (str, std::strlen(str));
	std::copy (std::istream_iterator <T> (is), std::istream_iterator <T> (), iterator); //istream_iterator separate by whitespace by defaults.
}

void Graph::buildEdge ()
{
	std::stringstream buf;
	std::map <std::string, unsigned int> tmp;

	unsigned int id = 0;
	for (unsigned int from = 0; from < size (); ++from) { // from each vertex
		for (Vertex::edge_iterator it = (*this)[from].edge.begin (); // edge from each vertex
				it != (*this)[from].edge.end (); ++it)
		{
			if (directed || from <= it->to) {
				buf << from << " " << it->to << " " << it->elabel << " ";
			}
			else {
				buf << it->to << " " << from << " " << it->elabel << " ";
			}

			// Assign unique id's to the edges.
			if (tmp.find (buf.str()) == tmp.end()) {
				it->id = id;
				tmp[buf.str()] = id;
				++id;
			} else {
				it->id = tmp[buf.str()];
			}
			buf.str(std::string());
		}
	}

	edge_size_ = id;
}

void Graph::_setEdgesize(unsigned int size) {
	edge_size_ = size;
}

std::istream &Graph::read (std::istream &is)
{   // read 
	std::vector <std::string> result;
	char line[1024];

	clear ();

	while (true) {

		unsigned int pos = is.tellg ();
		if (!is.getline(line, 1024)){
			break;
		}
		result.clear ();
		
		tokenize<std::string>(line, std::back_inserter (result));

		

		if (result.empty()) {
			// do nothing
		} else if (result[0] == "t") {
			if (! empty()) { // use as delimiter
				is.seekg (pos, std::ios_base::beg);
				break;
			} else {
				/*
				 * y = atoi (result[3].c_str());
				 */
			}
		} else if (result[0] == "v" && result.size() >= 3) {
			unsigned int id    = atoi (result[1].c_str());
			this->resize (id + 1);
			(*this)[id].label = atoi (result[2].c_str());
		} else if (result[0] == "e" && result.size() >= 4) {
			int from   = atoi (result[1].c_str());
			int to     = atoi (result[2].c_str());
			int elabel = atoi (result[3].c_str());

			if ((int)size () <= from || (int)size () <= to) {
				std::cerr << "Format Error:  define vertex lists before edges" << std::endl;
				exit (-1);
			}

			(*this)[from].push (from, to, elabel);
			if (directed == false)
				(*this)[to].push (to, from, elabel);
		} else if (result[0] == "y" && result.size() >= 2){
			this->y = atof(result[1].c_str());
			}
		else {
			std::cerr << "Format Error:  Unrecognized Format: " << line << std::endl;
			exit(-1);
		}
	}
	// assign unique id to edges
	buildEdge ();

	return is;
}

std::ostream &Graph::write (std::ostream &os)
{
	std::stringstream edge_buf;
	std::set <std::string> tmp;
	for (unsigned int from = 0; from < size(); ++from) {

		os << "v " << from << " " << (*this)[from].label << std::endl;

		for (Vertex::edge_iterator it = (*this)[from].edge.begin();
			it != (*this)[from].edge.end(); ++it) {
			if (directed || from <= it->to) {
				edge_buf << from << " " << it->to << " " << it->elabel << " ";
			}
			else {
				edge_buf << it->to << " " << from << " " << it->elabel << " ";
			}
			tmp.insert(edge_buf.str());
			edge_buf.str(std::string());
		}
	}
	for (std::set<std::string>::iterator it = tmp.begin(); it != tmp.end(); ++it) {
		os << "e " << *it << std::endl;
	}
	return os;
}

void Graph::check (void)
{
	/* Check all indices
	 */
	for (unsigned int from = 0 ; from < size () ; ++from) {
		//mexPrintf ("check vertex %d, label %d\n", from, (*this)[from].label);

		for (Vertex::edge_iterator it = (*this)[from].edge.begin ();
			it != (*this)[from].edge.end (); ++it)
		{
			//mexPrintf ("   check edge from %d to %d, label %d\n", it->from, it->to, it->elabel);
			assert (it->from >= 0 && it->from < size ());
			assert (it->to >= 0 && it->to < size ());
		}
	}
}

}

