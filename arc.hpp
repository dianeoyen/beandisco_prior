/*
 *  BEANDiscoPrior: Arc header file
 *
 *  Copyright 2014 Diane Oyen <doyen at cs.unm.edu>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>

#include "common.hpp"
#include "stacksubset.hpp"

#ifndef ARC_HPP
#define ARC_HPP

struct Arc {
	int tail;
	int head;
	
	void setFirst() {
		head = 0;
		tail = 1;
	}
	
	bool next(int nNodes) {
		++tail;
		if (head == tail)
			++tail;
		if (tail >= nNodes) {
			tail = 0;
			++head;
			if (head >= nNodes) {
				head = 0;
				return false;
			}
		}
		return true;
	}

	bool holds(int v, const StackSubset& pa) {
		return (head != v || pa.contains(tail));
	}
};

std::ostream& operator<<(std::ostream& os, const Arc& arc) {
	return os << arc.tail << " -> " << arc.head;
}

const Arc NullArc = { -1, -1 };



/**
 * A templated map data structure with Arc as index type.
 */
template <class T>
class ArcMap {
private:
	int nNodes_;
	T* data_;
	
	ArcMap(const ArcMap&); // disable copy constructor
	ArcMap& operator=(const ArcMap&); // disable copying
	
public:
	ArcMap(int nNodes) {
		nNodes_ = nNodes;
		data_ = new T[nNodes * nNodes];
	}
	

	~ArcMap() {
		delete[] data_;
	}
	
	void setAll(T value) {
		for (int i = 0; i < nNodes_ * nNodes_; ++i)
			data_[i] = value;
	}
	
	T& operator[] (Arc arc) {
		return data_[arc.head + arc.tail * nNodes_];
	}

	T operator[] (Arc arc) const {
		return data_[arc.head + arc.tail * nNodes_];
	}

  int getNnodes() {
    return nNodes_;
  }

};

#endif
