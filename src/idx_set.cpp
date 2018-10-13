#include "../include/q3c1_bits/idx_set.hpp"

#include <iostream>
#include <sstream>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/****************************************
	Index set
	****************************************/

	// Constructors
	IdxSet::IdxSet(int no_idxs) {
		_no_idxs = no_idxs;
		_idxs = new int[_no_idxs];
		std::fill_n(_idxs,_no_idxs,0);
	};
	IdxSet::IdxSet(int no_idxs, int* idxs) {
		_no_idxs = no_idxs;
		_idxs = new int[_no_idxs];
		std::copy(idxs,idxs+_no_idxs,_idxs);
	};

	IdxSet::IdxSet(std::vector<int> idxs) {
		_no_idxs = idxs.size();
		_idxs = new int[_no_idxs];
		for (auto i=0; i<_no_idxs; i++) {
			_idxs[i] = idxs[i];
		};
	};
	IdxSet::IdxSet(const IdxSet& other) {
		_copy(other);
	};
	IdxSet::IdxSet(IdxSet&& other) {
		_move(other);
	};
    IdxSet& IdxSet::operator=(const IdxSet& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    IdxSet& IdxSet::operator=(IdxSet&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	IdxSet::~IdxSet()
	{
		_clean_up();
	};

	// Helpers
	void IdxSet::_clean_up()
	{
		if (_idxs) {
			delete[] _idxs;
			_idxs = nullptr;
		};
	};
	void IdxSet::_copy(const IdxSet& other)
	{
		_no_idxs = other._no_idxs;
		_idxs = new int[_no_idxs];
		std::copy(other._idxs,other._idxs+_no_idxs,_idxs);
	};
	void IdxSet::_move(IdxSet& other)
	{
		_no_idxs = other._no_idxs;
		_idxs = other._idxs;

		// Reset other
		other._no_idxs = 0;
		other._idxs = nullptr;
	};

	// Accessors
	int IdxSet::operator [](int idx) const {
		return _idxs[idx];
	};
	int & IdxSet::operator [](int idx) {
		return _idxs[idx];
	};

	// Size
	int IdxSet::size() const {
		return _no_idxs;
	};

 	// Print
 	std::string IdxSet::print() const {
 		std::ostringstream stream;
		stream << "(";
		for (auto i=0; i<_no_idxs; i++) {
			stream << _idxs[i];
			if (i != _no_idxs-1) {
				stream << " ";
			};
		};
		stream << ")";
	    return stream.str();
 	};

	// Printing
	std::ostream& operator<<(std::ostream& stream, const IdxSet& idxs) {
		stream << idxs.print();
		return stream;
	};

    // Comparator
    bool operator==(const IdxSet &lhs, const IdxSet &rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		if (lhs[i] != rhs[i]) {
    			return false;
    		};
    	};
    	return true;
    };

	// Math
    IdxSet operator+(IdxSet lhs, const IdxSet& rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] += rhs[i];
    	};
    	return lhs;
    };
    IdxSet operator-(IdxSet lhs, const IdxSet& rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] -= rhs[i];
    	};
    	return lhs;
    };
    IdxSet operator+(IdxSet lhs, int rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] += rhs;
    	};
    	return lhs;
    };
    IdxSet operator-(IdxSet lhs, int rhs) {
    	for (auto i=0; i<lhs.size(); i++) {
    		lhs[i] -= rhs;
    	};
    	return lhs;
    };
};