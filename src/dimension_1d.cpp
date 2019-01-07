#include "../include/q3c1_bits/dimension_1d.hpp"

#include <iostream>
#include <cmath>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/****************************************
	Dimension1D
	****************************************/

	Dimension1D::Dimension1D(double delta, double zero) {
		_zero = zero;
		_delta = delta;
	};
	Dimension1D::Dimension1D(const Dimension1D& other) {
		_copy(other);
	};
	Dimension1D::Dimension1D(Dimension1D&& other) {
		_move(other);
	};
	Dimension1D& Dimension1D::operator=(const Dimension1D& other) {
		if (this != &other)
		{
			_clean_up();
			_copy(other);
		};
		return *this;
	};
	Dimension1D& Dimension1D::operator=(Dimension1D&& other) {
		if (this != &other)
		{
			_clean_up();
			_move(other);
		};
		return *this;
	};
	Dimension1D::~Dimension1D() {
		_clean_up();
	};
	void Dimension1D::_copy(const Dimension1D& other)
	{
		_zero = other._zero;
		_delta = other._delta;
	};
	void Dimension1D::_move(Dimension1D& other)
	{
		_zero = other._zero;
		_delta = other._delta;

		other._zero = 0.0;
		other._delta = 0.0;
	};
	void Dimension1D::_clean_up() {
	};

	// Accessors
	double Dimension1D::get_delta() const {
		return _delta;
	};
	double Dimension1D::get_zero() const {
		return _zero;
	};	

	// Get by idx
	double Dimension1D::get_pt_by_idx(int idx) const {
		return _zero + idx * _delta;
	};

	// Get closest index
	int Dimension1D::get_closest_idx(double x) const {
		int i = get_idxs_surrounding_pt(x);
		if (abs(x - get_pt_by_idx(i)) < abs(x - get_pt_by_idx(i+1))) {
			return i;
		} else {
			return i+1;
		};
	};

	// Get indexes surrounding a point
	// ie point is between i and i+1 where i is returned
	int Dimension1D::get_idxs_surrounding_pt(double x) const {
		int i = (x - _zero) / _delta;
		return i;
	};

	// Get fraction of a point between successive points
	double Dimension1D::get_frac_between(double x) const {
		return get_frac_between(x,get_idxs_surrounding_pt(x));
	};
	// Second optional specification: the return of the surrounding idxs
	double Dimension1D::get_frac_between(double x, int i) const {
		return (x - get_pt_by_idx(i)) / _delta;
	};
};