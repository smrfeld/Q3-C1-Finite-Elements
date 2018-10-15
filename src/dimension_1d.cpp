#include "../include/dcubic_bits/dimension_1d.hpp"

#include <iostream>
#include <cmath>

/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Dimension1D
	****************************************/

	Dimension1D::Dimension1D(double min_pt_pt, double max_pt_pt, int no_pts) {
		_min_pt = min_pt_pt;
		_max_pt = max_pt_pt;
		_no_pts = no_pts;
		_delta = (_max_pt - _min_pt) / (_no_pts - 1.0);
	};
	Dimension1D::Dimension1D(const Dimension1D& other) {
		_copy(other);
	};
	Dimension1D::Dimension1D(Dimension1D&& other) {
		_copy(other);
		other._reset();
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
			_copy(other);
			other._reset();
		};
		return *this;
	};
	Dimension1D::~Dimension1D() {
		_clean_up();
	};
	void Dimension1D::_copy(const Dimension1D& other)
	{
		_min_pt = other._min_pt;
		_max_pt = other._max_pt;
		_delta = other._delta;
		_no_pts = other._no_pts;
	};
	void Dimension1D::_reset()
	{
		_min_pt = 0.0;
		_max_pt = 0.0;
		_delta = 0.0;
		_no_pts = 0;
	};
	void Dimension1D::_clean_up() {
	};

	// Accessors
	int Dimension1D::get_no_pts() const {
		return _no_pts;
	};
	double Dimension1D::get_min_pt() const {
		return _min_pt;
	};
	double Dimension1D::get_max_pt() const {
		return _max_pt;
	};
	double Dimension1D::get_delta() const {
		return _delta;
	};

	// Get by idx
	double Dimension1D::get_pt_by_idx(int idx, bool start_at_one) const {
		if (start_at_one) {
			return _min_pt + (idx-1) * _delta;
		} else {
			return _min_pt + idx * _delta;
		};
	};

	// Check if point is in domain
	bool Dimension1D::check_if_pt_is_inside_domain(double x) const {
		if (x < _min_pt || x > _max_pt) { 
			return false; 
		} else {
			return true;
		};
	};

	// Get closest index
	int Dimension1D::get_closest_idx(double x, bool start_at_one) const {
		int i = get_idxs_surrounding_pt(x,start_at_one);
		if (abs(x - get_pt_by_idx(i,start_at_one)) < abs(x - get_pt_by_idx(i+1,start_at_one))) {
			return i;
		} else {
			return i+1;
		};
	};

	// Get indexes surrounding a point
	// ie point is between i and i+1 where i is returned
	int Dimension1D::get_idxs_surrounding_pt(double x, bool start_at_one) const {
		int i = (x - _min_pt) / _delta;
		if (i==_no_pts-1) {i--;};
		if (start_at_one) {
			return i+1;
		} else {
			return i;
		};
	};

	// Get fraction of a point between successive points
	double Dimension1D::get_frac_between(double x) const {
		return get_frac_between(x,get_idxs_surrounding_pt(x));
	};
	// Second optional specification: the return of the surrounding idxs
	double Dimension1D::get_frac_between(double x, int i, bool start_at_one) const {
		return (x - get_pt_by_idx(i,start_at_one)) / _delta;
	};
};