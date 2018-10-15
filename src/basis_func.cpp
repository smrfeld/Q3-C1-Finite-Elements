#include "../include/q3c1_bits/basis_func.hpp"

// Other headers
#include "../include/q3c1_bits/vertex.hpp"
#include "../include/q3c1_bits/general.hpp"

#include <iostream>
#include <sstream>
#include <cmath>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/****************************************
	BasisFunc
	****************************************/

	// Constructors
	BasisFunc::BasisFunc(Vertex* vertex, std::vector<DimType> dim_types) {
		if (!vertex) {
			show_error("BasisFunc","BasisFunc","Vertex must not be nullptr");
			exit(EXIT_FAILURE);
		};
		if (dim_types.size() != vertex->get_no_dims()) {
			show_error("BasisFunc","BasisFunc","# Dim types do not match vertex dim");
			exit(EXIT_FAILURE);
		};

		_no_dims = vertex->get_no_dims();
		_dim_types = dim_types;
		_vertex = vertex;
		_coeff = 0.0;
	};
	BasisFunc::BasisFunc(const BasisFunc& other) {
		_copy(other);
	};
	BasisFunc::BasisFunc(BasisFunc&& other) {
		_move(other);
	};
    BasisFunc& BasisFunc::operator=(const BasisFunc& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    BasisFunc& BasisFunc::operator=(BasisFunc&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	BasisFunc::~BasisFunc()
	{
		_clean_up();
	};

	// Helpers
	void BasisFunc::_clean_up()
	{
		// Nothing
	};
	void BasisFunc::_copy(const BasisFunc& other)
	{
		_no_dims = other._no_dims;
		_dim_types = other._dim_types;
		_vertex = new Vertex(*other._vertex);
		_coeff = other._coeff;
	};
	void BasisFunc::_move(BasisFunc& other)
	{
		_no_dims = other._no_dims;
		_dim_types = other._dim_types;
		_vertex = other._vertex;
		_coeff = other._coeff;

		// Reset other
		other._no_dims = 0;
		other._dim_types.clear();
		other._vertex = nullptr;
		other._coeff = 0.0;
	};

	/********************
	Dim types
	********************/

	DimType BasisFunc::get_dim_type(int dim) const {
		return _dim_types[dim];
	};
	const std::vector<DimType>& BasisFunc::get_all_dim_types() const {
		return _dim_types;
	};

	/********************
	Vertex
	********************/

	Vertex* BasisFunc::get_vertex() const {
		return _vertex;
	};

	/********************
	Coeff
	********************/

	double BasisFunc::get_coeff() const {
		return _coeff;
	};
	void BasisFunc::set_coeff(double val) {
		_coeff = val;
	};
	void BasisFunc::increment_coeff(double inc) {
		_coeff += inc;
	};

	/********************
	Get val/deriv based on local idx
	********************/

	double BasisFunc::get_bf_val(const IdxSet& local_idxs, const std::vector<double>& x_frac) const {
		double val=1.0;
		for (auto dim=0; dim<_no_dims; dim++) {
			// safety check
			/*
			if (_dim_types[dim] != 0 && _dim_types[dim] != 1) {
				show_error("BasisFunc","get_bf_val","idx should be 0,1");
				exit(EXIT_FAILURE);
			};
			*/

			if (_dim_types[dim] == DimType::VAL) {
				// VAL
				if (local_idxs[dim] == 0) {
					val *= 1.0 - 3.0*pow(x_frac[dim],2) + 2.0*pow(x_frac[dim],3);
				} else if (local_idxs[dim] == 1) {
					val *= 3.0*pow(x_frac[dim],2) - 2.0*pow(x_frac[dim],3);
				};
			} else {
				// DERIV
				if (local_idxs[dim] == 0) {
					val *= x_frac[dim] - 2.0*pow(x_frac[dim],2) + pow(x_frac[dim],3);
				} else if (local_idxs[dim] == 1) {
					val *= pow(x_frac[dim],3) - pow(x_frac[dim],2);
				};
			};
		};
		return val;
	};

	double BasisFunc::get_bf_deriv(const IdxSet& local_idxs, const std::vector<double>& x_frac, int deriv_dim) const {
		double val=1.0;
		for (auto dim=0; dim<_no_dims; dim++) {
			// safety check
			/*
			if (_dim_types[dim] != 0 && _dim_types[dim] != 1) {
				show_error("BasisFunc","get_bf_val","idx should be 0,1");
				exit(EXIT_FAILURE);
			};
			*/
			if (dim != deriv_dim) {

				// NO SPATIAL DERIV

				if (_dim_types[dim] == DimType::VAL) {
					// VAL
					if (local_idxs[dim] == 0) {
						val *= 1.0 - 3.0*pow(x_frac[dim],2) + 2.0*pow(x_frac[dim],3);
					} else if (local_idxs[dim] == 1) {
						val *= 3.0*pow(x_frac[dim],2) - 2.0*pow(x_frac[dim],3);
					};
				} else {
					// DERIV
					if (local_idxs[dim] == 0) {
						val *= x_frac[dim] - 2.0*pow(x_frac[dim],2) + pow(x_frac[dim],3);
					} else if (local_idxs[dim] == 1) {
						val *= pow(x_frac[dim],3) - pow(x_frac[dim],2);
					};
				};

			} else {


				// SPATIAL DERIV


				if (_dim_types[dim] == DimType::VAL) {
					// VAL
					if (local_idxs[dim] == 0) {
						val *= - 6.0*x_frac[dim] + 6.0*pow(x_frac[dim],2);
					} else if (local_idxs[dim] == 1) {
						val *= 6.0*x_frac[dim] - 6.0*pow(x_frac[dim],2);
					};
				} else {
					// DERIV
					if (local_idxs[dim] == 0) {
						val *= 1.0 - 4.0*x_frac[dim] + 3.0*pow(x_frac[dim],2);
					} else if (local_idxs[dim] == 1) {
						val *= 3.0*pow(x_frac[dim],2) - 2.0*x_frac[dim];
					};
				};
			};
		};

		return val;
	};
};