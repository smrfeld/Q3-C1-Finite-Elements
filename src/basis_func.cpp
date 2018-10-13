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
	BasisFunc::BasisFunc(Vertex* vertex) {
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
		_vertex = new Vertex(*other._vertex);
		_coeff = other._coeff;
	};
	void BasisFunc::_move(BasisFunc& other)
	{
		_vertex = other._vertex;
		_coeff = other._coeff;

		// Reset other
		other._vertex = nullptr;
		other._coeff = 0.0;
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



	/****************************************
	BasisFuncVal
	--- ABSTRACT BASE ---	
	****************************************/

	BasisFuncVal::BasisFuncVal(Vertex* vertex) : BasisFunc(vertex) {};
	BasisFuncVal::BasisFuncVal(const BasisFuncVal& other) : BasisFunc(other) {};
	BasisFuncVal::BasisFuncVal(BasisFuncVal&& other) : BasisFunc(std::move(other)) {};
	BasisFuncVal& BasisFuncVal::operator=(const BasisFuncVal &other) {
		if (this != &other) {
			BasisFunc::operator=(other);
		};
		return *this;
	};
	BasisFuncVal& BasisFuncVal::operator=(BasisFuncVal &&other) {
		if (this != &other) {
			BasisFunc::operator=(std::move(other));
		};
		return *this;
	};
	BasisFuncVal::~BasisFuncVal() {
		BasisFunc::~BasisFunc();
	};



	/****************************************
	BasisFuncDeriv
	--- ABSTRACT BASE ---	
	****************************************/

	BasisFuncDeriv::BasisFuncDeriv(Vertex* vertex) : BasisFunc(vertex) {};
	BasisFuncDeriv::BasisFuncDeriv(const BasisFuncDeriv& other) : BasisFunc(other) {};
	BasisFuncDeriv::BasisFuncDeriv(BasisFuncDeriv&& other) : BasisFunc(std::move(other)) {};
	BasisFuncDeriv& BasisFuncDeriv::operator=(const BasisFuncDeriv &other) {
		if (this != &other) {
			BasisFunc::operator=(other);
		};
		return *this;
	};
	BasisFuncDeriv& BasisFuncDeriv::operator=(BasisFuncDeriv &&other) {
		if (this != &other) {
			BasisFunc::operator=(std::move(other));
		};
		return *this;
	};
	BasisFuncDeriv::~BasisFuncDeriv() {
		BasisFunc::~BasisFunc();
	};


	/****************************************
	BasisFuncVal1D
	****************************************/

	double BasisFuncVal1D::get_val(IdxSet local_idxs, double x_frac) {
		if (local_idxs[0] == 0) {
			return 1.0 - 3.0*pow(x_frac,2) + 2.0*pow(x_frac,3);
		} else if (local_idxs[0] == 1) {
			return 3.0*pow(x_frac,2) - 2.0*pow(x_frac,3);
		} else {
			show_error("BasisFuncVal1D","get_val","idx should be 0,1");
			exit(EXIT_FAILURE);
		};
	};
	double BasisFuncVal1D::get_deriv(IdxSet local_idxs, double x_frac) {
		if (local_idxs[0] == 0) {
			return - 6.0*x_frac + 6.0*pow(x_frac,2);
		} else if (local_idxs[0] == 1) {
			return 6.0*x_frac - 6.0*pow(x_frac,2);
		} else {
			show_error("BasisFuncVal1D","get_deriv","idx should be 0,1");
			exit(EXIT_FAILURE);
		};
	};



	/****************************************
	BasisFuncDeriv1D
	****************************************/

	double BasisFuncDeriv1D::get_val(IdxSet local_idxs, double x_frac) {
		if (local_idxs[0] == 0) {
			return x_frac - 2.0*pow(x_frac,2) + pow(x_frac,3);
		} else if (local_idxs[0] == 1) {
			return pow(x_frac,3) - pow(x_frac,2);
		} else {
			show_error("BasisFuncDeriv1D","get_val","idx should be 0,1");
			exit(EXIT_FAILURE);
		};
	};
	double BasisFuncDeriv1D::get_deriv(IdxSet local_idxs, double x_frac) {
		if (local_idxs[0] == 0) {
			return 1.0 - 4.0*x_frac + 3.0*pow(x_frac,2);
		} else if (local_idxs[0] == 1) {
			return 3.0*pow(x_frac,2) - 2.0*x_frac;
		} else {
			show_error("BasisFuncDeriv1D","get_deriv","idx should be 0,1");
			exit(EXIT_FAILURE);
		};
	};

};