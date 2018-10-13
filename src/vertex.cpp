#include "../include/q3c1_bits/vertex.hpp"

// Other headers
#include "../include/q3c1_bits/general.hpp"
#include "../include/q3c1_bits/basis_func.hpp"
#include "../include/q3c1_bits/cell.hpp"

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
	Vertex::Vertex(IdxSet idxs, std::vector<double> abscissas) : _global_idxs(idxs) {
		// Check size
		if (idxs.size() != abscissas.size()) {
			show_error("Vertex","Vertex","idx size must match abscissa size");
			exit(EXIT_FAILURE);
		};

		_no_dims = abscissas.size();
		_abscissas = abscissas;

		_bf_val = nullptr;
		_bf_deriv = nullptr;
	};
	Vertex::Vertex(const Vertex& other) : _global_idxs(other._global_idxs) {
		_copy(other);
	};
	Vertex::Vertex(Vertex&& other) : _global_idxs(std::move(other._global_idxs)) {
		_move(other);
	};
    Vertex& Vertex::operator=(const Vertex& other) {
		if (this != &other) {
			_clean_up();
			_global_idxs = other._global_idxs;
			_copy(other);
		};
		return *this;
    };
    Vertex& Vertex::operator=(Vertex&& other) {
		if (this != &other) {
			_clean_up();
			_global_idxs = std::move(other._global_idxs);
			_move(other);
			other._global_idxs = IdxSet({});
		};
		return *this;
    };
	Vertex::~Vertex()
	{
		_clean_up();
	};

	// Helpers
	void Vertex::_clean_up()
	{
		if (_bf_val) { 
			delete _bf_val; 
			_bf_val = nullptr;
		};
		if (_bf_deriv) { 
			delete _bf_deriv; 
			_bf_deriv = nullptr;
		};
	};
	void Vertex::_copy(const Vertex& other)
	{
		_no_dims = other._no_dims;
		_abscissas = other._abscissas;
		_cells = other._cells;
		if (_no_dims == 1) {
			// Try to cast
			BasisFuncVal1D *bf_val_1d = dynamic_cast<BasisFuncVal1D*>(other._bf_val);
			// If cast fails, returns null
			if (bf_val_1d) {
				// Copy
				_bf_val = new BasisFuncVal1D(*bf_val_1d);
			};

			BasisFuncDeriv1D *bf_deriv_1d = dynamic_cast<BasisFuncDeriv1D*>(other._bf_deriv);
			if (bf_deriv_1d) {
				_bf_deriv = new BasisFuncDeriv1D(*bf_deriv_1d);
			};
		};
	};
	void Vertex::_move(Vertex& other)
	{
		_no_dims = other._no_dims;
		_abscissas = other._abscissas;
		_cells = other._cells;
		_bf_val = other._bf_val;
		_bf_deriv = other._bf_deriv;

		// Reset other
		other._no_dims = 0;
		other._abscissas.clear();
		other._cells.clear();
		other._bf_val = nullptr;
		other._bf_deriv = nullptr;
	};

	/********************
	Location
	********************/

	IdxSet Vertex::get_global_idxs() const {
		return _global_idxs;
	};
	std::vector<double> Vertex::get_abscissas() const {
		return _abscissas;
	};

	/********************
	Basis funcs
	********************/

	BasisFuncVal* Vertex::get_bf_val() const {
		return _bf_val;
	};
	BasisFuncDeriv* Vertex::get_bf_deriv() const {
		return _bf_deriv;
	};
};