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


		// Make bfs
		if (_no_dims == 1) {
			// Value
			_bfs.push_back(new BasisFunc(this,{DimType::VAL}));
			// Deriv
			_bfs.push_back(new BasisFunc(this,{DimType::DERIV}));
		} else if (_no_dims == 2) {
			// Value-value
			_bfs.push_back(new BasisFunc(this,{DimType::VAL,DimType::VAL}));
			// Value-deriv
			_bfs.push_back(new BasisFunc(this,{DimType::VAL,DimType::DERIV}));
			// Deriv-value
			_bfs.push_back(new BasisFunc(this,{DimType::DERIV,DimType::VAL}));
			// Deriv-deriv
			_bfs.push_back(new BasisFunc(this,{DimType::DERIV,DimType::DERIV}));
		} else if (_no_dims == 3) {
			// Value-value-value
			_bfs.push_back(new BasisFunc(this,{DimType::VAL,DimType::VAL,DimType::VAL}));
			// Value-value-deriv
			_bfs.push_back(new BasisFunc(this,{DimType::VAL,DimType::VAL,DimType::DERIV}));
			// Value-deriv-value
			_bfs.push_back(new BasisFunc(this,{DimType::VAL,DimType::DERIV,DimType::VAL}));
			// Deriv-value-value
			_bfs.push_back(new BasisFunc(this,{DimType::DERIV,DimType::VAL,DimType::VAL}));
			// Value-deriv-deriv
			_bfs.push_back(new BasisFunc(this,{DimType::VAL,DimType::DERIV,DimType::DERIV}));
			// Deriv-value-deriv
			_bfs.push_back(new BasisFunc(this,{DimType::DERIV,DimType::VAL,DimType::DERIV}));
			// Deriv-deriv-value
			_bfs.push_back(new BasisFunc(this,{DimType::DERIV,DimType::DERIV,DimType::VAL}));
			// Deriv-deriv-deriv
			_bfs.push_back(new BasisFunc(this,{DimType::DERIV,DimType::DERIV,DimType::DERIV}));
		} else {
			show_error("Vertex","Vertex","only no dims 1,2,3 are supported");
			exit(EXIT_FAILURE);
		};
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
		for (auto &bf: _bfs) {
			if (bf) {
				delete bf;
				bf = nullptr;
			};
		};
		_bfs.clear();
	};
	void Vertex::_copy(const Vertex& other)
	{
		_no_dims = other._no_dims;
		_abscissas = other._abscissas;
		_cells = other._cells;
		for (auto &bf: other._bfs) {
			_bfs.push_back(new BasisFunc(*bf));
		};
	};
	void Vertex::_move(Vertex& other)
	{
		_no_dims = other._no_dims;
		_abscissas = other._abscissas;
		_cells = other._cells;
		_bfs = other._bfs;

		// Reset other
		other._no_dims = 0;
		other._abscissas.clear();
		other._cells.clear();
		other._bfs.clear();
	};

	/********************
	Location
	********************/

	int Vertex::get_no_dims() const {
		return _no_dims;
	};

	int Vertex::get_global_idx(int dim) const {
		return _global_idxs[dim];
	};
	IdxSet Vertex::get_global_idxs() const {
		return _global_idxs;
	};

	double Vertex::get_abscissa(int dim) const {
		return _abscissas[dim];
	};
	std::vector<double> Vertex::get_abscissas() const {
		return _abscissas;
	};

	/********************
	Basis funcs
	********************/

	const std::vector<BasisFunc*>& Vertex::get_bfs() const {
		return _bfs;
	};
	BasisFunc* Vertex::get_bf(const std::vector<DimType> &dim_types) const {
		if (_no_dims == 1) {
			if (dim_types[0] == DimType::VAL) {
				return _bfs[0];
			} else {
				return _bfs[1];
			};
		} else if (_no_dims == 2) {
			if (dim_types[0] == DimType::VAL && dim_types[1] == DimType::VAL) {
				return _bfs[0];
			} else if (dim_types[0] == DimType::VAL && dim_types[1] == DimType::DERIV) {
				return _bfs[1];
			} else if (dim_types[0] == DimType::DERIV && dim_types[1] == DimType::VAL) {
				return _bfs[2];
			} else {
				return _bfs[3];
			};
		} else if (_no_dims == 3) {
			if (dim_types[0] == DimType::VAL && dim_types[1] == DimType::VAL && dim_types[2] == DimType::VAL) {
				return _bfs[0];
			} else if (dim_types[0] == DimType::VAL && dim_types[1] == DimType::VAL && dim_types[2] == DimType::DERIV) {
				return _bfs[1];
			} else if (dim_types[0] == DimType::VAL && dim_types[1] == DimType::DERIV && dim_types[2] == DimType::VAL) {
				return _bfs[2];
			} else if (dim_types[0] == DimType::DERIV && dim_types[1] == DimType::VAL && dim_types[2] == DimType::VAL) {
				return _bfs[3];
			} else if (dim_types[0] == DimType::VAL && dim_types[1] == DimType::DERIV && dim_types[2] == DimType::DERIV) {
				return _bfs[4];
			} else if (dim_types[0] == DimType::DERIV && dim_types[1] == DimType::VAL && dim_types[2] == DimType::DERIV) {
				return _bfs[5];
			} else if (dim_types[0] == DimType::DERIV && dim_types[1] == DimType::DERIV && dim_types[2] == DimType::VAL) {
				return _bfs[6];
			} else {
				return _bfs[7];
			};		
		} else {
			show_error("Vertex","get_bf","only dims 1,2,3 supported");
			exit(EXIT_FAILURE);
		};
	};
};