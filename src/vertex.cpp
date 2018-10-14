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
};