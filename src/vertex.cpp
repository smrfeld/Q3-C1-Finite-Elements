#include "../include/q3c1_bits/vertex.hpp"

// Other headers
#include "../include/q3c1_bits/general.hpp"
#include "../include/q3c1_bits/basis_func.hpp"
#include "../include/q3c1_bits/cell.hpp"

#include <iostream>
#include <sstream>
#include <cmath>

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
        std::vector<DimType> dim_types_possible({DimType::VAL,DimType::DERIV});
		if (_no_dims == 1) {
            
            for (auto const &dim0: dim_types_possible) {
                _bfs.push_back(new BasisFunc(this,{dim0}));
            };
            
		} else if (_no_dims == 2) {

            for (auto const &dim0: dim_types_possible) {
                for (auto const &dim1: dim_types_possible) {
                    _bfs.push_back(new BasisFunc(this,{dim0,dim1}));
                };
            };

        } else if (_no_dims == 3) {

            for (auto const &dim0: dim_types_possible) {
                for (auto const &dim1: dim_types_possible) {
                    for (auto const &dim2: dim_types_possible) {
                        _bfs.push_back(new BasisFunc(this,{dim0,dim1,dim2}));
                    };
                };
            };

        } else if (_no_dims == 4) {
            
            for (auto const &dim0: dim_types_possible) {
                for (auto const &dim1: dim_types_possible) {
                    for (auto const &dim2: dim_types_possible) {
                        for (auto const &dim3: dim_types_possible) {
                            _bfs.push_back(new BasisFunc(this,{dim0,dim1,dim2,dim3}));
                        };
                    };
                };
            };
            
        } else if (_no_dims == 5) {
            
            for (auto const &dim0: dim_types_possible) {
                for (auto const &dim1: dim_types_possible) {
                    for (auto const &dim2: dim_types_possible) {
                        for (auto const &dim3: dim_types_possible) {
                            for (auto const &dim4: dim_types_possible) {
                                _bfs.push_back(new BasisFunc(this,{dim0,dim1,dim2,dim3,dim4}));
                            };
                        };
                    };
                };
            };
            
        } else if (_no_dims == 6) {
            
            for (auto const &dim0: dim_types_possible) {
                for (auto const &dim1: dim_types_possible) {
                    for (auto const &dim2: dim_types_possible) {
                        for (auto const &dim3: dim_types_possible) {
                            for (auto const &dim4: dim_types_possible) {
                                for (auto const &dim5: dim_types_possible) {
                                    _bfs.push_back(new BasisFunc(this,{dim0,dim1,dim2,dim3,dim4,dim5}));
                                };
                            };
                        };
                    };
                };
            };

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
        int idx=0;
        for (auto dim=0; dim<_no_dims; dim++) {
            if (dim_types[dim] == DimType::DERIV) {
                idx += pow(2,_no_dims-dim-1);
                // _no_dims = 1
                //     VAL: 0
                //     DERIV: 0 + 2^(1-0-1) = 0 + 2^0 = 1
                // _no_dims = 2
                //     VAL VAL: 0
                //     VAL DERIV: 0 + 2^(2-1-1) = 0 + 2^0 = 1
                //     DERIV VAL: 2^(2-0-1) + 0 = 2^1 + 0 = 2
                //     DERIV DERIV: 2^(2-0-1) + 2^(2-1-1) = 2^1 + 2^0 = 3
                // _no_dims = 3
                //     VAL VAL VAL: 0
                //     VAL VAL DERIV: 0 + 0 + 2^(3-2-1) = 0 + 0 + 2^0 = 1
                //     VAL DERIV VAL: 0 + 2^(3-1-1) + 0 = 0 + 2^1 + 0 = 2
                //     VAL DERIV DERIV: 0 + 2^(3-1-1) + 2^(3-2-1) = 0 + 2^1 + 2^0 = 3
                //     DERIV VAL VAL: 2^(3-0-1) + 0 + 0 = 2^2 + 0 + 0 = 4
                //     DERIV VAL DERIV: 2^(3-0-1) + 0 + 2^(3-2-1) = 2^2 + 0 + 2^0 = 5
                //     DERIV DERIV VAL: 2^(3-0-1) + 2^(3-1-1) + 0 = 2^2 + 2^1 + 0 = 6
                //     DERIV DERIV DERIV: 2^(3-0-1) + 2^(3-1-1) + 2^(3-2-1) = 2^2 + 2^1 + 2^0 = 7
            };
        };
        return _bfs[idx];
	};
};
