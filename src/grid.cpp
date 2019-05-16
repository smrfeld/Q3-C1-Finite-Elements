#include "../include/q3c1_bits/grid.hpp"

#include "../include/q3c1_bits/vertex.hpp"
#include "../include/q3c1_bits/dimension_1d.hpp"
#include "../include/q3c1_bits/cell.hpp"
#include "../include/q3c1_bits/general.hpp"
#include "../include/q3c1_bits/basis_func.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/****************************************
	Grid
	****************************************/

	// Constructors
    Grid::Grid(std::vector<Dimension1D*> dims) : _idxs_get_cell(dims.size()) {
		if (dims.size() > 6 || dims.size() < 1) {
			std::cerr << ">>> Grid::Grid <<< Error: only 1,2,3,4,5,6 dims supported" << std::endl;
			exit(EXIT_FAILURE);
		};

		_no_dims = dims.size();
		_dims = dims;
        _frac_get_cell= std::vector<double>(_no_dims);
	};
    Grid::Grid(const Grid& other) : _idxs_get_cell(other._idxs_get_cell) {
		_copy(other);
	};
    Grid::Grid(Grid&& other) : _idxs_get_cell(other._idxs_get_cell)  {
		_move(other);
	};
    Grid& Grid::operator=(const Grid& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    Grid& Grid::operator=(Grid&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	Grid::~Grid()
	{
		_clean_up();
	};

	// Helpers
	void Grid::_clean_up()
	{
		for (auto &c: _cells) {
			if (c.second) {
				delete c.second;
				_cells[c.first] = nullptr;
			};
		};
		_cells.clear();
		for (auto &v: _verts) {
			if (v.second) {
				delete v.second;
				_verts[v.first] = nullptr;
			};
		};
		_verts.clear();
        _frac_get_cell.clear();
	};
	void Grid::_copy(const Grid& other)
	{
		_no_dims = other._no_dims;
		_dims = other._dims;
		for (auto &v: other._verts) {
			_verts[v.first] = new Vertex(*v.second);
		};
		for (auto &c: other._cells) {
			_cells[c.first] = new Cell(*c.second);
		};
        _frac_get_cell = other._frac_get_cell;
        _idxs_get_cell = other._idxs_get_cell;
	};
	void Grid::_move(Grid& other)
	{
		_no_dims = other._no_dims;
		_dims = other._dims;
		_verts = other._verts;
		_cells = other._cells;
        _frac_get_cell = other._frac_get_cell;
        _idxs_get_cell = other._idxs_get_cell;

		// Reset other
		other._no_dims = 0;
		other._dims.clear();
		other._verts.clear();
		other._cells.clear();
        other._frac_get_cell.clear();
	};

	/********************
	Make a cell
	********************/

	Vertex* Grid::_make_vertex(IdxSet idxs) const {		
		std::vector<double> abscissas(_no_dims);
		for (auto dim=0; dim<_no_dims; dim++) {
			abscissas[dim] = _dims[dim]->get_pt_by_idx(idxs[dim]);
		};

        // std::cout << "Making vertex: " << abscissas[0] << " " << abscissas[1] << " " << abscissas[2] << std::endl;
        
		Vertex* vnew = new Vertex(idxs, abscissas);
		_verts[idxs] = vnew;
		return _verts[idxs];
	};
	Vertex* Grid::_get_or_make_vertex(IdxSet idxs) const {
		// Similar to the get_vertex function, but this one either gives the vertex or directly makes it (no cell)
		auto it = _verts.find(idxs);
		if (it != _verts.end()) {
			return it->second;
		} else {
			return _make_vertex(idxs);
		};
	};

	Cell* Grid::_make_cell(IdxSet idxs) const {

		// std::cout << "Making cell: " << idxs << std::endl;

		// Collect needed verts
		IdxSet idxs_vert(_no_dims), idxs_local(_no_dims);
		std::vector<std::pair<IdxSet,Vertex*>> verts_local;
        
		// Go through all dims
		if (_no_dims == 1) {

            for (idxs_local[0]=0; idxs_local[0]<=1; idxs_local[0]++) {
                for (auto dim=0; dim<_no_dims; dim++) {
                    idxs_vert[dim] = idxs[dim] + idxs_local[dim];
                };

                verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));
            };
            
		} else if (_no_dims == 2) {

            for (idxs_local[0]=0; idxs_local[0]<=1; idxs_local[0]++) {
                for (idxs_local[1]=0; idxs_local[1]<=1; idxs_local[1]++) {
                    for (auto dim=0; dim<_no_dims; dim++) {
                        idxs_vert[dim] = idxs[dim] + idxs_local[dim];
                    };

                    verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));
                };
            };

		} else if (_no_dims == 3) {

            for (idxs_local[0]=0; idxs_local[0]<=1; idxs_local[0]++) {
                for (idxs_local[1]=0; idxs_local[1]<=1; idxs_local[1]++) {
                    for (idxs_local[2]=0; idxs_local[2]<=1; idxs_local[2]++) {
                        for (auto dim=0; dim<_no_dims; dim++) {
                            idxs_vert[dim] = idxs[dim] + idxs_local[dim];
                        };
                    
                        verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));
                    };
                };
            };

        } else if (_no_dims == 4) {
            
            for (idxs_local[0]=0; idxs_local[0]<=1; idxs_local[0]++) {
                for (idxs_local[1]=0; idxs_local[1]<=1; idxs_local[1]++) {
                    for (idxs_local[2]=0; idxs_local[2]<=1; idxs_local[2]++) {
                        for (idxs_local[3]=0; idxs_local[3]<=1; idxs_local[3]++) {
                            for (auto dim=0; dim<_no_dims; dim++) {
                                idxs_vert[dim] = idxs[dim] + idxs_local[dim];
                            };
                        
                            verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));
                        };
                    };
                };
            };
            
        } else if (_no_dims == 5) {
            
            for (idxs_local[0]=0; idxs_local[0]<=1; idxs_local[0]++) {
                for (idxs_local[1]=0; idxs_local[1]<=1; idxs_local[1]++) {
                    for (idxs_local[2]=0; idxs_local[2]<=1; idxs_local[2]++) {
                        for (idxs_local[3]=0; idxs_local[3]<=1; idxs_local[3]++) {
                            for (idxs_local[4]=0; idxs_local[4]<=1; idxs_local[4]++) {
                                for (auto dim=0; dim<_no_dims; dim++) {
                                    idxs_vert[dim] = idxs[dim] + idxs_local[dim];
                                };
                                
                                verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));
                            };
                        };
                    };
                };
            };
            
        } else if (_no_dims == 6) {
            
            for (idxs_local[0]=0; idxs_local[0]<=1; idxs_local[0]++) {
                for (idxs_local[1]=0; idxs_local[1]<=1; idxs_local[1]++) {
                    for (idxs_local[2]=0; idxs_local[2]<=1; idxs_local[2]++) {
                        for (idxs_local[3]=0; idxs_local[3]<=1; idxs_local[3]++) {
                            for (idxs_local[4]=0; idxs_local[4]<=1; idxs_local[4]++) {
                                for (idxs_local[5]=0; idxs_local[5]<=1; idxs_local[5]++) {
                                    for (auto dim=0; dim<_no_dims; dim++) {
                                        idxs_vert[dim] = idxs[dim] + idxs_local[dim];
                                    };
                                    
                                    verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));
                                };
                            };
                        };
                    };
                };
            };
            
        };

		// Make
		_cells[idxs] = new Cell(idxs, verts_local);
		return _cells[idxs];	
	};

	/********************
	Print
	********************/

	void Grid::print(const std::vector<DimType>& dim_types) const {
		for (auto &vert: _verts) {
			print(vert.second,dim_types);
		};
	};
	void Grid::print(Vertex* vert, const std::vector<DimType>& dim_types) const {
        for (auto dim=0; dim<_no_dims; dim++) {
            std::cout << vert->get_abscissa(dim) << std::flush;
            if (dim != _no_dims-1) {
                std::cout << " " << std::flush;
            } else {
                std::cout << ": " << std::flush;
            };
        };
		std::cout << vert->get_bf(dim_types)->get_coeff() << std::endl;
	};

	/********************
	Get dims
	********************/

	int Grid::get_no_dims() const {
		return _no_dims;
	};
	Dimension1D* Grid::get_dim(int dim) const {
		return _dims[dim];
	};
	const std::vector<Dimension1D*>& Grid::get_all_dims() const {
		return _dims;
	};

	/********************
	Get vertices
	********************/

	Vertex* Grid::get_vertex(IdxSet idxs) const {
		auto it = _verts.find(idxs);
		if (it != _verts.end()) {
			return it->second;
		} else {
			// Vertices should not exist without cells
			// Cell at this idx does not exist b/c vertex does not exist
			// Call make cell will make vertex!
			_make_cell(idxs);
			// Now it must exist; try again!
			return get_vertex(idxs);
		};
	};

	const std::map<IdxSet,Vertex*>& Grid::get_all_vertices() const {
		return _verts;
	};

	/********************
	Get cell
	********************/

	Cell* Grid::get_cell(const IdxSet& idxs) const {
        // std::cout << "Get cell: " << idxs << std::endl;
		auto it = _cells.find(idxs);
		if (it != _cells.end()) {
			return it->second;
		} else {
			return _make_cell(idxs);
		};
	};

	std::pair<Cell*,const std::vector<double>&> Grid::get_cell(const std::vector<double>& abscissas) const {
        // std::cout << "Get cell: " << abscissas[0] << " " << abscissas[1] << " " << abscissas[2] << std::endl;
		for (auto dim=0; dim<_no_dims; dim++) {
			_idxs_get_cell[dim] = _dims[dim]->get_idxs_surrounding_pt(abscissas[dim]);
			_frac_get_cell[dim] = _dims[dim]->get_frac_between(abscissas[dim],_idxs_get_cell[dim]);
		};
        // std::cout << "Surrounding idxs: " << idxs << std::endl;
		return std::make_pair(get_cell(_idxs_get_cell),_frac_get_cell);
	};

	const std::map<IdxSet,Cell*>& Grid::get_all_cells() const {
		return _cells;
	};

	/********************
	Get vals
	********************/

	double Grid::get_val(const std::vector<double>& abscissas) const {

		// Get cell and fraction of this abscissa in the cell
		std::pair<Cell*,const std::vector<double>&> pr = get_cell(abscissas);

		double val = 0.0;

		// Run through all verts of the cell
		for (auto const &v_pr: pr.first->get_all_vertices()) {
			// Run through all bfs defined on this vertex
			for (auto const &bf: v_pr.second->get_bfs()) {
				val += bf->get_coeff() * bf->get_bf_val(v_pr.first, pr.second);
			};
		};

		return val;
	};
	double Grid::get_deriv_wrt_abscissa(const std::vector<double>& abscissas, int deriv_dim) const {

		// Get cell and fraction of this abscissa in the cell
		std::pair<Cell*,const std::vector<double>&> pr = get_cell(abscissas);

		double val = 0.0;

		// Run through all verts of the cell
		for (auto const &v_pr: pr.first->get_all_vertices()) {
			// Run through all bfs defined on this vertex
			for (auto const &bf: v_pr.second->get_bfs()) {
				val += bf->get_coeff() * bf->get_bf_deriv(v_pr.first, pr.second, deriv_dim);
			};
		};

		return val;

	};
	double Grid::get_deriv_wrt_coeff(const std::vector<double>& abscissas, const IdxSet& global_vertex_idxs, const std::vector<DimType>& dim_types) const {

		// Get cell and fraction of this abscissa in the cell
		std::pair<Cell*,const std::vector<double>&> pr = get_cell(abscissas);

		// Get the vertex
		Vertex* vert = get_vertex(global_vertex_idxs);

		// What is the local idx of this vert?
		IdxSet idxs_local = pr.first->get_local_idxs_of_vertex(vert);

		// Get the basis func
		BasisFunc* bf = vert->get_bf(dim_types);

		return bf->get_bf_val(idxs_local,pr.second);
	};
	std::map<Vertex*,std::vector<double>> Grid::get_deriv_wrt_coeffs_for_all_surrounding_verts(const std::vector<double>& abscissas) const {

		// Get cell and fraction of this abscissa in the cell
		std::pair<Cell*,const std::vector<double>&> pr_cell = get_cell(abscissas);

		// Returned
		std::map<Vertex*,std::vector<double>> ret;

		// Go through all verts of the cell
        // First are the local idxs; second is the vertex itself
        std::vector<DimType> dim_types_possible({DimType::VAL,DimType::DERIV});
		if (_no_dims == 1) {
            
			for (auto &pr_v: pr_cell.first->get_all_vertices()) {
                for (auto const &dim0: dim_types_possible) {
                    ret[pr_v.second].push_back(pr_v.second->get_bf({dim0})->get_bf_val(pr_v.first,pr_cell.second));
                };
			};
            
		} else if (_no_dims == 2) {
            
			for (auto &pr_v: pr_cell.first->get_all_vertices()) {
                for (auto const &dim0: dim_types_possible) {
                    for (auto const &dim1: dim_types_possible) {
                        ret[pr_v.second].push_back(pr_v.second->get_bf({dim0,dim1})->get_bf_val(pr_v.first,pr_cell.second));
                    };
                };
            };
            
		} else if (_no_dims == 3) {
            
            for (auto &pr_v: pr_cell.first->get_all_vertices()) {
                for (auto const &dim0: dim_types_possible) {
                    for (auto const &dim1: dim_types_possible) {
                        for (auto const &dim2: dim_types_possible) {
                            ret[pr_v.second].push_back(pr_v.second->get_bf({dim0,dim1,dim2})->get_bf_val(pr_v.first,pr_cell.second));
                        };
                    };
                };
            };

        } else if (_no_dims == 4) {
            
            for (auto &pr_v: pr_cell.first->get_all_vertices()) {
                for (auto const &dim0: dim_types_possible) {
                    for (auto const &dim1: dim_types_possible) {
                        for (auto const &dim2: dim_types_possible) {
                            for (auto const &dim3: dim_types_possible) {
                                ret[pr_v.second].push_back(pr_v.second->get_bf({dim0,dim1,dim2,dim3})->get_bf_val(pr_v.first,pr_cell.second));
                            };
                        };
                    };
                };
            };

        } else if (_no_dims == 5) {
            
            for (auto &pr_v: pr_cell.first->get_all_vertices()) {
                for (auto const &dim0: dim_types_possible) {
                    for (auto const &dim1: dim_types_possible) {
                        for (auto const &dim2: dim_types_possible) {
                            for (auto const &dim3: dim_types_possible) {
                                for (auto const &dim4: dim_types_possible) {
                                ret[pr_v.second].push_back(pr_v.second->get_bf({dim0,dim1,dim2,dim3,dim4})->get_bf_val(pr_v.first,pr_cell.second));
                                };
                            };
                        };
                    };
                };
            };

        } else if (_no_dims == 6) {
            
            for (auto &pr_v: pr_cell.first->get_all_vertices()) {
                for (auto const &dim0: dim_types_possible) {
                    for (auto const &dim1: dim_types_possible) {
                        for (auto const &dim2: dim_types_possible) {
                            for (auto const &dim3: dim_types_possible) {
                                for (auto const &dim4: dim_types_possible) {
                                    for (auto const &dim5: dim_types_possible) {
                                        ret[pr_v.second].push_back(pr_v.second->get_bf({dim0,dim1,dim2,dim3,dim4,dim5})->get_bf_val(pr_v.first,pr_cell.second));
                                    };
                                };
                            };
                        };
                    };
                };
            };
        };
        
		return ret;
	};


	/********************
	Read/write grid
	********************/

	void Grid::read_from_file(std::string fname) {
		std::ifstream f;

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			show_error("Grid","read_from_file","Could not open file");
			exit(EXIT_FAILURE);
		};

		// Abscissas
		std::string abscissa_str;
		double abscissa;
		IdxSet vertex_idxs(_no_dims);
		Vertex *vert;

		// Ordinates
		std::string ordinate_str="";
		double ordinate;

        // Possible dims
        std::vector<DimType> dim_types_possible({DimType::VAL, DimType::DERIV});

		// Read
		std::string line;
		std::istringstream iss;
		while (getline(f,line)) {
			// Skip empty lines
			if (line == "") { continue; };
			// Line
			iss = std::istringstream(line);

			// Abscissas
			for (auto dim=0; dim<_no_dims; dim++) {
				abscissa_str = "";
				iss >> abscissa_str;
				abscissa = atof(abscissa_str.c_str());
				vertex_idxs[dim] = _dims[dim]->get_closest_idx(abscissa);
			};

			// Get vertex
			vert = get_vertex(vertex_idxs);

			// Ordinate
			if (_no_dims == 1) {
                
                for (auto const &dim0: dim_types_possible) {
                    ordinate_str = "";
                    iss >> ordinate_str;
                    ordinate = atof(ordinate_str.c_str());
                    vert->get_bf({dim0})->set_coeff(ordinate);
                };
                
			} else if (_no_dims == 2) {

                for (auto const &dim0: dim_types_possible) {
                    for (auto const &dim1: dim_types_possible) {
                        ordinate_str = "";
                        iss >> ordinate_str;
                        ordinate = atof(ordinate_str.c_str());
                        vert->get_bf({dim0,dim1})->set_coeff(ordinate);
                    };
                };

            } else if (_no_dims == 3) {
                
                for (auto const &dim0: dim_types_possible) {
                    for (auto const &dim1: dim_types_possible) {
                        for (auto const &dim2: dim_types_possible) {
                            ordinate_str = "";
                            iss >> ordinate_str;
                            ordinate = atof(ordinate_str.c_str());
                            vert->get_bf({dim0,dim1,dim2})->set_coeff(ordinate);
                        };
                    };
                };
                
            } else if (_no_dims == 4) {
                
                for (auto const &dim0: dim_types_possible) {
                    for (auto const &dim1: dim_types_possible) {
                        for (auto const &dim2: dim_types_possible) {
                            for (auto const &dim3: dim_types_possible) {
                                ordinate_str = "";
                                iss >> ordinate_str;
                                ordinate = atof(ordinate_str.c_str());
                                vert->get_bf({dim0,dim1,dim2,dim3})->set_coeff(ordinate);
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
                                    ordinate_str = "";
                                    iss >> ordinate_str;
                                    ordinate = atof(ordinate_str.c_str());
                                    vert->get_bf({dim0,dim1,dim2,dim3,dim4})->set_coeff(ordinate);
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
                                        ordinate_str = "";
                                        iss >> ordinate_str;
                                        ordinate = atof(ordinate_str.c_str());
                                        vert->get_bf({dim0,dim1,dim2,dim3,dim4,dim5})->set_coeff(ordinate);
                                    };
                                };
                            };
                        };
                    };
                };
            };
		};

		// Close
		f.close();
	};
	void Grid::write_to_file(std::string fname) const {
		std::ofstream f;
        f.precision(16);

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			show_error("Grid","write_to_file","Could not open file");
			exit(EXIT_FAILURE);
		};

        // Dims possible
        std::vector<DimType> dim_types_possible({DimType::VAL, DimType::DERIV});
        int space_ctr = 0;
        
		// Go through all grid pts
		bool first_line=true;
		for (auto &vert: _verts) {

			if (!first_line) {
				f << "\n";
			} else {
				first_line = false;
			};

			// Write abscissa
			for (auto dim=0; dim<_no_dims; dim++) {
				f << vert.second->get_abscissa(dim) << " ";
			};

			// Write vals for each bf
            space_ctr = 0;
			if (_no_dims == 1) {
                
                for (auto dim0: dim_types_possible) {
                    f << vert.second->get_bf({dim0})->get_coeff();
                    space_ctr++;
                    if (space_ctr != pow(2,_no_dims)-1) {
                        f << " ";
                    };
                };

            } else if (_no_dims == 2) {
                
                for (auto dim0: dim_types_possible) {
                    for (auto dim1: dim_types_possible) {
                        f << vert.second->get_bf({dim0,dim1})->get_coeff();
                        space_ctr++;
                        if (space_ctr != pow(2,_no_dims)-1) {
                            f << " ";
                        };
                    };
                };

            } else if (_no_dims == 3) {

                for (auto dim0: dim_types_possible) {
                    for (auto dim1: dim_types_possible) {
                        for (auto dim2: dim_types_possible) {
                            f << vert.second->get_bf({dim0,dim1,dim2})->get_coeff();
                            space_ctr++;
                            if (space_ctr != pow(2,_no_dims)-1) {
                                f << " ";
                            };
                        };
                    };
                };

            } else if (_no_dims == 4) {
                
                for (auto dim0: dim_types_possible) {
                    for (auto dim1: dim_types_possible) {
                        for (auto dim2: dim_types_possible) {
                            for (auto dim3: dim_types_possible) {
                                f << vert.second->get_bf({dim0,dim1,dim2,dim3})->get_coeff();
                                space_ctr++;
                                if (space_ctr != pow(2,_no_dims)-1) {
                                    f << " ";
                                };
                            };
                        };
                    };
                };
                
            } else if (_no_dims == 5) {
                
                for (auto dim0: dim_types_possible) {
                    for (auto dim1: dim_types_possible) {
                        for (auto dim2: dim_types_possible) {
                            for (auto dim3: dim_types_possible) {
                                for (auto dim4: dim_types_possible) {
                                    f << vert.second->get_bf({dim0,dim1,dim2,dim3,dim4})->get_coeff();
                                    space_ctr++;
                                    if (space_ctr != pow(2,_no_dims)-1) {
                                        f << " ";
                                    };
                                };
                            };
                        };
                    };
                };
                
            } else if (_no_dims == 6) {
                
                for (auto dim0: dim_types_possible) {
                    for (auto dim1: dim_types_possible) {
                        for (auto dim2: dim_types_possible) {
                            for (auto dim3: dim_types_possible) {
                                for (auto dim4: dim_types_possible) {
                                    for (auto dim5: dim_types_possible) {
                                        f << vert.second->get_bf({dim0,dim1,dim2,dim3,dim4,dim5})->get_coeff();
                                        space_ctr++;
                                        if (space_ctr != pow(2,_no_dims)-1) {
                                            f << " ";
                                        };
                                    };
                                };
                            };
                        };
                    };
                };
            };
		};

		// Close
		f.close();
	};

}
