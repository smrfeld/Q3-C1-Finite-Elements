#include "../include/q3c1_bits/grid.hpp"

#include "../include/q3c1_bits/vertex.hpp"
#include "../include/q3c1_bits/dimension_1d.hpp"
#include "../include/q3c1_bits/cell.hpp"
#include "../include/q3c1_bits/general.hpp"
#include "../include/q3c1_bits/basis_func.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/****************************************
	Grid
	****************************************/

	// Constructors
	Grid::Grid(std::vector<Dimension1D*> dims) {
		if (dims.size() > 3 || dims.size() < 1) {
			std::cerr << ">>> Grid::Grid <<< Error: only 1,2,3 dims supported" << std::endl;
			exit(EXIT_FAILURE);
		};

		_no_dims = dims.size();
		_dims = dims;
	};
	Grid::Grid(const Grid& other) {
		_copy(other);
	};
	Grid::Grid(Grid&& other) {
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
	};
	void Grid::_move(Grid& other)
	{
		_no_dims = other._no_dims;
		_dims = other._dims;
		_verts = other._verts;
		_cells = other._cells;

		// Reset other
		other._no_dims = 0;
		other._dims.clear();
		other._verts.clear();
		other._cells.clear();
	};

	/********************
	Make a cell
	********************/

	Vertex* Grid::_make_vertex(IdxSet idxs) const {		
		std::vector<double> abscissas(_no_dims);
		for (auto dim=0; dim<_no_dims; dim++) {
			abscissas[dim] = _dims[dim]->get_pt_by_idx(idxs[dim]);
		};

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

			// 0
			idxs_local[0] = 0;
			idxs_vert[0] = idxs[0];
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 1
			idxs_local[0] = 1;
			idxs_vert[0] = idxs[0]+1;
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

		} else if (_no_dims == 2) {

			// 0,0
			idxs_local[0] = 0;
			idxs_local[1] = 0;
			idxs_vert[0] = idxs[0];
			idxs_vert[1] = idxs[1];
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 0,1
			idxs_local[0] = 0;
			idxs_local[1] = 1;
			idxs_vert[0] = idxs[0];
			idxs_vert[1] = idxs[1]+1;
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 1,0
			idxs_local[0] = 1;
			idxs_local[1] = 0;
			idxs_vert[0] = idxs[0]+1;
			idxs_vert[1] = idxs[1];
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 1,1
			idxs_local[0] = 1;
			idxs_local[1] = 1;
			idxs_vert[0] = idxs[0]+1;
			idxs_vert[1] = idxs[1]+1;
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

		} else if (_no_dims == 3) {

			// 0,0,0
			idxs_local[0] = 0;
			idxs_local[1] = 0;
			idxs_local[2] = 0;
			idxs_vert[0] = idxs[0];
			idxs_vert[1] = idxs[1];
			idxs_vert[2] = idxs[2];
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 0,0,1
			idxs_local[0] = 0;
			idxs_local[1] = 0;
			idxs_local[2] = 1;
			idxs_vert[0] = idxs[0];
			idxs_vert[1] = idxs[1];
			idxs_vert[2] = idxs[2]+1;
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 0,1,0
			idxs_local[0] = 0;
			idxs_local[1] = 1;
			idxs_local[2] = 0;
			idxs_vert[0] = idxs[0];
			idxs_vert[1] = idxs[1]+1;
			idxs_vert[2] = idxs[2];
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 1,0,0
			idxs_local[0] = 1;
			idxs_local[1] = 0;
			idxs_local[2] = 0;
			idxs_vert[0] = idxs[0]+1;
			idxs_vert[1] = idxs[1];
			idxs_vert[2] = idxs[2];
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 0,1,1
			idxs_local[0] = 0;
			idxs_local[1] = 1;
			idxs_local[2] = 1;
			idxs_vert[0] = idxs[0];
			idxs_vert[1] = idxs[1]+1;
			idxs_vert[2] = idxs[2]+1;
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 1,0,1
			idxs_local[0] = 1;
			idxs_local[1] = 0;
			idxs_local[2] = 1;
			idxs_vert[0] = idxs[0]+1;
			idxs_vert[1] = idxs[1];
			idxs_vert[2] = idxs[2]+1;
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 1,1,0
			idxs_local[0] = 1;
			idxs_local[1] = 1;
			idxs_local[2] = 0;
			idxs_vert[0] = idxs[0]+1;
			idxs_vert[1] = idxs[1]+1;
			idxs_vert[2] = idxs[2];
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));

			// 1,1,1
			idxs_local[0] = 1;
			idxs_local[1] = 1;
			idxs_local[2] = 1;
			idxs_vert[0] = idxs[0]+1;
			idxs_vert[1] = idxs[1]+1;
			idxs_vert[2] = idxs[2]+1;
			verts_local.push_back(std::make_pair(idxs_local,_get_or_make_vertex(idxs_vert)));			
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
		if (_no_dims == 1) {
			std::cout << vert->get_abscissa(0) << ": " << std::flush;
		} else if (_no_dims == 2) {
			std::cout << vert->get_abscissa(0) << " " << vert->get_abscissa(1) << ": "<< std::flush;
		} else if (_no_dims == 3) {
			std::cout << vert->get_abscissa(0) << " " << vert->get_abscissa(1) << " "<< vert->get_abscissa(2) << ": " << std::flush;
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

	Cell* Grid::get_cell(IdxSet idxs) const {
		auto it = _cells.find(idxs);
		if (it != _cells.end()) {
			return it->second;
		} else {
			return _make_cell(idxs);
		};
	};

	std::pair<Cell*,std::vector<double>> Grid::get_cell(const std::vector<double>& abscissas) const {
		std::vector<double> frac;
		IdxSet idxs(_no_dims);
		for (auto dim=0; dim<_no_dims; dim++) {
			idxs[dim] = _dims[dim]->get_idxs_surrounding_pt(abscissas[dim]);
			frac.push_back(_dims[dim]->get_frac_between(abscissas[dim],idxs[dim]));
		};
		return std::make_pair(get_cell(idxs),frac);
	};

	const std::map<IdxSet,Cell*>& Grid::get_all_cells() const {
		return _cells;
	};

	/********************
	Get vals
	********************/

	double Grid::get_val(const std::vector<double>& abscissas) const {

		// Get cell and fraction of this abscissa in the cell
		std::pair<Cell*,std::vector<double>> pr = get_cell(abscissas);

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
		std::pair<Cell*,std::vector<double>> pr = get_cell(abscissas);

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
		std::pair<Cell*,std::vector<double>> pr = get_cell(abscissas);

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
		std::pair<Cell*,std::vector<double>> pr_cell = get_cell(abscissas);

		// Returned
		std::map<Vertex*,std::vector<double>> ret;

		// Go through all verts of the cell
		if (_no_dims == 1) {
			for (auto &pr_v: pr_cell.first->get_all_vertices()) {
				// First are the local idxs; second is the vertex itself
				// Val
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::VAL})->get_bf_val(pr_v.first,pr_cell.second));
				// Deriv
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::DERIV})->get_bf_val(pr_v.first,pr_cell.second));
			};
		} else if (_no_dims == 2) {
			for (auto &pr_v: pr_cell.first->get_all_vertices()) {
				// First are the local idxs; second is the vertex itself
				// Val-val
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::VAL,DimType::VAL})->get_bf_val(pr_v.first,pr_cell.second));
				// Val-deriv
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::VAL,DimType::DERIV})->get_bf_val(pr_v.first,pr_cell.second));
				// Deriv-val
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::DERIV,DimType::VAL})->get_bf_val(pr_v.first,pr_cell.second));
				// Deriv-deriv
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::DERIV,DimType::DERIV})->get_bf_val(pr_v.first,pr_cell.second));
			};
		} else if (_no_dims == 3) {
			for (auto &pr_v: pr_cell.first->get_all_vertices()) {
				// First are the local idxs; second is the vertex itself
				// Val-val-val
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::VAL,DimType::VAL,DimType::VAL})->get_bf_val(pr_v.first,pr_cell.second));
				// Val-val-deriv
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::VAL,DimType::VAL,DimType::DERIV})->get_bf_val(pr_v.first,pr_cell.second));
				// Val-deriv-val
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::VAL,DimType::DERIV,DimType::VAL})->get_bf_val(pr_v.first,pr_cell.second));
				// Deriv-val-val
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::DERIV,DimType::VAL,DimType::VAL})->get_bf_val(pr_v.first,pr_cell.second));
				// Val-deriv-deriv
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::VAL,DimType::DERIV,DimType::DERIV})->get_bf_val(pr_v.first,pr_cell.second));
				// Deriv-val-deriv
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::DERIV,DimType::VAL,DimType::DERIV})->get_bf_val(pr_v.first,pr_cell.second));
				// Deriv-deriv-val
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::DERIV,DimType::DERIV,DimType::VAL})->get_bf_val(pr_v.first,pr_cell.second));
				// Deriv-deriv-deriv
				ret[pr_v.second].push_back(pr_v.second->get_bf({DimType::DERIV,DimType::DERIV,DimType::DERIV})->get_bf_val(pr_v.first,pr_cell.second));
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
				// Val
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::VAL})->set_coeff(ordinate);

				// Deriv
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::DERIV})->set_coeff(ordinate);
			} else if (_no_dims == 2) {
				// Val-val
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::VAL,DimType::VAL})->set_coeff(ordinate);

				// Val-deriv
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::VAL,DimType::DERIV})->set_coeff(ordinate);

				// Deriv-val
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::DERIV,DimType::VAL})->set_coeff(ordinate);

				// Deriv-deriv
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::DERIV,DimType::DERIV})->set_coeff(ordinate);
			} else if (_no_dims == 3) {
				// Val-val-val
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::VAL,DimType::VAL,DimType::VAL})->set_coeff(ordinate);

				// Val-val-deriv
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::VAL,DimType::VAL,DimType::DERIV})->set_coeff(ordinate);

				// Val-deriv-val
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::VAL,DimType::DERIV,DimType::VAL})->set_coeff(ordinate);

				// Deriv-val-val
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::DERIV,DimType::VAL,DimType::VAL})->set_coeff(ordinate);

				// Val-deriv-deriv
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::VAL,DimType::DERIV,DimType::DERIV})->set_coeff(ordinate);

				// Deriv-val-deriv
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::DERIV,DimType::VAL,DimType::DERIV})->set_coeff(ordinate);

				// Deriv-deriv-val
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::DERIV,DimType::DERIV,DimType::VAL})->set_coeff(ordinate);

				// Deriv-deriv-deriv
				ordinate_str = "";
				iss >> ordinate_str;
				ordinate = atof(ordinate_str.c_str());
				vert->get_bf({DimType::DERIV,DimType::DERIV,DimType::DERIV})->set_coeff(ordinate);
			};			
		};

		// Close
		f.close();
	};
	void Grid::write_to_file(std::string fname) const {
		std::ofstream f;

		// Open
		f.open(fname);

		// Make sure we found it
		if (!f.is_open()) {
			show_error("Grid","write_to_file","Could not open file");
			exit(EXIT_FAILURE);
		};

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
			if (_no_dims == 1) {
				// Val
				f << vert.second->get_bf({DimType::VAL})->get_coeff() << " ";
				// Deriv
				f << vert.second->get_bf({DimType::DERIV})->get_coeff();
			} else if (_no_dims == 2) {
				// Val-val
				f << vert.second->get_bf({DimType::VAL,DimType::VAL})->get_coeff() << " ";
				// Val-deriv
				f << vert.second->get_bf({DimType::VAL,DimType::DERIV})->get_coeff() << " ";
				// Deriv-val
				f << vert.second->get_bf({DimType::DERIV,DimType::VAL})->get_coeff() << " ";
				// Deriv-deriv
				f << vert.second->get_bf({DimType::DERIV,DimType::DERIV})->get_coeff();
			} else if (_no_dims == 3) {
				// Val-val-val
				f << vert.second->get_bf({DimType::VAL,DimType::VAL,DimType::VAL})->get_coeff() << " ";
				// Val-val-deriv
				f << vert.second->get_bf({DimType::VAL,DimType::VAL,DimType::DERIV})->get_coeff() << " ";
				// Val-deriv-val
				f << vert.second->get_bf({DimType::VAL,DimType::DERIV,DimType::VAL})->get_coeff() << " ";
				// Deriv-val-val
				f << vert.second->get_bf({DimType::DERIV,DimType::VAL,DimType::VAL})->get_coeff() << " ";
				// Val-deriv-deriv
				f << vert.second->get_bf({DimType::VAL,DimType::DERIV,DimType::DERIV})->get_coeff() << " ";
				// Deriv-val-deriv
				f << vert.second->get_bf({DimType::DERIV,DimType::VAL,DimType::DERIV})->get_coeff() << " ";
				// Deriv-deriv-val
				f << vert.second->get_bf({DimType::DERIV,DimType::DERIV,DimType::VAL})->get_coeff() << " ";
				// Deriv-deriv-deriv
				f << vert.second->get_bf({DimType::DERIV,DimType::DERIV,DimType::DERIV})->get_coeff();
			};		
		};

		// Close
		f.close();
	};

}
