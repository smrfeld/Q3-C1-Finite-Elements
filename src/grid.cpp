#include "../include/q3c1_bits/grid.hpp"

#include "../include/q3c1_bits/vertex.hpp"
#include "../include/q3c1_bits/dimension_1d.hpp"
#include "../include/q3c1_bits/cell.hpp"
#include "../include/q3c1_bits/general.hpp"
#include "../include/q3c1_bits/basis_func.hpp"

#include <iostream>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/****************************************
	Grid
	****************************************/

	// Constructors
	Grid::Grid(std::vector<Dimension1D*> dims) {
		if (dims.size() > 3 || dims.size() == 0) {
			show_error("Grid","Grid","Only dimensions 1,2,3 supported");
			exit(EXIT_FAILURE);
		};

		_no_dims = dims.size();
		_dims = dims;
		for (auto const &dim: dims) {
			_no_verts_in_each_dim.push_back(dim->get_no_pts());
			_no_cells_in_each_dim.push_back(dim->get_no_pts()-1);
		};

		// Make all the verts, cells
		std::vector<double> abscissas;
		std::vector<std::pair<IdxSet,Vertex*>> verts_local;
		if (_no_dims == 1) {
			// Dim = 1

			// Verts
			IdxSet idxs(1);
			abscissas.push_back(0.0);
			for (idxs[0]=0; idxs[0]<_dims[0]->get_no_pts(); idxs[0]++) {
				abscissas[0] = _dims[0]->get_pt_by_idx(idxs[0]);
				_verts.push_back(new Vertex(idxs, abscissas));
			};

			// Cells
			IdxSet idxs_local(1);
			IdxSet idxs_find(1);
			for (idxs[0]=0; idxs[0]<_dims[0]->get_no_pts()-1; idxs[0]++) {

				// Surrounding verts
				verts_local.clear();
				idxs_find = idxs;

				idxs_find[0] = idxs[0];
				idxs_local[0] = 0;
				verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
				idxs_find[0] = idxs[0]+1;
				idxs_local[0] = 1;
				verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));

				// Make
				_cells.push_back(new Cell(idxs, verts_local));
			};


		} else if (_no_dims == 2) {
			// Dim = 2

			// Verts
			IdxSet idxs(2);
			abscissas.push_back(0.0);
			abscissas.push_back(0.0);
			for (idxs[0]=0; idxs[0]<_dims[0]->get_no_pts(); idxs[0]++) {
				for (idxs[1]=0; idxs[1]<_dims[1]->get_no_pts(); idxs[1]++) {
					abscissas[0] = _dims[0]->get_pt_by_idx(idxs[0]);
					abscissas[1] = _dims[1]->get_pt_by_idx(idxs[1]);
					_verts.push_back(new Vertex(idxs, abscissas));
				};
			};

			// Cells
			IdxSet idxs_local(2);
			IdxSet idxs_find(2);
			for (idxs[0]=0; idxs[0]<_dims[0]->get_no_pts()-1; idxs[0]++) {
				for (idxs[1]=0; idxs[1]<_dims[1]->get_no_pts()-1; idxs[1]++) {

					// Surrounding verts
					verts_local.clear();
					idxs_find = idxs;

					idxs_find[0] = idxs[0];
					idxs_find[1] = idxs[1];
					idxs_local[0] = 0;
					idxs_local[1] = 0;
					verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
					idxs_find[0] = idxs[0]+1;
					idxs_find[1] = idxs[1];
					idxs_local[0] = 1;
					idxs_local[1] = 0;
					verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
					idxs_find[0] = idxs[0];
					idxs_find[1] = idxs[1]+1;
					idxs_local[0] = 0;
					idxs_local[1] = 1;
					verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
					idxs_find[0] = idxs[0]+1;
					idxs_find[1] = idxs[1]+1;
					idxs_local[0] = 1;
					idxs_local[1] = 1;
					verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));

					// Make
					_cells.push_back(new Cell(idxs, verts_local));
				};
			};

		} else if (_no_dims == 3) {
			// Dim = 3

			// Verts
			IdxSet idxs(3);
			abscissas.push_back(0.0);
			abscissas.push_back(0.0);
			abscissas.push_back(0.0);
			for (idxs[0]=0; idxs[0]<_dims[0]->get_no_pts(); idxs[0]++) {
				for (idxs[1]=0; idxs[1]<_dims[1]->get_no_pts(); idxs[1]++) {
					for (idxs[2]=0; idxs[2]<_dims[2]->get_no_pts(); idxs[2]++) {
						abscissas[0] = _dims[0]->get_pt_by_idx(idxs[0]);
						abscissas[1] = _dims[1]->get_pt_by_idx(idxs[1]);
						abscissas[2] = _dims[2]->get_pt_by_idx(idxs[2]);
						_verts.push_back(new Vertex(idxs, abscissas));
					};
				};
			};

			// Cells
			IdxSet idxs_local(3);
			IdxSet idxs_find(3);
			for (idxs[0]=0; idxs[0]<_dims[0]->get_no_pts()-1; idxs[0]++) {
				for (idxs[1]=0; idxs[1]<_dims[1]->get_no_pts()-1; idxs[1]++) {
					for (idxs[2]=0; idxs[2]<_dims[2]->get_no_pts()-1; idxs[2]++) {

						// Surrounding verts
						verts_local.clear();
						idxs_find = idxs;

						idxs_find[0] = idxs[0];
						idxs_find[1] = idxs[1];
						idxs_find[2] = idxs[2];
						idxs_local[0] = 0;
						idxs_local[1] = 0;
						idxs_local[2] = 0;
						verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
						idxs_find[0] = idxs[0]+1;
						idxs_find[1] = idxs[1];
						idxs_find[2] = idxs[2];
						idxs_local[0] = 1;
						idxs_local[1] = 0;
						idxs_local[2] = 0;
						verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
						idxs_find[0] = idxs[0];
						idxs_find[1] = idxs[1]+1;
						idxs_find[2] = idxs[2];
						idxs_local[0] = 0;
						idxs_local[1] = 1;
						idxs_local[2] = 0;
						verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
						idxs_find[0] = idxs[0];
						idxs_find[1] = idxs[1];
						idxs_find[2] = idxs[2]+1;
						idxs_local[0] = 0;
						idxs_local[1] = 0;
						idxs_local[2] = 1;
						verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
						idxs_find[0] = idxs[0]+1;
						idxs_find[1] = idxs[1]+1;
						idxs_find[2] = idxs[2];
						idxs_local[0] = 1;
						idxs_local[1] = 1;
						idxs_local[2] = 0;
						verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
						idxs_find[0] = idxs[0]+1;
						idxs_find[1] = idxs[1];
						idxs_find[2] = idxs[2]+1;
						idxs_local[0] = 1;
						idxs_local[1] = 0;
						idxs_local[2] = 1;
						verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
						idxs_find[0] = idxs[0];
						idxs_find[1] = idxs[1]+1;
						idxs_find[2] = idxs[2]+1;
						idxs_local[0] = 0;
						idxs_local[1] = 1;
						idxs_local[2] = 1;
						verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));
						idxs_find[0] = idxs[0]+1;
						idxs_find[1] = idxs[1]+1;
						idxs_find[2] = idxs[2]+1;
						idxs_local[0] = 1;
						idxs_local[1] = 1;
						idxs_local[2] = 1;
						verts_local.push_back(std::make_pair(idxs_local,get_vertex(idxs_find)));

						// Make
						_cells.push_back(new Cell(idxs, verts_local));
					};
				};
			};
		};
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
		// Nothing
	};
	void Grid::_copy(const Grid& other)
	{
		_no_dims = other._no_dims;
		_dims = other._dims;
		_no_verts_in_each_dim = other._no_verts_in_each_dim;
		_no_cells_in_each_dim = other._no_cells_in_each_dim;
		for (auto &v: other._verts) {
			_verts.push_back(new Vertex(*v));
		};
		for (auto &c: other._cells) {
			_cells.push_back(new Cell(*c));
		};
	};
	void Grid::_move(Grid& other)
	{
		_no_dims = other._no_dims;
		_dims = other._dims;
		_no_verts_in_each_dim = other._no_verts_in_each_dim;
		_no_cells_in_each_dim = other._no_cells_in_each_dim;
		_verts = other._verts;
		_cells = other._cells;

		// Reset other
		other._no_dims = 0;
		other._dims.clear();
		other._no_verts_in_each_dim.clear();
		other._no_cells_in_each_dim.clear();
		other._verts.clear();
		other._cells.clear();
	};

	/********************
	Print
	********************/

	void Grid::print() const {
		// ..
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
		int idx=0;
		int term;
		for (auto dim=0; dim<_no_dims; dim++) {
			
			term = idxs[dim];
			for (auto dim2=dim+1; dim2<_no_dims; dim2++) {
				term *= _no_verts_in_each_dim[dim2];
			};

			idx += term;
		};
		return _verts[idx];
	};
	
	/********************
	Get cell
	********************/

	Cell* Grid::get_cell(IdxSet idxs) const {
		int idx=0;
		int term;
		for (auto dim=0; dim<_no_dims; dim++) {
			
			term = idxs[dim];
			for (auto dim2=dim+1; dim2<_no_dims; dim2++) {
				term *= _no_cells_in_each_dim[dim2];
			};

			idx += term;
		};
		return _cells[idx];
	};

	std::pair<Cell*,std::vector<double>> Grid::get_cell(const std::vector<double>& abscissas) const {
		std::vector<double> frac;
		if (_no_dims == 1) {
			IdxSet idxs(1);
			idxs[0] = _dims[0]->get_idxs_surrounding_pt(abscissas[0]);
			frac.push_back(_dims[0]->get_frac_between(abscissas[0],idxs[0]));
			return std::make_pair(get_cell(idxs),frac);
		} else if (_no_dims == 2) {
			IdxSet idxs(2);
			idxs[0] = _dims[0]->get_idxs_surrounding_pt(abscissas[0]);
			idxs[1] = _dims[1]->get_idxs_surrounding_pt(abscissas[1]);
			frac.push_back(_dims[0]->get_frac_between(abscissas[0],idxs[0]));
			frac.push_back(_dims[1]->get_frac_between(abscissas[1],idxs[1]));
			return std::make_pair(get_cell(idxs),frac);
		} else if (_no_dims == 3) {
			IdxSet idxs(3);
			idxs[0] = _dims[0]->get_idxs_surrounding_pt(abscissas[0]);
			idxs[1] = _dims[1]->get_idxs_surrounding_pt(abscissas[1]);
			idxs[2] = _dims[2]->get_idxs_surrounding_pt(abscissas[2]);
			frac.push_back(_dims[0]->get_frac_between(abscissas[0],idxs[0]));
			frac.push_back(_dims[1]->get_frac_between(abscissas[1],idxs[1]));
			frac.push_back(_dims[2]->get_frac_between(abscissas[2],idxs[2]));
			return std::make_pair(get_cell(idxs),frac);
		} else {
			return std::make_pair(nullptr,frac);
		};
	};


	/********************
	Get vals
	********************/

	double Grid::get_val(const std::vector<double>& abscissas) const {

		// Get cell and fraction of this abscissa in the cell
		std::pair<Cell*,std::vector<double>> pr = get_cell(abscissas);

		if (pr.first == nullptr) {
			show_error("Grid","get_val","Cell not found");
			exit(EXIT_FAILURE);
		};

		double val = 0.0;

		// Run through all verts of the cell
		for (auto const &v_pr: pr.first->get_all_vertices()) {
			// Run through all bfs defined on this vertex
			for (auto const &bf: v_pr.second->get_bfs()) {
				// val += bf->get_coeff() * bf->
				val += bf->get_bf_val(v_pr.first, pr.second);
			};
		};

		return val;
	};
	double Grid::get_deriv_wrt_abscissa(const std::vector<double>& abscissas, int deriv_dim) {


		// Get cell and fraction of this abscissa in the cell
		std::pair<Cell*,std::vector<double>> pr = get_cell(abscissas);

		if (pr.first == nullptr) {
			show_error("Grid","get_val","Cell not found");
			exit(EXIT_FAILURE);
		};

		double val = 0.0;

		// Run through all verts of the cell
		for (auto const &v_pr: pr.first->get_all_vertices()) {
			// Run through all bfs defined on this vertex
			for (auto const &bf: v_pr.second->get_bfs()) {
				// val += bf->get_coeff() * bf->
				val += bf->get_bf_deriv(v_pr.first, pr.second, deriv_dim);
			};
		};

		return val;

	};
	double Grid::get_deriv_wrt_coeff(const IdxSet& vertex_idxs, const std::vector<DimType>& dim_types) {

		
		
	};

}