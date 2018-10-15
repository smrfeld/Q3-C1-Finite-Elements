#include "../include/q3c1_bits/cell.hpp"

#include "../include/q3c1_bits/vertex.hpp"
#include "../include/q3c1_bits/general.hpp"

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/****************************************
	Cell
	****************************************/

	// Constructors
	Cell::Cell(IdxSet idxs, std::vector<std::pair<IdxSet,Vertex*>> verts) : _idxs(idxs) {
		if (_idxs.size() > 3 || _idxs.size() == 0) {
			show_error("Cell","Cell","Only dimensions 1,2,3 supported");
			exit(EXIT_FAILURE);		
		};

		_no_dims = _idxs.size();
		_verts = verts;
	};
	Cell::Cell(const Cell& other) : _idxs(other._idxs) {
		_copy(other);
	};
	Cell::Cell(Cell&& other) : _idxs(std::move(other._idxs)) {
		_move(other);
	};
    Cell& Cell::operator=(const Cell& other) {
		if (this != &other) {
			_clean_up();
			_copy(other);
		};
		return *this;
    };
    Cell& Cell::operator=(Cell&& other) {
		if (this != &other) {
			_clean_up();
			_move(other);
		};
		return *this;
    };
	Cell::~Cell()
	{
		_clean_up();
	};

	// Helpers
	void Cell::_clean_up()
	{
		// Nothing
	};
	void Cell::_copy(const Cell& other)
	{
		_no_dims = other._no_dims;
		_idxs = other._idxs;
		_verts = other._verts;
	};
	void Cell::_move(Cell& other)
	{
		_no_dims = other._no_dims;
		_idxs = other._idxs;
		_verts = other._verts;

		// Reset other
		other._no_dims = 0;
		other._idxs = IdxSet(0);
		other._verts.clear();
	};

	/********************
	Vertices
	********************/

	Vertex* Cell::get_vertex(IdxSet local_idxs) const {
		int idx=0;
		int term;
		for (auto dim=0; dim<_no_dims; dim++) {
			
			term = local_idxs[dim];
			for (auto dim2=dim+1; dim2<_no_dims; dim2++) {
				term *= 2;
			};

			idx += term;
		};
		return _verts[idx].second;
	};

	IdxSet Cell::get_local_idxs_of_vertex(Vertex *vert) const {
		for (auto &pr: _verts) {
			if (pr.second == vert) {
				return pr.first;
			};
		};

		// Never get here
		show_error("Cell","get_local_idxs_of_vertex","Could not find vertex");
		exit(EXIT_FAILURE);
	};

	const std::vector<std::pair<IdxSet,Vertex*>>& Cell::get_all_vertices() const {
		return _verts;
	};


};