#include <unordered_map>
#include <vector>

#ifndef IDX_SET_H
#define IDX_SET_H
#include "idx_set.hpp"
#endif

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	// Forwards
	class IdxSet;
	class Vertex;

	/****************************************
	Cell
	****************************************/

	class Cell {

	private:

		// Helpers
		void _clean_up();
		void _copy(const Cell& other);
		void _move(Cell& other);
		
		// No dims
		int _no_dims;

		// Idxs
		IdxSet _idxs;

		// Vertices
		std::vector<std::pair<IdxSet,Vertex*>> _verts;

	public:

		/********************
		Constructor
		********************/

		Cell(IdxSet idxs, std::vector<std::pair<IdxSet,Vertex*>>vert_dict);
		Cell(const Cell& other);
		Cell(Cell&& other);
		Cell& operator=(const Cell &other);
		Cell& operator=(Cell &&other);
		~Cell();

		/********************
		Vertices
		********************/

		Vertex* get_vertex(IdxSet local_idxs) const;

		const std::vector<std::pair<IdxSet,Vertex*>>& get_all_vertices() const;
	};
};

