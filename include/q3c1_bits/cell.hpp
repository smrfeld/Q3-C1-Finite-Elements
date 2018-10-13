#include <unordered_map>
#include <vector>

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

		// Vertices
		std::vector<std::pair<IdxSet,Vertex*>> _vertices;

	public:

		/********************
		Constructor
		********************/

		Cell();
		Cell(std::unordered_map<IdxSet,Vertex*> vertex_dic);
		Cell(const Cell& other);
		Cell(Cell&& other);
		Cell& operator=(const Cell &other);
		Cell& operator=(Cell &&other);
		~Cell();

		/********************
		Vertices
		********************/

		void add_vertex(IdxSet local_idxs, Vertex* vertex);

		Vertex* get_vertex(IdxSet local_idxs) const;

		const std::vector<std::pair<IdxSet,Vertex*>>& get_all_vertices() const;
	};
};

