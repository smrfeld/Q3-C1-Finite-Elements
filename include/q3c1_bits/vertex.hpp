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
	class Cell;
	class BasisFunc;
	enum class DimType: unsigned int;
	
	/****************************************
	Vertex
	****************************************/

	class Vertex {

	private:

		// Helpers
		void _clean_up();
		void _copy(const Vertex& other);
		void _move(Vertex& other);
		
		// No dims
		int _no_dims;

		// Idx set
		IdxSet _global_idxs;

		// Abscisses
		std::vector<double> _abscissas;

		// Cells bordered
		std::vector<Cell*> _cells;

		// Basis funcs defined on this vertex
		std::vector<BasisFunc*> _bfs;

	public:

		/********************
		Constructor
		********************/

		Vertex(IdxSet idxs, std::vector<double> abscissas);
		Vertex(const Vertex& other);
		Vertex(Vertex&& other);
		Vertex& operator=(const Vertex &other);
		Vertex& operator=(Vertex &&other);
		~Vertex();

		/********************
		Location
		********************/

		int get_no_dims() const;

		int get_global_idx(int dim) const;
		IdxSet get_global_idxs() const;

		double get_abscissa(int dim) const;
		std::vector<double> get_abscissas() const;

		/********************
		Basis funcs
		********************/

		const std::vector<BasisFunc*>& get_bfs() const;
		BasisFunc* get_bf(const std::vector<DimType> &dim_types) const;
	};
};

