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
	class BasisFuncVal;
	class BasisFuncDeriv;

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
		BasisFuncVal *_bf_val;
		BasisFuncDeriv *_bf_deriv;

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

		IdxSet get_global_idxs() const;
		std::vector<double> get_abscissas() const;

		/********************
		Basis funcs
		********************/

		BasisFuncVal* get_bf_val() const;
		BasisFuncDeriv* get_bf_deriv() const;
	};
};

