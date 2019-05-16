#include <vector>
#include <map>

#ifndef IDX_SET_H
#define IDX_SET_H
#include "idx_set.hpp"
#endif

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	// Forwards
	class Dimension1D;
	class Vertex;
	class Cell;
	enum class DimType: unsigned int;

	/****************************************
	Grid
	****************************************/

	class Grid {

	private:

		// Helpers
		void _clean_up();
		void _copy(const Grid& other);
		void _move(Grid& other);
		
		// No dims = 1,2,3
		int _no_dims;

		// Dimensions
		std::vector<Dimension1D*> _dims;

		// Vertices
		mutable std::map<IdxSet,Vertex*> _verts;

		// Grids
		mutable std::map<IdxSet,Cell*> _cells;

        // Frac point is between cell boundaries (faster for returning ref)
        mutable IdxSet _idxs_get_cell;
        mutable std::vector<double> _frac_get_cell;
        mutable std::map<int, std::map<DimType, std::map<int, double>>> _vals_get_val;
        
		// Make a cell/vertex (presumes that the cell does NOT exist!)
		Vertex* _make_vertex(IdxSet idxs) const;
		Vertex* _get_or_make_vertex(IdxSet idxs) const;
		Cell* _make_cell(IdxSet idxs) const;

	public:

		/********************
		Constructor
		********************/

		Grid(std::vector<Dimension1D*> dims);
		Grid(const Grid& other);
		Grid(Grid&& other);
		Grid& operator=(const Grid &other);
		Grid& operator=(Grid &&other);
		~Grid();

		/********************
		Print
		********************/

		void print(const std::vector<DimType>& dim_types) const;
		void print(Vertex* vert, const std::vector<DimType>& dim_types) const;

		/********************
		Get dims
		********************/

		int get_no_dims() const;
		Dimension1D* get_dim(int dim) const;
		const std::vector<Dimension1D*>& get_all_dims() const;

		/********************
		Get vertices
		********************/

		// Will make the vertex if needed!
		Vertex* get_vertex(IdxSet idxs) const;
		const std::map<IdxSet,Vertex*>& get_all_vertices() const;

		/********************
		Get cell
		********************/

		// Will make the cell if needed!
		Cell* get_cell(const IdxSet& idxs) const;
		std::pair<Cell*,const std::vector<double>*> get_cell(const std::vector<double>& abscissas) const;
		const std::map<IdxSet,Cell*>& get_all_cells() const;

		/********************
		Get vals
		********************/

		double get_val(const std::vector<double>& abscissas) const;
		double get_deriv_wrt_abscissa(const std::vector<double>& abscissas, int deriv_dim) const;
		double get_deriv_wrt_coeff(const std::vector<double>& abscissas, const IdxSet& global_vertex_idxs, const std::vector<DimType>& dim_types) const;
		std::map<Vertex*,std::vector<double>> get_deriv_wrt_coeffs_for_all_surrounding_verts(const std::vector<double>& abscissas) const;

		/********************
		Read/write grid
		********************/

		void read_from_file(std::string fname);
		void write_to_file(std::string fname) const;
	};
};

