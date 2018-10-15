#include <vector>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	// Forwards
	class Dimension1D;
	class IdxSet;
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

		// No x in each dim
		std::vector<int> _no_verts_in_each_dim;
		std::vector<int> _no_cells_in_each_dim;

		// Vertices
		std::vector<Vertex*> _verts;

		// Grids
		std::vector<Cell*> _cells;

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

		void print() const;

		/********************
		Get dims
		********************/

		int get_no_dims() const;
		Dimension1D* get_dim(int dim) const;
		const std::vector<Dimension1D*>& get_all_dims() const;

		/********************
		Get vertices
		********************/

		Vertex* get_vertex(IdxSet idxs) const;
		
		/********************
		Get cell
		********************/

		Cell* get_cell(IdxSet idxs) const;
		std::pair<Cell*,std::vector<double>> get_cell(const std::vector<double>& abscissas) const;

		/********************
		Get vals
		********************/

		double get_val(const std::vector<double>& abscissas) const;
		double get_deriv_wrt_abscissa(const std::vector<double>& abscissas, int deriv_dim);
		double get_deriv_wrt_coeff(const std::vector<double>& abscissas, const IdxSet& vertex_idxs, const std::vector<DimType>& dim_types);

	};
};

