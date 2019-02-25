#include <vector>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	enum class DimType: unsigned int {VAL, DERIV};

	// Forwards
	class IdxSet;
	class Vertex;

	/****************************************
	BasisFunc
	****************************************/

	class BasisFunc {

	private:

		// Helpers
		void _clean_up();
		void _copy(const BasisFunc& other);
		void _move(BasisFunc& other);
		
		// Dim types
		std::vector<DimType> _dim_types;

		// No dims
		int _no_dims;

		// Vertex
		Vertex* _vertex;

		// Coeff
		double _coeff;

        // Fixed vals?
        bool _are_vals_fixed;

	public:

		/********************
		Constructor
		********************/

		BasisFunc(Vertex* vertex, std::vector<DimType> dim_types);
		BasisFunc(const BasisFunc& other);
		BasisFunc(BasisFunc&& other);
		BasisFunc& operator=(const BasisFunc &other);
		BasisFunc& operator=(BasisFunc &&other);
		virtual ~BasisFunc();

		/********************
		Dim types
		********************/

		DimType get_dim_type(int dim) const;
		const std::vector<DimType>& get_all_dim_types() const;

		/********************
		Vertex
		********************/

		Vertex* get_vertex() const;

		/********************
		Coeff
		********************/

		double get_coeff() const;
		void set_coeff(double val);
		void increment_coeff(double inc);

        // ***************
        // MARK: - Fix vals
        // ***************
        
        void set_is_val_fixed(bool fixed);
        bool get_is_val_fixed() const;
        
		/********************
		Get val/deriv based on local idx
		********************/
	
		double get_bf_val(const IdxSet& local_idxs, const std::vector<double>& x_frac) const;
		double get_bf_deriv(const IdxSet& local_idxs, const std::vector<double>& x_frac, int deriv_dim) const;

	};

};

