/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

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
		
		// Vertex
		Vertex* _vertex;

		// Coeff
		double _coeff;

	public:

		/********************
		Constructor
		********************/

		BasisFunc(Vertex* vertex);
		BasisFunc(const BasisFunc& other);
		BasisFunc(BasisFunc&& other);
		BasisFunc& operator=(const BasisFunc &other);
		BasisFunc& operator=(BasisFunc &&other);
		~BasisFunc();

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
	};



	/****************************************
	BasisFuncVal
	--- ABSTRACT BASE ---	
	****************************************/

	class BasisFuncVal : public BasisFunc {

	public:

		/********************
		Constructor
		********************/

		BasisFuncVal(Vertex* vertex);
		BasisFuncVal(const BasisFuncVal& other);
		BasisFuncVal(BasisFuncVal&& other);
		BasisFuncVal& operator=(const BasisFuncVal &other);
		BasisFuncVal& operator=(BasisFuncVal &&other);
		virtual ~BasisFuncVal();

		/********************
		Get val/deriv based on local idx
		--- CAUTION: PURE VIRTUAL ---
		********************/

		virtual double get_val(IdxSet local_idxs, double x_frac) = 0;
		virtual double get_deriv(IdxSet local_idxs, double x_frac) = 0;

	};



	/****************************************
	BasisFuncDeriv
	--- ABSTRACT BASE ---	
	****************************************/

	class BasisFuncDeriv : public BasisFunc {

	public:

		/********************
		Constructor
		********************/

		BasisFuncDeriv(Vertex* vertex);
		BasisFuncDeriv(const BasisFuncDeriv& other);
		BasisFuncDeriv(BasisFuncDeriv&& other);
		BasisFuncDeriv& operator=(const BasisFuncDeriv &other);
		BasisFuncDeriv& operator=(BasisFuncDeriv &&other);
		virtual ~BasisFuncDeriv();
		
		/********************
		Get val/deriv based on local idx
		--- CAUTION: PURE VIRTUAL ---
		********************/

		virtual double get_val(IdxSet local_idxs, double x_frac) = 0;
		virtual double get_deriv(IdxSet local_idxs, double x_frac) = 0;

	};


	/****************************************
	1D
	****************************************/


	/****************************************
	BasisFuncVal1D
	****************************************/

	class BasisFuncVal1D : public BasisFuncVal {

	public:

		/********************
		Constructor
		********************/

		using BasisFuncVal::BasisFuncVal;

		/********************
		Get val/deriv based on local idx
		********************/

		double get_val(IdxSet local_idxs, double x_frac);
		double get_deriv(IdxSet local_idxs, double x_frac);

	};


	/****************************************
	BasisFuncDeriv1D
	****************************************/

	class BasisFuncDeriv1D : public BasisFuncDeriv {

	public:

		/********************
		Constructor
		********************/

		using BasisFuncDeriv::BasisFuncDeriv;

		/********************
		Get val/deriv based on local idx
		********************/

		double get_val(IdxSet local_idxs, double x_frac);
		double get_deriv(IdxSet local_idxs, double x_frac);

	};
};

