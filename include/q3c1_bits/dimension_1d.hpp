/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/****************************************
	Dimension1D
	****************************************/

	class Dimension1D {

	private:

		double _delta;
		double _zero;

		// Internal copy func/clean up
		void _clean_up();
		void _move(Dimension1D& other);
		void _copy(const Dimension1D& other);

	public:

		Dimension1D(double delta, double zero);
		Dimension1D(const Dimension1D& other);
		Dimension1D(Dimension1D&& other);
		Dimension1D& operator=(const Dimension1D& other);
		Dimension1D& operator=(Dimension1D&& other);
		~Dimension1D();

		// Accessors
		double get_delta() const;
		double get_zero() const;

		// Get by idx
		double get_pt_by_idx(int idx) const;

		// Get closest index
		int get_closest_idx(double x) const;

		// Get indexes surrounding a point
		// ie point is between i and i+1 where i is returned
		int get_idxs_surrounding_pt(double x) const; 

		// Get fraction of a point between successive points
		double get_frac_between(double x) const;
		// Second optional specification: the return of the surrounding idxs
		double get_frac_between(double x, int i) const;
	};
};