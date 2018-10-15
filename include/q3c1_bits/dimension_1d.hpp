/************************************
* Namespace for dcu
************************************/

namespace dcu {

	/****************************************
	Dimension1D
	****************************************/

	class Dimension1D {

	private:

		double _min_pt;
		double _max_pt;
		double _delta;
		int _no_pts;

		// Internal copy func/clean up
		void _clean_up();
		void _reset();
		void _copy(const Dimension1D& other);

	public:

		Dimension1D(double min_pt_pt, double max_pt_pt, int no_pts);
		Dimension1D(const Dimension1D& other);
		Dimension1D(Dimension1D&& other);
		Dimension1D& operator=(const Dimension1D& other);
		Dimension1D& operator=(Dimension1D&& other);
		~Dimension1D();

		// Accessors
		int get_no_pts() const;
		double get_min_pt() const;
		double get_max_pt() const;
		double get_delta() const;

		// Get by idx
		double get_pt_by_idx(int idx, bool start_at_one=false) const;

		// Check if point is in domain
		bool check_if_pt_is_inside_domain(double x) const;

		// Get closest index
		int get_closest_idx(double x, bool start_at_one=false) const;

		// Get indexes surrounding a point
		// ie point is between i and i+1 where i is returned
		int get_idxs_surrounding_pt(double x, bool start_at_one=false) const; 

		// Get fraction of a point between successive points
		double get_frac_between(double x) const;
		// Second optional specification: the return of the surrounding idxs
		double get_frac_between(double x, int i, bool start_at_one=false) const;
	};
};