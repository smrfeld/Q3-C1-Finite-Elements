#include <q3c1>

#include <iostream>
#include <vector>

using namespace std;
using namespace q3c1;

int main() {

	// Make a dim
	Dimension1D dim(0.0,3.0,5);

	// Make a grid
	Grid grid({&dim});

	// Randomly populate
	IdxSet idx(1);
	std::vector<DimType> dim_types;
	dim_types.push_back(DimType::VAL);
	for (idx[0]=0; idx[0]<5; idx[0]++) {
		dim_types[0] = DimType::VAL;
		grid.get_vertex(idx)->get_bf(dim_types)->set_coeff(randD(-1.0,1.0));
		dim_types[0] = DimType::DERIV;
		grid.get_vertex(idx)->get_bf(dim_types)->set_coeff(randD(-1.0,1.0));
	};

	// Write grid
	grid.write_to_file("test_1d.txt");

	// Try a point
	std::vector<double> pt;
	pt.push_back(2.22);
	std::cout << "Val @ " << pt[0] << " : " << grid.get_val(pt) << std::endl;
	std::cout << "Deriv @ " << pt[0] << " wrt x: " << grid.get_deriv_wrt_abscissa(pt,0) << std::endl;
	// std::cout << "Deriv @ " << pt[0] << " wrt val coeff: " << grid.get_deriv_wrt_coeff(pt,0) << std::endl;

	return 0;
};