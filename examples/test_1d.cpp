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

	// Try a point
	std::vector<double> pt;
	pt.push_back(2.22);
	std::cout << "Val @ " << pt[0] << " : " << grid.get_val(pt) << std::endl;
	std::cout << "Deriv @ " << pt[0] << " wrt x: " << grid.get_deriv_wrt_abscissa(pt,0) << std::endl;
	// std::cout << "Deriv @ " << pt[0] << " wrt val coeff: " << grid.get_deriv_wrt_coeff(pt,0) << std::endl;

	return 0;
};