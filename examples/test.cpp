#include <q3c1>

#include <iostream>
#include <vector>

using namespace std;
using namespace q3c1;

int main() {

	// Make vertexes
	IdxSet idxs(1);
	idxs[0] = 0;
	Vertex v0(idxs,{0.0});
	idxs[0] = 1;
	Vertex v1(idxs,{0.1});

	// Make basis funcs for each
	BasisFuncVal1D bf0_val(&v0);
	BasisFuncVal1D bf1_val(&v1);
	BasisFuncDeriv1D bf0_deriv(&v0);
	BasisFuncDeriv1D bf1_deriv(&v1);

	// Evaluate
	IdxSet local(1);
	local[0] = 0;
	for (int i=0; i<=10; i++) {
		std::cout << "val bf @ vertex: " << local[0] << " : abscissa: " << i/10.0 << ": " << bf0_val.get_val(local,i/10.0) << std::endl;
	};
	for (int i=0; i<=10; i++) {
		std::cout << "deriv bf @ vertex: " << local[0] << " : abscissa: " << i/10.0 << ": " << bf0_deriv.get_val(local,i/10.0) << std::endl;
	};

	local[0] = 1;
	for (int i=0; i<=10; i++) {
		std::cout << "val bf @ vertex: " << local[0] << " : abscissa: " << i/10.0 << ": " << bf1_val.get_val(local,i/10.0) << std::endl;
	};
	for (int i=0; i<=10; i++) {
		std::cout << "deriv bf @ vertex: " << local[0] << " : abscissa: " << i/10.0 << ": " << bf1_deriv.get_val(local,i/10.0) << std::endl;
	};


	return 0;
};