#include <q3c1>

#include <iostream>
#include <vector>

using namespace std;
using namespace q3c1;

int main() {

	// Make vertexes
	std::vector<Vertex> vertices;
	IdxSet idxs(1);
	for (idxs[0]=0; idxs[0]<2; idxs[0]++) {
		vertices.push_back(Vertex(idxs,{idxs[0]*0.1}));
	};

	// Make basis funcs for each
	std::vector<BasisFunc> bfs;
	for (auto &v: vertices) {
		// Value
		bfs.push_back(BasisFunc(&v,{DimType::VAL}));
		// Deriv
		bfs.push_back(BasisFunc(&v,{DimType::DERIV}));
	};

	// Evaluate
	IdxSet local(1);

	for (auto const &bf: bfs) {
		std::cout << "--- Basis function on vertex: " << bf.get_vertex()->get_abscissa(0) << " of type: ";
		if (bf.get_dim_type(0) == DimType::VAL) {
			std::cout << "VAL ---" << std::endl;
		} else {
			std::cout << "DERIV ---" << std::endl;
		};
		for (double x=0.0; x<=1.0; x+=0.1) {
			if (bf.get_vertex()->get_global_idx(0) == 0) { 
				local[0] = 0;
			} else {
				local[0] = 1;
			};
			std::cout << bf.get_val(local,{x}) << std::endl;
		};
	};

	return 0;
};