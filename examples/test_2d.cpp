#include <q3c1>

#include <iostream>
#include <vector>

using namespace std;
using namespace q3c1;

int main() {

	// Make vertexes
	std::vector<Vertex> vertices;
	IdxSet idxs(2);
	for (idxs[0]=0; idxs[0]<2; idxs[0]++) {
		for (idxs[1]=0; idxs[1]<2; idxs[1]++) {
			vertices.push_back(Vertex(idxs,{idxs[0]*0.1,idxs[1]*0.1}));
		};
	};

	// Make basis funcs for each
	std::vector<BasisFunc> bfs;
	for (auto &v: vertices) {
		// Value-value
		bfs.push_back(BasisFunc(&v,{DimType::VAL,DimType::VAL}));
		// Value-deriv
		bfs.push_back(BasisFunc(&v,{DimType::VAL,DimType::DERIV}));
		// Deriv-value
		bfs.push_back(BasisFunc(&v,{DimType::DERIV,DimType::VAL}));
		// Deriv-deriv
		bfs.push_back(BasisFunc(&v,{DimType::DERIV,DimType::DERIV}));
	};

	// Evaluate
	IdxSet local(2);
	for (auto const &bf: bfs) {
		std::cout << "--- Basis function on vertex: (" << bf.get_vertex()->get_abscissa(0) << "," << bf.get_vertex()->get_abscissa(1) << ") of type: ";
		if (bf.get_dim_type(0) == DimType::VAL && bf.get_dim_type(1) == DimType::VAL) {
			std::cout << "VAL-VAL ---" << std::endl;
		} else if (bf.get_dim_type(0) == DimType::VAL && bf.get_dim_type(1) == DimType::DERIV) {
			std::cout << "VAL-DERIV ---" << std::endl;
		} else if (bf.get_dim_type(0) == DimType::DERIV && bf.get_dim_type(1) == DimType::VAL) {
			std::cout << "DERIV-VAL ---" << std::endl;
		} else {
			std::cout << "DERIV-DERIV ---" << std::endl;
		};

		for (double x=0.0; x<=1.0; x+=0.1) {
			for (double y=0.0; y<=1.0; y+=0.1) {
				if (bf.get_vertex()->get_global_idx(0) == 0) { 
					local[0] = 0;
				} else {
					local[0] = 1;
				};
				if (bf.get_vertex()->get_global_idx(1) == 0) { 
					local[1] = 0;
				} else {
					local[1] = 1;
				};
				std::cout << "(" << x << "," << y << "): " << bf.get_val(local,{x,y}) << std::endl;
			};
		};
	};

	return 0;
};