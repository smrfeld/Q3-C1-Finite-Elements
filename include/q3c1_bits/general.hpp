#include <string>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/************************************
	* Show error
	************************************/

	void show_error(std::string class_name, std::string func, std::string message);

	/********************
	Zero pad a string
	********************/

	std::string pad_str(int i, int n_zeros);
	
	/********************
	Random numbers
	********************/

	double randD(double dMin, double dMax); // inclusive
	int randI(int iMin, int iMax); // inclusive

};