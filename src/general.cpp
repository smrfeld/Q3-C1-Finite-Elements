#include "../include/q3c1_bits/general.hpp"

#include <iostream>
#include <sstream>
#include <iomanip>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/************************************
	* Show error
	************************************/

	void show_error(std::string class_name, std::string func, std::string message) {

		std::cerr << ">>> Error: " << class_name << "::" << func << " <<<" << std::endl;
		std::cerr << "   " << message << std::endl;
		std::cerr << ">>> Quiting <<<" << std::endl;
	};

	/********************
	Random numbers
	********************/

	double randD(double dMin, double dMax)
	{
	    return dMin + ((double)rand() / RAND_MAX) * (dMax - dMin);
	};

	int randI(int iMin, int iMax)
	{
		return iMin + rand() % (iMax - iMin + 1);
	};

};