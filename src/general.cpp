#include "../include/q3c1_bits/idx_set.hpp"

#include <iostream>
#include <sstream>

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

}