#include <vector>
#include <string>

/************************************
* Namespace for q3c1
************************************/

namespace q3c1 {

	/****************************************
	Index set
	****************************************/

	class IdxSet {

	private:

		// Helpers
		void _clean_up();
		void _copy(const IdxSet& other);
		void _move(IdxSet& other);

	protected:
		
		// Idxs
		int _no_idxs;
		int* _idxs;

	public:

		/********************
		Constructor
		********************/

		IdxSet(int no_idxs);
		IdxSet(int no_idxs, int* idxs);
		IdxSet(std::vector<int> idxs);
		IdxSet(const IdxSet& other);
		IdxSet(IdxSet&& other);
		IdxSet& operator=(const IdxSet &other);
		IdxSet& operator=(IdxSet &&other);
		~IdxSet();

		/********************
		Accessors
		********************/

		int operator [](int idx) const;
		int & operator [](int idx);

		int size() const;

		std::string print() const;
	};

	// Printing
    std::ostream& operator<< (std::ostream& stream, const IdxSet& idxs);

    // Comparator
    bool operator==(const IdxSet &lhs, const IdxSet &rhs);
    bool operator<(const IdxSet &lhs, const IdxSet &rhs);

    // Math
    IdxSet operator+(IdxSet lhs, const IdxSet& rhs);
    IdxSet operator-(IdxSet lhs, const IdxSet& rhs);
    IdxSet operator+(IdxSet lhs, int rhs);
    IdxSet operator-(IdxSet lhs, int rhs);
};

