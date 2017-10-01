#pragma once

#include <vector>
#include <fstream>
#include <sstream>

class DataStream {
public:
	DataStream(std::string file_name, size_t column, size_t number,
			size_t to_skip);
	size_t size();
	const std::vector<double> &get();

	const double operator[](size_t n) const;
	double & operator[](size_t n);

	void fill_first_elements(size_t n, double v);
private:
	std::vector<double> data;
	void parseFile(std::ifstream& file, size_t column, size_t number,
			size_t to_skip);
};
