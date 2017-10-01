#include <iostream>
#include "data_preprosessor.hpp"

DataStream::DataStream(std::string file_name, size_t column, size_t number,
		size_t to_skip) {
	std::ifstream data_stream(file_name);
	if (data_stream)
		parseFile(data_stream, column, number, to_skip);
	else
		std::cout << "Cannot open the " << file_name << " file" << std::endl;
}
size_t DataStream::size() {
	return data.size();
}
const std::vector<double> &DataStream::get() {
	return data;
}
const double DataStream::operator[](size_t n) const {
	return data[n];
}
double & DataStream::operator[](size_t n) {
	return data[n];
}

void DataStream::parseFile(std::ifstream& file, size_t column, size_t number,
		size_t to_skip) {
	size_t read_number = 0;
	size_t skipped = 0;
	for (std::string line; std::getline(file, line);) {
		if (read_number >= number)
			break;
		std::istringstream in(line);
		double piece; // of data
		size_t counter = 0;
		while (in >> piece) {
			if (counter == (column)) {
				if (skipped < to_skip)
					++skipped;
				else {
					data.push_back(piece);
					read_number++;
				}
			}
			counter++;
		}
	}
}

void DataStream::fill_first_elements(size_t n, double v) {
	std::fill(data.begin(), data.begin() + n, v);
}
