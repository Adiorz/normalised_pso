#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#ifndef ADDITIONAL_H_
#define ADDITIONAL_H_

void donothing(std::vector<double> *realdata, std::vector<double> &fitnesses,
		std::vector<std::vector<double> > &X);

static bool abs_compare(double a, double b) {
	return (std::abs(a) > std::abs(b));
}

void init_minmax(std::vector<double> &xmin, std::vector<double> &xmax,
		size_t numofdims, std::vector<double> &data, size_t fs = -1);
void get_frequencies(std::vector<double> &freq, size_t numofsamples, double fs);
void get_times(std::vector<double> &times, size_t numofsamples, double ts);

void calc_response(std::vector<std::vector<double>> results,
		size_t numofsamples, double ts, std::vector<double> &response);
void approximate_amp(std::vector< std::vector<double> > factors,
		std::vector<double> &time_real, size_t numofsamples, double *out);

size_t find_max_02_id(double max, std::vector<double> *realdata);

double min(double *values, size_t size);

void findMinimas(std::vector<double> &v, size_t start, size_t end,
		std::vector<size_t> &idx);

void findMinima(std::vector<double> &v, size_t maxIDx, size_t &idxL,
		size_t &idxR);

std::vector<double> gaussian_filter(std::vector<double> &input, size_t sigma);

bool should_skip(size_t f, std::vector<std::pair<size_t, size_t>> &to_skip);

std::string doubleToText(const double & d);

std::vector<size_t> get_left_middle_right(std::vector<double>& A_gauss, size_t numofsamples_2, size_t f);

double get_signal_energy(std::vector<double> signal, size_t start = 1, size_t end = 0);

void prepare_log_file_for_visualisation(std::string log_file_name,
		size_t num_workers, std::vector< std::vector<double> > factors,
		std::set<std::pair<size_t, size_t> > freq_ranges,
		std::vector<double> &time, std::vector<double> &amp,
		std::vector<double> &amp_gauss, size_t numofsamples);

#endif /* ADDITIONAL_H_ */

