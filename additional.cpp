//#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "additional.hpp"
#include "fft.hpp"
#include "pso.hpp"

void init_minmax(std::vector<double> &xmin, std::vector<double> &xmax,
		size_t numofdims, std::vector<double> &data, size_t fs) {

	double coeff = 2.0;
	double max_dump = 0.05;

	xmin = std::vector<double>(numofdims);
	xmax = std::vector<double>(numofdims);
	//find max
	// nie działa!
	std::vector<double> data_copy(data.size());
	std::copy(data.begin(), data.end(), data_copy.begin());
	std::sort(data_copy.begin(), data_copy.end(), &abs_compare);
	double max_abs = data_copy[0];
	xmin[1] = 20; //freq
	xmax[1] = fs / 2;
	xmin[2] = 0.000001; //damping
	xmax[2] = max_dump;

	max_abs = std::abs(
			max_abs
					/ (exp(-2 * M_PI * xmax[2])
							* sin(2 * M_PI * sqrt(1 - xmax[2] * xmax[2]))));
	max_abs *= coeff;
	std::cout << "max_abs: " << max_abs << std::endl;
	xmin[0] = 0.0; //amp
	xmax[0] = max_abs;
}

void get_frequencies(std::vector<double> &freq, size_t numofsamples,
		double fs) {
	size_t numofsamples_2 = (size_t) (numofsamples / 2) + 1;

	freq = std::vector<double>(numofsamples_2);

	for (size_t i = 0; i < freq.size(); ++i)
		freq[i] = i * fs / (numofsamples_2 - 1) / 2;
}

void get_times(std::vector<double> &times, size_t numofsamples, double ts) {
	times = std::vector<double>(numofsamples);
	std::generate(times.begin(), times.begin() + numofsamples,
			[ts] () mutable ->double {static double value; return value += ts;});
}

void calc_response(std::vector<std::vector<double>> results,
		size_t numofsamples, double ts, std::vector<double> &response) {
	response = std::vector<double>(numofsamples, 0.);
	for (size_t i = 0; i < results.size(); ++i) {
		double amp = results[i][0];
		double freq = results[i][1];
		double damp = results[i][2];
		for (size_t j = 0; j < numofsamples; ++j) {
			double t = j * ts;
			response[j] += amp
					* sin(2 * M_PI * freq * sqrt(1 - damp * damp) * t)
					* exp(-2 * M_PI * freq * damp * t);
		}
	}
}

void approximate_amp(std::vector<std::vector<double>> factors,
		std::vector<double> &time_real, size_t numofsamples, double *out) {
	std::vector<double> appr(numofsamples, 0.0);
	for (size_t i = 0; i < factors.size(); ++i) {
		double amp = factors[i][0];
		double freq = factors[i][1];
		double damp = factors[i][2];
		for (size_t j = 0; j < numofsamples; ++j) {
			appr[j] += amp * exp(-2 * M_PI * freq * damp * time_real[j])
					* sin(
							2 * M_PI * freq * sqrt(1 - damp * damp)
									* time_real[j]);
		}
	}

	std::vector<double> A_appr;
	std::vector<double> P_appr;
	fft(appr, A_appr, P_appr);
	size_t numofsamples_2 = (size_t) (numofsamples / 2) + 1;
	for (size_t i = 0; i < numofsamples_2; ++i)
		out[i] = A_appr[i];
}

double min(double *values, size_t size) {
	size_t counter = 0;
	double min = values[0];
	while (std::isinf(min)) {
		min = values[counter];
		++counter;
	}
	for (size_t i = counter; i < size; ++i) {
		if (values[i] < min and !std::isinf(values[i]))
			min = values[i];
	}
	return min;
}

void findMinimas(std::vector<double> &v, size_t start, size_t end,
		std::vector<size_t> &idx) {
	idx = std::vector<size_t>();
	for (unsigned int i = start + 1; i < end; ++i) {
		if (v[i] < v[i - 1]) {
			unsigned int j = i;
			while (v[j] == v[j + 1]) {
				++j;
			}
			++j;
			if (v[j] > v[i]) {
				idx.push_back(i);
				i = j;
			}
		}
	}
	if (start == 0 && idx.size() == 0)
		idx.push_back(0);
}

std::vector<double> gaussian_filter(std::vector<double> &input, size_t sigma) {
	size_t n = 6 * sigma;
	std::vector<double> kernel(n);
	double half_n = (n - 1) / 2.0;
	double sum = 0.0;
	for (size_t i = 0; i < n; ++i) {
		double x = i - half_n;
		kernel[i] = 1.0 / sqrtf(2 * M_PI * sigma * sigma)
				* exp(-x * x / (2 * sigma * sigma));
		sum += kernel[i];
	}
	std::vector<double> output = std::vector<double>(input.size());
	size_t cols = input.size();
	size_t kCols = kernel.size();
	size_t kCenterX = kernel.size() / 2;
	for (size_t j = 0; j < cols; ++j)          // columns
			{
		for (size_t n = 0; n < kCols; ++n) // kernel columns
				{
			size_t nn = kCols - 1 - n; // column void prepare_log_file_for_visualisation()index of flipped kernel
			// index of input signal, used for checking boundary
			size_t jj = j + (n - kCenterX);

			// ignore input samples which are out of bound
			if (jj >= 0 && jj < cols)
				output[j] += input[jj] * kernel[nn];
		}
	}
	return output;
}

bool should_skip(size_t f, std::vector<std::pair<size_t, size_t>> &to_skip) {
	for (size_t i = 0; i < to_skip.size(); ++i)
		if (f >= to_skip[i].first && f <= to_skip[i].second)
			return true;
	return false;
}

std::string doubleToText(const double & d) {
	std::stringstream ss;
	ss << std::setprecision(std::numeric_limits<double>::max_digits10);
	ss << d;
	return ss.str();
}

std::vector<size_t> get_left_middle_right(std::vector<double>& A_gauss,
		size_t numofsamples_2, size_t f) {
	std::vector<size_t> idxL;
	std::vector<size_t> idxR;
	findMinimas(A_gauss, 0, f, idxL);
	findMinimas(A_gauss, f, numofsamples_2 - 1, idxR);
	if (idxR.size() == 0)
		idxR.push_back(numofsamples_2 - 1);
	return std::vector<size_t> { idxL[idxL.size() - 1] + 1, f, idxR[0] };
}

double get_signal_energy(std::vector<double> signal, size_t start, size_t end) {
	if (start == 1 && end == 0) {
		start = 0;
		end = signal.size();
	}
	//signal energy
//	return std::inner_product(signal.begin() + start,
//			signal.end() + (end - start), signal.begin() + start, 0.0);
	return std::accumulate(signal.begin(), signal.begin() + (end - start), 0.0, // start with first element
			[](double x, double y) {return x + y*y;});

	//RMS
//	return std::sqrt(
//			(std::inner_product(signal.begin() + start,
//					signal.end() + (end - start), signal.begin() + start, 0.0))
//					/ static_cast<double>(end - start));

}

void prepare_log_file_for_visualisation(std::string log_file_name,
		size_t num_workers, std::vector<std::vector<double> > factors,
		std::set<std::pair<size_t, size_t> > freq_ranges,
		std::vector<double> &time, std::vector<double> &amp,
		std::vector<double> &amp_gauss, size_t numofsamples) {

	size_t numofsamples_2 = size_t(numofsamples / 2) + 1;

	std::ofstream myfile;

	std::vector<std::vector<double>> results;
	std::vector<double> response(numofsamples);

	std::vector<double *> individual_responses(num_workers);
	for (size_t i = 0; i < num_workers; ++i) {
		individual_responses[i] = new double[numofsamples_2];
		approximate_amp(std::vector<std::vector<double> > { factors[i] }, time,
				numofsamples, individual_responses[i]);
	}

	for (size_t i = 0; i < num_workers; ++i) {
//		results.push_back({factors[i][0], factors[i][1], factors[i][2]});
		results.push_back(factors[i]);
	}
	for (size_t j = 0; j < numofsamples; ++j) {
		double t = time[j];
		response[j] += PSO::calc_response(results, t);
	}
	std::vector<double> A_found;
	std::vector<double> P;

	fft(response, A_found, P);

	myfile.open(log_file_name, std::ios::out);

	size_t num_dims = factors[0].size();

	myfile << num_workers << " " << num_dims << std::endl;

	std::set<std::pair<size_t, size_t>>::iterator it = freq_ranges.begin();
//	size_t w = 0;
	for (std::vector<double> factor : factors) {
		for (double x : factor) {
			myfile << doubleToText(x) << " ";
		}
		myfile << doubleToText((*it).first) << " " << doubleToText((*it).second)
				<< " ";
		myfile << std::endl;
		std::advance(it, 1);
//		++w;
	}

	for (size_t i = 0; i < numofsamples_2; ++i) {
		double real_A = amp[i];
		double found = A_found[i];
		myfile << i << ";" << real_A << ";" << found;
		for (size_t h = 0; h < num_workers; ++h)
			myfile << ";" << individual_responses[h][i];
		myfile << ";" << amp_gauss[i];
		myfile << std::endl;
	}
	myfile.close();

	for (auto r : individual_responses)
		delete[] r;
}

