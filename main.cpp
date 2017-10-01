#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>
#include <set>

#include "data_preprosessor.hpp"
#include "scheduler.hpp"
#include "worker.hpp"
#include "pso.hpp"
#include "fft.hpp"
#include "additional.hpp"
#include "barrier.hpp"

#define num_dims 3

// Declarations of static members
Scheduler Worker::scheduler = Scheduler(0);
size_t Worker::num_workers;
size_t Worker::num_iterations;
std::vector<std::vector<std::vector<double> > > Worker::data;
std::vector<size_t> PSO::ids;
std::vector<std::vector<double>> PSO::helper_founds;

int main(int argc, char *argv[]) {
	// get input arguments
	std::string file_name(argv[1]);
	std::cout << "Filename: " << file_name << std::endl;
	std::string log_file_name(argv[2]);
	std::cout << "Log filename: " << log_file_name << std::endl;
	size_t channel = atoi(argv[3]);
	std::cout << "Channel: " << channel << std::endl;

	size_t numofsamples_2 = atoi(argv[4]);
	std::cout << "numofsamples_2: " << numofsamples_2 << std::endl;

	int W = atoi(argv[5]);
	std::cout << "Number of workers (swarms): " << W << std::endl; //modes number
	int I = atoi(argv[6]);
	std::cout << "Number of iterations: " << I << std::endl; //number of iterations for each worker (swarm)
	int P = atoi(argv[7]);
	std::cout << "Number of particles per swarm: " << P << std::endl;


	// pre-process data
	// 0 freqs, 1 hammer, 2-16 15 channels
	//size_t numofsamples_2 = 501;
	size_t numofsamples = (numofsamples_2 - 1) * 2;
	DataStream freq_stream(file_name, 0, numofsamples_2, 0);
	DataStream channel_stream(file_name, channel, numofsamples_2, 0);

	// store captured frequencies and corresponding amplitudes in double vectors
	std::vector<double> freq(freq_stream.get());
	channel_stream.fill_first_elements(10, 0);    // zero first 10 values


	double fs = numofsamples; //5000 is real for the data
	std::vector<double> time(numofsamples);

	get_times(time, numofsamples, 1 / fs);

	std::vector<double> amp(channel_stream.get());
	std::vector<double> amp_gauss(gaussian_filter(amp, 2));

	/////////////////////// simulated signal - star

//    std::vector<std::vector<double>> sim {{10., 200, 0., 0.005},{3, 350, 0., 0.01},{5, 750, 0., 0.02},{1, 700, 0., 0.01}};

    std::vector<std::vector<double>> sim {
		{0.1, 450, 0.015},
		{0.08, 150, 0.009},
		{0.05, 250, 0.006},
		{0.15, 550, 0.03},
		{0.5, 270, 0.03}
	};

    for (size_t i = 0; i < sim.size(); ++i)
    	for (size_t j = 0; j < num_dims; ++j)
    		std::cout << std::fixed << sim[i][j] << std::endl;
    std::vector<double> simulated(numofsamples, 0);
    for (size_t j = 0; j < numofsamples; ++j) {
		double t = time[j];
		simulated[j] += PSO::calc_response(sim, t);
	}

    std::vector<double> sim_A;
	std::vector<double> sim_P;
	fft(simulated, sim_A, sim_P);

	amp = std::vector<double>(sim_A);
	amp_gauss = std::vector<double>(gaussian_filter(amp, 2));
	/////////////////////// simulated signal - end

	double total_e = get_signal_energy(amp);
	std::cout << "Signal e: " << total_e << std::endl;

	// init min/max values
	std::vector<double> xmin, xmax;
	init_minmax(xmin, xmax, num_dims, amp, fs);

	Worker::scheduler = Scheduler(W);
	Worker::num_workers = W;
	Worker::num_iterations = I;
	Worker::data = std::vector<std::vector<std::vector<double> > >(W);
	PSO::helper_founds = std::vector<std::vector<double>>(W);

	for (size_t i = 0; i < W; ++i) {
		Worker::data[i] = std::vector<std::vector<double> >(i + 1);
		for (size_t j = 0; j < i + 1; ++j) {
			Worker::data[i][j] = std::vector<double>(num_dims);
		}
	}

	for (size_t i = 0; i < W; ++i)
		std::cout << Worker::data[i].size() << std::endl;

    std::vector<PSO*> workers = std::vector<PSO*>();
	std::vector<std::thread> t_workers(W);

	for (size_t i = 0; i < W; ++i)
		workers.push_back(
				new PSO(P, num_dims, xmin, xmax, &time, &amp,
						numofsamples, i));

	for (size_t i = 0; i < W; ++i)
		t_workers[i] = std::thread(&PSO::run, std::ref(*workers[i]));

	for (size_t i = 0; i < W; ++i)
		t_workers[i].join();

	std::cout << "Found parameters" << std::endl;
	std::cout << "----" << std::endl;
	for (size_t i = 0; i < workers.size(); ++i) {
		std::cout << workers[i]->get_id() << " ";
		for (double x: workers[i]->getgbest()) {
			std::cout << x << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "----" << std::endl;


	std::vector<std::vector<double> > founds(Worker::data[W - 1]);


	 std::vector<std::vector<double> > xmin_h;
	 std::vector<std::vector<double> > xmax_h;

	//find max
	std::vector<std::vector<double> > data_copy(W);

	size_t f;
	size_t f_l;
	size_t f_r;
	std::vector<size_t> idxL;
	std::vector<size_t> idxR;
	std::vector<std::vector<size_t> > found_ranges;

	std::set<std::pair<size_t, size_t>> myset;
	std::set<std::pair<size_t, size_t>>::iterator it;
	std::pair<std::set<std::pair<size_t, size_t>>::iterator,bool> ret;

	std::vector< std::vector<double> > founds_filtered;

	for (size_t m = 0; m < W; ++m) {
		// determine frequency ranges for each of modes
		f = founds[m][1] + freq[1]; // TODO: moved all results by one freq sample

		std::vector<size_t> l_m_r = get_left_middle_right(amp_gauss, numofsamples_2, f);
		f_l = l_m_r[0];
		f_r = l_m_r[2];
		found_ranges.push_back({f_l, f_r});

		std::cout << "helper " << m << ": <" << freq[f_l] << ", " << f << ", " << freq[f_r] << ">" << std::endl;

		size_t f_range = f_r - f_l;
		std::vector<double>::const_iterator first = amp.begin() + f_l;
		std::vector<double>::const_iterator last = amp.begin() + f_r;

		data_copy[m] = std::vector<double>(first, last);
		std::sort(data_copy[m].begin(), data_copy[m].end(), &abs_compare);
		double max_abs = data_copy[m][0];

		double coeff = 0.2;
		double max_damp = 0.05;
		max_abs = std::abs(max_abs / (exp(-2 * M_PI * max_damp) * sin(2 * M_PI * sqrt(1 - max_damp * max_damp))));
		max_abs *= coeff;
		std::cout << "max: " << max_abs << std::endl;


		ret = myset.insert(std::pair<size_t, size_t>(f_l, f_r));
		if (ret.second) {
			std::cout << "inserting" << std::endl;
			founds_filtered.push_back(founds[m]);

			xmin_h.push_back(std::vector<double> ({0.0, freq[f_l], 0.0001}));
			xmax_h.push_back(std::vector<double> ({max_abs, freq[f_r], max_damp}));
		}
	}

	for (std::pair<size_t, size_t> range: myset)
		std::cout << "<" << range.first << ", " << range.second << ">" << std::endl;

	// Prepare data log for visualisation
//	prepare_log_file_for_visualisation(log_file_name+"_main", W, Worker::data[W - 1], found_ranges, time, amp, amp_gauss, numofsamples);
	prepare_log_file_for_visualisation(log_file_name+"_main", founds_filtered.size(), founds_filtered, myset, time, amp, amp_gauss, numofsamples);
	for (auto w: workers)
		delete w;

	return 0;

	double used = 0.0;
	for (std::pair<size_t, size_t> range: myset) {
		double energy = get_signal_energy(amp, range.first, range.second);
		std::cout << "Energy from: " << range.first << " to: " << range.second  << " is: " << doubleToText(energy) << std::endl;
		used += energy;
	}
	std::cout << "Used energy: " << doubleToText(used) << std::endl;
	std::cout << "Used energy: " << 100.0*used/total_e << "[%]" << std::endl;

	size_t H = founds_filtered.size();
	Worker::scheduler = Scheduler(H);
	Worker::num_workers = H;
	Worker::data = std::vector<std::vector<std::vector<double> > >(Worker::data.begin(), Worker::data.begin()+H);
	Barrier barrier(H);
	PSO::ids = std::vector<size_t> (H);
	for (size_t i = 0; i < H; ++i)
		PSO::ids[i] = i;
//	PSO::helper_founds = std::vector<std::vector<double>>(H);

	std::cout << "Number of helpers: " << H << std::endl;

    std::vector<PSO*> helpers = std::vector<PSO*>();
	std::vector<std::thread> t_helpers(H);

	for (size_t i = 0; i < H; ++i)
		helpers.push_back(
				new PSO(P, num_dims, xmin_h[i], xmax_h[i], &time, &amp,
						numofsamples, i, &barrier, true));

	for (size_t i = 0; i < H; ++i)
		t_helpers[i] = std::thread(&PSO::run, std::ref(*helpers[i]));

	for (size_t i = 0; i < H; ++i)
		t_helpers[i].join();

//	PSO pso(2*P, num_dims, xmin_h[0], xmax_h[0], &time, &amp,
//							numofsamples, 0, true);
//	pso.run();
	std::vector<std::vector<double> > helper_founds;
	for (auto h: helpers) {
		helper_founds.push_back(h->getgbest());
		for (double x: h->getgbest()) {
			std::cout << x << " ";
		}
		std::cout << std::endl;
	}

//	prepare_log_file_for_visualisation(log_file_name, W, std::vector<std::vector<double> > {{pso.getgbest()[0] , pso.getgbest()[1], pso.getgbest()[2]}}, time, amp, amp_gauss, numofsamples);
	prepare_log_file_for_visualisation(log_file_name, H, Worker::data[H - 1], myset, time, amp, amp_gauss, numofsamples);
//	for (auto w: workers)
//		delete w;
	for (auto h: helpers)
		delete h;
	return 0;
}
