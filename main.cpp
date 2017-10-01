#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>

#include "data_preprosessor.hpp"
#include "scheduler.hpp"
#include "worker.hpp"
#include "pso.hpp"
#include "fft.hpp"
#include "additional.hpp"

#define num_dims 3
#define numofparticles 50

// Declarations of static members
Scheduler Worker::scheduler = Scheduler(0);
size_t Worker::num_workers;
size_t Worker::num_iterations;
std::vector<std::vector<std::vector<double> > > Worker::data;

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

	// pre-process data
	// 0 freqs, 1 hammer, 2-16 15 channels
	//size_t numofsamples_2 = 501;
	size_t numofsamples = (numofsamples_2 - 1) * 2;
	DataStream freq_stream(file_name, 0, numofsamples_2, 0);
	DataStream channel_stream(file_name, channel, numofsamples_2, 0);

	// store captured frequencies and corresponding amplitudes in double vectors
	std::vector<double> freq(freq_stream.get());
	channel_stream.fill_first_elements(10, 0);    // zero first 10 values
	std::vector<double> amp(channel_stream.get());
	std::vector<double> amp_gauss(gaussian_filter(amp, 1));

	double fs = numofsamples; //5000 is real for the data
	std::vector<double> time(numofsamples);

	get_times(time, numofsamples, 1 / fs);

	// init min/max values
	std::vector<double> xmin, xmax;
	init_minmax(xmin, xmax, num_dims, amp, fs);

	Worker::scheduler = Scheduler(W);
	Worker::num_workers = W;
	Worker::num_iterations = I;
	Worker::data = std::vector<std::vector<std::vector<double> > >(W);

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
				new PSO(numofparticles, num_dims, xmin, xmax, &time, &amp,
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

	// Prepare data log for visualisation
	prepare_log_file_for_visualisation(log_file_name, W, Worker::data[W - 1], time, amp, amp_gauss, numofsamples);

	for (auto w: workers)
		delete w;
	return 0;
}
