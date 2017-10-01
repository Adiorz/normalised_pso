#pragma once

#include <vector>
#include "basic_worker.hpp"
#include "scheduler.hpp"

class Worker: public BasicWorker {
	//size_t id;
	size_t counter = 0;

protected:
	static void write_data(size_t i, std::vector<double> d);
	static void read_data(size_t i);

public:
	static size_t num_workers;
	static size_t num_iterations;
	static std::vector<std::vector<std::vector<double> > > data;
	static Scheduler scheduler;

	Worker(size_t id);
	size_t get_id() const {
		return id;
	}
	void write(std::vector<double>);
	void read();
	void run();
};
