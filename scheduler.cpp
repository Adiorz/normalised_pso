#include <chrono>
#include <iostream>

#include "worker.hpp"
#include "scheduler.hpp"

using namespace std::chrono_literals;

Scheduler::Scheduler(size_t num_schedulees) { //number of processes
	this->num_schedulees = num_schedulees;
	ms = std::vector < std::mutex > (num_schedulees);
	cvs = std::vector < std::condition_variable > (num_schedulees);
	readys = std::vector<bool>(num_schedulees, false);
	processeds = std::vector<bool>(num_schedulees, false);
}

Scheduler::Scheduler(): Scheduler(0) {}

void Scheduler::send_data_to_worker(BasicWorker *basic_worker, std::vector<double> d) {
	size_t i = basic_worker->get_id();
	{
		std::lock_guard < std::mutex > lk(ms[i]); //send my data

		//write data to bufer
		((Worker*)basic_worker)->write(d);
		//after writing

		readys[i] = true; //my data
	}
	cvs[i].notify_one(); //notify worker//i+1
}

void Scheduler::wait_for_worker_to_read(BasicWorker *basic_worker) { //i<T-1
	{
		size_t i = basic_worker->get_id();
		std::unique_lock < std::mutex > lk(ms[i]); // my data
		if (!cvs[i].wait_for(lk, 15*1000ms, [this, i] {return processeds[i];})) //I (i) am waiting for my data to be read
			std::cerr << i << "timed out: processeds[" << i << "] = " << processeds[i] << std::endl;
		processeds[i] = false;
	}
}

void Scheduler::wait_for_data_and_read(BasicWorker *basic_worker) { //for i>0
	size_t i = basic_worker->get_id();
	std::unique_lock < std::mutex > lk(ms[i - 1]); // data held by to i-1
	cvs[i - 1].wait(lk, [this, i] {return readys[i - 1];}); //I (i) am waiting for data held by to i-1
	readys[i - 1] = false;

	//read data
	((Worker*)basic_worker)->read();
	//after read

	processeds[i - 1] = true; //read data belonging to i-1

	lk.unlock();
	cvs[i - 1].notify_one(); //notify main i-1
}

