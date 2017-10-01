#include <algorithm>
#include <random>
#include <chrono>
#include <complex>
#include <map>
#include <iostream>

#include <fstream>

#include <thread>
#include <iomanip>

#include "pso.hpp"
#include "fft.hpp"
#include "additional.hpp"
#include "scheduler.hpp"

#define rand_01 ((double)rand() / (double)RAND_MAX)

std::mutex m;

double PSO::init_param(size_t j) {
	if (j == 1)
		return PSO::distribution(generator) * PSO::distribution(generator)
				* (Xmax.at(j) - Xmin.at(j)) + Xmin.at(j);
	if (j == 2)
		return PSO::distribution(generator) * PSO::distribution(generator)
				* PSO::distribution(generator) * (Xmax.at(j) - Xmin.at(j))
				+ Xmin.at(j);
	return PSO::distribution(generator) * (Xmax.at(j) - Xmin.at(j)) + Xmin.at(j);
}

void PSO::init_max_velocities() {
	for (size_t i = 0; i < numofdims; i++) {
		Vmax.at(i) = 0.2 * (Xmax.at(i) - Xmin.at(i));
		Vmin.at(i) = -Vmax.at(i);
	}
}

void PSO::initpopulation() {
	for (size_t i = 0; i < numofparticles; i++) {
		X.at(i) = std::vector<double>(numofdims);
		V.at(i) = std::vector<double>(numofdims);
		pbests.at(i) = std::vector<double>(numofdims);
		for (size_t j = 0; j < numofdims; j++) {
			X.at(i).at(j) = init_param(j);
		}
	}
}

void PSO::addparticle() {
	numofparticles++;
	numofsamples_2 = size_t(numofsamples / 2) + 1;
	std::vector<double> new_x(numofdims);
	std::vector<double> new_v(numofdims, 0);
	for (size_t i = 0; i < numofdims; ++i) {
		new_x[i] = init_param(i);
	}
	X.push_back(new_x);
	V.push_back(new_v);
	pbests.push_back(new_x);
	double fitness = fitnessfunc_singleparticle(X[numofparticles - 1]);
	fitnesses.push_back(fitness);
	pbestfits.push_back(fitness);
}

void PSO::reposition_particle(size_t i) {
	std::vector<double> new_x(numofdims);
	std::vector<double> new_v(numofdims, 0);
	for (size_t i = 0; i < numofdims; ++i) {
		new_x[i] = init_param(i);
	}
	X[i] = std::vector<double> (new_x);
	V[i] = std::vector<double> (new_v);

	double fitness = fitnessfunc_singleparticle(X[i]);
	fitnesses[i] = fitness;
	if (fitness < pbestfits[i]) {
		pbestfits[i] = fitness;
		pbests[i] = std::vector<double> (new_x);
	}
}

double PSO::fitnessfunc_singleparticle(std::vector<double> &p) {
	std::vector<double> response(numofsamples, 0);

	//Calc response

	for (size_t j = 0; j < numofsamples; ++j) {
		double t = time->at(j);
		response[j] += calc_response(std::vector<std::vector<double>> { p }, t);
	}

	std::vector<double> A;
	std::vector<double> P;
	fft(response, A, P);
	double total_e = 0.;
	for (size_t i = 0; i < this->A->size(); ++i)
		total_e += this->A->at(i) * this->A->at(i);

	//Calc fitness
	//TODO: work on this first
	double fitness = 0.;
	for (size_t j = 0; j < numofsamples_2; ++j) {
		double residue = this->A->at(j) - A[j];
		residue /= total_e;
		double temp_fit = residue * residue;// * residue * residue;
		if (id > 0) {
			size_t l_skip, m_skip, h_skip;
			bool b_skip = false;
			for (size_t i = 0; i < to_skip.size(); ++i) {
				l_skip = to_skip[i][0];
				h_skip = to_skip[i][2];
				if (j >= l_skip && j <= h_skip && p[1] >= l_skip
						&& p[1] <= h_skip) {
					b_skip = true;
					m_skip = to_skip[i][1];
					break;
				}
			}
			if (b_skip) {
				double penalty = 0.;
				if (j < m_skip) {
					penalty = (j - l_skip) / (double) (m_skip - l_skip);
				} else {
					penalty = (h_skip - j) / (double) (h_skip - m_skip);
				}
				double r_temp = this->A->at(j) / total_e * (1.f + penalty);
				temp_fit += r_temp * r_temp;// * r_temp * r_temp;
			}
		}
		fitness += temp_fit;
	}
	return fitness;
}

double PSO::calc_response(std::vector<std::vector<double>> results, double t,
		size_t id) {
	double response = 0.;
	for (size_t i = 0; i < results.size(); ++i) {
		double amp = results[i][0];
		double freq = results[i][1];
		double damp = results[i][2];

		response += amp * exp(-2 * M_PI * freq * damp * t)
				* sin(2 * M_PI * freq * sqrt(1 - damp * damp) * t);
	}
	return response;
}

void PSO::fitnessfunc() {
	for (size_t p = 0; p < numofparticles; p++) {
		fitnesses[p] = fitnessfunc_singleparticle(X[p]);
	}
}

void PSO::calcgbest(bool first) {
	sorted_indices = std::vector<size_t>(fitnesses.size()); // indices of sorted fitnesses (from min to max)
	std::size_t n(0);
	std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return n++; });
	std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](double d1, double d2) { return fitnesses[d1] < fitnesses[d2]; } );

	size_t minfitidx = sorted_indices[0];
	double minfit = fitnesses[minfitidx];

	// recalculate the best particle of lower rank swarms using new info from higher rank swarms
	if (id > 0) {
		gbestfit = fitnessfunc_singleparticle(gbest);
	}
	// did minfit improved?
	if (first || minfit < gbestfit) {
		gbestfit = minfit;
		//TODO: change for fast vector copying
		for (size_t i = 0; i < numofdims; ++i)
			gbest.at(i) = X.at(minfitidx).at(i);
	}
}

void PSO::update() {
	for (size_t i = 0; i < numofparticles; i++) {
		for (size_t j = 0; j < numofdims; j++) {
			// update velocity
			V.at(i).at(j) = update_velocity(w, X.at(i).at(j), V.at(i).at(j),
					Vmin[j], Vmax[j], gbest.at(j), pbests.at(i).at(j), c1, c2);
			// update position
			X.at(i).at(j) = update_position(X.at(i).at(j), V.at(i).at(j),
					Xmin[j], Xmax[j]);
		}
	}
}

double PSO::update_velocity(double w, double X, double V, double V_min,
		double V_max, double gbest, double pbests, double c1, double c2) {
	return std::min(
			std::max(
					(w * V + rand_01 * c1 * (pbests - X)
							+ rand_01 * c2 * (gbest - X)), V_min), V_max);
}

double PSO::update_position(double X, double V, double X_min, double X_max) {
	return std::min(std::max((X + V), X_min), X_max);
}

PSO::PSO(size_t numofparticles, size_t numofdims, std::vector<double> Xmin,
		std::vector<double> Xmax, std::vector<double> *time,
		std::vector<double> *A, size_t numofsamples, size_t id, double c1, double c2) :
		Worker(id) {

	this->id = id;
	this->numofparticles = numofparticles;
	this->numofdims = numofdims;
	this->numofsamples = numofsamples;
	this->numofsamples_2 = size_t(numofsamples / 2) + 1;
	this->c1 = c1;
	this->c2 = c2;

	V = std::vector<std::vector<double> >(numofparticles);
	X = std::vector<std::vector<double> >(numofparticles);
	Vmax = std::vector<double>(numofdims);
	Vmin = std::vector<double>(numofdims);
	pbests = std::vector<std::vector<double> >(numofparticles);
	pbestfits = std::vector<double>(numofparticles);
	fitnesses = std::vector<double>(numofparticles);

	gbest = std::vector<double>(numofdims);
	bests = std::vector<double>(num_iterations);

	this->time = time;
	this->A = A;

	to_skip = std::vector<std::vector<size_t>>();

	this->Xmax = Xmax;
	this->Xmin = Xmin;

	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();
	myclock::duration d = myclock::now() - beginning;
	unsigned seed = d.count();

	generator = std::default_random_engine(seed);
	distribution = std::uniform_real_distribution<double>(0.0, 1.0);

	A_gauss = gaussian_filter(*A, 1);

	init_max_velocities();
	initpopulation();
}

void PSO::run() {
	if (!initialized) {
		fitnessfunc();
		calcgbest(true);
		initialized = true;
	}
	std::map<size_t, bool> prog;
	for (size_t t = 0; t < num_iterations; t++) {
		size_t progress = t * 100 / num_iterations;
		if ((progress % 10 == 0) && (!prog[progress])) {
			std::cout << "id: " << id << ": " << progress << ": " << std::fixed
					<< std::setprecision(17) << gbestfit << std::setprecision(5)
					<< std::endl;
			prog[progress] = true;
		}
		if (num_iterations)
			w = 0.9 - 0.7 * t / num_iterations;

		for (size_t i = 0; i < numofparticles; i++) {
			if (fitnesses[i] < pbestfits[i]) {
				pbestfits[i] = fitnesses[i];
				for (size_t j = 0; j < numofdims; j++)
					pbests[i][j] = X[i][j];
			}
		}
		update();

		to_skip.clear();
		if (id > 0) {
			scheduler.wait_for_data_and_read(this);

			for (size_t i = 0; i < id; ++i) {
				size_t f = data[id][i][1]; // 1 is freq
				std::vector<size_t> idxL;
				std::vector<size_t> idxR;
				findMinimas(A_gauss, 0, f, idxL);
				findMinimas(A_gauss, f, numofsamples_2 - 1, idxR);
				if (idxR.size() == 0)
					idxR.push_back(numofsamples_2 - 1);
				to_skip.push_back(std::vector<size_t> { idxL[idxL.size() - 1],
						f, idxR[0] });
			}
		}

		fitnessfunc();
		calcgbest();
		bests[t] = gbestfit;


		if (progress % 10 == 0 && t != 0) {
			size_t num_worst = size_t(numofparticles/(100/10)); // 10%
			std::vector<size_t> worst_indices(&sorted_indices[numofparticles-1-num_worst], &sorted_indices[numofparticles-1]);
			for (size_t i: worst_indices) {
				reposition_particle(i);
			}
		}

		scheduler.send_data_to_worker(this, std::vector<double>(getgbest()));
		if (id < num_workers - 1)
			scheduler.wait_for_worker_to_read(this);
	}
}

std::vector<double> PSO::getgbest() const {
	return gbest;
}

double PSO::getgbestfit() {
	return gbestfit;
}
