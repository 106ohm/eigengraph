/*
 * File: matrix_sketo.cpp
 * ----------------------
 * A matrix addition performed with Sketo.
 */

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sketo/matrix_skeletons.h>

struct gen_t : public sketo::functions::base<double (sketo::mindex)> {
	double operator()(sketo::mindex) const {
		return ((double) rand() / (double) RAND_MAX);
	}
} gen;

int sketo::main(int argc, char *argv[]) {
	size_t rows, cols;

	if (argc != 3) {
		// Using sketo::cout since Sketo has no sketo::cerr.
		sketo::cout << "Usage: " << argv[0] << " rows cols" << std::endl;
		return -1;
	} else {
		rows = atoi(argv[1]);
		cols = atoi(argv[2]);
	}	

	dist_matrix<double> A = sketo::matrix_skeletons::generate(msize(rows, cols), gen);
	dist_matrix<double> B = sketo::matrix_skeletons::generate(msize(rows, cols), gen);

	std::clock_t begin = std::clock();
	dist_matrix<double> C = sketo::matrix_skeletons::zipwith(std::plus<double>(), A, B);
	std::clock_t end = std::clock();

	sketo::cout << "Time elapsed with Sketo: " << ((double) (end - begin)) / CLOCKS_PER_SEC << "s" << std::endl;

	return 0;
}

