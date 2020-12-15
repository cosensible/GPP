#include <iostream>
#include <vector>
#include <ctime>

#include "Solver.h"

using namespace std;

int main() {
	for (int i = 0; i < 10; ++i) {
		Solver solver;
		solver.input();

		solver.init();
		solver.coarsen();
		solver.uncoarsen();
		solver.output(i);
	}
}