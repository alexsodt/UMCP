#define __globalsc__
#include "globals.h"
#include "parallel.h"
#include "input.h"
#include <stdlib.h>

parameterBlock *global_block = NULL;
double kT = 0.592;
int global_debug_mode = 0;
parallel_info par_info; 
double lipid_DC = 1e10;
double solution_DC = 1e10;
