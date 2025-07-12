#pragma once

#include "global_variables.h"

void build_MILP_min_sum(instance* inst);

void solve_MILP_min_sum(instance* inst);

void clean_MILP_min_sum(instance* inst);
