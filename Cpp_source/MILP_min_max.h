#pragma once

#include "global_variables.h"

void build_MILP_min_max(instance* inst);

void solve_MILP_min_max(instance* inst);

void clean_MILP_min_max(instance* inst);