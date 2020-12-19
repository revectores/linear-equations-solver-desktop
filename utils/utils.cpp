#include <cmath>
#include "utils.h"

int max(int a, int b){
    return a > b ? a : b;
}

bool is_zero(long double x){
    return (abs(x) < 1e-10);
}
