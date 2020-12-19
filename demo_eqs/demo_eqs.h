#pragma once
#include <vector>
#include "matrix.h"


Matrix A1({
    {   4,   2,  -3,  -1,   2,   1,   0,   0,   0,   0},
    {   8,   6,  -5,  -3,   6,   5,   0,   1,   0,   0},
    {   4,   2,  -2,  -1,   3,   2,  -1,   0,   3,   1},
    {   0,  -2,   1,   5,  -1,   3,  -1,   1,   9,   4},
    {  -4,   2,   6,  -1,   6,   7,  -3,   3,   2,   3},
    {   8,   6,  -8,   5,   7,  17,   2,   6,  -3,   5},
    {   0,   2,  -1,   3,  -4,   2,   5,   3,   0,   1},
    {  16,  10, -11,  -9,  17,  34,   2,  -1,   2,   2},
    {   4,   6,   2,  -7,  13,   9,   2,   0,  12,   4},
    {   0,   0,  -1,   8,  -3, -24,  -8,   6,   3,  -1}
});
Matrix b1({{5}, {12}, {3}, {2}, {3}, {46}, {13}, {38}, {19}, {-21}});
Equation eq1(A1, b1);
Matrix s1({{1}, {-1}, {0}, {1}, {2}, {0}, {3}, {1}, {-1}, {2}});

Matrix A2({
    {   4,   2,  -4,   0,   2,   4,   0,   0},
    {   2,   2,  -1,  -2,   1,   3,   2,   0},
    {  -4,  -1,  14,   1,  -8,  -3,   5,   6},
    {   0,  -2,   1,   6,  -1,  -4,  -3,   3},
    {   2,   1,  -8,  -1,  22,   4, -10,  -3},
    {   4,   3,  -3,  -4,   4,  11,   1,  -4},
    {   0,   2,   5,  -3, -10,   1,  14,   2},
    {   0,   0,   6,   3,  -3,  -4,   2,  19},
});
Matrix b2({{0}, {-6}, {6}, {23}, {11}, {-22}, {-15}, {45}});
Equation eq2(A2, b2);
Matrix s2({{1}, {-1}, {0}, {2}, {1}, {-1}, {0}, {2}});

Matrix A3({
    {   4,  -1,   0,   0,   0,   0,   0,   0,   0,   0},
    {  -1,   4,  -1,   0,   0,   0,   0,   0,   0,   0},
    {   0,  -1,   4,  -1,   0,   0,   0,   0,   0,   0},
    {   0,   0,  -1,   4,  -1,   0,   0,   0,   0,   0},
    {   0,   0,   0,  -1,   4,  -1,   0,   0,   0,   0},
    {   0,   0,   0,   0,  -1,   4,  -1,   0,   0,   0},
    {   0,   0,   0,   0,   0,  -1,   4,  -1,   0,   0},
    {   0,   0,   0,   0,   0,   0,  -1,   4,  -1,   0},
    {   0,   0,   0,   0,   0,   0,   0,  -1,   4,  -1},
    {   0,   0,   0,   0,   0,   0,   0,   0,  -1,   4},
});
Matrix b3({{7}, {5}, {-13}, {2}, {6}, {-12}, {14}, {-4}, {5}, {-5}});
Equation eq3(A3, b3);
Matrix s3({{2}, {1}, {-3}, {0}, {1}, {-2}, {3}, {0}, {1}, {-1}});

std::vector<Equation> DEMO_EQUATIONS = {eq1, eq2, eq3};
