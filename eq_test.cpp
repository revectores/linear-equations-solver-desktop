#include "eq_test.h"


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

Matrix A4({
    {2, -1, -1},
    {1,  5, -1},
    {1,  1, 10}
});
Matrix b4 = Matrix(matrix_t({{-5}, {8}, {11}}));
Equation eq4(A4, b4);
Matrix s4(matrix_t({{-1}, {2}, {1}}));

Matrix A5({
    {4, 3, 0},
    {3, 4, -1},
    {0, -1, 4}
});
Matrix b5(matrix_t({{24}, {30}, {-24}}));
Equation eq5(A5, b5);
Matrix s5(matrix_t({{3}, {4}, {-5}}));


std::vector<Equation> TEST_EQUATIONS = {eq1, eq2, eq3};



void test_Gaussian_elimination() {
    Equation eq1_copy(eq1);

    Matrix r1 = eq1.Gaussian_elimination();
    Matrix r2 = eq2.Gaussian_elimination();
    Matrix r3 = eq3.Gaussian_elimination();
    Matrix rr1 = eq1.Gaussian_elimination_with_column_pivot();
    Matrix rr2 = eq2.Gaussian_elimination_with_column_pivot();
    Matrix rr3 = eq3.Gaussian_elimination_with_column_pivot();

    assert(eq1 == eq1_copy);
    assert(r1 == s1);
    assert(r2 == s2);
    assert(r3 == s3);
    assert(rr1 == s1);
    assert(rr2 == s2);
    assert(rr3 == s3);

    // std::cout << s1 << std::endl;
    // std::cout << s4 << std::endl << std::endl;
    // std::cout << s2 << std::endl;
    // std::cout << s5 << std::endl << std::endl;
    // std::cout << s3 << std::endl;
    // std::cout << s6 << std::endl << std::endl;
}


void test_Doolittle_decompose(){
    std::vector<Matrix> LU1 = A1.Doolittle_decompose();
    std::vector<Matrix> LU2 = A2.Doolittle_decompose();
    std::vector<Matrix> LU3 = A3.Doolittle_decompose();
    Equations eqs1(LU1, b1);
    Equations eqs2(LU2, b2);
    Equations eqs3(LU3, b3);

    assert(LU1[0].is_lower() && LU1[1].is_upper());
    assert(LU2[0].is_lower() && LU2[1].is_upper());
    assert(LU3[0].is_lower() && LU3[1].is_upper());
    assert(LU1[0] * LU1[1] == A1);
    assert(LU2[0] * LU2[1] == A2);
    assert(LU3[0] * LU3[1] == A3);
    assert(eqs1.solve() == s1);
    assert(eqs2.solve() == s2);
    assert(eqs3.solve() == s3);
}


void test_Cholesky_decompose(){
    std::vector<Matrix> LLt2 = A2.Cholesky_decompose();
    Equations eqs2(LLt2, b2);

    assert(LLt2[0].is_lower() && LLt2[1].is_upper());
    assert(LLt2[0].get_transpose() == LLt2[1]);
    assert(LLt2[0] * LLt2[1] == A2);
    assert(eqs2.solve() == s2);
}


void test_refined_Cholesky_decompose(){
    std::vector<Matrix> LDLt2 = A2.refined_Cholesky_decompose();
    Equations eqs2(LDLt2, b2);

    assert(LDLt2[0].is_lower() && LDLt2[2].is_upper());
    assert(LDLt2[0].get_transpose() == LDLt2[2]);
    assert(LDLt2[1].is_diag());
    assert(LDLt2[0] * LDLt2[1] * LDLt2[2] == A2);
    assert(eqs2.solve() == s2);
}


void test_Crout_decompose(){
    std::vector<Matrix> LU1 = A1.Crout_decompose();
    std::vector<Matrix> LU2 = A2.Crout_decompose();
    std::vector<Matrix> LU3 = A3.Crout_decompose();
    Equations eqs1(LU1, b1);
    Equations eqs2(LU2, b2);
    Equations eqs3(LU3, b3);

    assert(LU1[0].is_lower() && LU1[1].is_upper());
    assert(LU2[0].is_lower() && LU2[1].is_upper());
    assert(LU3[0].is_lower() && LU3[1].is_upper());
    assert(LU1[0] * LU1[1] == A1);
    assert(LU2[0] * LU2[1] == A2);
    assert(LU3[0] * LU3[1] == A3);
    assert(eqs1.solve() == s1);
    assert(eqs2.solve() == s2);
    assert(eqs3.solve() == s3);
}


void test_get_augment(){
    Matrix M({{4, 7}, {2, 6}});
    Matrix MM = M.get_augment(M);
    // std::cout << MM << std::endl;
}


void test_get_col_slice(){
    Matrix M({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    Matrix expected_slice_2_3((matrix_t){{3}, {6}, {9}});
    Matrix slice_2_3 = M.get_col_slice(2, 3);
    
    Matrix MM = M.get_augment(M);
    Matrix M_new = MM.get_col_slice(M.cols, M.cols * 2);
    
    assert(slice_2_3 == expected_slice_2_3);
    assert(M_new == M);
}


void test_get_inverse(){
    Matrix M1({{4, 7}, {2, 6}});
    Matrix I({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    Matrix expected_M1_inv({{0.6, -0.7}, {-0.2, 0.4}});

    assert(M1.get_inverse() == expected_M1_inv);
    assert(I.get_inverse() == I);
}


void test_Jacobi_iteration(){
    // Matrix r1 = eq1.Jacobi_iteration();
    // Matrix r2 = eq2.Jacobi_iteration();
    Matrix r2 = eq2.Jacobi_iteration();
    Matrix r3 = eq3.Jacobi_iteration();
    Matrix r4 = eq4.Jacobi_iteration();
    Matrix r5 = eq5.Jacobi_iteration();

    Matrix init_vector = Matrix(matrix_t({{2}, {2}, {2}}));
    Matrix r5_2 = eq5.Jacobi_iteration(init_vector);

    // std::cout << r1 << std::endl << std::endl;
    // std::cout << s1 << std::endl << std::endl;
    // assert(r2 == s2);
    // std::cout << r3 << std::endl;
    // std::cout << s3 << std::endl;
    // std::cout << r2 << std::endl;
    assert((r3 - s3).inf_norm() < 1e-3);
    assert((r4 - s4).inf_norm() < 1e-3);
    assert((r5 - s5).inf_norm() < 1e-3);
    assert((r5_2 - s5).inf_norm() < 1e-3);
}


void test_Gauss_Seidel_iteration(){
    // Matrix r1 = eq1.Gauss_Seidol_iteration();
    Matrix r2 = eq2.Gauss_Seidol_iteration();
    Matrix r3 = eq3.Gauss_Seidol_iteration();
    Matrix r4 = eq4.Gauss_Seidol_iteration();
    Matrix r5 = eq5.Gauss_Seidol_iteration();

    // std::cout << r1 << std::endl << std::endl;
    // std::cout << s1 << std::endl << std::endl;
    // assert(r1 == s1);
    // std::cout << r2 << std::endl;
    // std::cout << s2 << std::endl;
    // std::cout << r5 << std::endl;
    // std::cout << r4 << std::endl;
    assert((r2 - s2).inf_norm() < 1e-3);
    assert((r3 - s3).inf_norm() < 1e-5);
    assert((r4 - s4).inf_norm() < 1e-5);
    assert((r5 - s5).inf_norm() < 1e-5);
}


void test_SOR(){
    Matrix r1 = eq1.SOR(1.25);
    Matrix r2 = eq2.SOR(1.25);
    Matrix r3 = eq3.SOR(1.25);
    Matrix r4 = eq4.SOR(1.25);
    Matrix r5 = eq5.SOR(1.25);
    // std::cout << r1 << std::endl;
    // std::cout << r2 << std::endl;
    // std::cout << r3 << std::endl;
    // std::cout << r4 << std::endl;
    // std::cout << r5 << std::endl;
    // assert((r1 - s1).inf_norm() < 1e-3);
    assert((r2 - s2).inf_norm() < 1e-3);
    assert((r3 - s3).inf_norm() < 1e-3);
    assert((r4 - s4).inf_norm() < 1e-3);
    assert((r5 - s5).inf_norm() < 1e-3);
}


// int main(){
//     test_Gaussian_elimination();
//     test_Doolittle_decompose();
//     test_Cholesky_decompose();
//     test_refined_Cholesky_decompose();
//     test_Crout_decompose();

//     test_get_augment();
//     test_get_col_slice();
//     test_get_inverse();

//     test_Jacobi_iteration();
//     test_Gauss_Seidel_iteration();
//     test_SOR();
// }
