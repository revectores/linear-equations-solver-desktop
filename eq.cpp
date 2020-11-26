#include "eq.h"
#include "utils.h"


bool Matrix::operator==(Matrix other){
    if ((rows != other.rows) || (cols != other.cols)){
        return false;
    }
    for (int r = 0; r < rows; r++){
        for (int c = 0; c < cols; c++){
            if (abs(m[r][c] - other[r][c]) > 1e-10) return false;
            // if (m[r][c] != other[r][c]) return false;
        }
    }
    return true;
}


Matrix Matrix::operator*(long double coef){
    Matrix prod(*this);
    for (int r = 0; r < prod.rows; r++) {
        for (int c = 0; c < prod.cols; c++) {
            prod[r][c] *= coef;
        }
    }
    return prod;
}


Matrix Matrix::operator*(Matrix other){
    if (cols != other.rows) {
        throw Matrix::dimension_not_fit();
    }
    Matrix prod(rows, other.cols);

    for (int r = 0; r < rows; r++){
        for (int c = 0; c < other.cols; c++){
            for (int i = 0; i < cols; i++){
                prod[r][c] += m[r][i] * other[i][c];
            }
        }
    }
    return prod;
}


std::ostream& operator<< (std::ostream& os, Matrix &m){
    if (m.cols == 1) {
        printf("(");
        for (int r = 0; r < m.rows; r++){
            if (abs(m[r][0]) < 1e-10){
                printf("           0");
            } else {
                printf("%12.6Lf", m[r][0]);
            }
        }
        printf(")^T");
    } else {
        for (int r = 0; r < m.rows; r++){
            for (int c = 0; c < m.cols; c++){
                if (abs(m[r][c]) < 1e-10){
                    printf("           0");
                } else {
                    printf("%12.6Lf", m[r][c]);
                }
            }
            printf("\n");
        }
    }
    return os;
}



Matrix Matrix::col(int c){
    if (c >= cols) {
        throw std::out_of_range("column range exceeded");
    }
    std::vector<std::vector<long double>> col_vector;
    for (int r = 0; r < rows; r++){
        col_vector.push_back(std::vector<long double>{m[r][c]});
    }
    return Matrix(col_vector);
}


Matrix Matrix::last_col(){
    return col(cols-1);
}


Matrix Matrix::get_transpose(){
    Matrix T(cols, rows);
    for (int r = 0; r < rows; r++) for (int c = 0; c < cols; c++) T[c][r] = m[r][c];
    return T;
}


bool Matrix::is_symmetric(){
    return (*this) == get_transpose();
}


bool Matrix::is_upper(){
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < r; c++) {
            if (!is_zero(m[r][c])) return false;
        }
    }
    return true;
}


bool Matrix::is_lower(){
    return get_transpose().is_upper();
}


bool Matrix::is_diag(){
    return is_upper() && is_lower();
}


Matrix Matrix::swap(int r1, int r2){
    if ((r1 >= rows) || (r2 >= rows)) throw std::out_of_range("row range exceeded");
    // std::vector<double> temp = m[r2];
    auto temp = m[r2];
    m[r2] = m[r1];
    m[r1] = temp;
    return *this;
}


Matrix Matrix::row_add(int r1, int r2) {
    if ((r1 >= rows) || (r2 >= rows)) throw std::out_of_range("row range exceeded");
    for (int c = 0; c < rows; c++) m[r1][c] = m[r1][c] + m[r2][c];
    return *this;
}


Matrix Matrix::row_sub(int r1, int r2) {
    if ((r1 >= rows) || (r2 >= rows)) throw std::out_of_range("row range exceeded");
    for (int c = 0; c < rows; c++) m[r1][c] = m[r1][c] - m[r2][c];
    return *this;
}


Matrix Matrix::line_elimination(int r1, int r2){
    if ((r1 >= rows) || (r2 >= rows)) throw std::out_of_range("row range exceeded");
    return *this;
}





Matrix Matrix::to_upper() {
    for (int base_row = 0; base_row < rows; base_row++) {
        for (int row = base_row + 1; row < rows; row++) {
            long double coef = m[row][base_row] / m[base_row][base_row];
            for (int c = 0; c < cols; c++){
                m[row][c] -= coef * m[base_row][c];
            }
        }
    }
    return *this;
}


Matrix Matrix::to_upper_with_column_pivot() {
    for (int base_row = 0; base_row < rows; base_row++) {
        int pivot_row = base_row;
        for (int r = base_row; r < rows; r++){
            if (m[r][base_row] > m[pivot_row][base_row]) {
                pivot_row = r;
            }
        }
        this->swap(base_row, pivot_row);

        for (int r = base_row + 1; r < rows; r++) {
            long double coef = m[r][base_row] / m[base_row][base_row];
            for (int c = 0; c < cols; c++){
                m[r][c] -= coef * m[base_row][c];
            }
        }
    }
    return *this;
}


Matrix Matrix::upper2diag(){
    for (int base_row = rows - 1; base_row >= 1; base_row--) {
        for (int row = base_row - 1; row >= 0; row--) {
            long double coef = m[row][base_row] / m[base_row][base_row];
            for (int c = 0; c < cols; c++) {
                m[row][c] -= coef * m[base_row][c];
            }
        }
    }

    return *this;
}


Matrix Matrix::diag2identity(){
    for (int r = 0; r < rows; r++){
        m[r][cols - 1] /= m[r][r];
        m[r][r] = 1;
    }
    // for (int r = 0; r < rows; r++) {
    //     for (int c = 0; c < cols; c++){
    //         m[r][c] /= m[r][r];
    //     }
    // }
    return *this;
}


Matrix Matrix::get_augment(Matrix b){
    if (b.rows != rows){
        throw Matrix::dimension_not_fit();
    }
    Matrix Mb(*this);
    for (int row = 0; row < rows; row++){
        Mb.m[row].push_back(b[row][0]);
    }
    Mb.cols += 1;
    return Mb;
}


std::vector<Matrix> Matrix::Doolittle_decompose(){
    if (rows != cols) throw Matrix::not_square();

    Matrix L(rows, cols);
    Matrix U(rows, cols);

    for (int k = 0; k < rows; k++) {
        for (int j = k; j < rows; j++){
            long double sum = 0;
            for (int r = 0; r <= k-1; r++) sum += L[k][r] * U[r][j];
            U[k][j] = m[k][j] - sum;
        }
        for (int i = k; i < rows; i++){
            long double sum = 0;
            for (int r = 0; r <= k-1; r++) sum += L[i][r] * U[r][k];
            L[i][k] = (m[i][k] - sum) / U[k][k];
        }
    }
    return {L, U};
}


std::vector<Matrix> Matrix::Crout_decompose(){
    if (rows != cols) throw Matrix::not_square();

    Matrix L(rows, cols);
    Matrix U(rows, cols);

    for (int k = 0; k < rows; k++){
        for (int i = 0; i < rows; i++){
            long double sum = 0;
            for(int r = 0; r <= k-1; r++) sum += L[i][r] * U[r][k];
            L[i][k] = m[i][k] - sum;
        }

        for (int j = k; j < rows; j++){
            long double sum = 0;
            for(int r = 0; r <= k-1; r++) sum += L[k][r] * U[r][j];
            U[k][j] = (m[k][j] - sum) / L[k][k];
        }
    }
    return {L, U};
}


std::vector<Matrix> Matrix::Cholesky_decompose(){
    if (rows != cols)    throw Matrix::not_square();
    if (!is_symmetric()) throw Matrix::not_symmetric_positive_definite();

    Matrix L(rows, cols);
    for (int j = 0; j < rows; j++){
        long double sum = 0;
        for (int k = 0; k <= j-1; k++) sum += pow(L[j][k], 2);
        L[j][j] = sqrt(m[j][j] - sum);

        for (int i = j+1; i < rows; i++){
            long double sum = 0;
            for (int k = 0; k <= j-1; k++) sum += L[i][k] * L[j][k];
            L[i][j] = (m[i][j] - sum) / L[j][j];
        }
    }
    return {L, L.get_transpose()};
}


std::vector<Matrix> Matrix::refined_Cholesky_decompose(){
    if (rows != cols)    throw Matrix::not_square();
    if (!is_symmetric()) throw Matrix::not_symmetric_positive_definite();

    Matrix L(rows, cols);
    Matrix D(rows, cols);

    for (int j = 0; j < rows; j++){
        long double sum = 0;
        for(int k = 0; k <= j-1; k++) sum += pow(L[j][k], 2) * D[k][k];
        D[j][j] = m[j][j] - sum;

        for (int i = j; i < rows; i++){
            long double sum = 0;
            for(int k = 0; k <= j-1; k++) sum += L[i][k] * D[k][k] * L[j][k];
            L[i][j] = (m[i][j] - sum) / D[j][j];
        }
    }
    return {L, D, L.get_transpose()};
}



bool Equation::operator==(Equation other){
    return (A == other.A) && (b == other.b);
}


std::ostream& operator<< (std::ostream& os, Equation &eq){
    os << eq.A << std::endl;
    os << eq.b << std::endl;
    return os;
}


Matrix Equation::Gaussian_elimination() {
    Matrix Ab = A.get_augment(b);
    return Ab.to_upper().upper2diag().diag2identity().last_col();
}


Matrix Equation::Gaussian_elimination_with_column_pivot() {
    Matrix Ab = A.get_augment(b);
    return Ab.to_upper_with_column_pivot().upper2diag().diag2identity().last_col();
}


Equations Equation::Doolittle_decompose(){
    return Equations(A.Doolittle_decompose(), b);
}


Equations Equation::Crout_decompose(){
    return Equations(A.Crout_decompose(), b);
}


Equations Equation::Cholesky_decompose(){
    return Equations(A.Cholesky_decompose(), b);
}


Equations Equation::refined_Cholesky_decompose(){
    return Equations(A.refined_Cholesky_decompose(), b);
}


Matrix Equations::solve(){
    for(auto matrix_it = As.begin(); matrix_it < As.end(); matrix_it++) {
        Equation eq(*matrix_it, b);
        b = eq.Gaussian_elimination();
    }
    return b;
}
