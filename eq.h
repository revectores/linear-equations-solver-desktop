#pragma once

#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

typedef std::vector<std::vector<long double>> matrix_t;

class Matrix {
public:
    int rows;
    int cols;
    matrix_t m;

    class Row {
        friend class Matrix;
    private:
        Matrix& parent;
        int row;
        Row(Matrix &parent_, int row_): parent(parent_), row(row_) {}
    public:
        long double& operator[](int col){
            return parent.m[row][col];
        }
        // Row& operator=(const Row& r){
        //     std::cout << parent.cols << " " << r.parent.cols << std::endl;
        //     assert(parent.cols == r.parent.cols);
        //     return *this;
        // }
    };

    class not_matrix: public std::exception {
    public:
        const char* what() const throw() {
            return "The two dimensional vector cannot be parsed as a matrix";
        }
    };

    class not_square: public std::exception {
    public:  
        const char* what() const throw() {
            return "The matrix is not a square matrix";
        }
    };  

    class dimension_not_fit: public std::exception {
    public:  
        const char* what() const throw() {
            return "Matrix dimension not fit";
        }
    };

    class not_symmetric_positive_definite: public std::exception {
    public:
        const char* what() const throw()  {
            return "The matrix is not a symmetric positive definite matrix";
        }
    };

    class matrix_slice_error: public std::exception {
    public:
        const char* what() const throw()  {
            return "The slice end index must be larger than the slice start";
        }
    };

    Matrix(): rows(0), cols(0), m({}) {}
    Matrix(int rows_, int cols_): rows(rows_), cols(cols_) {
        m = matrix_t(rows, std::vector<long double>(cols, 0));
    }
    Matrix(int rows_, int cols_, long double value): rows(rows_), cols(cols_) {
        m = matrix_t(rows, std::vector<long double>(cols, value));
    }
    Matrix(matrix_t m_){
        if (m_.size() > 0){
            size_t prev_col = m_[0].size();
            for (size_t r = 0; r < m_.size(); r++){
                // printf("prev_col = %d, m[r].size() = %d\n", prev_col, m_[r].size()); fflush(stdout);
                // if (prev_col != m_[r].size()) throw not_matrix();
                prev_col = m_[r].size();
            }
        }

        m = m_;
        rows = m_.size();
        cols = rows > 0 ? m_[0].size() : 0;
    }

    bool operator==(Matrix other) const;
    Matrix operator+(Matrix other) const;
    Matrix operator-() const;
    Matrix operator-(Matrix other) const;
    Matrix operator*(long double coef) const;
    Matrix operator*(Matrix other) const;
    Matrix get_inverse() const;
    Row operator[](int row) {
        return Row(*this, row);
    }

    bool is_empty() const;
    bool is_symmetric() const;
    bool is_upper() const;
    bool is_lower() const;
    bool is_diag() const;
    long double inf_norm() const;
    Matrix col(int c) const;
    Matrix get_col_slice(int col_st, int col_ed) const;
    Matrix last_col() const;
    Matrix get_transpose() const;
    Matrix get_diag() const;
    Matrix get_upper() const;
    Matrix get_lower() const;
    Matrix get_identity() const;

    Matrix swap(int r1, int r2);
    Matrix row_add(int r1, int r2);
    Matrix row_sub(int r1, int r2);
    Matrix line_elimination(int r1, int r2);
    Matrix to_upper();
    Matrix to_upper_with_column_pivot();
    Matrix upper2diag();
    Matrix diag2identity();

    Matrix get_augment(Matrix A) const;
    std::vector<Matrix> Doolittle_decompose() const;
    std::vector<Matrix> Crout_decompose() const;
    std::vector<Matrix> Cholesky_decompose() const;
    std::vector<Matrix> refined_Cholesky_decompose() const;

    friend std::ostream& operator<< (std::ostream& os, Matrix &m);
};


class Equations {
public:
    std::vector<Matrix> As;
    Matrix b;

    Equations(): As({}), b(){}
    Equations(std::vector<Matrix> As_, Matrix b_): As(As_), b(b_){}
    Matrix solve();
};


class Equation {
public:
    Matrix A;
    Matrix b;

    class row_not_fit: public std::exception {
    public:
        const char* what() const throw() {
            return "Coefficient matrix should have same rows with the constant vector";
        }
    };

    Equation(Matrix A_, Matrix b_): A(A_), b(b_) {
        // std::cout << A_.rows << " " << b_.rows << std::endl;
        // if (A_.rows != b_.rows) throw row_not_fit();
    }

    bool operator==(Equation other);
    Matrix augment();
    Matrix Gaussian_elimination();
    Matrix Gaussian_elimination_with_column_pivot();

    Matrix Jacobi_iteration(Matrix init_vector=Matrix(), long double precision=1e-6, size_t* iteration_counter=nullptr);
    Matrix Gauss_Seidol_iteration(Matrix init_vector=Matrix(), long double precision=1e-6, size_t* iteration_counter=nullptr);
    Matrix SOR(long double w, Matrix init_vector=Matrix(), long double precision=1e-6, size_t* iteration_counter=nullptr);

    Equations Doolittle_decompose() const;
    Equations Crout_decompose() const;
    Equations Cholesky_decompose() const;
    Equations refined_Cholesky_decompose() const;


    friend std::ostream& operator<< (std::ostream& os, Equation &eq);
};
