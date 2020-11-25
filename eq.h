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

    // class RangeOverFlowException: public std::exception {
    //     virtual const char* msg() const throw() {
    //         return "Range OverFlow";
    //     }
    // };

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

    class not_symmetric_positive_definite: public std::exception {
    public:
        const char* what() const throw()  {
            return "The matrix is not a symmetric positive definite matrix";
        }
    };

    Matrix(): rows(0), cols(0), m({}) {}
    Matrix(int rows_, int cols_): rows(rows_), cols(cols_) {
        m = matrix_t(rows, std::vector<long double>(cols, 0));
    }
    Matrix(matrix_t m_){
        size_t prev_col = m_[0].size();
        for (size_t r = 0; r < m_.size(); r++){
            // printf("prev_col = %d, m[r].size() = %d\n", prev_col, m_[r].size()); fflush(stdout);
            if (prev_col != m_[r].size()) throw not_matrix();
            prev_col = m_[r].size();
        }

        m = m_;
        rows = m_.size();
        cols = rows > 0 ? m_[0].size() : 0;
    }

    bool operator==(Matrix other);
    Matrix operator*(long double coef);
    Matrix operator*(Matrix other);
    Row operator[](int row){
        return Row(*this, row);
    }

    Matrix col(int c);
    Matrix last_col();
    Matrix get_transpose();
    bool is_symmetric();
    bool is_upper();
    bool is_lower();
    bool is_diag();
    Matrix swap(int r1, int r2);
    Matrix row_add(int r1, int r2);
    Matrix row_sub(int r1, int r2);
    Matrix line_elimination(int r1, int r2);
    Matrix to_upper();
    Matrix to_upper_with_column_pivot();
    Matrix upper2diag();
    Matrix get_augment(Matrix b);

    std::vector<Matrix> Doolittle_decompose();
    std::vector<Matrix> Crout_decompose();
    std::vector<Matrix> Cholesky_decompose();
    std::vector<Matrix> refined_Cholesky_decompose();

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

    Equation(Matrix A_, Matrix b_): A(A_), b(b_) {
        std::cout << A_.rows << " " << b_.rows << std::endl;
        if (A_.rows != b_.rows) throw row_not_fit();
    }

    class row_not_fit: public std::exception {
    public:
        const char* what() const throw() {
            return "Coefficient matrix should have same rows with the constant vector";
        }
    };

    bool operator==(Equation other);
    Matrix augment();
    Matrix Gaussian_elimination();
    Matrix Gaussian_elimination_with_column_pivot();
    Equations Doolittle_decompose();
    Equations Crout_decompose();
    Equations Cholesky_decompose();
    Equations refined_Cholesky_decompose();

    friend std::ostream& operator<< (std::ostream& os, Equation &eq);
};


