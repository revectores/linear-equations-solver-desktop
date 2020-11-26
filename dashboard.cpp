#include "dashboard.h"
#include <cstdio>



const int METHOD_COUNT = 6;
enum method {
    GAUSSIAN_ELIMINATOIN,
    GAUSSIAN_ELIMINATOIN_WITH_COLUMN_PIVOT,
    DOOLITTLE_DECOMPOSE,
    CROUT_DECOMPOSE,
    CHOLESKY_DECOMPOSE,
    REFINED_CHOLESKY_DECOMPOSE
};

std::string RADIO_STRINGS[6] = {
    "Gaussian Elimination",
    "Gaussian Elimination with Column Pivot",
    "Doolittle Decompose",
    "Crout Decompose",
    "Cholesky Decompose",
    "Refined Cholesky Decompose",
};



Dashboard::Dashboard(QWidget *parent) : QWidget(parent) {
    QGridLayout *layout = new QGridLayout(this);
    DecomposeWindow decompose_window;
    decompose_window.show();

    matrix_table_init();
    constant_table_init();
    solution_table_init();
    method_group_init();
    action_group_init();
    info_label_init();

    layout->addWidget(matrix_table, 0, 0, 2, 1);
    layout->addWidget(solution_table, 0, 1, 2, 1);
    layout->addWidget(constant_table, 0, 2, 2, 1);
    layout->addWidget(method_group, 0, 3);
    layout->addWidget(action_group, 1, 3);
    layout->addWidget(info_label, 2, 0, 1, 4);

    cur_demo_id = 1;
    clear();
}


void Dashboard::matrix_table_init(){

    matrix_table = new QTableWidget(10, 10);
    matrix_table->horizontalHeader()->setDefaultSectionSize(40);
    matrix_table->verticalHeader()->setDefaultSectionSize(20);
    matrix_table->horizontalHeader()->hide();
    matrix_table->verticalHeader()->hide();

    matrix_table->setFixedWidth(402);
    matrix_table->setFixedHeight(212);
}

void Dashboard::constant_table_init(){
    constant_table = new QTableWidget(10, 1);

    constant_table->horizontalHeader()->setDefaultSectionSize(40);
    constant_table->verticalHeader()->setDefaultSectionSize(20);
    constant_table->horizontalHeader()->hide();
    constant_table->verticalHeader()->hide();
    constant_table->setFixedWidth(42);
    constant_table->setFixedHeight(212);
}

void Dashboard::solution_table_init(){
    solution_table = new QTableWidget(10, 1);
    solution_table->horizontalHeader()->setDefaultSectionSize(40);
    solution_table->verticalHeader()->setDefaultSectionSize(20);
    solution_table->horizontalHeader()->hide();
    solution_table->verticalHeader()->hide();
    solution_table->setFixedWidth(42);
    solution_table->setFixedHeight(212);
}

void Dashboard::method_group_init(){
    method_group = new QGroupBox;
    QGridLayout *layout = new QGridLayout;

    for (int i = 0; i < METHOD_COUNT; i++){
        method_radios[i] = new QRadioButton(QString::fromStdString(RADIO_STRINGS[i]));
        layout->addWidget(method_radios[i], i, 0);
    }
    method_radios[0]->setChecked(true);

    method_group->setFixedHeight(140);
    method_group->setLayout(layout);
}

void Dashboard::action_group_init(){
    action_group = new QGroupBox();
    QGridLayout *layout = new QGridLayout;

    QPushButton *solve_btn = new QPushButton("Solve");
    QPushButton *demo_btns[TEST_EQUATIONS.size()];

    QPushButton *clear_btn = new QPushButton("Clear");
    QPushButton *exit_btn = new QPushButton("Exit");

    connect(solve_btn, &QPushButton::clicked, this, &Dashboard::solve);
    connect(clear_btn, &QPushButton::clicked, this, &Dashboard::clear);

    for (int i = 0; i < TEST_EQUATIONS.size(); i++){
        demo_btns[i] = new QPushButton(QString("Equation ") + QString::number(i + 1));
        connect(demo_btns[i], &QPushButton::clicked, this, [=]{demo(i);});
        layout->addWidget(demo_btns[i], 1, i);
    }

    layout->addWidget(solve_btn, 0, 0, 1, 4);
    layout->addWidget(clear_btn, 1, 3);

    action_group->setFixedHeight(212-140);
    action_group->setLayout(layout);
}


void Dashboard::info_label_init(){
    info_label = new QLabel;
    info_label->setText("The execution info and error message will be shown here.");
    // info_label->setStyleSheet("QLabel {color: blue;}");
}


void Dashboard::info_msg(std::string msg){
    info_label->clear();
    // info_label->setStyleSheet("QLabel {color: blue;}");
    info_label->setText(QString::fromStdString("SUCCESS: " + msg));
}


void Dashboard::err_msg(std::string msg){
    info_label->clear();
    // info_label->setStyleSheet("QLabel {color: yellow;}");
    info_label->setText(QString::fromStdString("ERROR: " + msg));
}




void Dashboard::demo(int demo_id){
    Equation eq = TEST_EQUATIONS[demo_id];
    Matrix A = eq.A;
    Matrix b = eq.b;

    clear();

    for (int r = 0; r < A.rows; r++){
        for (int c = 0; c < A.cols; c++){
            QTableWidgetItem* e = new QTableWidgetItem(QString::number((double)A[r][c]));
            matrix_table->setItem(r, c, e);
        }
    }

    for (int r = 0; r < b.rows; r++){
        QTableWidgetItem* e = new QTableWidgetItem(QString::number((double)b[r][0]));
        constant_table->setItem(r, 0, e);
    }
}


void Dashboard::clear(){
    for (int r = 0; r < matrix_table->rowCount(); r++){
        for (int c = 0; c < matrix_table->columnCount(); c++){
            QTableWidgetItem* empty = new QTableWidgetItem("");
            matrix_table->setItem(r, c, empty);
        }
    }

    for (int r = 0; r < constant_table->rowCount(); r++){
        QTableWidgetItem* empty = new QTableWidgetItem("");
        constant_table->setItem(r, 0, empty);
    }

    for (int r = 0; r < solution_table->rowCount(); r++){
        QTableWidgetItem* empty = new QTableWidgetItem("");
        empty->setFlags(empty->flags() & ~Qt::ItemIsEditable);
        solution_table->setItem(r, 0, empty);
    }

    decompose_window.clear();
}


Equation Dashboard::read_equation(){
    size_t max_row = matrix_table->rowCount();
    size_t max_col = matrix_table->columnCount();
    size_t r = 0;
    size_t c = 0;
    matrix_t Am;

    while(r < max_row) {
        QTableWidgetItem* e = matrix_table->item(r, c);
        bool is_empty = (!e) || (e->text().isEmpty());
        // printf("(%d, %d)\n", r, c);

        if (c == 0) {
            if (is_empty) break;
            Am.push_back({});
        }

        if (c < max_col && !is_empty) {
            long double e_value = e->text().toDouble();
            Am.back().push_back(e_value);
            c++;
        } else {
            c = 0;
            r++;
        }
    }
    // printf("1111");

    max_row = constant_table->rowCount();
    matrix_t bm;
    r = 0;
    while(r < max_row) {
        QTableWidgetItem* e = constant_table->item(r, 0);
        bool is_empty = (!e) || (e->text().isEmpty());

        if (is_empty) break;
        long double e_value = e->text().toDouble();
        bm.push_back({e_value});
        r++;
    }
    // for (size_t r = 0; r < m.size(); r++){
    //     for (size_t c = 0; c < m[r].size(); c++){
    //         std::cout << m[r][c] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    try {
        Matrix A(Am);
        Matrix b(bm);
        Equation eq(A, b);
        return eq;
    } catch (Matrix::not_matrix e) {
        err_msg(e.what());
    } catch (Equation::row_not_fit e) {
        err_msg(e.what());
    }
}


void Dashboard::fill_solution(Matrix solution){
    for (int r = 0; r < solution.rows; r++) {
        double e_value = solution[r][0];
        if (is_zero(e_value)) e_value = 0;
        QTableWidgetItem* e = new QTableWidgetItem(QString::number(e_value));
        // e->setFlags(e->flags() & ~Qt::ItemIsEditable);
        solution_table->setItem(r, 0, e);
    }
}


void Dashboard::solve(){
    Equation eq = read_equation();

    int method_id;
    for (method_id = 0; method_id < METHOD_COUNT; method_id++){
        if (method_radios[method_id]->isChecked()) break;
    }

    Equations eqs;
    Matrix solution;

    try {
        switch (method_id){
        case GAUSSIAN_ELIMINATOIN:
            solution = eq.Gaussian_elimination();
            break;

        case GAUSSIAN_ELIMINATOIN_WITH_COLUMN_PIVOT:
            solution = eq.Gaussian_elimination_with_column_pivot();
            break;

        case DOOLITTLE_DECOMPOSE:
            eqs = eq.Doolittle_decompose();
            solution = eqs.solve();
            break;

        case CROUT_DECOMPOSE:
            eqs = eq.Crout_decompose();
            solution = eqs.solve();
            break;

        case CHOLESKY_DECOMPOSE:
            try {
                eqs = eq.Cholesky_decompose();
                solution = eqs.solve();
            } catch (Matrix::not_symmetric_positive_definite e) {
                err_msg(e.what());
            }
            break;

        case REFINED_CHOLESKY_DECOMPOSE:
            try {
                eqs = eq.refined_Cholesky_decompose();
                solution = eqs.solve();
            } catch (Matrix::not_symmetric_positive_definite e) {
                err_msg(e.what());
            }
            break;
        }
    } catch (Matrix::dimension_not_fit e) {
        err_msg(e.what());
        return;
    }

    // std::cout << solution << std::endl; fflush(stdout);
    fill_solution(solution);
    info_msg("The equation has been successfully solved by " + RADIO_STRINGS[method_id]);

    decompose_window.clear();
    if (method_id >= DOOLITTLE_DECOMPOSE){
        for (int i = 0; i < eqs.As.size(); i++) {
            decompose_window.fill(i, eqs.As[i]);
        }
        decompose_window.show();
    }

}
