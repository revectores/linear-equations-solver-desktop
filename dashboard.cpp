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

char* RADIO_STRINGS[6] = {
    "Gaussian Elimination",
    "Gaussian Elimination with Column Pivot",
    "Doolittle Decompose",
    "Crout Decompose",
    "Cholesky Decompose",
    "Refined Cholesky Decompose",
};


Dashboard::Dashboard(QWidget *parent) : QWidget(parent) {
    QGridLayout *layout = new QGridLayout(this);

    matrix_table_init();
    constant_table_init();
    solution_table_init();
    method_group_init();
    action_group_init();
    decompose_label_init();

    layout->addWidget(matrix_table, 0, 0, 2, 1);
    layout->addWidget(solution_table, 0, 1, 2, 1);
    layout->addWidget(constant_table, 0, 2, 2, 1);
    layout->addWidget(method_group, 0, 3);
    layout->addWidget(action_group, 1, 3);
    layout->addWidget(decompose_label, 2, 0, 1, 4);

    cur_demo_id = 0;
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
        method_radios[i] = new QRadioButton(RADIO_STRINGS[i]);
        layout->addWidget(method_radios[i], i, 0);
    }

    method_group->setFixedHeight(140);
    method_group->setLayout(layout);
}

void Dashboard::action_group_init(){
    action_group = new QGroupBox();
    QGridLayout *layout = new QGridLayout;

    QPushButton *solve_btn = new QPushButton("Solve");
    QPushButton *demo_btn = new QPushButton("Demo 1");
    QPushButton *clear_btn = new QPushButton("Clear");
    QPushButton *exit_btn = new QPushButton("Exit");

    connect(solve_btn, &QPushButton::clicked, this, &Dashboard::solve);
    connect(demo_btn, &QPushButton::clicked, this, &Dashboard::demo);
    connect(clear_btn, &QPushButton::clicked, this, &Dashboard::clear);

    layout->addWidget(solve_btn, 0, 0, 1, 2);
    layout->addWidget(demo_btn, 1, 0);
    layout->addWidget(clear_btn, 1, 1);

    action_group->setFixedHeight(212-140);
    action_group->setLayout(layout);
}


void Dashboard::decompose_label_init(){
    decompose_label = new QLabel;
}


void Dashboard::demo(){
    Equation eq = TEST_EQUATIONS[cur_demo_id];
    Matrix A = eq.A;
    Matrix b = eq.b;

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
}


void Dashboard::read_equation(){
    int max_r = -1;
    int max_c = -1;

    for (int r = 0; r < matrix_table->rowCount(); r++){
        for (int c = 0; c < matrix_table->columnCount(); c++){
            QTableWidgetItem* e = matrix_table->item(r, c);
            if (e && !e->text().isEmpty()) {
                max_r = max(max_r, r);
                max_c = max(max_c, c);
            }
        }
    }

    // if ((max_r == -1) && (max_c == -1)) return;

    // Matrix coef_matrix;
    // Matrix constant_matrix;
}



void Dashboard::solve(){
    read_equation();

    int method_id;

    for (method_id = 0; method_id < METHOD_COUNT; method_id++){
        if (method_radios[method_id]->isChecked()) break;
    }
    switch (method_id){
    case GAUSSIAN_ELIMINATOIN:

    break;
    case DOOLITTLE_DECOMPOSE:

    break;

    }
}






