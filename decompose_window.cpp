#include "decompose_window.h"

DecomposeWindow::DecomposeWindow(QWidget *parent) : QWidget(parent) {
    resize(1080, 220);
    setWindowTitle("Matrix Decomposition");

    QGridLayout *layout = new QGridLayout(this);

    tables_init();

    for (int i = 0; i < TABLE_COUNT; i++){
        layout->addWidget(matrix_tables[i], 0, i);
    }
}


void DecomposeWindow::tables_init(){
    QGridLayout *layout = new QGridLayout;

    for (int i = 0; i < TABLE_COUNT; i++){
        matrix_tables[i] = new QTableWidget(10, 10);
        matrix_tables[i]->horizontalHeader()->setDefaultSectionSize(40);
        matrix_tables[i]->verticalHeader()->setDefaultSectionSize(20);
        matrix_tables[i]->horizontalHeader()->hide();
        matrix_tables[i]->verticalHeader()->hide();

        matrix_tables[i]->setFixedWidth(402);
        matrix_tables[i]->setFixedHeight(212);

    }
}


void DecomposeWindow::clear(){
    for (int i = 0; i < TABLE_COUNT; i++){
        for (int r = 0; r < matrix_tables[i]->rowCount(); r++){
            for (int c = 0; c < matrix_tables[i]->columnCount(); c++){
                QTableWidgetItem* empty = new QTableWidgetItem("");
                matrix_tables[i]->setItem(r, c, empty);
            }
        }
    }
}


void DecomposeWindow::fill(int i, Matrix matrix){
    if (i < 0 || i >= TABLE_COUNT) return;

    for (int r = 0; r < matrix.rows; r++) {
        for (int c = 0; c < matrix.cols; c++){
            double e_value = matrix[r][c];
            if (is_zero(e_value)) e_value = 0;
            QTableWidgetItem* e = new QTableWidgetItem(QString::number(e_value));
            matrix_tables[i]->setItem(r, c, e);
        }
    }
}
