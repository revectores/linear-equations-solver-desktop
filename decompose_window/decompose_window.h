#include <QApplication>
#include <QWidget>
#include <QGridLayout>
#include <QGroupBox>

#include <QLabel>
#include <QRadioButton>
#include <QCheckBox>
#include <QPushButton>
#include <QTableWidget>
#include <QtCharts>

#include <QMessageBox>
#include "matrix.h"
#include "utils.h"


const int TABLE_COUNT = 3;


class DecomposeWindow : public QWidget {
private:
    void tables_init();
    QTableWidget* matrix_tables[TABLE_COUNT];

public:
    DecomposeWindow(QWidget *parent = 0);
    void clear();
    void fill(int i, Matrix matrix);
};

