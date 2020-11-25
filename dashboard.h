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
#include "utils.h"
#include "eq_test.h"


extern std::vector<Equation> TEST_EQUATIONS;


class Dashboard : public QWidget {
private:
    QTableWidget* matrix_table;
    QTableWidget* constant_table;
    QTableWidget* solution_table;
    QGroupBox*    method_group;
    QRadioButton* method_radios[6];
    QGroupBox*    action_group;
    QLabel*       decompose_label;
    int           cur_demo_id;

    void matrix_table_init();
    void constant_table_init();
    void solution_table_init();
    void method_group_init();
    void action_group_init();
    void decompose_label_init();
    void read_equation();

    void demo();
    void clear();
    void solve();


public:
    Dashboard(QWidget *parent = 0);
};

