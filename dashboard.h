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

#include <string>
#include "utils.h"
#include "eq_test.h"
#include "decompose_window.h"



extern std::vector<Equation> TEST_EQUATIONS;


class Dashboard : public QWidget {
private:
    DecomposeWindow decompose_window;
    QTableWidget* matrix_table;
    QTableWidget* constant_table;
    QTableWidget* solution_table;
    QGroupBox*    method_group;
    QRadioButton* method_radios[6];
    QGroupBox*    action_group;
    QLabel*       info_label;
    int           cur_demo_id;

    void matrix_table_init();
    void constant_table_init();
    void solution_table_init();
    void method_group_init();
    void action_group_init();
    void info_label_init();
    Equation read_equation();
    void fill_solution(Matrix solution);

    void demo(int demo_id);
    void clear();
    void solve();
    void err_msg(std::string msg);
    void info_msg(std::string msg);

public:
    Dashboard(QWidget *parent = 0);
};

