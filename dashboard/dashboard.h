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
#include "decompose_window.h"


// extern std::vector<Equation> DEMO_EQUATIONS;
const int METHOD_COUNT = 9;


class Dashboard : public QWidget {
private:
    DecomposeWindow decompose_window;
    QTableWidget* matrix_table;
    QTableWidget* constant_table;
    QTableWidget* solution_table;
    QGroupBox*    method_group;
    QGroupBox*    parameter_group;
    QLineEdit*    precision_edit;
    QLineEdit*    relaxation_factor_edit;
    QRadioButton* method_radios[METHOD_COUNT];
    QGroupBox*    action_group;
    QLabel*       info_label;
    int           cur_demo_id;
    int           cur_method_id;

    void matrix_table_init();
    void constant_table_init();
    void solution_table_init();
    void method_group_init();
    void parameter_group_init();
    void action_group_init();
    void info_label_init();
    Equation read_equation();
    Matrix read_init_vector();
    long double read_precision();
    long double read_relaxation_factor();
    void fill_solution(Matrix solution);

    int get_cur_method_type();
    void method_radio_checked(int method_radio_id);
    void demo(int demo_id);
    void clear_matrix_table();
    void clear_solution_table();
    void clear_constant_table();
    void clear();
    void solve();
    void err_msg(std::string msg);
    void info_msg(std::string msg);

public:
    Dashboard(QWidget *parent = 0);
};

