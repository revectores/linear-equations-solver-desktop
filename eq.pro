QT += widgets
QT += charts
QT += core

LIBS += -lmatrix -L/usr/local/lib
INCLUDEPATH += /usr/local/include
INCLUDEPATH += "include"
INCLUDEPATH += "demo_eqs"
INCLUDEPATH += "utils"
INCLUDEPATH += "dashboard"
INCLUDEPATH += "decompose_window"

SOURCES += "main.cpp"
SOURCES += "demo_eqs/demo_eqs.h"
SOURCES += "utils/utils.h"
SOURCES += "utils/utils.cpp"
SOURCES += "dashboard/dashboard.h"
SOURCES += "dashboard/dashboard.cpp"
SOURCES += "decompose_window/decompose_window.h"
SOURCES += "decompose_window/decompose_window.cpp"

