#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.setWindowTitle("Первая краевая задача. Применение метода баланса.");
    w.show();

    return a.exec();
}
