#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "numbersolution.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <QMessageBox>
#include <QStandardItemModel>
#include <QScrollArea>
#include <ctime>

int show_status = 0;

using namespace std;

//U"xx + U"yy = -f(x, y)
//a <= x <= b, c <= y <= d
//U(a, y) = M1(y), U(b, y) = M2(y)
//U(x, c) = M3(x), U(x, d) = M3(x)

int MainWindow::main(void)
{
    setlocale(LC_ALL, "rus");

    int		Nmax = 10000;
    int		S = 0;
    double	eps = 0.000000001;
    double	epsMax = 0;
    int		n = 0, m = 0;
    double	**V = NULL;
    double	**F = NULL;
    double	a, b, c, d;
    double start, finish;
    int Exit = 1, Show = 0;

    a = 0;
    b = 1;
    c = 0;
    d = 3;

    n = 0;
    m = 0;
    Nmax = 0;

    while (n <= 0)
    {
        n = ui->spinBox_3->value();
    }
    while (m <= 0)
    {
        m = ui->spinBox_4->value();
    }
    while (Nmax <= 0)
    {
        Nmax = ui->spinBox_5->value();
    }

    start = clock();
    V = MemoryAllocator(n + 1, m + 1);
    QStandardItemModel *model2 = new QStandardItemModel;
    QStandardItem *item2;

    //Заголовки столбцов
    QStringList horizontalHeader2;
    //Заголовки строк
    QStringList verticalHeader2;
    model2->setHorizontalHeaderLabels(horizontalHeader2);
    model2->setVerticalHeaderLabels(verticalHeader2);

    FillStartSolution(V, n, m, a, b, c, d);

    {
        if (show_status == 1 || show_status == 3)
        {
            for (int j = 0; j <= m; j++)
            {
                for (int i = 0; i <= n; i++)
                {
                    item2 = new QStandardItem(QString::number(V[i][j])); // ZeidelsMethod
                    model2->setItem(i, j, item2);
                }
            }
        }
    }
    ui->tableView->setModel(model2);
    ui->tableView->resizeRowsToContents();
    ui->tableView->resizeColumnsToContents();

    ZeidelsMethod(V, n, m, a, b, c, d, eps, Nmax, epsMax, S);
    finish = clock();

    QStandardItemModel *model1 = new QStandardItemModel;
    QStandardItem *item;

    //Заголовки столбцов
    QStringList horizontalHeader;
    //Заголовки строк
    QStringList verticalHeader;

    model1->setHorizontalHeaderLabels(horizontalHeader);
    model1->setVerticalHeaderLabels(verticalHeader);

    {
        if (show_status == 2 || show_status == 3)
        {
            for (int j = 0; j <= m; j++)
            {
                for (int i = 0; i <= n; i++)
                {
                    item = new QStandardItem(QString::number(V[i][j])); // ZeidelsMethod
                    model1->setItem(i, j, item);
                }
            }
        }
    }
    ui->tableView_2->setModel(model1);
    ui->tableView_2->resizeRowsToContents();
    ui->tableView_2->resizeColumnsToContents();

    MemoryCleaner(V, n);
    return 0;
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QStandardItemModel *model1 = new QStandardItemModel;
    QStandardItemModel *model2 = new QStandardItemModel;
    QStandardItem *item;

    //Заголовки столбцов
    QStringList horizontalHeader;

    QStringList horizontalHeader2;

    //Заголовки строк
    QStringList verticalHeader, verticalHeader2;

    model1->setHorizontalHeaderLabels(horizontalHeader);
    model2->setHorizontalHeaderLabels(horizontalHeader2);
    model1->setVerticalHeaderLabels(verticalHeader);
    model2->setVerticalHeaderLabels(verticalHeader2);

    for (int i = 0; i < 3; i++)
     for (int j = 0; j < 4; j++) {
        item = new QStandardItem("");
        model2->setItem(i, j, item);
    }

    for (int i = 0; i < 3; i++)
     for (int j = 0; j < 4; j++) {
        item = new QStandardItem("");
        model1->setItem(i, j, item);
    }

    ui->tableView->setModel(model2);

    ui->tableView->resizeRowsToContents();
    ui->tableView->resizeColumnsToContents();

    ui->tableView_2->setModel(model1);

    ui->tableView_2->resizeRowsToContents();
    ui->tableView_2->resizeColumnsToContents();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    if (show_status == 2)
        show_status = 3; // test
    else show_status = 1;
    main();
}

void MainWindow::on_pushButton_2_clicked()
{
    if (show_status == 1)
        show_status = 3; // test
    else show_status = 2;
    main();
}
