#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <math.h>
#define _USE_MATH_DEFINES

double k1(double x) {
    return sqrt(2) * (sin(x) + 1);
}
double k2(double x) {
    return pow(cos(x), 2);
}
double q1(double x) {
    return 1.0;
}
double q2(double x) {
    return x * x;
}
double f1(double x) {
    return sin(2 * x);
}
double f2(double x) {
    return cos(x);
}
double fi1(double x) {
    return -cos(2 * x) / 2;
}
double fi2(double x) {
    return sin(x);
}
double d1(double x) {
    return x;
}
double d2(double x) {
    return x * x * x / 3;
}
double a_1(double x) {
    //return Math.Log(Math.Abs(Math.Tan(x * .5))) / sqrt(2);
    //изменил k1(x) = 2^.5*sin(x) на k1(x) = 2^.5*(sin(x)+1)
    return -sqrt(2) / (tan(x / 2) + 1);
}
double a_2(double x) {
    return tan(x);
}

class Widget : public QWidget {
protected:
  void paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.drawPixmap(0, 0, QPixmap("/home/dariamityagina/ChM_lab_work_2/Numerical-methods/lab_2/srcspravka.jpeg"));
  }
};

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
       QStandardItemModel *model1 = new QStandardItemModel;
       QStandardItemModel *model2 = new QStandardItemModel;
       QStandardItem *item;
       ui->tabWidget->setStyleSheet("color: rgb(255, 255, 255);");

       //Заголовки столбцов
       QStringList horizontalHeader;
       horizontalHeader.append("       xi      ");
       horizontalHeader.append("      v(xi)    ");
       horizontalHeader.append("     v2(x2i)   ");
       horizontalHeader.append("v(xi) − v2(x2i)");

       QStringList horizontalHeader2;
       horizontalHeader2.append("       xi      ");
       horizontalHeader2.append("      u(xi)    ");
       horizontalHeader2.append("     v2(xi)    ");
       horizontalHeader2.append(" u(xi) − v(xi) ");

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

void MainWindow::on_pushButton_5_clicked()
{
   QWidget* Form = new QWidget;
   Form->setWindowTitle("Справка");
   Form->setAttribute(Qt::WA_DeleteOnClose, true);

   QPixmap pix("/home/dariamityagina/ChM_lab_work_2/Numerical-methods/lab_2/src/new.jpeg");
   Form->resize(1020, 443);
   QPalette palette;
   palette.setBrush(QPalette::Background, pix);
   Form->setPalette(palette);
   Form->show();
}

void MainWindow::on_pushButton_2_clicked()
{
    if(ui->tab->isActiveWindow() && ui->tab_2->isHidden())
        CalculateTestProblem();
    else
        CalculateMainProblem();
}

void MainWindow::CalculateMainProblem()
{
    double ksi = M_PI / 4;
    int n = ui->doubleSpinBox_3->value();
//    double step_main = ui->doubleSpinBox_5->value();
//    double h = step_main / (double)n;
    double h = 1.0 / (double)n;
    double* x = new double[n + 1];
    double mu1 = 1.0, mu2 =0.0;
    x[0] = 0.0;
    x[n] = 1.0;
    for (int i = 1; i < n; i++) {
        x[i] = (double)i * h;
    }
    double* xs = new double[n];
    for (int i = 0; i < n; i++) {
        xs[i] = ((double)i + 0.5) * h;
    }

    double* a = new double[n];
    for (int i = 0; i < n; i++) {
        if (x[i + 1] < ksi && x[i] < ksi) {
            a[i] = h / (a_1(x[i + 1]) - a_1(x[i]));
        }
        else if (x[i + 1] > ksi && x[i] < ksi) {
            a[i] = h / (a_2(x[i + 1]) - a_2(ksi) + a_1(ksi) - a_1(x[i]));
        }
        else if (x[i + 1] > ksi && x[i] > ksi) {
            a[i] = h / (a_2(x[i + 1]) - a_2(x[i]));
        }
    }

    double* fi = new double[n - 1];
    for (int i = 0; i < n - 1; i++) {
        if (xs[i + 1] < ksi && xs[i] < ksi) {
            fi[i] = (fi1(xs[i + 1]) - fi1(xs[i])) / h;
        }
        else if (xs[i + 1] > ksi && xs[i] < ksi) {
            fi[i] = (fi2(xs[i + 1]) - fi2(ksi) + fi1(ksi) - fi1(xs[i])) / h;
        }
        else if (xs[i + 1] > ksi && xs[i] > ksi) {
            fi[i] = (fi2(xs[i + 1]) - fi2(xs[i])) / h;
        }
    }

    double* d = new double[n - 1];
    for (int i = 0; i < n - 1; i++) {
        if (xs[i + 1] < ksi && xs[i] < ksi) {
            d[i] = (d1(xs[i + 1]) - d1(xs[i])) / h;
        }
        else if (xs[i + 1] > ksi && xs[i] < ksi) {
            d[i] = (d2(xs[i + 1]) - d2(ksi) + d1(ksi) - d1(xs[i])) / h;
        }
        else if (xs[i + 1] > ksi && xs[i] > ksi) {
            d[i] = (d2(xs[i + 1]) - d2(xs[i])) / h;
        }
    }

    //прямой ход прогонки
    double* V = new double[n + 1];
    double* alpha = new double[n], *betta = new double[n];
    double A, B, C;
    alpha[0] = 0.0;
    betta[0] = mu1;
    for (int i = 1; i < n; i++) {
        A = a[i - 1] / (h * h);
        C = (a[i - 1] + a[i]) / (h * h) + d[i - 1];
        B = a[i] / (h * h);
        betta[i] = (fi[i - 1] + A * betta[i - 1]) / (C - alpha[i - 1] * A);
        alpha[i] = B / (C - alpha[i - 1] * A);
    }

    //обратный ход
    V[n] = mu2;
    for (int i = n - 1; i >= 0; i--) {
        V[i] = alpha[i] * V[i + 1] + betta[i];
    }

    QPen pen1, pen2;
    pen1.setColor(Qt::blue);
    pen2.setColor(Qt::red);
    ui->widget_2->addGraph();
    ui->widget_2->graph(0)->setPen(pen1);
    ui->widget_2->graph(0)->setLineStyle((QCPGraph::LineStyle)QCPGraph::lsLine);
    ui->widget_2->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 1));
    QVector<double> fir(n + 1), sec(n + 1), thir(n + 1);

    //отрисовка числ решения ~ blue
    for (int i = 0; i < n + 1; i++) {
        fir.push_back(x[i]);
        sec.push_back(V[i]);
    }
    ui->widget_2->graph(0)->setData(fir, sec);
    ui->widget_2->xAxis->setLabel("y");
    ui->widget_2->yAxis->setLabel("x");
    ui->widget_2->xAxis->setRange(0.01, 1.5);
    ui->widget_2->yAxis->setRange(0.01, 1.5);
    ui->widget_2->replot();

    //шаг в 2 раза мельче
    int n2 = 2 * n;
    double h2 = 1 / (double)n2;
    double* x2 = new double[n2 + 1];
    x2[0] =0.0;
    x2[n2] = 1.0;
    for (int i = 1; i < n2; i++) {
        x2[i] = (double)i * h2;
    }
    double* xs2 = new double[n2];
    for (int i = 0; i < n2; i++) {
        xs2[i] = ((double)i + .5) * h2;
    }

    a = new double[n2];
    for (int i = 0; i < n2; i++) {
        if (x2[i + 1] < ksi && x2[i] < ksi) {
            a[i] = h2 / (a_1(x2[i + 1]) - a_1(x2[i]));
        }
        else if (x2[i + 1] > ksi && x2[i] < ksi) {
            a[i] = h2 / (a_2(x2[i + 1]) - a_2(ksi) + a_1(ksi) - a_1(x2[i]));
        }
        else if (x2[i + 1] > ksi && x2[i] > ksi) {
            a[i] = h2 / (a_2(x2[i + 1]) - a_2(x2[i]));
        }
    }

    fi = new double[n2 - 1];
    for (int i = 0; i < n2 - 1; i++) {
        if (xs2[i + 1] < ksi && xs2[i] < ksi) {
            fi[i] = (fi1(xs2[i + 1]) - fi1(xs2[i])) / h2;
        }
        else if (xs2[i + 1] > ksi && xs2[i] < ksi) {
            fi[i] = (fi2(xs2[i + 1]) - fi2(ksi) + fi1(ksi) - fi1(xs2[i])) / h2;
        }
        else if (xs2[i + 1] > ksi && xs2[i] > ksi) {
            fi[i] = (fi2(xs2[i + 1]) - fi2(xs2[i])) / h2;
        }
    }

    d = new double[n2 - 1];
    for (int i = 0; i < n2 - 1; i++) {
        if (xs2[i + 1] < ksi && xs2[i] < ksi) {
            d[i] = (d1(xs2[i + 1]) - d1(xs2[i])) / h2;
        }
        else if (xs2[i + 1] > ksi && xs2[i] < ksi) {
            d[i] = (d2(xs2[i + 1]) - d2(ksi) + d1(ksi) - d1(xs2[i])) / h2;
        }
        else if (xs2[i + 1] > ksi && xs2[i] > ksi) {
            d[i] = (d2(xs2[i + 1]) - d2(xs2[i])) / h2;
        }
    }

    //прямой ход прогонки
    double* V2 = new double[n2 + 1];
    alpha = new double[n2 + 1];
    betta = new double[n2 + 1];
    // alpha[0] =0.0;
    // betta[0] = mu1;
    alpha[1] =0.0;
    betta[1] = mu1;

    for (int i = 1; i < n2; i++) {
        // A = a[i - 1] / (h2 * h2);
        // C = (a[i - 1] + a[i]) / (h2 * h2) + d[i - 1];
        // B = a[i] / (h2 * h2);
        // betta[i] = (fi[i - 1] + A * betta[i - 1]) / (C - alpha[i - 1] * A);
        // alpha[i] = B / (C - alpha[i - 1] * A);
        A = a[i] / (h * h);
        B = a[i + 1] / (h * h);
        C = (a[i] + a[i + 1]) / (h * h) + d[i];

        alpha[i + 1] = B / (C - alpha[i] * A);
        betta[i] = (fi[i] + A * betta[i]) / (C - alpha[i] * A);
    }

    //обратный ход
    V2[n2] = mu2;
    for (int i = n2 - 1; i >= 0; i--) {
        V2[i] = alpha[i + 1] * V2[i + 1] + betta[i + 1];
    }

    //отрисовка числ решения с шагом в 2 раза мельче ~ red
    // for (int i = 0; i < n2 + 1; i++) {
    //     fir.push_back(x[i]);
    //     thir.push_back(V2[i]);
    // }
    // ui->widget_2->addGraph();
    // ui->widget_2->graph(1)->setPen(pen2);
    // ui->widget_2->graph(1)->setLineStyle((QCPGraph::LineStyle)QCPGraph::lsLine);
    // ui->widget_2->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 1));
    // ui->widget_2->graph(1)->setData(fir, thir);
    // ui->widget_2->xAxis->setLabel("y");
    // ui->widget_2->yAxis->setLabel("x");
    // ui->widget_2->xAxis->setRange(0.01, 1.5);
    // ui->widget_2->yAxis->setRange(0.01, 1.5);
    // ui->widget_2->replot();

    ui->widget_2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
            QCP::iSelectLegend | QCP::iSelectPlottables);

    QStandardItemModel *model1 = new QStandardItemModel;
    QStandardItem *item;
    ui->tabWidget->setStyleSheet("color: rgb(255, 255, 255);");

    //Заголовки столбцов
    QStringList horizontalHeader;
    horizontalHeader.append("        xi       ");
    horizontalHeader.append("       v(xi)     ");
    horizontalHeader.append("      v2(x2i)    ");
    horizontalHeader.append("|v(xi) − v2(x2i)|");

    QStringList horizontalHeader2;
    horizontalHeader2.append("        xi       ");
    horizontalHeader2.append("       u(xi)     ");
    horizontalHeader2.append("      v2(xi)     ");
    horizontalHeader2.append(" |u(xi) − v(xi)| ");

    //Заголовки строк
    QStringList verticalHeader;

    model1->setHorizontalHeaderLabels(horizontalHeader);
    model1->setVerticalHeaderLabels(verticalHeader);
    //заполняем dataGrid
    double* VV = new double[n + 1];
    double max =0.0;
    int jj = 0;
    double x_needed;
    for (int i = 0; i < n + 1; i++) {
        VV[i] = abs(V2[2 * i] - V[i]);
        if (VV[i] > max) {
            max = VV[i];
            x_needed = x[i];
            jj = i;
        }
    }
    for (int i = 0; i < n + 1; i++)
        for (int j = 0; j < 4; j++) {
           if (j == 0)
               item = new QStandardItem(QString::number(x[i]));
           if (j == 1)
               item = new QStandardItem(QString::number(V2[2 * i]));
           if (j == 2)
               item = new QStandardItem(QString::number(V[2 * i]));
           if (j == 3)
               item = new QStandardItem(QString::number(VV[i]));
           model1->setItem(i, j, item);
        }
    ui->tableView_2->setModel(model1);
    ui->tableView_2->resizeRowsToContents();
    ui->tableView_2->resizeColumnsToContents();

    delete []a;
    delete []fi;
    delete []d;
    delete []x;
    delete []x2;
    delete []xs;
    delete []V;
    delete []alpha;
    delete []betta;
    delete []V2;
    delete []VV;
    ui->plainTextEdit->document()->setPlainText("Основная задача : \n\n - Для решения задачи использована равномерная сетка с числом разбиений n = " + QString::number(n) + "; \n - Задача должна быть решена с погрешностью не более E=0.5*10^-6; \n - Задача решена с точностью E2 = " + QString::number(max) + "; \n - Максимальная  разность численных решений в общих узлах сетки наблюдается в точке x = " + QString::number(x_needed) + " - что соответствует узлу " + QString::number(jj + 1));
}

void MainWindow::CalculateTestProblem()
{
    double ksi = M_PI / 4;
    int n = ui->doubleSpinBox_3->value();

    double h = 1.00 / (double)n;
    double* x = new double[n + 1];
    double mu1 = 1.0, mu2 = 0.0;
    // x[0] = 1;
    x[n] = 1.0;
    for (int i = 0; i < n; i++) {
        x[i] = (double)i * h;
    }
    double* xs = new double[n];
    for (int i = 0; i < n; i++) {
        xs[i] = ((double)i + 0.5) * h;
    }
    double* fi = new double[n - 1];
    for (int i = 0; i < n - 1; i++) {
        if (xs[i + 1] < ksi && xs[i] < ksi) {
            fi[i] = (xs[i + 1] - xs[i]) / h;
        }
        else if (xs[i + 1] > ksi && xs[i] < ksi) {
            fi[i] = (sqrt(2.0) / 2.0 * (xs[i + 1] - ksi) + (ksi - xs[i])) / h;
        }
        else if (xs[i + 1] > ksi && xs[i] > ksi) {
            fi[i] = (sqrt(2.0) / 2.0 * xs[i + 1] - sqrt(2.0) / 2.0 * xs[i]) / h;
        }
    }
    double* a = new double[n];
    for (int i = 0; i < n; i++) {
        if (x[i + 1] < ksi && x[i] < ksi) {
            a[i] = h / (x[i + 1] - x[i]);
        }
        else if (x[i + 1] > ksi && x[i] < ksi) {
            a[i] = h / (2.0 * (x[i + 1] - ksi) + (ksi - x[i]));
        }
        else if (x[i + 1] > ksi && x[i] > ksi) {
            a[i] = h / (2.0 * x[i + 1] - 2.0 * x[i]);
        }
    }
    double* d = new double[n - 1];
    for (int i = 0; i < n - 1; i++) {
        if (xs[i + 1] < ksi && xs[i] < ksi) {
            d[i] = (xs[i + 1] - xs[i]) / h;
        }
        else if (xs[i + 1] > ksi && xs[i] < ksi) {
            d[i] = ((M_PI  * M_PI  / 16.0) * (xs[i + 1] - ksi) + (ksi - xs[i])) / h;
        }
        else if (xs[i + 1] > ksi && xs[i] > ksi) {
            d[i] = ((M_PI  * M_PI  / 16.0) * xs[i + 1] - (M_PI  * M_PI  / 16.0) * xs[i]) / h;
        }
    }
    double* V = new double[n + 1];
    double* alpha = new double[n], *betta = new double[n];

    //прямой ход прогонки
    double A, B, C;
    alpha[0] =0.0;
    betta[0] = mu1;
    for (int i = 1; i < n; i++) {
        A = a[i - 1] / (h * h);
        C = (a[i - 1] + a[i]) / (h * h) + d[i - 1];
        B = a[i] / (h * h);
        betta[i] = (fi[i - 1] + A * betta[i - 1]) / (C - alpha[i - 1] * A);
        alpha[i] = B / (C - alpha[i - 1] * A);
    }

    //обратный ход прогонки
    V[n] = mu2;
    for (int i = n - 1; i >= 0; i--) {
        V[i] = alpha[i] * V[i + 1] + betta[i];
    }

    QStandardItemModel *model2 = new QStandardItemModel;
    QStandardItem *item;
    ui->tabWidget->setStyleSheet("color: rgb(255, 255, 255);");

    //Заголовки столбцов
    QStringList horizontalHeader2;
    horizontalHeader2.append("        xi       ");
    horizontalHeader2.append("       u(xi)     ");
    horizontalHeader2.append("      v2(xi)     ");
    horizontalHeader2.append(" |u(xi) − v(xi)| ");

    //Заголовки строк
    QStringList verticalHeader2;

    model2->setHorizontalHeaderLabels(horizontalHeader2);
    model2->setVerticalHeaderLabels(verticalHeader2);

    QPen pen1, pen2;
    pen1.setColor(Qt::blue);
    pen2.setColor(Qt::red);
    ui->widget->addGraph();
    ui->widget->graph(0)->setPen(pen1);
    ui->widget->graph(0)->setLineStyle((QCPGraph::LineStyle)QCPGraph::lsLine);
    ui->widget->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 1));
    QVector<double> fir(n + 1), sec(n + 1), thir(n + 1);

    //отрисовка численного решения ~ blue
    for (int i = 0; i < n + 1; i++) {
        fir.push_back(x[i]);
        sec.push_back(V[i]);
    }
    ui->widget->graph(0)->setData(fir, sec);
    ui->widget->xAxis->setLabel("y");
    ui->widget->yAxis->setLabel("x");
    ui->widget->xAxis->setRange(0.01, 1.5);
    ui->widget->yAxis->setRange(0.01, 1.5);
    ui->widget->replot();

    //счет точного решения
    double c1, c2, c3, c4, tmp_a, tmp_b;
    double e4 = exp((M_PI / 4));
    double e_4 = exp(-(M_PI / 4));
//    double e4 = exp((M_PI / 4));
//    double e_4 = exp(-(M_PI / 4));
    double e_8 = exp(-(M_PI) * (M_PI) / (8 * sqrt(2)));
    double e8_2 = exp((M_PI) * (M_PI) / (8 * sqrt(2)) - (M_PI) / sqrt(2));
    double e8_22 = exp((M_PI) * (M_PI) / (8 * sqrt(2)) - (M_PI) / (2 * sqrt(2)));
    double pi_4 = (M_PI) / (4 * sqrt(2));
    double pi_8 = 8 * sqrt(2) / ((M_PI) * (M_PI));
    tmp_a = -e_8 + e8_2 - pi_4 * (e_8 + e8_2) / (e4 + e_4) * (-e_4 + e4);
    tmp_b = pi_8 - 1 - pi_8 * e8_22 + 2 / (M_PI) * e8_22 / (e4 + e_4) * (-e_4 + e4);
    c1 = -2 / (M_PI) * e8_22 / (e4 + e_4) - tmp_b / tmp_a * pi_4 * (e_8 + e8_2) / (e4 + e_4);
    c2 = 2 / (M_PI) * e8_22 / (e4 + e_4) + tmp_b / tmp_a * pi_4 * (e_8 + e8_2) / (e4 + e_4);
    c3 = -pi_8 * exp(-(M_PI) / (2 * sqrt(2))) - tmp_b / tmp_a * exp(-(M_PI) / sqrt(2));
    c4 = tmp_b / tmp_a;
    double* u = new double[n + 1];
    for (int i = 0; i < n + 1; i++) {
        if (x[i] < ksi) {
            u[i] = c1 * exp(x[i]) + c2 * exp(-x[i]) + 1;
        }
        else if (x[i] > ksi) {
            u[i] = c3 * exp((M_PI) / sqrt(8) * x[i]) + c4 * exp(-(M_PI) / sqrt(8) * x[i]) + pi_8;
        }
    }

    //отрисовка точного решения ~ red
    for (int i = 0; i < n + 1; i++) {
        fir.push_back(x[i]);
        thir.push_back(u[i]);
    }
    ui->widget->addGraph();
    ui->widget->graph(1)->setPen(pen2);
    ui->widget->graph(1)->setLineStyle((QCPGraph::LineStyle)QCPGraph::lsLine);
    ui->widget->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 1));
    ui->widget->graph(1)->setData(fir, thir);
    ui->widget->xAxis->setLabel("y");
    ui->widget->yAxis->setLabel("x");
    ui->widget->xAxis->setRange(0.01, 1.5);
    ui->widget->yAxis->setRange(0.01, 1.5);
    ui->widget->replot();

    ui->widget->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
            QCP::iSelectLegend | QCP::iSelectPlottables);
    //заполняем dataGrid
    double* uV = new double[n + 1];
    double max =0.0;
    int jj = 0;
    double x_needed;
    for(int i = 0; i < n+1; i++) {
        uV[i] = abs(u[i] - V[i]);
        if (uV[i] > max) {
            max = uV[i];
            x_needed = x[i];
            jj = i;
        }
    }

    for (int i = 0; i < n + 1; i++)
        for (int j = 0; j < 4; j++) {
             if (j == 0)
                item = new QStandardItem(QString::number(x[i]));
             if (j == 1)
                item = new QStandardItem(QString::number(u[i]));
             if (j == 2)
                item = new QStandardItem(QString::number(V[i]));
             if (j == 3)
                item = new QStandardItem(QString::number(uV[i]));
             model2->setItem(i, j, item);
        }
    ui->tableView->setModel(model2);
    ui->tableView->resizeRowsToContents();
    ui->tableView->resizeColumnsToContents();

    delete []a;
    delete []fi;
    delete []d;
    delete []x;
    delete []xs;
    delete []V;
    delete []alpha;
    delete []betta;
    delete []uV;
    delete []u;

    ui->plainTextEdit->document()->setPlainText("Тестовая задача : \n\n - Для решения задачи использована равномерная сетка с числом разбиений n = " + QString::number(n) + "; \n - Задача должна быть решена с погрешностью не более E=0.5*10^-6; \n - Задача решена с точностью E1 = " + QString::number(max) + "; \n - Максимальная  разность численных решений в общих узлах сетки наблюдается в точке x = " + QString::number(x_needed) + " - что соответствует узлу " + QString::number(jj + 1));
}
