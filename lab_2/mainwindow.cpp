#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "all_functions.h"
#include <QMessageBox>
//#include <QPixmap>

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
//    if(ui->tab_2->isHidden())
    else
        CalculateMainProblem();
}

void MainWindow::CalculateMainProblem()
{
    int n = ui->doubleSpinBox_3->value();
    double min_step_main = ui->doubleSpinBox_10->value();
    double max_step_main = ui->doubleSpinBox_9->value();
    double step_main = ui->doubleSpinBox_5->value();
    double mult_const_main = ui->doubleSpinBox_4->value();
    double integral_points_main = ui->doubleSpinBox_11->value();

    double h = step_main / n;
    //    double h = 1.0 / n;
    double s = h / 2;
    double *x;
    double *x2;
    double *fi;
    double *d;
    double *a;
    double *v;
    double *v2;
    double *alpha, *beta;
    double Epsilon2 = 0.0;
    double xEpsilon;
    double *ind;
    x = new double[n + 1];
    fi = new double[n];
    d = new double[n];
    a = new double[n + 1];
    v = new double[n + 1];
    alpha = new double[n + 1];
    beta = new double[n + 1];
    v[0] = 1.0;
    v[n] = 0.0;
    a[0] = -1.0;
    d[0] = -1.0;
    fi[0] = -1.0;
    alpha[0] = -1.0;
    beta[0] = -1.0;
    alpha[1] = 0.0;
    beta[1] = 1.0;

    for (int i = 0; i <= n; i++)
    {
        x[i] = i * h;
    }

    for (int i = 1; i <= n - 1; i++)
    {
        if (x[i] + s <= (M_PI / 4))
        {
            fi[i] = f1(x[i]);
            d[i] = q1(x[i]);
        }
        else
        {
            if (x[i] - s >= (M_PI / 4))
            {
                fi[i] = f2(x[i]);
                d[i] = q2(x[i]);
            }
            else
            {
                fi[i] = (1 / h) * (f1(((M_PI / 4) + x[i] - s) / 2.0) * ((M_PI / 4) - x[i] + s) + f2((x[i] + s + (M_PI / 4)) / 2.0) * (x[i] + s - (M_PI / 4)));
                d[i] = (1 / h) * (q1(((M_PI / 4) + x[i] - s) / 2.0) * ((M_PI / 4) - x[i] + s) + q2(((M_PI / 4) + x[i] + s) / 2.0) * (x[i] + s - (M_PI / 4)));
            }
        }
    }

    for (int i = 1; i <= n; i++)
    {
        if (x[i] <= (M_PI / 4))
        {
            a[i] = k1(x[i] - s);
        }
        else
        {
            if (x[i - 1] >= (M_PI / 4))
            {
                a[i] = k2(x[i] - s);
            }
            else
            {
                a[i] = pow((1 / h) * ((1 / k1(((M_PI / 4) + x[i - 1]) / 2.0)) * ((M_PI / 4) - x[i - 1]) + (1 / k2(((M_PI / 4) + x[i]) / 2.0)) * (x[i] - (M_PI / 4))), -1);
            }
        }
    }

    for (int i = 1; i <= n - 1; i++)
    {
        double Ai = a[i] / (h * h);
        double Bi = a[i + 1] / (h * h);
        double Ci = ((a[i] + a[i + 1]) / (h * h)) + d[i];

        alpha[i + 1] = Bi / (Ci - alpha[i] * Ai);
        beta[i + 1] = (fi[i] + Ai * beta[i]) / (Ci - alpha[i] * Ai);
    }

    for (int i = n - 1; i >= 1; i--)
    {
        v[i] = alpha[i + 1] * v[i + 1] + beta[i + 1];
    }

    delete[]fi;
    delete[]a;
    delete[]d;
    delete[]alpha;
    delete[]beta;

    n = 2 * n;
    h = 1.0 / n;
    s = h / 2;

    ind = new double[n + 1];
    x2 = new double[n + 1];
    v2 = new double[n + 1];
    fi = new double[n];
    d = new double[n];
    a = new double[n + 1];
    alpha = new double[n + 1];
    beta = new double[n + 1];

    a[0] = -1.0;
    d[0] = -1.0;
    fi[0] = -1.0;
    alpha[0] = -1.0;
    beta[0] = -1.0;
    alpha[1] = 0.0;
    beta[1] = 1.0;
    v2[0] = 1.0;
    v2[n] = 0.0;

    for (int i = 0; i <= n; i++)
    {
        x2[i] = i * h;
    }

    for (int i = 1; i <= n - 1; i++)
    {
        if (x2[i] + s <= (M_PI / 4))
        {
            fi[i] = f1(x2[i]);
            d[i] = q1(x2[i]);
        }
        else
        {
            if (x2[i] - s >= (M_PI / 4))
            {
                fi[i] = f2(x2[i]);
                d[i] = q2(x2[i]);
            }
            else
            {
                fi[i] = (1 / h) * (f1(((M_PI / 4) + x2[i] - s) / 2.0) * ((M_PI / 4) - x2[i] + s) + f2((x2[i] + s + (M_PI / 4)) / 2.0) * (x2[i] + s - (M_PI / 4)));
                d[i] = (1 / h) * (q1(((M_PI / 4) + x2[i] - s) / 2.0) * ((M_PI / 4) - x2[i] + s) + q2(((M_PI / 4) + x2[i] + s) / 2.0) * (x2[i] + s - (M_PI / 4)));
            }
        }
    }

    for (int i = 1; i <= n; i++)
    {
        if (x2[i] <= (M_PI / 4))
        {
            a[i] = k1(x2[i] - s);
        }
        else
        {
            if (x2[i - 1] >= (M_PI / 4))
            {
                a[i] = k2(x2[i] - s);
            }
            else
            {
                a[i] = pow((1 / h) * ((1 / k1(((M_PI / 4) + x2[i - 1]) / 2.0)) * ((M_PI / 4) - x2[i - 1]) + (1 / k2(((M_PI / 4) + x2[i]) / 2.0)) * (x2[i] - (M_PI / 4))), -1);
            }
        }
    }

    for (int i = 1; i <= n - 1; i++)
    {
        double Ai = a[i] / (h * h);
        double Bi = a[i + 1] / (h * h);
        double Ci = ((a[i] + a[i + 1]) / (h * h)) + d[i];

        alpha[i + 1] = Bi / (Ci - alpha[i] * Ai);
        beta[i + 1] = (fi[i] + Ai * beta[i]) / (Ci - alpha[i] * Ai);
    }


    for (int i = n - 1; i >= 1; i--)
    {
        v2[i] = alpha[i + 1] * v2[i + 1] + beta[i + 1];
    }

    n = n / 2;
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

    QPen pen1, pen2;
    pen1.setColor(Qt::blue);
    pen2.setColor(Qt::red);
    ui->widget_2->addGraph();
    ui->widget_2->graph(0)->setPen(pen1);
//    ui->widget_2->graph(0)->removeData(0.0, 0.0);
    ui->widget_2->graph(0)->setLineStyle((QCPGraph::LineStyle)QCPGraph::lsLine);
    ui->widget_2->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    QVector<double> fir(n + 1), sec(n + 1), thir(n + 1);
    for (int i = 0; i <= n; i++)
    {
        if (abs(v[i] - v2[2 * i]) > Epsilon2)
        {
            Epsilon2 = abs(v[i] - v2[2 * i]);
            xEpsilon = x[i];
        }

        for (int j = 0; j < 4; j++) {
           if (j == 0)
               item = new QStandardItem(QString::number(x[i]));
           if (j == 1)
               item = new QStandardItem(QString::number(v[i]));
           if (j == 2)
               item = new QStandardItem(QString::number(v2[2 * i]));
           if (j == 3)
               item = new QStandardItem(QString::number(abs(v[i] - v2[2 * i])));
           model1->setItem(i, j, item);
        }
        ind[i] = i;
        fir.push_back(i);
        if (i == 0)
        {
            if (v2[i] > 0 && v[i] > 0)
            {
                sec.push_back(v[i]);
                thir.push_back(v2[i]);
            }
        }
        else
        {
            sec.push_back(v[i]);
            thir.push_back(v2[i]);
        }
    }
    ui->tableView_2->setModel(model1);

    ui->tableView_2->resizeRowsToContents();
    ui->tableView_2->resizeColumnsToContents();

    ui->widget_2->graph(0)->setData(fir, sec);
//    ui->widget_2->graph(0)->removeData(0.0, 0.0);
//    ui->widget_2->graph(0)->removeDataBefore(0.1);
    ui->widget_2->xAxis->setLabel("y");
    ui->widget_2->yAxis->setLabel("x");
    ui->widget_2->xAxis->setRange(0, n + 1);
    ui->widget_2->yAxis->setRange(0, n + 1);
    ui->widget_2->replot();

    ui->widget_2->addGraph();
    ui->widget_2->graph(1)->setPen(pen2);
//    ui->widget_2->graph(1)->removeData(0.0, 0.0);
    ui->widget_2->graph(1)->setLineStyle((QCPGraph::LineStyle)QCPGraph::lsLine);
    ui->widget_2->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->widget_2->graph(1)->setData(fir, thir);
//    ui->widget_2->graph(1)->removeData(0.0, 0.0);
    ui->widget_2->xAxis->setLabel("y");
    ui->widget_2->yAxis->setLabel("x");
    ui->widget_2->xAxis->setRange(-1, 5);
    ui->widget_2->yAxis->setRange(-1, 5);
    ui->widget_2->replot();

    ui->widget_2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
            QCP::iSelectLegend | QCP::iSelectPlottables);

    ui->plainTextEdit->document()->setPlainText("Основная задача : \n\n - Для решения задачи использована равномерная сетка с числом разбиений n = " + QString::number(n) + "; \n - Задача должна быть решена с погрешностью не более E=0.5*10^-6; \n - Задача решена с точностью E2 = " + QString::number(Epsilon2) + "; \n - Максимальная  разность численных решений в общих узлах сетки наблюдается в точке x = " + QString::number(xEpsilon));

    delete[]a;
    delete[]ind;
    delete[]fi;
    delete[]d;
    delete[]x;
    delete[]x2;
    delete[]v;
    delete[]v2;
    delete[]alpha;
    delete[]beta;
}

void MainWindow::CalculateTestProblem()
{
    int n = ui->doubleSpinBox_3->value();

    // Получение параметров из GUI
    double min_step_test = ui->doubleSpinBox_20->value();
    double max_step_test = ui->doubleSpinBox_21->value();
    double step_test = ui->doubleSpinBox_18->value();
    double mult_const_test = ui->doubleSpinBox_19->value();
    double integral_points_test = ui->doubleSpinBox_17->value();

//    int min_step_main = ui->doubleSpinBox_10->value();
//    int max_step_main = ui->doubleSpinBox_9->value();
//    int step_main = ui->doubleSpinBox_5->value();
//    int mult_const_main = ui->doubleSpinBox_4->value();
//    int integral_points_main = ui->doubleSpinBox_11->value();
    // Получение параметров из GUI

    double h = step_test / n;
//    double h = 1.0 / n;
    double s = h / 2;
    double *x;
    double *fi;
    double *d;
    double *a;
    double *u;
    double *v;
    double *alpha, *beta;
    double Epsilon1 = 0.0;
    double xEpsilon;
    x = new double[n + 1];
    fi = new double[n];
    d = new double[n];
    a = new double[n + 1];
    v = new double[n + 1];
    u = new double[n + 1];
    alpha = new double[n + 1];
    beta = new double[n + 1];
    v[0] = 1.0;
    v[n] = 0.0;
    a[0] = -1.0;
    d[0] = -1.0;
    fi[0] = -1.0;
    alpha[0] = -1.0;
    beta[0] = -1.0;
    alpha[1] = 0.0;
    beta[1] = 1.0;

    for (int i = 0; i <= n; i++)
    {
      x[i] = i * h;
    }

    for (int i = 1; i <= n - 1; i++)
    {
        if (x[i] + s <= (M_PI / 4))
        {
         fi[i] = 1;
         d[i] = (M_PI / 4);
        }
        else
        {
            if (x[i] - s >= (M_PI / 4))
            {
                fi[i] = sin((M_PI / 4)*M_PI);
                d[i] = 0.09;
            }
            else
            {
                fi[i] = (1 / h) * (((M_PI / 4) - x[i] + s) + sin((M_PI / 4) * M_PI) * (x[i] + s - (M_PI / 4)));
                d[i] = (1 / h) * ((M_PI / 4) * ((M_PI / 4) - x[i] + s) + 0.09 * (x[i] + s - (M_PI / 4)));
            }
        }
    }

    for (int i = 1; i <= n; i++)
    {
        if (x[i] <= (M_PI / 4))
        {
            a[i] = 2.09;
        }
        else
        {
            if (x[i - 1] >= (M_PI / 4))
            {
                a[i] = 0.09;
            }
            else
            {
                a[i] = pow((1 / h)*((1 / 2.09)*((M_PI / 4) - x[i - 1]) + (1 / 0.09)*(x[i] - (M_PI / 4))), -1);
            }
        }
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

    for (int i = 1; i <= n - 1; i++)
    {
        double Ai = a[i] / (h*h);
        double Bi = a[i + 1] / (h*h);
        double Ci = ((a[i] + a[i + 1]) / (h*h)) + d[i];

        alpha[i + 1] = Bi / (Ci - alpha[i] * Ai);
        beta[i + 1] = (fi[i] + Ai*beta[i]) / (Ci - alpha[i] * Ai);
    }

    for (int i = n - 1; i >= 1; i--)
    {
        v[i] = alpha[i + 1] * v[i + 1] + beta[i + 1];
    }

    QPen pen1, pen2;
    pen1.setColor(Qt::blue);
    pen2.setColor(Qt::red);
    ui->widget->addGraph();
    ui->widget->graph(0)->setPen(pen1);
    ui->widget->graph(0)->setLineStyle((QCPGraph::LineStyle)QCPGraph::lsLine);
    ui->widget->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    QVector<double> fir(n + 1), sec(n + 1), thir(n + 1);

    for (int i = 0; i <= n ; i++)
    {
        if (x[i] <= (M_PI / 4))
        {
            u[i] = u1(x[i]);
        }
        else
        {
            u[i] = u2(x[i]);
        }

        if (abs(u[i] - v[i]) > Epsilon1)
        {
            Epsilon1 = abs(u[i] - v[i]);
            xEpsilon = x[i];
        }

         for (int j = 0; j < 4; j++) {
             if (j == 0)
                item = new QStandardItem(QString::number(x[i]));
             if (j == 1)
                item = new QStandardItem(QString::number(u[i]));
             if (j == 2)
                item = new QStandardItem(QString::number(v[i]));
             if (j == 3)
                item = new QStandardItem(QString::number(abs(u[i] - v[i])));
            model2->setItem(i, j, item);
        }
        fir.push_back(i);
        if ((u[i] > 0 && v[i] > 0 && i == 0)||(i > 0))
        if (i == 0)
        {
            if (u[i] > 0 && v[i] > 0)
            {
                sec.push_back(u[i]);
                thir.push_back(v[i]);
            }
        }
        else
        {
            sec.push_back(u[i]);
            thir.push_back(v[i]);
        }
    }
    ui->tableView->setModel(model2);

    ui->tableView->resizeRowsToContents();
    ui->tableView->resizeColumnsToContents();
//    ui->widget->graph(0)->removeData(0, 0);
    ui->widget->graph(0)->setData(fir, sec);
//    ui->widget->graph(0)->removeData(0, 0);
    ui->widget->xAxis->setLabel("y");
    ui->widget->yAxis->setLabel("x");
    ui->widget->xAxis->setRange(0, n + 1);
    ui->widget->yAxis->setRange(0, n + 1);
    ui->widget->replot();

    ui->widget->addGraph();
    ui->widget->graph(1)->setPen(pen2);
    ui->widget->graph(1)->setLineStyle((QCPGraph::LineStyle)QCPGraph::lsLine);
    ui->widget->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
//    ui->widget_2->graph(1)->removeData(0, 0);
    ui->widget->graph(1)->setData(fir, thir);
//    ui->widget_2->graph(1)->removeData(0, 0);
    ui->widget->xAxis->setLabel("y");
    ui->widget->yAxis->setLabel("x");
    ui->widget->xAxis->setRange(-1, 5);
    ui->widget->yAxis->setRange(-1, 5);
    ui->widget->replot();

    ui->widget->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
            QCP::iSelectLegend | QCP::iSelectPlottables);

    ui->plainTextEdit->document()->setPlainText("Тестовая задача : \n\n - Для решения задачи использована равномерная сетка с числом разбиений n = " + QString::number(n) + "; \n - Задача должна быть решена с погрешностью не более E=0.5*10^-6; \n - Задача решена с точностью E1 = " + QString::number(Epsilon1) + "; \n - Максимальная  разность численных решений в общих узлах сетки наблюдается в точке x = " + QString::number(xEpsilon));

    delete[]x;
    delete[]fi;
    delete[]d;
    delete[]a;
    delete[]u;
    delete[]v;
    delete[]alpha;
    delete[]beta;
}
