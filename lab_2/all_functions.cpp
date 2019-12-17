//#include "all_functions.h"

//double u1(double x)
//{
//    // (M_PI / 4)
//    double a1 = exp(sqrt((M_PI / 4) / 2.09) * (M_PI / 4));
//    double a2 = exp(-sqrt((M_PI / 4) / 2.09) * (M_PI / 4));
//    double a3 = 1 / (M_PI / 4);
//    double a4 = sin((M_PI / 4) * M_PI) / 0.09;
//    double a5 = exp((M_PI / 4));
//    double a6 = exp(-(M_PI / 4));
//    double a7 = sqrt(2.09 * (M_PI / 4));
//    double c4 = ((a2 - a1) * (a7 * a1 * (a3 - 1) - 0.09 * a5 * a4 / exp(1)) + a7 * (a1 * (a4 * (1 - (a5 / exp(1))) + a1 * (a3 - 1) - a3) + a2 * (a4 * (1 - (a5 / exp(1))) + a1 * (a3 - 1) - a3))) / (0.09 * (a2 - a1) * ((a5 / exp(2)) + a6) - a7 * (a6 - (a5 / exp(2))) * (a1 + a2));
//    double c2 = (c4 * (a6 - a5 / exp(2)) + a4 * (1 - a5 / exp(1)) + a1 * (a3 - 1) - a3) / (a2 - a1);
//    double c1 = 1 - c2 - a3;

//    return c1 * exp(sqrt((M_PI / 4) / 2.09) * x) + c2 * exp(-sqrt((M_PI / 4) / 2.09) * x) + a3;
//}

//double u2(double x)
//{
//    double a1 = exp(sqrt((M_PI / 4) / 2.09) * (M_PI / 4));
//    double a2 = exp(-sqrt((M_PI / 4) / 2.09) * (M_PI / 4));
//    double a3 = 1 / (M_PI / 4);
//    double a4 = sin((M_PI / 4) * M_PI) / 0.09;
//    double a5 = exp((M_PI / 4));
//    double a6 = exp(-(M_PI / 4));
//    double a7 = sqrt(2.09 * (M_PI / 4));
//    double c4 = ((a2 - a1) * (a7 * a1 * (a3 - 1) - 0.09 * a5 * a4 / exp(1)) + a7 * (a1 * (a4 * (1 - (a5 / exp(1))) + a1 * (a3 - 1) - a3) + a2 * (a4 * (1 - (a5 / exp(1))) + a1 * (a3 - 1) - a3))) / (0.09 * (a2 - a1) * ((a5 / exp(2)) + a6) - a7 * (a6 - (a5 / exp(2))) * (a1 + a2));
//    double c3 = -c4 / exp(2) - a4 / exp(1);

//    return c3 * exp(x) + c4 * exp(-x) + a4;
//}

//double k1(double x)
//{
//    return sqrt(2) * sin(x);
//}

//double k2(double x)
//{
//    return pow(cos(x), 2);
//}

//double q1(double x)
//{
//    return 1;
//}

//double q2(double x)
//{
//    return x*x;
//}

//double f1(double x)
//{
//    return sin(2 * x);
//}

//double f2(double x)
//{
//    return cos(x);
//}
