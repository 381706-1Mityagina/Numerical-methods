#ifndef NUMSOL
#define NUMSOL

#include <iomanip>
#include <iostream>
using namespace std;

//-------------------------------------------------------------------------��������-������--------------------------------------------------------------------------------------
double Uxy(double x, double y)
{
    return 1 - (x * x) - (y * y);
}
double f(double x, double y)
{
    return 4;
}
double M1(double y)
{
    return 1 - y * y;
}
double M2(double y)
{
    return -y * y;
}
double M3(double x)
{
    return 1 - x * x;
}
double M4(double x)
{
    return -8 - x * x;
}
//-------------------------------------------------------------------------���������������-�������------------------------------------------------------------------------------
#define WIDTH_FIELD 6
#define PRECISION_AFTER_POINT 1

double** MemoryAllocator(int n, int m)
{
    double **Matrix = NULL;

    Matrix = new double*[n];
    for (int i = 0; i < n; i++)
        Matrix[i] = new double[m];

    return Matrix;
}
void MemoryCleaner(double **arr, int n)
{
    for (int i = 0; i < n; i++)
        delete[] arr[i];

    delete[] arr;
    arr = NULL;
}
void ShowSolution(double** V, int n, int m)
{
    for (int j = 0; j <= m; j++)
    {
        for (int i = 0; i <= n; i++)
            printf("\nV[%d][%d] = %5.2lf", i, j, V[i][j]);
        printf("\n");
    }
}

void FillRightSide(double **F, int n, int m, double a, double c, double h, double k)
{
    for (int j = 1; j < m; j++)
        for (int i = 1; i < n; i++)
        {
            double Xi, Yj, sum = 0;
            Xi = a + i * h;
            Yj = c + j * k;

            if (j == 1)
                sum += (1 / (k * k)) * M3(Xi);
            else
                if (j == m - 1)
                    sum += (1 / (k * k)) * M4(Xi);
            if (i == 1)
                sum += (1 / (h * h)) * M1(Yj);
            else
                if (i == n - 1)
                    sum += (1 / (h * h)) * M2(Yj);

            F[i][j] = -f(Xi, Yj) - sum;
        }
}
void FillStartSolution(double **V, int n, int m, double a, double b, double c, double d)
{
    double h, k;
    h = (b - a) / n;
    k = (d - c) / m;

    for (int j = 0; j <= m; j++)
        for (int i = 0; i <= n; i++)
        {
            if (i == 0 || j == 0 || i == n || j == m)
            {
                double Xi, Yj, sum = 0;
                Xi = a + i * h;
                Yj = c + j * k;
                if (j == 0)
                    V[i][j] = M3(Xi);
                else
                    if (j == m)
                        V[i][j] = M4(Xi);
                if (i == 0)
                    V[i][j] = M1(Yj);
                else
                    if (i == n)
                        V[i][j] = M2(Yj);
            }
            else
                V[i][j] = 0;
        }
}

void ZeidelsMethod(double **V, int n, int m, double a, double b, double c, double d, double eps, int Nmax, double &epsMax, int &S)
{
    double	epsCur = 0;
    double	a2, k2, h2;
    double	v_old;
    double	v_new;

    h2 = -((n / (b - a)) * (n / (b - a)));
    k2 = -((m / (d - c)) * (m / (d - c)));
    a2 = -2 * (h2 + k2);

    while (true)
    {
        epsMax = 0;
        for (int j = 1; j < m; j++)
            for (int i = 1; i < n; i++)
            {
                double Xi, Yj;
                Xi = a + i * ((b - a) / n);
                Yj = c + j * ((d - c) / m);

                v_old = V[i][j];
                v_new = -(h2*(V[i + 1][j] + V[i - 1][j]) + k2*(V[i][j + 1] + V[i][j - 1]));
                v_new = v_new + f(Xi, Yj);
                v_new = v_new / a2;

                epsCur = abs(v_old - v_new);
                if (epsCur > epsMax)
                    epsMax = epsCur;

                V[i][j] = v_new;
            }
        ++S;

        if ((epsMax < eps) || (S >= Nmax))
            break;
    }
}
double DiscrepancyOfSolution(double **V/*, double **F,*/, int n, int m, double a, double b, double c, double d)
{
    double	a2, k2, h2;
    double  h, k;
    double** F;
    double rs = 0;

    h = (b - a) / n;
    k = (d - c) / m;

    h2 = ((n / (b - a)) * (n / (b - a)));
    k2 = ((m / (d - c)) * (m / (d - c)));
    a2 = -2 * (h2 + k2);

    F = MemoryAllocator(n + 1, m + 1);
    FillRightSide(F, n, m, a, c, h, k);

    for (int j = 1; j < m; j++)
    {
        for (int i = 1; i < n; i++)
        {
            double r;
            double mult;

            if (j != 1 && j != m - 1)
            {
                if (i != 1 && i != n - 1)
                    mult = k2 * V[i][j - 1] + h2 * V[i - 1][j] + a2 * V[i][j] + h2 * V[i + 1][j] + k2 * V[i][j + 1];
                else
                    if (i == 1)
                        mult = k2 * V[i][j - 1] + a2 * V[i][j] + h2 * V[i + 1][j] + k2 * V[i][j + 1];
                    else
                        if (i == n - 1)
                            mult = k2 * V[i][j - 1] + h2 * V[i - 1][j] + a2 * V[i][j] + k2 * V[i][j + 1];
            }
            else
                if (j == 1)
                {
                    if (i == 1)
                        mult = a2 * V[i][j] + h2 * V[i + 1][j] + k2 * V[i][j + 1];
                    else
                        if (i != n - 1)
                            mult = h2 * V[i - 1][j] + a2 * V[i][j] + h2 * V[i + 1][j] + k2 * V[i][j + 1];
                        else
                            if (i == n - 1)
                                mult = h2 * V[i - 1][j] + a2 * V[i][j] + k2 * V[i][j + 1];
                }
                else
                    if (j == m - 1)
                    {
                        if (i == 1)
                            mult = k2 * V[i][j - 1] + a2 * V[i][j] + h2 * V[i + 1][j];
                        else
                            if (i != n - 1)
                                mult = k2 * V[i][j - 1] + h2 * V[i - 1][j] + a2 * V[i][j] + h2 * V[i + 1][j];
                            else
                                if (i == n - 1)
                                    mult = k2 * V[i][j - 1] + h2 * V[i - 1][j] + a2 * V[i][j];
                    }

            r = abs(mult - F[i][j]);

            if (r > rs)
                rs = r;
        }
    }
    MemoryCleaner(F, n);

    return rs;
}

double CheckNumSolution(double** V, int n, int m, double a, double b, double c, double d)
{
    double **U = MemoryAllocator(n + 1, m + 1);
    double h, k;
    double zs = 0;

    h = (b - a) / n;
    k = (d - c) / m;

    for (int j = 0; j <= m; j++)
        for (int i = 0; i <= n; i++)
        {
            double Xi, Yj;
            Xi = a + i * h;
            Yj = c + j * k;

            U[i][j] = Uxy(Xi, Yj);
        }
    for (int j = 1; j < m; j++)
        for (int i = 1; i < n; i++)
        {
            double z = abs(U[i][j] - V[i][j]);

            if (z > zs)
                zs = z;
        }
    MemoryCleaner(U, n);

    return zs;
}

#endif
