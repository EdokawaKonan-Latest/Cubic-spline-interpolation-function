#include "eigen.h"
#include <QString>
eigen::eigen()
{

}

void eigen::getans()
{
    const double min = 0;
    const double max = 1.0;
    int size = sizeof(x) / sizeof(double);
    xx.clear();
    yy.clear();
    double step = 0.01;
    double value = x[0];
    int num = (max - min) / step;
    for (double i = 0; i <= num; i++) {
        xx.push_back(value);
        value = value + step;
    }
    int size_xx = xx.size();

    vector<double> dx;
    vector<double> dy;
    for (int i = 0; i < size - 1; i++) {
        double temp_dx = x[i + 1] - x[i];
        dx.push_back(temp_dx);
        double temp_dy = y[i + 1] - y[i];
        dy.push_back(temp_dy);
    }

    MatrixXd H = MatrixXd::Random(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            H(i, j) = 0;
        }
    }
    VectorXd Y(size);
    for (int i = 0; i < size; i++) {
        Y(i) = 0;
    }
    VectorXd M(size);
    for (int i = 0; i < size; i++) {
        M(i) = 0;
    }

    H(0, 0) = 1;
    H(size - 1, size - 1) = 1;
    for (int i = 1; i < size - 1; i++) {
        H(i, i - 1) = dx[i - 1];
        H(i, i) = 2 * (dx[i - 1] + dx[i]);
        H(i, i + 1) = dx[i];
        Y(i) = 3 * (dy[i] / dx[i] - dy[i - 1] / dx[i - 1]);
    }

    M = H.colPivHouseholderQr().solve(Y);

    vector<double> ai, bi, ci, di;
    for (int i = 0; i < size - 1; i++) {
        ai.push_back(y[i]);
        di.push_back((M(i + 1) - M(i)) / (3 * dx[i]));
        bi.push_back(dy[i] / dx[i] - dx[i] * (2 * M(i) + M(i + 1)) / 3);
        ci.push_back(M(i));
    }

    vector<int> x_, xx_;
    for (int i = 0; i < size; i++) {
        int temp = x[i] / 0.01;
        x_.push_back(temp);
    }
    for (int i = 0; i < size_xx; i++) {
        int temp = xx[i] / 0.01;
        xx_.push_back(temp);
    }

    for (int i = 0; i < size_xx; i++) {
        int k = 0;
        for (int j = 0; j < size - 1; j++) {
            if (xx_[i] >= x_[j] && xx_[i] < x_[j + 1]) {
                k = j;
                break;
            }
            else if (xx[i] == x[size - 1]) {
                k = size - 1;
            }
        }
        //yy(i) = y[i] + bi(k) * (xx[i] - x[k]) + 1 / 2.0 * M(i) * pow((xx[i] - x[k]) , 2) + di(k) * pow((xx[i] - x[k]),3);
        double temp = ai[k] + bi[k] * (xx[i] - x[k]) + M(k) * pow((xx[i] - x[k]), 2) + di[k] * pow((xx[i] - x[k]), 3);
        yy.push_back(temp);
    }

    std::ofstream output;
    output.open("Spline.txt");
    for (int i = 0; i < size_xx; i++) {
        output << xx[i] << '\t' << yy[i] << std::endl;
        //cout << xx[i] << '\t' << yy[i] << std::endl;
    }
    emit datechange();
    output.close();
}

void eigen::changey(int x_id, double y_value)
{
    y[x_id] = y_value;
    getans();
}

double eigen::getyVal(int id)
{
    return yy[id];
}

double eigen::getxVal(int id)
{
    return (xx[id]);
}
