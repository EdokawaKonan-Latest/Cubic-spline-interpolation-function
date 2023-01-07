#ifndef EIGEN_H
#define EIGEN_H

#include <QObject>
#include<vector>
#include <Eigen/SVD>
#include <Eigen/Core>
#include<iostream>
#include<cmath>
#include <fstream>
#include <qDebug>
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;
class eigen : public QObject
{
    Q_OBJECT
public:
    eigen();
    Q_INVOKABLE void getans();
    Q_INVOKABLE void changey(int x_id, double y_value);
    Q_INVOKABLE double getyVal(int id);
    Q_INVOKABLE double getxVal(int id);
    double x[6] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    double y[6] = {0.0, 0.8, 0.2, 0.3, 0.4, 0.5};
    vector<double> xx, yy;
signals:
    void datechange();
};

#endif // EIGEN_H
