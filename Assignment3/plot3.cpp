#define _USE_MATH_DEFINES
#include<iostream>
#include<cmath>
#include<vector>

#include"matplot/matplot.h"

void plotresult1c_(std::vector<double> chargeX, std::vector<double>chargeY, std::vector<double>chargeZ) {

    const int resolution = 50;
    std::vector<std::vector<double>> x(resolution, std::vector<double>(resolution));
    std::vector<std::vector<double>> y(resolution, std::vector<double>(resolution));
    std::vector<std::vector<double>> z(resolution, std::vector<double>(resolution));

    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            double theta = M_PI * i / (resolution - 1);
            double phi = 2 * M_PI * j / (resolution - 1);

            x[i][j] = sin(theta) * cos(phi);
            y[i][j] = sin(theta) * sin(phi);
            z[i][j] = cos(theta);
        }
    }

    auto s = matplot::surf(x, y, z);
    s->face_alpha(0);
    s->edge_color("none");
    s->contour_surface(matplot::off);
    matplot::hold(matplot::on);
    auto p = matplot::plot3(chargeX, chargeY, chargeZ, "r.");
    p->marker_size(14);
    matplot::xlabel("X");
    matplot::ylabel("Y");
    matplot::zlabel("Z");
    matplot::axis(matplot::equal);
    matplot::hold(matplot::off);
    matplot::view(3);
    matplot::show();

}

void plotresult3b_(std::vector<std::vector<double>> ri, std::vector<std::vector<double>> lagrangePoints, double mu) {

    std::vector<double> Rx = { -mu, 1-mu };
    std::vector<double> Ry = { 0,0 };


    std::vector<double> xi;
    std::vector<double> yi;
    std::vector<double> zi;
    for (const auto& pos : ri) {
        xi.push_back(pos[0]);
        yi.push_back(pos[1]);
        zi.push_back(pos[2]);
    }

    std::vector<double> Lxi;
    std::vector<double> Lyi;
    for (const auto& pos : lagrangePoints) {
        Lxi.push_back(pos[0]);
        Lyi.push_back(pos[1]);
    }

    matplot::plot(Lxi, Lyi, "g.");
    matplot::hold(matplot::on);
    matplot::plot(xi, yi);
    matplot::plot(Rx, Ry, "r.");

    std::vector<double> startpointX = { xi[0] };
    std::vector<double> startpointY = { yi[0] };

    matplot::plot(startpointX, startpointY, "b.");
    matplot::text(xi[0], yi[0], "Startpoint");

    matplot::text(Rx[0], Ry[0], "M1");
    matplot::text(Rx[1], Ry[1], "M2");

    matplot::text(lagrangePoints[0][0], lagrangePoints[0][1], "L1");
    matplot::text(lagrangePoints[1][0], lagrangePoints[1][1], "L2");
    matplot::text(lagrangePoints[2][0], lagrangePoints[2][1], "L3");
    matplot::text(lagrangePoints[3][0], lagrangePoints[3][1], "L4");
    matplot::text(lagrangePoints[4][0], lagrangePoints[4][1], "L5");

    matplot::xlabel("X / 1");
    matplot::ylabel("Y / 1");
    matplot::hold(matplot::off);
    matplot::grid(matplot::on);
    matplot::show();
}