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
    s->face_alpha(0.1);
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