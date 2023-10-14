//Testing equations
#include "boyer.h"
#include "rk4.h"
#include "rk4_adaptive.h"
#include "rk45_dormand_prince.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>


int main(){
    double a = 1.0;
    double M = 1.0;
    double rH = M + sqrt(pow(M,2) - pow(a, 2));

    boyer_lindquist_metric metric(a, M);
    
    auto ut = [&metric](const std::vector<double> &y){
        metric.compute_metric(y[0], y[1]);
        return (sqrt((metric.gamma11*pow(y[3],2) + metric.gamma22*pow(y[4],2) + metric.gamma33*pow(y[5],2))))/metric.alpha;

    };
    auto utsub = [&metric,ut](const std::vector<double> &y){
        metric.compute_metric(y[0], y[1]);
        return -pow(metric.alpha, 2)*ut(y) + y[5]*metric.beta3;

    };

    auto H = [&metric](const std::vector<double> &y){
        metric.compute_metric(y[0], y[1]);
        return (sqrt((metric.gamma11*y[3]*y[3] + metric.gamma22*y[4]*y[4] + metric.gamma33*y[5]*y[5])))*metric.alpha-y[5]*metric.beta3;

    };
    


    auto dydx = [&metric,ut,utsub,H] (double x, const std::vector<double> &y) {
        metric.compute_metric(y[0], y[1]); // y[0] is r, y[1] is theta // compute the right hand sides of equations (17)-(22)
        std::vector<double> result1(6);
        // dr/dt = gamma11(y[3]/ut)
        // dtheta/dt = gamma22(y[4]/ut)
        // dphi/dt = gamma33(y[5]/ut) - beta3
        // dur/dt =alpha*ut*d_alpha_dr ...
        // dutheta/dt  = ...
        // duphi/dt = 0
        result1[0] = metric.gamma11*(y[3]/ut(y));
        result1[1] = metric.gamma22*(y[4]/ut(y));
        result1[2] = metric.gamma33*(y[5]/ut(y)) -metric.beta3;

        result1[3] = -metric.alpha*ut(y)*metric.d_alpha_dr + y[5]*metric.d_beta3_dr - (1.0/(2.0*ut(y)))*(pow(y[3],2)*metric.d_gamma11_dr + pow(y[4], 2)*metric.d_gamma22_dr
        + pow(y[5], 2)*metric.d_gamma33_dr);

        result1[4] = -metric.alpha*ut(y)*metric.d_alpha_dth + y[5]*metric.d_beta3_dth - (1.0/(2.0*ut(y)))*(pow(y[3],2)*metric.d_gamma11_dth + pow(y[4], 2)*metric.d_gamma22_dth
        + pow(y[5], 2)*metric.d_gamma33_dth);

        result1[5] = 0.0; 

        
        std::cout << "r: " << y[0] << std::endl;
        std::cout << "phi: " <<y[2] << std::endl;
        
        return result1;

    };

    // Define our stop condition
   

    // //A stop condition
    // auto stopA = [](double x, const std::vector<double> &y) { return y[0] > 3.0; };

    // //B stop condition
    // auto stopB = [](double x, const std::vector<double> &y) { return y[0] > 1.0; };

    // //C stop condition
    // auto stopC = [](double x, const std::vector<double> &y) { return y[0] > 1.0; };

    // //D stop condition
    // auto stopD = [](double x, const std::vector<double> &y) { return y[0] > 1.0; };

    // //E stop condition
    // auto stopE = [](double x, const std::vector<double> &y) { return y[0] > 1.0; };
    
    
    

    //Initial conditions
    double theta0 = M_PI/2.0;
    double phi0 = 0.0;
    double ur0 = 0.0;
    double h = 1e-3;
    double x0 = 0.0;

    //Case A
    double r0A = 1.0 + sqrt(2.0);
    double utheta0A = sqrt(11.0 + 8.0*sqrt(2.0));
    double uphi0A = 0.0;

    //Case B
    double r0B = 1.0 + sqrt(3.0);
    double utheta0B = sqrt(12 + 8*sqrt(3.0));
    double uphi0B = -1.0;

    //Case C
    double r0C = 3.0;
    double utheta0C = sqrt(27.0);
    double uphi0C = -2.0;

    //Case D
    double r0D = 1.0 + 2.0*sqrt(2.0);
    double utheta0D = sqrt(-13.0 + 16.0*sqrt(2.0));
    double uphi0D = -6.0;

    //Case E
    double r0E = 2.0;
    double utheta0E = sqrt(16.0);
    double uphi0E = 1.0;

      //A stop condition
    auto stopA = [r0A,rH,h](double x, const std::vector<double> &y) { return ((y[0] < r0A-h )|| (y[0]> r0A)) ; };

    //B stop condition
    auto stopB = [r0B,rH,h](double x, const std::vector<double> &y) { return ((y[0] < rH + h) || (y[0]> r0B)); };

    //C stop condition
    auto stopC = [r0C,rH,h](double x, const std::vector<double> &y) { return ((y[0] < rH + h) || (y[0]> r0C)); };

    //D stop condition
    auto stopD = [r0D,rH,h](double x, const std::vector<double> &y) { return ((y[0] < rH + h) || (y[0]> r0D)); };

    //E stop condition
    auto stopE = [r0E,rH,h](double x, const std::vector<double> &y) { return ((y[0] < rH + h) || (y[0]> r0E)); };



    //No, I think this intializes the class?
    //runge_kutta_4 rk4(6);
    //rk4_adaptive rk4adapt(6, 1e-13, 1e-13);
    rk45_dormand_prince rkdpA(6, 1e-17, 1e-17);
    rk45_dormand_prince rkdpB(6, 1e-10, 1e-10);
    rk45_dormand_prince rkdpC(6, 1e-10, 1e-10);
    rk45_dormand_prince rkdpD(6, 1e-10, 1e-10);
    rk45_dormand_prince rkdpE(6, 1e-10, 1e-10);
    

    //Rest of initial conditions
    
    //Wait.. doesn't that mean y0 has to have six variables...
    // YES! Because it's dr/dt, dtheta/dt, dphi/dt, dur/dt, dutheta/dt, duphi/dt
    std::vector<double> y0A = {r0A, theta0,phi0 ,ur0 ,utheta0A ,uphi0A};
    std::vector<double> y0B = {r0B, theta0,phi0 ,ur0 ,utheta0B ,uphi0B};
    std::vector<double> y0C = {r0C, theta0,phi0 ,ur0 ,utheta0C ,uphi0C};
    std::vector<double> y0D = {r0D, theta0,phi0 ,ur0 ,utheta0D ,uphi0D};
    std::vector<double> y0E = {r0E, theta0,phi0 ,ur0 ,utheta0E ,uphi0E};

    //intialize y, we never did
    //auto y =rk4.integrate(dydx, stop, h, x0, y0);
    //auto y =rk4adapt.integrate(dydx, stop, h, x0, y0);
    auto y =rkdpA.integrate(dydx, stopA, h, x0, y0A);
    y =rkdpB.integrate(dydx, stopB, h, x0, y0B);
    y =rkdpC.integrate(dydx, stopC, h, x0, y0C);
    y =rkdpD.integrate(dydx, stopD, h, x0, y0D);
    y =rkdpE.integrate(dydx, stopE, h, x0, y0E);

    // Output the results

    // std::ofstream output_file("output_einstein_ring_RK4.csv");
    // for (int i = 0; i < rk4.xs.size(); i++) {
    //     output_file << rk4.xs[i] << "," << rk4.result[i][0] << "," << rk4.result[i][1]<< "," << rk4.result[i][2]
    //      << "," << rk4.result[i][3] << "," << rk4.result[i][4]<< "," << rk4.result[i][5]<< std::endl;
    // }
    // output_file.close();

    // std::ofstream output_file_adapt("output_einstein_ring.csv");
    // for (int i = 0; i < rk4adapt.xs.size(); i++) {
    //     output_file_adapt << rk4adapt.xs[i] << "," << rk4adapt.result[i][0] << "," << rk4adapt.result[i][1]<< "," << rk4adapt.result[i][2]
    //     << "," << rk4adapt.result[i][3] << "," << rk4adapt.result[i][4]<< "," << rk4adapt.result[i][5]<< std::endl;
    // }
    // output_file_adapt.close();

    //First example
    std::ofstream output_file_A("output_case_A.csv");
    for (int i = 0; i < rkdpA.xs.size(); i++) {
        output_file_A << rkdpA.xs[i] << "," << rkdpA.result[i][0] << "," << rkdpA.result[i][1]<< "," << rkdpA.result[i][2]
         << "," << rkdpA.result[i][3] << "," << rkdpA.result[i][4]<< "," << rkdpA.result[i][5]<< std::endl;
    }
    output_file_A.close();

    std::ofstream output_file_B("output_case_B.csv");
    for (int i = 0; i < rkdpB.xs.size(); i++) {
        output_file_B << rkdpB.xs[i] << "," << rkdpB.result[i][0] << "," << rkdpB.result[i][1]<< "," << rkdpB.result[i][2]
         << "," << rkdpB.result[i][3] << "," << rkdpB.result[i][4]<< "," << rkdpB.result[i][5]<< std::endl;
    }
    output_file_B.close();

    std::ofstream output_file_C("output_case_C.csv");
    for (int i = 0; i < rkdpC.xs.size(); i++) {
        output_file_C << rkdpC.xs[i] << "," << rkdpC.result[i][0] << "," << rkdpC.result[i][1]<< "," << rkdpC.result[i][2]
         << "," << rkdpC.result[i][3] << "," << rkdpC.result[i][4]<< "," << rkdpC.result[i][5]<< std::endl;
    }
    output_file_C.close();

    std::ofstream output_file_D("output_case_D.csv");
    for (int i = 0; i < rkdpD.xs.size(); i++) {
        output_file_D << rkdpD.xs[i] << "," << rkdpD.result[i][0] << "," << rkdpD.result[i][1]<< "," << rkdpD.result[i][2]
         << "," << rkdpD.result[i][3] << "," << rkdpD.result[i][4]<< "," << rkdpD.result[i][5]<< std::endl;
    }
    output_file_D.close();

    std::ofstream output_file_E("output_case_E.csv");
    for (int i = 0; i < rkdpE.xs.size(); i++) {
        output_file_E << rkdpE.xs[i] << "," << rkdpE.result[i][0] << "," << rkdpE.result[i][1]<< "," << rkdpE.result[i][2]
         << "," << rkdpE.result[i][3] << "," << rkdpE.result[i][4]<< "," << rkdpE.result[i][5]<< std::endl;
    }
    output_file_E.close();

    return 0;
};