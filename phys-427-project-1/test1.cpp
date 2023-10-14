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
    double a = 0.0;
    double M = 1.0;

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

        
        std::cout << "utsub: " << utsub(y) << std::endl;
        std::cout << "H: " <<H(y) << std::endl;
        
        return result1;

    };

    // Define our stop condition
    //TO-DO: Figure out stop condition (DID IT)
    //This is the actual stop condtion, uncomment when you fix code

    //First stop condition
    //auto stop = [M](double x, const std::vector<double> &y) { return y[0] > 10.0*M; };

    //Second stop condition
    auto stop = [M](double x, const std::vector<double> &y) { return y[0] > 20.0*M; };
    
    
    

    //Initial conditions
    //TO-DO:Figure out initial conditions for r0 and theta0 (DID IT)

    //First case of r0
    //double r0 = 10.0*M;

    //Second case of r0
    double r0 = 20.0*M;

    double theta0 = M_PI/2.0;
    double phi0 = 0.0;

    //First case for xi
    //double xi = 0.48857874 ;

    //Second case for xi
    double xi = 0.24904964;


    //No, I think this intializes the class?
    //runge_kutta_4 rk4(6);
    //rk4_adaptive rk4adapt(6, 1e-13, 1e-13);
    rk45_dormand_prince rkdp(6, 1e-18, 1e-18);
    

    //Need this to compute other initial conditions
    boyer_lindquist_metric metriciv(a, M);
    metriciv.compute_metric(r0, theta0);

    //Rest of initial conditions
    double ur0 = -cos(xi)*sqrt(metriciv.g_11);
    double utheta0 = 0.0;
    double uphi0 = sin(xi)*sqrt(metriciv.g_33);

    double h = 1e-9;
    double x0 = 0.0;
    //Wait.. doesn't that mean y0 has to have six variables...
    // YES! Because it's dr/dt, dtheta/dt, dphi/dt, dur/dt, dutheta/dt, duphi/dt
    std::vector<double> y0 = {r0, theta0,phi0 ,ur0 ,utheta0 ,uphi0};

    //intialize y, we never did
    //auto y =rk4.integrate(dydx, stop, h, x0, y0);
    //auto y =rk4adapt.integrate(dydx, stop, h, x0, y0);
    auto y =rkdp.integrate(dydx, stop, h, x0, y0);

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
    // std::ofstream output_file("output_einstein_ring_RKDP_2_pi.csv");
    // for (int i = 0; i < rkdp.xs.size(); i++) {
    //     output_file << rkdp.xs[i] << "," << rkdp.result[i][0] << "," << rkdp.result[i][1]<< "," << rkdp.result[i][2]
    //      << "," << rkdp.result[i][3] << "," << rkdp.result[i][4]<< "," << rkdp.result[i][5]<< std::endl;
    // }
    // output_file.close();

    //Second example
    std::ofstream output_file("output_einstein_ring_RKDP_4_pi.csv");
    for (int i = 0; i < rkdp.xs.size(); i++) {
        output_file << rkdp.xs[i] << "," << rkdp.result[i][0] << "," << rkdp.result[i][1]<< "," << rkdp.result[i][2]
         << "," << rkdp.result[i][3] << "," << rkdp.result[i][4]<< "," << rkdp.result[i][5]<< std::endl;
    }
    output_file.close();

    return 0;
};