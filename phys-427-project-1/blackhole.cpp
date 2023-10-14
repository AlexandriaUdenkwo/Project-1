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
        metric.compute_metric(y[0], y[1]); 
        std::vector<double> result1(6);
       
        result1[0] = metric.gamma11*(y[3]/ut(y));
        result1[1] = metric.gamma22*(y[4]/ut(y));
        result1[2] = metric.gamma33*(y[5]/ut(y)) -metric.beta3;

        result1[3] = -metric.alpha*ut(y)*metric.d_alpha_dr + y[5]*metric.d_beta3_dr - (1.0/(2.0*ut(y)))*(pow(y[3],2)*metric.d_gamma11_dr + pow(y[4], 2)*metric.d_gamma22_dr
        + pow(y[5], 2)*metric.d_gamma33_dr);

        result1[4] = -metric.alpha*ut(y)*metric.d_alpha_dth + y[5]*metric.d_beta3_dth - (1.0/(2.0*ut(y)))*(pow(y[3],2)*metric.d_gamma11_dth + pow(y[4], 2)*metric.d_gamma22_dth
        + pow(y[5], 2)*metric.d_gamma33_dth);

        result1[5] = 0.0; 

        
        // std::cout << "r: " << y[0] << std::endl;
        // std::cout << "phi: " <<y[2] << std::endl;
        
        return result1;

    };

    
    

    //Initial conditions
    double theta0 = M_PI/2.0;
    double phi0 = 0.0;
    double ur0 = 0.0;
    double h = 1e-3;
    double tol = 1e-5;
    double x0 = 0.0;

    //Case (A)ccrection
    double r0A = 8.0*M;
    double utheta0A = sqrt(11.0 + 8.0*sqrt(2.0));
    double uphi0A = 0.0;

    //Case (B)lack Hole
    

    //Case (I)nfinity

   

    //A stop condition
    auto stopA = [r0A,tol](double x, const std::vector<double> &y) { return ((y[0] < r0A+tol )&& (y[0]> r0A-tol)) ; };

    //B stop conditions
    auto stopB = [r0A,tol](double x, const std::vector<double> &y) { return (y[0]< r0A-tol) ; };

    //I stop conditions
    auto stopI = [r0A,tol](double x, const std::vector<double> &y) { return (y[0]> r0A+tol) ; };

   



    //No, I think this intializes the class?
    //runge_kutta_4 rk4(6);
    //rk4_adaptive rk4adapt(6, 1e-13, 1e-13);
    rk45_dormand_prince rkdpA(6, 1e-17, 1e-17);
    
    

    //Rest of initial conditions
    
    //Wait.. doesn't that mean y0 has to have six variables...
    // YES! Because it's dr/dt, dtheta/dt, dphi/dt, dur/dt, dutheta/dt, duphi/dt
    std::vector<double> y0A = {r0A, theta0,phi0 ,ur0 ,utheta0A ,uphi0A};

    //intialize y, we never did
    //auto y =rk4.integrate(dydx, stop, h, x0, y0);
    //auto y =rk4adapt.integrate(dydx, stop, h, x0, y0);
    auto y =rkdpA.integrate(dydx, stopA, h, x0, y0A);

    
    //Make grid
    //alpha = beta = 0.01 rad = 1000 (radx10^5)
    double length = 1000.0; //radx10^5
    int points = 100;
    double spacing = length/points;

    //Centers
    double x1 = length/2.0;
    double y1 = 0;

    //Intialize grid
    double grid[points][points];

    //Calculate points
    for (int i =0; i < points; i++){


        
    }


    double num = rkdpA.result.size();
    std::cout << "num: " << num << std::endl;

    double last_element_r = rkdpA.result[num-1][1];
    std::cout << "Last element of r: " << last_element_r << std::endl;

    //double last_element = rkdpA.result[num-1][0];

  

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
    // std::ofstream output_file_A("output_case_A.csv");
    // for (int i = 0; i < rkdpA.xs.size(); i++) {
    //     output_file_A << rkdpA.xs[i] << "," << rkdpA.result[i][0] << "," << rkdpA.result[i][1]<< "," << rkdpA.result[i][2]
    //      << "," << rkdpA.result[i][3] << "," << rkdpA.result[i][4]<< "," << rkdpA.result[i][5]<< std::endl;
    // }
    // output_file_A.close();


    return 0;
};