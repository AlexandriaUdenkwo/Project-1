//#pragma once
//#include "rk4_adaptive.h"
#include <vector>
#include <cmath>
#include <iostream>




class boyer_lindquist_metric {
 public:
    boyer_lindquist_metric(double a0, double M0) { // Initialize the parameters a and M
        a = a0;
        M = M0;

    }
    void compute_metric(double r, double th) {
        //Mine
        // g_11 = (r*r + a*a*std::pow(std::cos(th),2.0))/(r*r + a*a - 2.0*M*r);
        // g_22 = r*r + a*a*std::pow(std::cos(th),2.0);
        // g_33 = (((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0))*std::pow(std::sin(th),2.0))/(r*r + a*a*std::pow(std::cos(th),2.0));
        // g_03 = (-2.0*M*a*r* std::pow(std::sin(th), 2.0))/(r*r + a*a*std::pow(std::cos(th),2.0));

        // gamma11 = 1.0/g_11;
        // gamma22 = 1.0/g_22;
        // gamma33 = 1.0/g_33;

        // alpha = std::sqrt(((r*r + a*a*std::pow(std::cos(th),2.0))*(r*r + a*a - 2.0*M*r))/((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0)));
        // beta3 = (-2.0*M*a*r)/((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0));

        // d_alpha_dr = ( (r* std::pow(r*r + a*a*std::pow(std::cos(th),2.0), -0.5)*std::pow(r*r + a*a - 2.0*M*r,0.5)*
        // std::pow((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0), -0.5)) 

        // + ((r - M)*std::pow(r*r + a*a - 2.0*M*r,-0.5)*std::pow(r*r + a*a*std::pow(std::cos(th),2.0), 0.5)*
        // std::pow((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0), -0.5))           

        // - ((2.0*r*(r*r + a*a) - a*a*(r - M)*std::pow(std::sin(th),2.0))*
        // std::pow((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0), -1.5)*std::pow(r*r + a*a*std::pow(std::cos(th),2.0), 0.5)*
        // std::pow(r*r + a*a - 2.0*M*r,0.5))

        // );

        // d_beta3_dr = -2.0*M*a*( std::pow((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0), -1.0)

        // - (r*(4.0*r*(r*r + a*a) - 2.0*a*a*(r - M)*std::pow(std::sin(th),2.0))*
        // std::pow((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0), -2.0))   

        // );

        // d_gamma11_dr = 2.0*(r-M)* std::pow((r*r + a*a*std::pow(std::cos(th),2.0)),-1.0) - (2.0*r*std::pow((r*r + a*a*std::pow(std::cos(th),2.0)),-2.0)*
        // (r*r + a*a - 2.0*M*r));

        // d_gamma22_dr = -2.0*r*(std::pow((r*r + a*a*std::pow(std::cos(th),2.0)),-2.0));

        // d_gamma33_dr =  2.0*r*(std::pow(((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0))*std::pow(sin(th), 2.0), -1.0))
        
        // -  ( 4.0*r*(r*r + a*a)*std::pow(sin(th),2.0) - 2.0*a*a*(r-M)*std::pow(sin(th),4))*
        // (std::pow(((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0))*std::pow(sin(th), 2.0), -2.0))*(r*r + a*a*std::pow(std::cos(th),2.0));

        // d_alpha_dth =( (-a*a*std::sin(th)*std::cos(th)*std::pow((r*r +a*a*std::pow(cos(th), 2.0)),-0.5))*      
        // std::pow(r*r + a*a - 2.0*M*r,0.5)*
        // std::pow((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0), -0.5))
        
        // + a*a*(r*r + a*a - 2.0*M*r)*std::sin(th)*std::cos(th)*( std::pow(r*r + a*a - 2.0*M*r,0.5))*
        // (std::pow((r*r +a*a*std::pow(cos(th), 2.0)),0.5))*(std::pow((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0), -1.5))
        // ;

        // d_beta3_dth =  -4.0*M*std::pow(a,3)*r*(r*r + a*a - 2.0*M*r)*std::sin(th)*std::cos(th)*
        // (std::pow((r*r + a*a)*(r*r + a*a) - a*a*(r*r + a*a - 2.0*M*r)*std::pow(std::sin(th),2.0), -2.0))
        // ;

        // d_gamma11_dth = -(2.0*a*a*std::cos(th)*std::sin(th))*std::pow((r*r + a*a*std::pow(std::cos(th), 2)), -2.0)*(r*r + a*a - 2.0*M*r);

        // d_gamma22_dth = -2.0*a*a*std::cos(th)* std::sin(th)* std::pow((r*r + a*a*std::pow(std::cos(th),2)) ,-2);

        // d_gamma33_dth = 2.0*a*a*std::sin(th)*std::cos(th)*std::pow((std::pow((r*r+a*a),2) - a*a*(r*r + a*a -2.0*M*r)*std::pow(std::sin(th), 2)), -1)*
        // std::pow(std::sin(th), -2.0)
        
        // + (2.0*a*a*(r*r + a*a - 2.0*M*r)*std::sin(th)*std::cos(th))*std::pow((std::pow((r*r+a*a),2) - a*a*(r*r + a*a -2.0*M*r)*std::pow(std::sin(th), 2)), -2)*
        // std::pow(std::sin(th), -2.0)*(r*r + a*a*std::pow(std::cos(th), 2))

        // + std::pow((std::pow((r*r+a*a),2) - a*a*(r*r + a*a -2.0*M*r)*std::pow(std::sin(th), 2)), -1)*(r*r + a*a*std::pow(std::cos(th), 2))* 
        // (-2.0*std::cos(th)*std::pow(std::sin(th), -3));

        //Mathematica's
        g_11 = (pow(r,2) + pow(a,2)*pow(cos(th),2))/(pow(a,2) - 2.0*M*r + pow(r,2)); 
        g_22 = pow(r,2) + pow(a,2)*pow(cos(th),2);
        g_33 = (pow(sin(th),2)*(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2)))/
                (pow(r,2) + pow(a,2)*pow(cos(th),2));
        
        gamma11 = 1.0/g_11;
        gamma22 = 1.0/g_22;
        gamma33 = 1.0/g_33;

        alpha = sqrt(((pow(a,2) - 2.0*M*r + pow(r,2))*(pow(r,2) + pow(a,2)*pow(cos(th),2)))/
                (pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2)));

        beta3 = (-2.0*a*M*r)/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2)); 

        d_alpha_dr = (-(((pow(a,2) - 2.0*M*r + pow(r,2))*(pow(r,2) + pow(a,2)*pow(cos(th),2))*(4.0*r*(pow(a,2) + pow(r,2)) - pow(a,2)*(-2.0*M + 2.0*r)*pow(sin(th),2)))/pow(pow(pow(a,2) + pow(r,2),2) - 
                        pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2),2)) + (2.0*r*(pow(a,2) - 2.0*M*r + pow(r,2)))/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2)) + ((-2.0*M + 2.0*r)*(pow(r,2) + pow(a,2)*pow(cos(th),2)))/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2)))/(2.0*sqrt(((pow(a,2) - 2.0*M*r + pow(r,2))*(pow(r,2) + pow(a,2)*pow(cos(th),2)))/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2))));


        d_beta3_dr = (2.0*a*M*r*(4.0*r*(pow(a,2) + pow(r,2)) - pow(a,2)*(-2.0*M + 2.0*r)*pow(sin(th),2)))/pow(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2),2) - (2.0*a*M)/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2));

        d_gamma11_dr = (-2.0*r*(pow(a,2) - 2.0*M*r + pow(r,2)))/pow(pow(r,2) + pow(a,2)*pow(cos(th),2),2) + (-2.0*M + 2.0*r)/(pow(r,2) + pow(a,2)*pow(cos(th),2));

        d_gamma22_dr = (-2.0*r)/pow(pow(r,2) + pow(a,2)*pow(cos(th),2),2);

        d_gamma33_dr = -(((pow(r,2) + pow(a,2)*pow(cos(th),2))*pow((1.0/sin(th)),2)*(4.0*r*(pow(a,2) + pow(r,2)) - pow(a,2)*(-2.0*M + 2.0*r)*pow(sin(th),2)))/pow(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2),2)) + (2.0*r*pow((1.0/sin(th)),2))/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2));
        
        d_alpha_dth = ((2.0*pow(a,2)*pow(pow(a,2) - 2.0*M*r + pow(r,2),2)*cos(th)*(pow(r,2) + pow(a,2)*pow(cos(th),2))*sin(th))/pow(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2),2) - (2.0*pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*cos(th)*sin(th))/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2)))/(2.0*sqrt(((pow(a,2) - 2.0*M*r + pow(r,2))*(pow(r,2) + pow(a,2)*pow(cos(th),2)))/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2))));
        
        d_beta3_dth = (-4.0*pow(a,3)*M*r*(pow(a,2) - 2.0*M*r + pow(r,2))*cos(th)*sin(th))/pow(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2),2);

        d_gamma11_dth = (2.0*pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*cos(th)*sin(th))/pow(pow(r,2) + pow(a,2)*pow(cos(th),2),2);

        d_gamma22_dth = (2.0*pow(a,2)*cos(th)*sin(th))/pow(pow(r,2) + pow(a,2)*pow(cos(th),2),2); 

        d_gamma33_dth = (2.0*pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*(pow(r,2) + pow(a,2)*pow(cos(th),2))*(1.0/tan(th)))/pow(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2),2) - (2.0*pow(a,2)*(1.0/tan(th)))/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2)) - (2.0*(pow(r,2) + pow(a,2)*pow(cos(th),2))*(1.0/tan(th))*pow((1.0/sin(th)),2))/(pow(pow(a,2) + pow(r,2),2) - pow(a,2)*(pow(a,2) - 2.0*M*r + pow(r,2))*pow(sin(th),2));
        

    // Compute all the required metric components in one go
    }


    double a, M;
    double alpha, beta3;
    double gamma11, gamma22, gamma33; // components of upper gamma^ij
    double g_00, g_03, g_11, g_22, g_33; // components of lower g_\mu\nu
    double d_alpha_dr, d_beta3_dr, d_gamma11_dr, d_gamma22_dr, d_gamma33_dr; 
    double d_alpha_dth, d_beta3_dth, d_gamma11_dth, d_gamma22_dth, d_gamma33_dth;
};
