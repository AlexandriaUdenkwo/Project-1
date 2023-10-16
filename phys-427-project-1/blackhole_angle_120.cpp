#include "boyer.h"
#include "rk45_dormand_prince.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>


int main(){
    double a = 0.99; //How fast blackhole is spinning
    double M = 1.0; //Masss of blackhole
    double rH = M + sqrt(pow(M,2) - pow(a, 2)); //event horizon, radius of blackhole
    double Dobs = 500; //Distance away observer is from blackhole

    double r0Accin = 5.0*(M);
    double r0Accout = 20.0*M; //Radius of light ring 
    double thetaiv = 120.0*(M_PI/180.0);

    //Precision tolerance 
    double tol = 0.5e-2;

    double x0 = 0.0; //time
    double h = 1e-3; //step size for differential equations

    boyer_lindquist_metric metric(a, M);
    boyer_lindquist_metric metriciv(a, M);
    boyer_lindquist_metric metricdop(a, M);
    
    auto ut = [&metric](const std::vector<double> &y){
        metric.compute_metric(y[0], y[1]);
        return (sqrt((metric.gamma11*pow(y[3],2) + metric.gamma22*pow(y[4],2) + metric.gamma33*pow(y[5],2))))/metric.alpha;

    };
    auto lastut = [&metricdop](const std::vector<double> &y){
        metricdop.compute_metric(y[0], y[1]);
        return (sqrt((metricdop.gamma11*pow(y[3],2) + metricdop.gamma22*pow(y[4],2) + metricdop.gamma33*pow(y[5],2))))/metricdop.alpha;

    };
    auto utsub = [&metric,ut](const std::vector<double> &y){
        metric.compute_metric(y[0], y[1]);
        return -pow(metric.alpha, 2)*ut(y) + y[5]*metric.beta3;

    };

    auto H = [&metric](const std::vector<double> &y){
        metric.compute_metric(y[0], y[1]);
        return (sqrt((metric.gamma11*y[3]*y[3] + metric.gamma22*y[4]*y[4] + metric.gamma33*y[5]*y[5])))*metric.alpha-y[5]*metric.beta3;

    };
    


    auto dydx = [&metric,ut,utsub,H,Dobs] (double x, const std::vector<double> &y) {
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

        
        
        return result1;

    };

    //Initialize class
    rk45_dormand_prince rkdpA(6, 1e-15, 1e-15);
    
    

    
    //Make grid
    double xrange = 44.0; 
    double yrange = 22.0; 
    double lengthx = xrange/Dobs; //alpha rads
    double lengthy = yrange/Dobs;
    int pointsx = 512;//64;//960;//1920; //64*1;
    int pointsy = 288;//36;//540;//1080; //36*1;
    double spacingx = lengthx/pointsx;
    double spacingy = lengthy/pointsy;


    //Intialize grid
    std::vector<std::vector<double>> grid(pointsy, std::vector<double>(pointsx));

    //Calculate points
    for (int i =-pointsy/2; i < pointsy/2; i++){
        double ygrid = Dobs*spacingy* i;
        for (int j =-pointsx/2; j < pointsx/2; j++){
            double xgrid = Dobs*spacingx* j;
            double r0 = sqrt(pow(Dobs,2) + pow(xgrid, 2) + pow(ygrid, 2));
           
        

            auto stop = [rH,tol, Dobs, r0Accout, r0Accin](double x, const std::vector<double> &y) { return ((((y[1] < (M_PI/2.0)+tol )&&(y[1] > (M_PI/2.0)-tol ))&&( (y[0]> +r0Accin-tol)&&(y[0] < r0Accout+tol )))||
            (y[0]< 1.01*rH)||(y[0]>Dobs*2.5));  };

            
            //Initial conditions
            double alpha = ygrid/Dobs;
            double beta = xgrid/Dobs;


            double theta0 = thetaiv - alpha;
            double phi0 = beta;

    
            
            metriciv.compute_metric(r0, theta0);
         

            double ur0 = -sqrt(metriciv.g_11)*cos(beta)*cos(alpha); 
            double utheta0= sqrt(metriciv.g_22)*sin(alpha);
            double uphi0 = sqrt(metriciv.g_33)*sin(beta)*cos(alpha);
            std::vector<double> y0 = {r0, theta0,phi0 ,ur0 ,utheta0 ,uphi0};
           
                

            auto y =rkdpA.integrate(dydx, stop, h, x0, y0);
            

            double num = rkdpA.result.size();
            
            
            double last_element_r = rkdpA.result[num-1][0];
            double last_element_theta = rkdpA.result[num-1][1];
            double last_element_phi = rkdpA.result[num-1][2];
            double last_element_ur = rkdpA.result[num-1][3];
            double last_element_utheta = rkdpA.result[num-1][4];
            double last_element_uphi = rkdpA.result[num-1][5];
            std::vector<double> lasty = {last_element_r, last_element_theta,last_element_phi ,last_element_ur ,last_element_utheta ,last_element_uphi};
            double time = rkdpA.xs[num-1];
            
            // std::cout << "xs: " << xgrid << std::endl;
            // std::cout << "ys: " << ygrid << std::endl; 
            // std::cout << "pointsx: " << i << std::endl;
            // std::cout << "pointsy: " << j << std::endl; 
            // std::cout << "first_element_r: " << r0 << std::endl;
            // std::cout << "first_element_theta: " << theta0 << std::endl;

            // std::cout << "last_element_r: " << last_element_r << std::endl;
            // std::cout << "last_element_theta: " << last_element_theta << std::endl;

            // std::cout << "time: " << time << std::endl;

         
           
            if(((last_element_theta < M_PI/2.0+tol )&&(last_element_theta > M_PI/2.0-tol ))&&( (last_element_r> +r0Accin-tol)&&(last_element_r < r0Accout+tol ))){

                ///Now include doppler effect
                metricdop.compute_metric(last_element_r, last_element_theta);
                double big_omega = 1.0/(a + pow(last_element_r, 1.5)/sqrt(M));
                double valz = (1.0 +big_omega*(last_element_uphi/lastut(lasty)))/sqrt(-metricdop.g_00 - pow(big_omega,2)*metricdop.g_33 - 2.0*big_omega*metricdop.g_03);

                double Intensity = 1/(pow(1 + valz,3));

                grid[i+pointsy/2][j+pointsx/2] = Intensity*1.0;
                //std::cout << "YEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH"  << std::endl;
                
            } else{
                grid[i+pointsy/2][j+pointsx/2] = 0.0;
            };
            

        }


        std::cout << "pointsx: " << i << std::endl;
    };

    std::cout << "done"  << std::endl;
   



  

    // Output the results

    std::ofstream output_file("output_blackhole_120.csv");
    for (int i = 0; i < pointsy; i++) {
        for (int j = 0; j < pointsx; j++){
       

        if (j < (pointsx - 1)) {
                output_file << grid[i][j] << ",";
            }
            else if (j == (pointsx - 1)) {
                output_file << grid[i][j] <<"\n" <<std::endl;
            }

        }

    }
    output_file.close();


    return 0;
};