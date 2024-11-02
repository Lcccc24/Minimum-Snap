#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    //VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);


    int index_col = 0;
    int index_row = 0;

    //M矩阵
    MatrixXd M = MatrixXd::Zero(m*d_order*2,p_num1d*m); 

    //循环迭代
    // for(int i = 0; i < m; ++i){
    //     for(int j = 0; j < d_order; ++j){
    //         for(int k = j; k < p_num1d; ++k){
    //             int index_row = d_order*2*i;
    //             int index_col = p_num1d*i;
    //             M(index_row+j+1,index_col+k+1) = Factorial(k)/Factorial(k-j)*std::pow(0,k-j);
    //             M(index_row-d_order+j+1,index_col+k+1) = Factorial(k)/Factorial(k-j)*std::pow(Time[i+1],k-j);
    //         }
    //     }
    // }

    //拼接M矩阵
    ROS_INFO("M ING");
    for(int i = 0; i < m; ++i)
    {
        int M_ROW = i * d_order*2;
        int M_COL = i * p_num1d;
        MatrixXd M_sub = MatrixXd::Zero(d_order*2,p_num1d); 
        for(int j = 0; j < d_order; ++j){
            for(int k = j; k < p_num1d; ++k){
                M_sub(j,k) = Factorial(k)/Factorial(k-j)*std::pow(0,k-j);
                M_sub(j+d_order,k) = Factorial(k)/Factorial(k-j)*std::pow(Time[i],k-j);
            }
        }

        M.block(M_ROW,M_COL,d_order*2,p_num1d) = M_sub;
    }
    ROS_INFO("M done");
    ROS_INFO_STREAM("Matrix M:\n" << M);

    //Q矩阵
    ROS_INFO("Q ING");
    MatrixXd Q = MatrixXd::Zero(p_num1d*m,p_num1d*m);
    for(int i = 0; i < m; ++i)
    {
        int Q_ROW_COL = i * p_num1d;
        MatrixXd Q_sub = MatrixXd::Zero(p_num1d,p_num1d);
        for(int j = d_order; j < p_num1d; ++j){
            for(int k = d_order; k < p_num1d; ++k){
                Q_sub(j,k) = Factorial(k)/Factorial(k-d_order)*Factorial(j)/Factorial(j-d_order)/(j+k-2*d_order+1)*pow(Time[i],j+k-2*d_order-1);
            }
        }

        Q.block(Q_ROW_COL,Q_ROW_COL,p_num1d,p_num1d) = Q_sub;
    }
    ROS_INFO("Q done");
    ROS_INFO_STREAM("Matrix Q:\n" << Q);

    //C矩阵
    ROS_INFO("C_T ING");
    MatrixXd C_T = MatrixXd::Zero(2*m*d_order,(m+1)*d_order);

    //起点  
    for(int i = 0; i < d_order; ++i)
        C_T(i,i) = 1;

    //中间点
    for(int i = 0; i < m-1; ++i){
        index_col = d_order+i; 
        index_row = d_order*(2*i+1); //d_order+2*d_order*i

        C_T(index_row,index_col) = 1;
        C_T(index_row+d_order,index_col) = 1;
    }

    //终点
    for(int i = 0; i < d_order; ++i){
        index_col = d_order+m-1;
        index_row = d_order*(2*m-1); //2m*d_order-d_order
        C_T(index_row+i,index_col+i) = 1;
    }
    
    //中间点微分
    for(int i = 0; i < m-1; ++i){
        index_col = d_order*(i+2) +m-i-1; //d_order*2 + m-1 + (d_order-1)*i
        index_row = d_order*(2*i+1); //d_order+2*d_order*i
        for(int j = 0; j < d_order-1; ++j){
            C_T(index_row+j+1,index_col+j) = 1;
            C_T(index_row+d_order+j+1,index_col+j) = 1;
        }
    }
    ROS_INFO("C_T done");
    ROS_INFO_STREAM("Matrix C_T:\n" << C_T);

    MatrixXd R = C_T.transpose() * M.transpose().inverse() * Q * M.inverse() * C_T;
    ROS_INFO("R done");

    std::vector<MatrixXd> PolyCoeff_xyz;

    for(int xyz = 0; xyz < 3; ++xyz){

        VectorXd dF = VectorXd::Zero(d_order*2+m-1);
        VectorXd dFP = VectorXd::Zero(d_order*(m+1));
        
        //起点pva j=0
        dF(0) = Path(0,xyz);
        dF(1) = Vel(0,xyz);
        dF(2) = Acc(0,xyz);

        //终点pva j=0
        dF(d_order+m-1) = Path(Path.rows()-1,xyz);
        dF(d_order+m-1 +1) = Vel(1,xyz);
        dF(d_order+m-1 +2) = Acc(1,xyz);

        //中间点P
        for(int i = 0; i < m-1; ++i)
            dF(d_order+i) = Path(i+1,xyz);

        //ROS_INFO("dF done");

        MatrixXd Rpp = R.block(2*d_order+m-1,2*d_order+m-1,(d_order-1)*(m-1),(d_order-1)*(m-1)); //从第几行第几列开始，插入多少行多少列
        MatrixXd Rfp = R.block(0,2*d_order+m-1,2*d_order+m-1,(d_order-1)*(m-1));
            
        MatrixXd dP = -Rpp.inverse()*Rfp.transpose()*dF;

        VectorXd d(dF.size()+dP.size());
        d << dF, dP;

        VectorXd P = M.inverse() * C_T * d;

        MatrixXd Coeff = MatrixXd::Zero(m, p_num1d);  
        for(int i = 0; i < m; ++i){
            Coeff.block(i,0,1,p_num1d) = P.block(p_num1d*i,0,p_num1d,1).transpose();
        }

        PolyCoeff_xyz.push_back(Coeff);
    }

    // 将xyz3个矩阵按列复制到 PolyCoeff 中
    for (int i = 0; i < PolyCoeff_xyz.size(); ++i) {
        PolyCoeff.block(0, i*p_num1d, m, p_num1d) = PolyCoeff_xyz[i];
    }

    ROS_INFO("PolyCoeff done");

    return PolyCoeff;
}
