#include "trajectory.h"

double PredictPitchXY::v0now;

PredictPitchXY PredictPitchXY::instance; // 一定要初始化静态变量

double PredictPitchXY::setBulletSpeed(float bulletSpeed)
{
    if (bulletSpeed != 0)
    {
        return bulletSpeed;
    }
    else
    {
        return 15.7;
    }
}

void PredictPitchXY::InitPredictPitch(float bulletSpeedNow, double distance, float hrRec)
{
    // 赋值初始量
    v0now = setBulletSpeed(bulletSpeedNow);
    X_r = distance;
    Y_r = hrRec;
}

// 8.27：应当以时间为单位作为步长，或者说我下面使用的步长正是时间而我理解为了距离（受到沈航代码先入为主的影响），所以才会出现角度变动如此之大和稳定的问题
double PredictPitchXY::dropshotRK45()
{
    double error = 0;
    theta = 15 * 3.14 / 180.0; // 初始化角度
    for (int i = 0; i < maxstep; i++)
    { // 设置最大步长
        // 初始化
        theta_d = theta;                   // 将theta赋给过程量
        X_d = 0.17 * cos(theta_d);         // 枪管水平长度
        Y_d = 0.424 + 0.17 * sin(theta_d); // 枪管垂直高度
        v0_d = v0now;
        // double x_step=(X_r-X_d)/step;//迭代步长
        double time_step = 0.0004; // 时间步长
        while (X_d < X_r)          // 迭代
        {
            // 计算各阶

            auto k1_v = (-k * pow(v0_d, 2) - g * sin(theta_d)) * time_step;
            auto k1_theta = (-g * cos(theta_d) / v0_d) * time_step; // 由于大弹丸很难产生马格努斯效应（横向摩擦轮怎么可能产生后向旋转），故忽略升力
            // 其实这里和公式是不完全相同的，未考虑t随步长的改变，因为t在推导过程中简化约去了

            auto k1_v_2 = v0_d + k1_v / 4.0;
            auto k1_theta_2 = theta_d + k1_theta / 4.0;

            auto k2_v = (-k * pow(k1_v_2, 2) - g * sin(k1_theta_2)) * time_step;
            auto k2_theta = (-g * cos(k1_theta_2) / k1_v_2) * time_step;
            auto k12_v_3 = v0_d + 3.0 / 32.0 * k1_v + 9.0 / 32.0 * k2_v;
            auto k12_theta_3 = theta_d + 3.0 / 32.0 * k1_theta + 9.0 / 32.0 * k2_theta;

            auto k3_v = (-k * pow(k12_v_3, 2) - g * sin(k12_theta_3)) * time_step;
            auto k3_theta = (-g * cos(k12_theta_3) / k12_v_3) * time_step;
            auto k123_v_4 = v0_d + 1932.0 / 2179.0 * k1_v - 7200.0 / 2179.0 * k2_v + 7296.0 / 2179.0 * k3_v;
            auto k123_theta_4 = theta_d + 1932.0 / 2179.0 * k1_theta - 7200.0 / 2179.0 * k2_theta + 7296.0 / 2179.0 * k3_theta;

            auto k4_v = (-k * pow(k123_v_4, 2) - g * sin(k123_theta_4)) * time_step;
            auto k4_theta = (-g * cos(k123_theta_4) / k123_v_4) * time_step;
            auto k1234_v_5 = v0_d + 439.0 / 216.0 * k1_v - 8.0 * k2_v + 3680.0 / 513.0 * k3_v - 845.0 / 4140.0 * k4_v;
            auto k1234_theta_5 = theta_d + 439.0 / 216.0 * k1_theta - 8.0 * k2_theta + 3680.0 / 513.0 * k3_theta - 845.0 / 4140.0 * k4_theta;

            auto k5_v = (-k * pow(k1234_v_5, 2) - g * sin(k1234_theta_5)) * time_step;
            auto k5_theta = (-g * cos(k1234_theta_5) / k1234_v_5) * time_step;
            auto k12345_v_6 = v0_d - 8.0 / 27.0 * k1_v + 2.0 * k2_v - 3544.0 / 2565.0 * k3_v + 1859.0 / 4104.0 * k4_v - 11.0 / 40.0 * k5_v;
            auto k12345_theta_6 = theta_d - 8.0 / 27.0 * k1_theta + 2.0 * k2_theta - 3544.0 / 2565.0 * k3_theta + 1859.0 / 4104.0 * k4_theta - 11.0 / 40.0 * k5_theta;

            auto k6_v = (-k * pow(k12345_v_6, 2) - g * sin(k12345_theta_6)) * time_step;
            auto k6_theta = (-g * cos(k12345_theta_6) / k12345_v_6) * time_step;

            // 计算近似解以及相对误差
            // 注：计算相对误差采用五阶精度与四阶精度系数结果之差
            auto vclass_4 = v0_d + 25.0 / 216.0 * k1_v + 1408.0 / 2565.0 * k3_v + 2197.0 / 4104.0 * k4_v - 1.0 / 5.0 * k5_v;
            auto thetaclass_4 = theta_d + 25.0 / 216.0 * k1_theta + 1408.0 / 2565.0 * k3_theta + 2197.0 / 4104.0 * k4_theta - 1.0 / 5.0 * k5_theta;

            auto vclass_5 = v0_d + 16.0 / 135.0 * k1_v + 6656.0 / 12825.0 * k3_v + 28561.0 / 56430.0 * k4_v - 9.0 / 50.0 * k5_v + 2.0 / 55.0 * k6_v;
            auto thetaclass_5 = theta_d + 16.0 / 135.0 * k1_theta + 6656.0 / 12825.0 * k3_theta + 28561.0 / 56430.0 * k4_theta - 9.0 / 50.0 * k5_theta + 2.0 / 55.0 * k6_theta;

            auto error_v = abs(vclass_5 - vclass_4);
            auto error_theta = abs(thetaclass_5 - thetaclass_4);
            auto error_y = tan(error_theta) * time_step;

            // cout << "error_v:" << error_v << endl;
            // cout << "error_theta:" << error_theta << endl;
            // cout << "error_y:" << error_y << endl;

            v0_d = vclass_5;
            theta_d = thetaclass_5;
            X_d += time_step * v0_d * cos(theta_d);
            Y_d += time_step * v0_d * sin(theta_d);
        }

        // 评估迭代结果，修正theta
        error = Y_r - Y_d; // error可以pid调节
        if (abs(error) < minerr_y)
        {
            std::cout << i << std::endl;
            std::cout << "..............................................." << std::endl;
            std::cout << "theta:" << theta * 180 / 3.1415 << std::endl;
            std::cout << "..............................................." << std::endl;
            return theta;
        } // 合适则输出本次迭代使用的theta
        else
        {
            theta += atan((error) / X_r);
            // std::cout << "...................................................." << std::endl;
            // std::cout << "error:" << error << std::endl;
            // std::cout << "theta_d:" << theta_d * 180 / 3.1415 << endl;
            // std::cout << "...................................................." << std::endl;
        }
    }

    // 迭代失败则输出上次的迭代值(暂定为0)
    return 0;
}
