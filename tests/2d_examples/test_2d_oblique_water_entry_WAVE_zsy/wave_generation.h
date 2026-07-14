#ifndef WAVE_GENERATION_H
#define WAVE_GENERATION_H

#include "sphinxsys.h"  // 包含 SPHinXsys 的类型定义，如 Real、Vec2d 等
#include <functional>
#include <vector>
#include <cmath>

using namespace SPH;

// 波形函数类型：输入时间 t，输出推板位移 disp 和速度 vel
using WaveFormFunc = std::function<void(Real t, Real &disp, Real &vel)>;


// statements only, definitions are in .cpp
Real solveDispersionEquation(Real omega, Real h, Real g, Real tol = 1e-6);
Real computePistonStroke(Real H, Real k, Real h);

// ---------- 1. 规则波（使用传递函数） ----------
// 参数: H - 目标波高 (m), T - 周期 (s), phase - 初始相位 (rad), h - 水深 (m), g - 重力加速度
inline WaveFormFunc createRegularWave(Real H, Real T, Real phase, Real h, Real g = 9.81)
{
    Real omega = 2.0 * Pi / T;
    Real k = solveDispersionEquation(omega, h, g);
    Real stroke = computePistonStroke(H, k, h);
    Real amplitude = 0.5 * stroke;   // 推板振幅
    return [amplitude, omega, phase](Real t, Real &disp, Real &vel) {
        disp = amplitude * sin(omega * t + phase);
        vel  = amplitude * omega * cos(omega * t + phase);
    };
}

// ---------- 2. 双色波（基于波高和周期） ----------
// 参数: H1, T1 - 第一成分波高、周期; H2, T2 - 第二成分波高、周期; delta_phi - 第二成分相对第一成分的相位差 (rad)
inline WaveFormFunc createBiChromaticWave(Real H1, Real T1, Real H2, Real T2,
                                          Real delta_phi, Real h, Real g = 9.81)
{
    Real omega1 = 2.0 * Pi / T1;
    Real omega2 = 2.0 * Pi / T2;
    Real k1 = solveDispersionEquation(omega1, h, g); // k1和k2是根据周期T1和T2通过色散关系求解得到的波数
    Real k2 = solveDispersionEquation(omega2, h, g);
    Real a1 = 0.5 * computePistonStroke(H1, k1, h);
    Real a2 = 0.5 * computePistonStroke(H2, k2, h);
    return [a1, omega1, a2, omega2, delta_phi](Real t, Real &disp, Real &vel) {
        disp = a1 * sin(omega1 * t) + a2 * sin(omega2 * t + delta_phi);
        vel  = a1 * omega1 * cos(omega1 * t) + a2 * omega2 * cos(omega2 * t + delta_phi);
    };
}

// ---------- 3. 双色波（基于波陡，可选） ----------
// 参数: f1 - 第一成分频率 (Hz), delta_f - 频差 (f2 = f1 - delta_f), epsilon1 - 第一成分波陡 (a1*k1),
//       ratio_epsilon - 第二成分波陡与第一成分的比值, delta_phi - 相位差
inline WaveFormFunc createBiChromaticWaveFromSteepness(Real f1, Real delta_f, Real epsilon1,
                                                       Real ratio_epsilon, Real delta_phi,
                                                       Real h, Real g = 9.81)
{
    Real omega1 = 2.0 * Pi * f1;
    Real omega2 = 2.0 * Pi * (f1 - delta_f);
    Real k1 = solveDispersionEquation(omega1, h, g);
    Real k2 = solveDispersionEquation(omega2, h, g);
    Real a1_target = epsilon1 / k1;               // 目标波浪振幅 (m)
    Real a2_target = ratio_epsilon * epsilon1 / k2;
    // 转换为推板振幅
    Real H1 = 2.0 * a1_target, H2 = 2.0 * a2_target;
    Real a1_push = 0.5 * computePistonStroke(H1, k1, h);
    Real a2_push = 0.5 * computePistonStroke(H2, k2, h);
    return [a1_push, omega1, a2_push, omega2, delta_phi](Real t, Real &disp, Real &vel) {
        disp = a1_push * sin(omega1 * t) + a2_push * sin(omega2 * t + delta_phi);
        vel  = a1_push * omega1 * cos(omega1 * t) + a2_push * omega2 * cos(omega2 * t + delta_phi);
    };
}

// ---------- 4. 聚焦波（高斯谱 + 聚焦条件） ----------
// 参数:
//   Af        - 谱峰处对应的目标波浪振幅 (m)（即波高的一半）
//   fp        - 谱峰频率 (Hz)
//   bandwidth - 频率带宽 (Hz)，频率范围为 [fp - bandwidth/2, fp + bandwidth/2]
//   Nf        - 离散频率数量（至少为2）
//   tf        - 聚焦时刻 (s)
//   xf        - 聚焦位置 (m)（通常取水槽中点）
//   h         - 水深 (m)
//   g         - 重力加速度 (m/s²)，默认 9.81
//
// 返回: WaveFormFunc，可直接用于 WaveMaking 类

inline WaveFormFunc createFocusedWave(Real Af, Real fp, Real bandwidth, int Nf,
                                      Real tf, Real xf, Real h, Real g = 9.81)
{
    std::cout << "Inside createFocusedWave: tf = " << tf << ", xf = " << xf << std::endl;
    if (Nf < 2)
        Nf = 2;
    Real omega_peak = 2.0 * Pi * fp;
    Real domega = 2.0 * Pi * (bandwidth / (Nf - 1));
    Real omega_min = omega_peak - 0.5 * (2.0 * Pi * bandwidth);
    Real sigma_omega = 2.0 * Pi * (bandwidth / 6.0);

    std::vector<Real> omegas(Nf), k(Nf), push_amps(Nf), phases(Nf);

    // 1. 计算频率数组
    for (int i = 0; i < Nf; ++i)
    {
        omegas[i] = omega_min + i * domega;
    }
    // 2. 计算高斯谱密度值（未归一化）
    std::vector<Real> S(Nf);
    Real sum_S = 0.0;
    for (int i = 0; i < Nf; ++i)
    {
        S[i] = std::exp(-0.5 * std::pow((omegas[i] - omega_peak) / sigma_omega, 2));
        sum_S += S[i] * domega;
    }
    // 3. 设计振幅（聚焦点最大波幅）
    Real A_design = Af; // Af 现在代表期望的聚焦点最大波幅

    // 4. 计算每个频率的推板振幅和相位
    for (int i = 0; i < Nf; ++i)
    {
        Real a_i = A_design * (S[i] * domega) / sum_S;
        k[i] = solveDispersionEquation(omegas[i], h, g);
        Real stroke = computePistonStroke(2.0 * a_i, k[i], h);
        push_amps[i] = 0.5 * stroke;
        phases[i] = -omegas[i] * tf + k[i] * xf;
    }

    if (Nf > 0)
    {
        Real c0 = omegas[0] / k[0];
        Real travel_time = xf / c0;
        std::cout << "Focus wave: central speed = " << c0 << " m/s, travel time to xf = " << travel_time << " s" << std::endl;
        Real t_peak = (Pi / 2.0 - phases[0]) / omegas[0];
        Real t_expected = tf - xf / (omegas[0] / k[0]);
        std::cout << "First component: omega=" << omegas[0] << ", k=" << k[0]
                  << ", phase=" << phases[0] << ", t_peak=" << t_peak
                  << ", expected t_peak=" << t_expected << std::endl;
    }

    return [push_amps, omegas, phases](Real t, Real &disp, Real &vel)
    {
        disp = 0.0;
        vel = 0.0;
        for (size_t i = 0; i < push_amps.size(); ++i)
        {
            disp += push_amps[i] * cos(omegas[i] * t + phases[i]);
            vel -= push_amps[i] * omegas[i] * sin(omegas[i] * t + phases[i]);//求导负数
        }

    };
}

// ---------- 辅助函数：计算聚焦波在任意时刻、任意位置的波面高度（用于初始条件）----------
// 参数与 createFocusedWave 完全相同，额外输入位置 x (m) 和时间 t (s)
// 返回波面高度 η (m)（注意：该函数不涉及推板传递函数，直接输出目标波浪振幅叠加结果）
//inline Real evaluateFocusedWaveElevation(Real Af, Real fp, Real bandwidth, int Nf,
//                                         Real tf, Real xf, Real h, Real g,
//                                         Real x, Real t)
//{
//    if (Nf < 2)
//        Nf = 2;
//    Real omega_peak = 2.0 * Pi * fp;
//    Real domega = 2.0 * Pi * (bandwidth / (Nf - 1));
//    Real omega_min = omega_peak - 0.5 * (2.0 * Pi * bandwidth);
//    Real sigma_omega = 2.0 * Pi * (bandwidth / 6.0);
//
//    std::vector<Real> target_amps(Nf), omegas(Nf), k(Nf), phases(Nf);
//    for (int i = 0; i < Nf; ++i)
//    {
//        omegas[i] = omega_min + i * domega;
//        Real env = std::exp(-0.5 * std::pow((omegas[i] - omega_peak) / sigma_omega, 2));
//        target_amps[i] = env;
//    }
//    Real max_env = *std::max_element(target_amps.begin(), target_amps.end());
//    for (int i = 0; i < Nf; ++i)
//    {
//        target_amps[i] = Af * (target_amps[i] / max_env);
//        k[i] = solveDispersionEquation(omegas[i], h, g);
//        // 与造波板相位保持一致
//        phases[i] = -omegas[i] * tf +  k[i] * xf;
//    }
//
//    Real eta = 0.0;
//    for (int i = 0; i < Nf; ++i)
//    {
//        //  波面 η = a * cos(ωt - kx + φ)
//        eta += target_amps[i] * cos(omegas[i] * t - k[i] * x + phases[i]);
//    }
//    return eta;
//}
inline Real evaluateFocusedWaveElevation(Real Af, Real fp, Real bandwidth, int Nf,
                                         Real tf, Real xf, Real h, Real g,
                                         Real x, Real t)
{
    if (Nf < 2)
        Nf = 2;
    Real omega_peak = 2.0 * Pi * fp;
    Real domega = 2.0 * Pi * (bandwidth / (Nf - 1));
    Real omega_min = omega_peak - 0.5 * (2.0 * Pi * bandwidth);
    Real sigma_omega = 2.0 * Pi * (bandwidth / 6.0);

    std::vector<Real> omegas(Nf), S(Nf), target_amps(Nf), k(Nf), phases(Nf);

    // 1. 计算频率和高斯谱密度 S(ω)
    for (int i = 0; i < Nf; ++i)
    {
        omegas[i] = omega_min + i * domega;
        S[i] = std::exp(-0.5 * std::pow((omegas[i] - omega_peak) / sigma_omega, 2));
    }

    // 2. 计算谱密度积分（能量归一化因子）
    Real sum_S = 0.0;
    for (int i = 0; i < Nf; ++i)
    {
        sum_S += S[i] * domega;
    }

    // 3. 分配各频率成分的目标波幅（与 createFocusedWave 一致）
    for (int i = 0; i < Nf; ++i)
    {
        target_amps[i] = Af * (S[i] * domega) / sum_S; // 与造波板振幅分配一致

        // 色散关系和相位（与造波板保持一致）
        k[i] = solveDispersionEquation(omegas[i], h, g);
        phases[i] = -omegas[i] * tf + k[i] * xf;
    }

    // 4. 计算给定位置 x、时间 t 的波面高度
    Real eta = 0.0;
    for (int i = 0; i < Nf; ++i)
    {
        eta += target_amps[i] * cos(omegas[i] * t - k[i] * x + phases[i]);
    }
    return eta;
}

#endif // WAVE_GENERATION_H