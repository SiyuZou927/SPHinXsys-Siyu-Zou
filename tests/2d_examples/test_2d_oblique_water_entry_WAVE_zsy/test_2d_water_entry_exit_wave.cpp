/**
 * @file	oblique water entry with wave.cpp
 * @brief	2D water entry and exit example with surface wetting considered.
 * @details	This is the one of FSI test cases, also one case for
 * 			understanding spatial temporal identification approach,
 *          especially when coupled with the wetting.
 *          Modified to include wave generation by a piston-type wavemaker.
 * @author  Shuoguo Zhang and Xiangyu Hu (original)
 * @modified Siyu Zou
 */
// #include "2d_flow_around_cylinder.h"
#include "sphinxsys.h" //SPHinXsys Library.
#include "wave_generation.h"
using namespace SPH;   // Namespace cite here.

//----------------------------------------------------------------------
//	核心物理参数（实际工程值）
//----------------------------------------------------------------------
// 几何参数
Real actual_cylinder_length = 0.4796; /**< 圆柱总长 (m) */
Real actual_cylinder_radius = 0.025;  /**< 圆柱半径 (m) */
// 质量与转动惯量（kg·m²）
Real actual_total_mass = 1.855; /**< 圆柱总质量 (kg) */
Real Ix = 0.0009;               /**< 绕圆柱轴线（x轴）转动惯量 */
Real Iy = 0.020;                /**< 绕y轴转动惯量 */
Real Iz = 0.020;                /**< 绕z轴转动惯量（2D模拟主惯量） */
Vec2d centroid(0.0, 0.0);       /**< 圆柱质心位置（相对于圆柱中心） */
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4;                         /**< Water tank length. */
Real DH = 5;                        /**< Water tank height. */
Real LH = 2;                         /**< Water column height. */
Real cavity_length = 1;             // leftover length for wave development
Real particle_spacing_ref = 0.02;    /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;    /**< Thickness of tank wall. */

//波浪控制参数
Real release_time = 5; // 造波时长
bool released = false; // 是否已释放

// output control parameters
int pre_output_count = release_time*10;// vtp output count before release
int post_output_count = 20;        // vtp output count after release
int post_froce_output_count = 500;  // force output count after release

// 圆柱初始位置（根据波浪动态调整）
Vec2d cylinder_center; /**< 实际将在 main 中计算 */

// 初始速度参数
Real initial_speed = 70;                        /**< Initial velocity magnitude (m/s). */
Real initial_angle = -20 * Pi / 180.0;          /**< Initial velocity angle (radians, negative = downward). */
Real initial_rotation_angle = -20 * Pi / 180.0; /**< Initial body rotation (radians). */
Real initial_angular_velocity = 0.0;            // 可选：初始角速度（rad/s）（如果需要旋转入水，可设置）
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Fluid density. */
Real rho0_s = 2.063;                     /**< Cylinder density. */
Real gravity_g = 9.81;                   /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max; /**< Reference sound speed. */
Real mu_f = 8.9e-7;      /**< Water dynamics viscosity. */
//----------------------------------------------------------------------
//	Wetting parameters
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";                     //  ϕ∗
Real diffusion_coeff = 100.0 * pow(particle_spacing_ref, 2); /**< Wetting coefficient. γ∗ */
Real fluid_moisture = 1.0;                                   /**< fluid moisture. */
Real cylinder_moisture = 0.0;                                /**< cylinder moisture. */
Real wall_moisture = 0.0;                                    /**< wall moisture. */

//----------------------------------------------------------------------
// Create the object shape
//----------------------------------------------------------------------
std::vector<Vecd> createObjectShape()
{
    std::vector<Vecd> object_shape = {
        Vecd(0.2780, 0.01055),
        Vecd(0.2780, -0.01055), 
        Vecd(0.0843, -0.02925),
        Vecd(-0.1249, -0.02925),
        Vecd(-0.1249, -0.0180),
        Vecd(-0.2016, -0.0180),
        Vecd(-0.2016, 0.0180), 
        Vecd(-0.1249, 0.0180),
        Vecd(-0.1249, 0.02925),
        Vecd(0.0843, 0.02925),
        Vecd(0.2780, 0.01055) 
    };
    return object_shape;
}

//----------------------------------------------------------------------
// Reset the front center observation point position after rotation and translation
//----------------------------------------------------------------------
Vecd resetFrontCenterObserverPosition()
{
    Vec2d front_center_initial(0.2780, -0.01055);
    Vec2d rotation_radius = front_center_initial - centroid;
    Rotation2d rotation(initial_rotation_angle);
    Vec2d front_center_rotated = rotation * rotation_radius;
    // here the cylinder center should be the relativa distance to the global center.
    Vec2d front_center_translated = front_center_rotated + cylinder_center;
    return front_center_translated;
}

//----------------------------------------------------------------------
// 坐标变换函数 - 全局力转换为圆柱自身坐标系力
//----------------------------------------------------------------------
Vec2d transformGlobalForceToLocal(const Vecd &global_force, Real rotation_angle)
{
    Real cos_theta = std::cos(rotation_angle);
    Real sin_theta = std::sin(rotation_angle);
    Real local_x = global_force[0] * cos_theta + global_force[1] * sin_theta;  // local-x force
    Real local_y = -global_force[0] * sin_theta + global_force[1] * cos_theta; // local-y force
    return Vec2d(local_x, local_y);
}
//----------------------------------------------------------------------
// 获取圆柱实时旋转角度（从Simbody State中提取）
//----------------------------------------------------------------------
Real getCylinderRotationAngle(SimTK::MobilizedBody::Planar &tethered_spot, const SimTK::State &state)
{
    SimTK::Rotation rot = tethered_spot.getBodyRotation(state);
    SimTK::Vec3 angles = rot.convertRotationToBodyFixedXYZ();
    Real angle_from_rot = angles[2]; // z-axis rotation
    return angle_from_rot;          
}

//----------------------------------------------------------------------
// 波浪相关辅助函数
//----------------------------------------------------------------------
// 求解线性波色散方程：omega^2 = g * k * tanh(k * h)
Real solveDispersionEquation(Real omega, Real h, Real g, Real tol)
{
    Real k0 = omega * omega / g; // 深水近似初值
    Real k = k0;
    for (int iter = 0; iter < 20; ++iter)
    {
        Real f = g * k * tanh(k * h) - omega * omega;
        Real df = g * (tanh(k * h) + k * h * (1.0 - tanh(k * h) * tanh(k * h)));
        Real k_new = k - f / df;
        if (fabs((k_new - k) / k_new) < tol)
            return k_new;
        k = k_new;
    }
    std::cout << "Warning: Dispersion equation not converged. Using last k = " << k << std::endl;
    return k;
}

// 计算推板冲程（线性Biesel传递函数）
Real computePistonStroke(Real H, Real k, Real h)
{
    Real kh = k * h;
    Real sinh_kh = sinh(kh);
    Real cosh_kh = cosh(kh);
    // 传递函数 T = H / S = 2 * (sinh^2(kh) / (sinh(kh)cosh(kh) + kh))
    Real transfer = 2.0 * sinh_kh * sinh_kh / (sinh_kh * cosh_kh + kh);
    return H / transfer;
}

// 定义一个辅助函数，生成一个竖直条状的探头形状（宽度 2*h，从底部到水面以上）
MultiPolygon createWaveProbeShape(Real x_center, Real h, Real water_depth, Real tank_height)
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(x_center - h, 0.0));
    pnts.push_back(Vecd(x_center - h, tank_height));
    pnts.push_back(Vecd(x_center + h, tank_height));
    pnts.push_back(Vecd(x_center + h, 0.0));
    pnts.push_back(Vecd(x_center - h, 0.0));
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}

// 具体探头位置（单位：米）
Real probe_x1 = 6;                       // 近造波板
Real probe_x2 = 11.54;                       // 水槽中部
Real probe_x3 = 13.08;                       // 近物体（物体初始位置 x = 0.2*DL = 1.6，可调整）
Real probe_h = 1.3 * particle_spacing_ref; // 探头宽度半高

//----------------------------------------------------------------------
// 造波板动力学类
//----------------------------------------------------------------------
class WaveMaking : public BodyPartMotionConstraint
{
    WaveFormFunc wave_func_;
    Real *physical_time_;
  public:
    static Real current_vel;
    static Real getCurrentVelocity() { return current_vel; }
    WaveMaking(BodyPartByParticle &body_part, WaveFormFunc func)
        : BodyPartMotionConstraint(body_part),
          wave_func_(func),
          physical_time_(body_part.getSPHBody().getSPHSystem().getSystemVariableDataByName<Real>("PhysicalTime")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        Real time = *physical_time_;
        Real disp = 0.0, vel = 0.0;
        wave_func_(time, disp, vel);
        pos_[index_i] = pos0_[index_i] + Vecd(disp, 0.0);
        vel_[index_i] = Vecd(vel, 0.0);
        current_vel = vel;
    }
};
Real WaveMaking::current_vel = 0.0;
//----------------------------------------------------------------------
// 消波区域形状
//----------------------------------------------------------------------
MultiPolygon createDampingBufferShape()
{
    std::vector<Vecd> pnts;
    Real damping_start = DL - 0.5; // 从距离右端2m 处开始阻尼
    pnts.push_back(Vecd(damping_start, 0.0));
    pnts.push_back(Vecd(damping_start, DH));
    pnts.push_back(Vecd(DL + BW, DH));
    pnts.push_back(Vecd(DL + BW, 0.0));
    pnts.push_back(Vecd(damping_start, 0.0));
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
// observer_location
//----------------------------------------------------------------------
StdVec<Vecd> front_center_observer_location = {Vecd(0, 0)}; // 占位，稍后更新
//----------------------------------------------------------------------
// Output summary file.
//----------------------------------------------------------------------
class SummaryOutput
{
  public:
    SummaryOutput(const std::string &filename)
        : output_file_(filename)
    {
        output_file_ << "Time[s]   "
                     << "ViscousForceGlobal_X ViscousForceGlobal_Y "
                     << "PressureForceGlobal_X PressureForceGlobal_Y "
                     << "ViscousForceLocalX PressureForceLocalX TotalForceLocalX "
                     << "ViscousForceLocalY PressureForceLocalY TotalForceLocalY "
                     << "FrontCenterObserver_Position_X FrontCenterObserver_Position_Y\n";
    }

    ~SummaryOutput()
    {
        output_file_.close();
    }
    void writeData(Real time,
                   const Vec2d &total_viscous_force_global,
                   const Vec2d &total_pressure_force_global,
                   Vec2d viscous_local,
                   Vec2d pressure_local,
                   const Vecd &front_center_position)

    {
        Real total_force_local_x = viscous_local[0] + pressure_local[0];
        Real total_force_local_y = viscous_local[1] + pressure_local[1];

        output_file_ << std::scientific << std::setprecision(9)
                     << time << " "
                     << total_viscous_force_global[0] << " " << total_viscous_force_global[1] << " "
                     << total_pressure_force_global[0] << " " << total_pressure_force_global[1] << " "
                     << viscous_local[0] << " " << pressure_local[0] << " " << total_force_local_x << " "
                     << viscous_local[1] << " " << pressure_local[1] << " " << total_force_local_y << " "
                     << front_center_position[0] << " " << front_center_position[1] << "\n";

        output_file_.flush();
    }

  protected:
    std::ofstream output_file_;
};
//----------------------------------------------------------------------
//	Definition for water body
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block;
    water_block.push_back(Vecd(0.0, 0.0));
    water_block.push_back(Vecd(0.0, LH));
    water_block.push_back(Vecd(DL, LH));
    water_block.push_back(Vecd(DL, 0.0));
    water_block.push_back(Vecd(0.0, 0.0));

    return water_block;
}
class WettingFluidBody : public MultiPolygonShape
{
  public:
    explicit WettingFluidBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};

class WettingFluidBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingFluidBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = fluid_moisture;
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	Definition for wall body
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall;
    outer_wall.push_back(Vecd(-cavity_length - BW, -BW));      
    outer_wall.push_back(Vecd(-cavity_length - BW, DH + BW)); 
    outer_wall.push_back(Vecd(DL + BW, DH + BW)); 
    outer_wall.push_back(Vecd(DL + BW, -BW));    
    outer_wall.push_back(Vecd(-cavity_length - BW, -BW));      
    return outer_wall;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall;
    inner_wall.push_back(Vecd(0.0, 0.0));
    inner_wall.push_back(Vecd(0.0, DH + BW)); // 顶部提升到外壁顶部
    inner_wall.push_back(Vecd(DL, DH + BW));  // 顶部提升到外壁顶部
    inner_wall.push_back(Vecd(DL, 0.0));
    inner_wall.push_back(Vecd(0.0, 0.0));
    return inner_wall;
}

std::vector<Vecd> createCavityInnerShape()
{
    std::vector<Vecd> inner;
    inner.push_back(Vecd(-cavity_length, 0.0));
    inner.push_back(Vecd(-cavity_length, DH + BW)); // 顶部提升到外壁顶部
    inner.push_back(Vecd(-BW, DH + BW));            // 顶部提升到外壁顶部
    inner.push_back(Vecd(-BW, 0.0));
    inner.push_back(Vecd(-cavity_length, 0.0));
    return inner;
}
class WettingWallBody : public MultiPolygonShape
{
  public:
    explicit WettingWallBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(createCavityInnerShape(), ShapeBooleanOps::sub);
    }
};
class WettingWallBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingWallBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = wall_moisture;
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	Definition for ObjectBody body
//----------------------------------------------------------------------
class ObjectBody : public MultiPolygonShape
{
  public:
    explicit ObjectBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vec2d> oringial_shape = createObjectShape();
        std::vector<Vec2d> transformed_shape;
        Rotation2d rotation(initial_rotation_angle);
        for (const auto &oringial_point : oringial_shape)
        {
            Vec2d rotated_point = rotation * (oringial_point - centroid) + centroid;
            // here the cylinder center should be the relativa distance to the global center.
            Vec2d final_point = rotated_point + cylinder_center;
            transformed_shape.push_back(final_point);
        }
        multi_polygon_.addAPolygon(transformed_shape, ShapeBooleanOps::add);
    }
};

class WettingCylinderBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingCylinderBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = cylinder_moisture;
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	The diffusion model of wetting
//----------------------------------------------------------------------
using CylinderFluidDiffusionDirichlet =
    DiffusionRelaxationRK2<DiffusionRelaxation<Dirichlet<KernelGradientContact>, IsotropicDiffusion>>;
//------------------------------------------------------------------------------
// Constrained part for Simbody
//------------------------------------------------------------------------------
MultiPolygon createSimbodyConstrainShape(SPHBody &sph_body)
{
    MultiPolygon multi_polygon;
    Rotation2d rotation(initial_rotation_angle);
    std::vector<Vec2d> oringial_shape = createObjectShape();
    std::vector<Vec2d> transformed_shape;
    for (const auto &oringial_point : oringial_shape)
    {
        Vec2d rotated_point = rotation * (oringial_point - centroid) + centroid;
        // here the cylinder center should be the relativa distance to the global center.
        Vec2d final_point = rotated_point + cylinder_center;
        transformed_shape.push_back(final_point);
    }
    multi_polygon.addAPolygon(transformed_shape, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
//	Set object initial velocity.
//----------------------------------------------------------------------
class CylinderInitialCondition : public LocalDynamics
{
  public:
    CylinderInitialCondition(SPHBody &sph_body) : LocalDynamics(sph_body),
                                                  vel_(particles_->getVariableDataByName<Vecd>("Velocity")) {}

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][0] = initial_speed * cos(initial_angle);
        vel_[index_i][1] = initial_speed * sin(initial_angle);
    }

  protected:
    Vecd *vel_;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{

    //----------------------------------------------------------------------
    // 规则波参数
    //----------------------------------------------------------------------
    //Real H = 0.15;     // 波高 (m)
    //Real T = 1.0;     // 周期 (s)
    //Real phase = 0.0; // 相位 (rad)
    //WaveFormFunc wave_func = createRegularWave(H, T, phase, LH, gravity_g);
    //// 计算初始波面高度和圆柱位置
    //Real cylinder_x = 0.2 * DL;
    //Real omega = 2.0 * Pi / T;
    //Real k = solveDispersionEquation(omega, LH, gravity_g);
    //Real eta0 = 0.5 * H * cos(k * cylinder_x + phase); // 规则波波面
    //Real cylinder_y = eta0 + LH + 0.2;
    //cylinder_center = Vecd(cylinder_x, cylinder_y);

    //----------------------------------------------------------------------
    // 双色波参数
    //----------------------------------------------------------------------
    //Real H1 = 0.0521, T1 = 1.007;
    //Real H2 = 0.0862, T2 = 1.2866;
    //Real delta_phi = Pi / 3.0; // 相位差
    //WaveFormFunc wave_func = createBiChromaticWave(H1, T1, H2, T2, delta_phi, LH, gravity_g);

    //// 计算初始波面高度
    //Real cylinder_x = 0.2 * DL;
    //Real omega1 = 2.0 * Pi / T1, omega2 = 2.0 * Pi / T2;
    //Real k1 = solveDispersionEquation(omega1, LH, gravity_g);
    //Real k2 = solveDispersionEquation(omega2, LH, gravity_g);
    //Real eta0 = 0.5 * H1 * cos(k1 * cylinder_x) + 0.5 * H2 * cos(k2 * cylinder_x + delta_phi);
    //Real cylinder_y = eta0 + LH + 0.3;
    //cylinder_center = Vecd(cylinder_x, cylinder_y);

    //----------------------------------------------------------------------
    // 聚焦波参数
    //----------------------------------------------------------------------
    Real Af = 0.09;       // 谱峰处目标波浪振幅 (m)
    Real fp = 0.83;        // 谱峰频率 (Hz) 能量集中的中心频率，决定波浪周期
    Real bandwidth = 0.6; // 带宽 (Hz)，频率范围 [0.5, 1.1] Hz 频率成分的分布范围，影响波群长度和聚焦程度
    int Nf = 29;          // 离散频率数量（奇数可得到对称谱）
    Real tf = 20; // 聚焦时刻 (s)
    Real xf = 11.54;   // 聚焦位置 (m) - 水槽中央
                     // DL = 22;  DH = 5; LH = 0.5; 
    WaveFormFunc wave_func = createFocusedWave(Af, fp, bandwidth, Nf, tf, xf, LH, gravity_g);
    std::cout << "=== Focusing wave: tf = " << tf << ", xf = " << xf << " m" << std::endl;
     //计算初始波面高度（t=0，x=cylinder_x 处）
    Real cylinder_x = 0.3 * DL;
    Real eta0 = evaluateFocusedWaveElevation(Af, fp, bandwidth, Nf, tf, xf, LH, gravity_g,
                                             cylinder_x, 0.0);
    Real cylinder_y = Af*10 + LH + 1;
    cylinder_center = Vecd(cylinder_x, cylinder_y);



    std::cout << "Cylinder initial position: (" << cylinder_center[0] << ", " << cylinder_center[1] << ")" << std::endl;
    front_center_observer_location = {resetFrontCenterObserverPosition()}; // 更新依赖 cylinder_center 的观测点位置
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------

    BoundingBoxd system_domain_bounds(Vec2d(-cavity_length - BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    //sph_system.setRunParticleRelaxation(true);
    //sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WettingFluidBody>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WettingWallBody>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody cylinder(sph_system, makeShared<ObjectBody>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 4.0);//Multi-resolution
    cylinder.defineBodyLevelSetShape();

    cylinder.defineClosure<Solid, IsotropicDiffusion>(
        rho0_s, ConstructArgs(diffusion_species_name, diffusion_coeff));
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        : cylinder.generateParticles<BaseParticles, Lattice>();

    BaseParticles &cylinder_particles = cylinder.getBaseParticles();
    cylinder_particles.registerStateVariableData<Vecd>("Velocity");
    cylinder_particles.registerStateVariableData<Real>("AngularVelocity");

    // 前段中心观测体
    ObserverBody front_center_observer(sph_system, "FrontCenterObserver");
    front_center_observer.generateParticles<ObserverParticles>(front_center_observer_location);

    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation cylinder_inner(cylinder);
    ContactRelation water_block_contact(water_block, {&wall_boundary, &cylinder});
    ContactRelation cylinder_contact(cylinder, {&water_block});
    ContactRelation wetting_observer_contact(front_center_observer, {&cylinder});
    ContactRelation front_center_observer_contact(front_center_observer, {&cylinder});

    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation cylinder_inner(cylinder);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(cylinder);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(cylinder);
        /** A  Physics relaxation step. */
        RelaxationStepInner relaxation_step_inner(cylinder_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }

    //----------------------------------------------------------------------
    //	Define the fluid dynamics used in the simulation.
    //----------------------------------------------------------------------
    GetDiffusionTimeStepSize get_thermal_time_step(cylinder);
    CylinderFluidDiffusionDirichlet cylinder_wetting(cylinder_contact);
    SimpleDynamics<WettingFluidBodyInitialCondition> wetting_water_initial_condition(water_block);
    SimpleDynamics<WettingWallBodyInitialCondition> wetting_wall_initial_condition(wall_boundary);
    SimpleDynamics<WettingCylinderBodyInitialCondition> wetting_cylinder_initial_condition(cylinder);

    // 创建圆柱初始速度设置的动力学对象

    SimpleDynamics<CylinderInitialCondition> cylinder_set_initial_velocity(cylinder);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);
    InteractionWithUpdate<WettingCoupledSpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_block_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);

    /** Kernel correction matrix and transport velocity formulation. */
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> kernel_correction_complex(DynamicsArgs(water_block_inner, 0.5), water_block_contact);
    // Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact); // with KGC correction

    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_time_step(water_block);

    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(fluid_density_relaxation)>> pressure_force_from_fluid(cylinder_contact);

    // 全局坐标系总粘性力/压力力的统计
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_viscous_force_global(cylinder, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force_global(cylinder, "PressureForceFromFluid");

    ReduceDynamics<QuantitySummation<Vecd>> calculate_cylinder_total_pressure_force(cylinder, "PressureForceFromFluid");
    ReduceDynamics<QuantitySummation<Vecd>> calculate_cylinder_total_viscous_force(cylinder, "ViscousForceFromFluid");

    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //  波浪相关动力学（新增）
    //----------------------------------------------------------------------
    // 定义造波区域：左侧 x∈[-BW, 0] 的固体壁面粒子
    std::vector<Vecd> wavemaker_shape_pnts;
    wavemaker_shape_pnts.push_back(Vecd(-BW, 0));
    wavemaker_shape_pnts.push_back(Vecd(-BW, DH + BW));
    wavemaker_shape_pnts.push_back(Vecd(0.0, DH + BW));
    wavemaker_shape_pnts.push_back(Vecd(0.0, 0));
    wavemaker_shape_pnts.push_back(Vecd(-BW, 0));
    MultiPolygon wavemaker_poly;
    wavemaker_poly.addAPolygon(wavemaker_shape_pnts, ShapeBooleanOps::add);
    BodyRegionByParticle wave_maker(wall_boundary, makeShared<MultiPolygonShape>(wavemaker_poly, "WaveMaker"));
    SimpleDynamics<WaveMaking> wave_making(wave_maker, wave_func);

    // 定义消波区域（流体右侧）
    BodyRegionByCell damping_buffer(water_block, makeShared<MultiPolygonShape>(createDampingBufferShape(), "DampingBuffer"));
    SimpleDynamics<fluid_dynamics::DampingBoundaryCondition> damping_wave(damping_buffer);
    //----------------------------------------------------------------------
    //	Building Simbody.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass properties of the fixed spot. */
    SimTK::Body::Rigid fixed_spot_info(SimTK::MassProperties(1.0, SimTKVec3(0), SimTK::UnitInertia(1)));
    SolidBodyPartForSimbody cylinder_constraint_area(cylinder, makeShared<MultiPolygonShape>(createSimbodyConstrainShape(cylinder), "cylinder"));
    /** Mass properties of the constrained spot. */
    Vecd tethering_point = cylinder_constraint_area.initial_mass_center_; // 使用实际计算的质心
    SimTK::MassProperties cylinder_mass_props(
        cylinder_constraint_area.body_part_mass_properties_->getMass(),       // 保留原质量
        SimTKVec3(tethering_point[0], tethering_point[1], 0.0),               // 强制质心为tethering_point
        cylinder_constraint_area.body_part_mass_properties_->getUnitInertia()); // 保留原转动惯量

    SimTK::Body::Rigid tethered_spot_info(cylinder_mass_props);
    /** Mobility of the fixed spot. */
    SimTK::MobilizedBody::Weld fixed_spot(matter.Ground(), SimTK::Transform(SimTKVec3(tethering_point[0], tethering_point[1], 0.0)),
                                          fixed_spot_info, SimTK::Transform(SimTKVec3(0)));
    /** Mobility of the tethered spot. */
    Vecd displacement0 = cylinder_constraint_area.initial_mass_center_ - tethering_point;

    SimTK::MobilizedBody::Planar tethered_spot(fixed_spot,
                                               SimTK::Transform(SimTKVec3(displacement0[0], displacement0[1], 0.0)),
                                               tethered_spot_info, SimTK::Transform(SimTKVec3(0)));
    //SimTK::MobilizedBody::Planar tethered_spot(matter.Ground(), // connect to ground, not the fixed spot
    //                                           SimTK::Transform(SimTKVec3(displacement0[0], displacement0[1], 0.0)),
    //                                           tethered_spot_info, SimTK::Transform(SimTKVec3(0)));
    // discrete forces acting on the bodies.
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, Real(-9.81), 0.0), 0.0);
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    fixed_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.02));
    tethered_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.4));
    SimTK::State state = MBsystem.realizeTopology();

    state.updQ()[0] = displacement0[0];       // x位移（从系留点到质心）
    state.updQ()[1] = displacement0[1];       // y位移（从系留点到质心）
    state.updQ()[2] = 0.0;

    SimTK::Vec3 mobilizer_vel(0.0, initial_speed * cos(initial_angle), initial_speed * sin(initial_angle)); 
    tethered_spot.setU(state, mobilizer_vel); // 设置初始速度（U）：通过Mobilizer的setU方法

    // 设置完Q/U后，让Simbody重新感知状态
    MBsystem.realize(state, SimTK::Stage::Velocity);
    MBsystem.realize(state, SimTK::Stage::Acceleration);
    MBsystem.realize(state, SimTK::Stage::Dynamics);
     
    /** Time stepping method for multibody system.*/
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH.
    //----------------------------------------------------------------------
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_tethered_spot(cylinder_constraint_area, MBsystem, tethered_spot, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_tethered_spot(cylinder_constraint_area, MBsystem, tethered_spot, integ);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");          // output for debug
    body_states_recording.addToWrite<Real>(water_block, "Density");           // output for debug
    body_states_recording.addToWrite<int>(water_block, "Indicator");          // output for debug
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection"); // output for debug
    RestartIO restart_io(sph_system);

    /** WaveProbes. */
    Real probe_width = 1.3 * particle_spacing_ref;
    BodyRegionByCell wave_probe_1(water_block, makeShared<MultiPolygonShape>(createWaveProbeShape(probe_x1, probe_width, LH, DH), "WaveProbe_1"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>> wave_probe_1_recorder(wave_probe_1, "FreeSurfaceHeight");

    BodyRegionByCell wave_probe_2(water_block, makeShared<MultiPolygonShape>(createWaveProbeShape(probe_x2, probe_width, LH, DH), "WaveProbe_2"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>> wave_probe_2_recorder(wave_probe_2, "FreeSurfaceHeight");

    BodyRegionByCell wave_probe_3(water_block, makeShared<MultiPolygonShape>(createWaveProbeShape(probe_x3, probe_width, LH, DH), "WaveProbe_3"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>> wave_probe_3_recorder(wave_probe_3, "FreeSurfaceHeight");

    ObservedQuantityRecording<Real> write_cylinder_wetting("Phi", wetting_observer_contact);
    ObservedQuantityRecording<Vecd> write_front_center_position("Position", front_center_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    cylinder_normal_direction.exec();
    wetting_water_initial_condition.exec();
    wetting_wall_initial_condition.exec();
    wetting_cylinder_initial_condition.exec();
    Real dt_thermal = get_thermal_time_step.exec();
    free_stream_surface_indicator.exec();
    constant_gravity.exec();
    cylinder_set_initial_velocity.exec(); // 执行圆柱初始速度设置
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100; 
    //int observation_sample_interval = screen_output_interval * 1;
    //int restart_output_interval = screen_output_interval * 2000;
    Real end_time = release_time+0.05;

// 计算释放前和释放后的输出间隔
    Real pre_interval = release_time / pre_output_count;
    Real output_interval_vtp = (end_time - release_time) / post_output_count;
    Real output_interval_force = (end_time - release_time) / post_froce_output_count;

    Real next_vtp_output = pre_interval;     // 释放前第一个 VTP 输出时刻
    Real next_force_output = end_time + 1.0; // 初始时力输出不启用（设为大值）

    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    wave_probe_1_recorder.writeToFile(0);
    wave_probe_2_recorder.writeToFile(0);
    wave_probe_3_recorder.writeToFile(0);
    SummaryOutput summary_output("./output/SummaryOutput.dat"); // 创建自定义的汇总输出对象
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------

    while (physical_time < end_time)
    {
        Real current_vtp_interval = released ? output_interval_vtp : pre_interval;
        Real next_event_time = end_time;
        next_event_time = std::min(next_event_time, next_vtp_output);
        if (released)
            next_event_time = std::min(next_event_time, next_force_output);
        Real target_time = std::max(0.0, next_event_time - physical_time);
        Real integration_time = 0.0;
        while (integration_time < target_time && physical_time < end_time)
        {
            time_instance = TickCount::now();
            Real Dt = fluid_advection_time_step.exec()*0.5;//
            fluid_density_by_summation.exec();
            viscous_force.exec();
            kernel_correction_complex.exec(); // with KGC correction
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real dt = 0.0;
            viscous_force_from_fluid.exec(); 

            while (relaxation_time < Dt && physical_time < end_time)
            {
                dt = SMIN(SMIN(dt_thermal, fluid_acoustic_time_step.exec()), Dt);
                fluid_pressure_relaxation.exec(dt);
                pressure_force_from_fluid.exec();
                fluid_density_relaxation.exec(dt);
                cylinder_wetting.exec(dt);
                wave_making.exec(dt);
                wall_boundary.updateCellLinkedList();
                water_block_contact.updateConfiguration();

                if (released)
                {
                    integ.stepBy(dt);
                    SimTK::State &state_for_update = integ.updAdvancedState();
                    force_on_bodies.clearAllBodyForces(state_for_update);
                    force_on_bodies.setOneBodyForce(state_for_update, tethered_spot, force_on_tethered_spot.exec());
                    constraint_tethered_spot.exec();
                }

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << "N=" << number_of_iterations << "  Time = " << physical_time
                          << "  Dt = " << Dt << "  dt = " << dt
                          << "  V_wavemaker = " << WaveMaking::getCurrentVelocity() << "\n";
            }
            number_of_iterations++;

            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                particle_sorting.exec();
            water_block.updateCellLinkedList();
            cylinder.updateCellLinkedList();
            water_block_inner.updateConfiguration();
            cylinder_inner.updateConfiguration();
            cylinder_contact.updateConfiguration();
            water_block_complex.updateConfiguration();
            front_center_observer_contact.updateConfiguration();
            free_stream_surface_indicator.exec();
            damping_wave.exec(Dt);
            interval_updating_configuration += TickCount::now() - time_instance;
        }

        if (physical_time >= next_vtp_output )
        {
            TickCount t4 = TickCount::now();
            body_states_recording.writeToFile();
            wave_probe_1_recorder.writeToFile();
            wave_probe_2_recorder.writeToFile();
            wave_probe_3_recorder.writeToFile();
            TickCount t5 = TickCount::now();
            interval += t5 - t4; 
            next_vtp_output += current_vtp_interval;
            if (next_vtp_output <= physical_time)
                next_vtp_output = physical_time + current_vtp_interval;
        }

        if (released && physical_time >= next_force_output)
        {
            TickCount t2 = TickCount::now();

            viscous_force_from_fluid.exec();
            pressure_force_from_fluid.exec();
            write_total_viscous_force_global.writeToFile(number_of_iterations);
            write_total_pressure_force_global.writeToFile(number_of_iterations);

            Vec2d total_viscous_force_g = calculate_cylinder_total_viscous_force.exec();
            Vec2d total_pressure_force_g = calculate_cylinder_total_pressure_force.exec();
            Real current_rotation = getCylinderRotationAngle(tethered_spot, integ.getAdvancedState());
            Real total_rotation_angle = initial_rotation_angle + current_rotation;
            Vec2d viscous_local = transformGlobalForceToLocal(total_viscous_force_g, total_rotation_angle);
            Vec2d pressure_local = transformGlobalForceToLocal(total_pressure_force_g, total_rotation_angle);

            auto &front_center_particles = front_center_observer.getBaseParticles();
            Vecd front_center_pos = front_center_particles.getVariableDataByName<Vecd>("Position")[0];
            summary_output.writeData(physical_time,
                                     total_viscous_force_g, total_pressure_force_g,
                                     viscous_local, pressure_local, front_center_pos);

            TickCount t3 = TickCount::now();
            interval += t3 - t2; 
            next_force_output += output_interval_force;
            if (next_force_output <= physical_time)
                next_force_output = physical_time + output_interval_force;
        }


        if (!released && physical_time >= release_time)
        {
            released = true;
            next_vtp_output = physical_time + output_interval_vtp;
            next_force_output = physical_time + output_interval_force;
            SimTK::State &state_for_update = integ.updAdvancedState();
            state_for_update.updU()[0] = initial_angular_velocity;
            state_for_update.updU()[1] = initial_speed * cos(initial_angle);
            state_for_update.updU()[2] = initial_speed * sin(initial_angle);
            MBsystem.realize(state_for_update, SimTK::Stage::Velocity);
        }
    }

        TickCount t6 = TickCount::now();
       
    TimeInterval tt;
    tt = t6 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
};



