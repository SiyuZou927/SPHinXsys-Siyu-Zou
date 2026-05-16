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
Real DL = 4;                     /**< Water tank length. */
Real DH = 5;                     /**< Water tank height. */
Real LH = 2;                     /**< Water column height. */
Real particle_spacing_ref = 0.02; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;                    /**< Thickness of tank wall. */

// 波浪控制参数
Real wave_height = 0.2;                  /**< 目标波高 (m) */
Real wave_period = 1.0;                   /**< 波浪周期 (s) */
Real wave_phase = 0.0;                    /**< 初始相位 (rad)，用于控制圆柱入水位置 */
Real wave_k = 0.0;                        /**< 波数（将自动计算） */
Real wave_stroke = 0.0;                   /**< 推板冲程（将自动计算） */
Real wave_omega = 2.0 * Pi / wave_period; /**< 角频率 */
Real release_time = 5.0;                  // 造波时长 间接控制结束时刻，建议设置为波浪周期的整数倍以获得完整波形
bool released = false; // 是否已释放
Real cavity_length = 0.5; // 左侧空腔长度，应大于最大推板位移

int pre_output_count = 50;  // 释放前 VTP 输出次数（不含 t=0）
int post_output_count = 20; // 释放后 VTP 输出次数
int post_froce_output_count = 500; // 释放后力输出次数

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
Real wall_moisture = 1.0;                                    /**< wall moisture. */


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
    Real angle_from_rot = angles[2]; // 绕Z轴的旋转
    return angle_from_rot; // 使用旋转矩阵的角度
}

//----------------------------------------------------------------------
// 波浪相关辅助函数（新增）
//----------------------------------------------------------------------
// 求解线性波色散方程：omega^2 = g * k * tanh(k * h)
Real solveDispersionEquation(Real omega, Real h, Real g, Real tol = 1e-6)
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
//----------------------------------------------------------------------
// 造波板动力学类（新增）
//----------------------------------------------------------------------
class WaveMaking : public BodyPartMotionConstraint
{
    Real amplitude_;      // 冲程的一半（0.5 * stroke）
    Real omega_;          // 角频率
    Real phase_;          // 初始相位（可设为零）
    Real *physical_time_; // 指向物理时间的指针
                          // Vecd *acc_;  // 若需要记录加速度可取消注释，并确保粒子注册了"Acceleration"

  public:
    WaveMaking(BodyPartByParticle &body_part)
        : BodyPartMotionConstraint(body_part),
          amplitude_(0.5 * wave_stroke),
          omega_(wave_omega),
          phase_(0.0),
          physical_time_(body_part.getSPHBody().getSPHSystem().getSystemVariableDataByName<Real>("PhysicalTime"))
    {
        // 如果需要加速度变量，取消下面注释，并确保粒子注册了"Acceleration"
        // acc_ = particles_->template getVariableDataByName<Vecd>("Acceleration");
    }

    virtual void update(size_t index_i, Real dt = 0.0)
    {
        Real time = *physical_time_;
        Real displacement = amplitude_ * sin(omega_ * time + phase_);
        pos_[index_i] = pos0_[index_i] + Vecd(displacement, 0.0);
        vel_[index_i] = Vecd(amplitude_ * omega_ * cos(omega_ * time + phase_), 0.0);
        // 如果不需要加速度，可以注释掉下面一行
        // acc_[index_i] = Vecd(-amplitude_ * omega_ * omega_ * sin(omega_ * time + phase_), 0.0);
    }
};
//----------------------------------------------------------------------
// 消波区域形状（新增）
//----------------------------------------------------------------------
MultiPolygon createDampingBufferShape()
{
    std::vector<Vecd> pnts;
    Real damping_start = DL - 0.5; // 从距离右端 0.5m 处开始阻尼
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
// 观测点位置（原有）
//----------------------------------------------------------------------

// 新增：前段中心观测点初始位置（注意 cylinder_center 尚未计算，将在 main 中重新赋值）
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
    inner_wall.push_back(Vecd(0.0, DH));
    inner_wall.push_back(Vecd(DL, DH));
    inner_wall.push_back(Vecd(DL, 0.0));
    inner_wall.push_back(Vecd(0.0, 0.0));

    return inner_wall;
}
std::vector<Vecd> createCavityInnerShape()
{
    std::vector<Vecd> inner;
    inner.push_back(Vecd(-cavity_length, 0.0)); 
    inner.push_back(Vecd(-cavity_length, DH));  
    inner.push_back(Vecd(-BW, DH));             
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
    //  1. 波浪参数预计算（新增）
    //----------------------------------------------------------------------
    wave_omega = 2.0 * Pi / wave_period;
    wave_k = solveDispersionEquation(wave_omega, LH, gravity_g);
    wave_stroke = computePistonStroke(wave_height, wave_k, LH);
    std::cout << "Wave number k = " << wave_k << ", stroke = " << wave_stroke << std::endl;

    // 根据波浪相位确定圆柱初始位置
    Real cylinder_x = 0.5 * DL;                                        
    Real eta0 = 0.5 * wave_height * cos(wave_k * cylinder_x + wave_phase); // t=0时的波面高度
    Real cylinder_y = eta0 + LH + 0.2; // 设置圆柱中心 y 坐标
    cylinder_center = Vecd(cylinder_x, cylinder_y);
    std::cout << "Cylinder initial position: (" << cylinder_center[0] << ", " << cylinder_center[1] << ")" << std::endl;

    // 更新依赖 cylinder_center 的观测点位置

    front_center_observer_location = {resetFrontCenterObserverPosition()};
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

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);
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
    // 插值更新前段中心观测点位置
    //InteractionDynamics<InterpolatingAQuantity<Vecd>>
    //    interpolation_front_center_position(front_center_observer_contact, "Position", "Position");
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
    SimpleDynamics<WaveMaking> wave_making(wave_maker);

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
    // SimTK::MobilizedBody::Planar tethered_spot(matter.Ground(), // connect to ground, not the fixed spot
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
    //state.updQ()[2] = initial_rotation_angle; // 相对于父体的旋转
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
    int screen_output_interval = 1; 
    int observation_sample_interval = screen_output_interval * 1;
    int restart_output_interval = screen_output_interval * 2000;
    Real end_time = release_time+0.02;

// 计算释放前和释放后的输出间隔
    Real pre_interval = release_time / pre_output_count;

    Real output_interval_vtp = (end_time - release_time) / post_output_count;
    Real output_interval_force = (end_time - release_time) / post_froce_output_count;
    //Real next_force_output = output_interval_force;
    //Real next_vtp_output = output_interval_vtp;
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
    //write_cylinder_displacement.writeToFile(number_of_iterations);
    //write_cylinder_wetting.writeToFile(number_of_iterations);
    //write_front_center_position.writeToFile(number_of_iterations);
    SummaryOutput summary_output("./output/SummaryOutput.dat"); // 创建自定义的汇总输出对象
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------

        /** Integrate time (loop) until the next output time. */
        //----------------------------------------------------------------------
        //	Main loop starts here.
        //----------------------------------------------------------------------
        // ========== 第一阶段：释放前（仅输出 VTP） ==========
        while (physical_time < release_time) // 微小量避免浮点误差
        {
            Real integration_time = 0.0;
            while (integration_time < pre_interval && physical_time < release_time )
            {
             
                time_instance = TickCount::now();
                Real Dt = fluid_advection_time_step.exec();

                fluid_density_by_summation.exec();
                viscous_force.exec();
                transport_velocity_correction.exec();
                interval_computing_time_step += TickCount::now() - time_instance;

                time_instance = TickCount::now();
                Real relaxation_time = 0.0;
                Real dt = 0.0;
                viscous_force_from_fluid.exec();

                while (relaxation_time < Dt)
                {
                    dt = SMIN(SMIN(dt_thermal, fluid_acoustic_time_step.exec()), Dt);
                    fluid_pressure_relaxation.exec(dt);
                    pressure_force_from_fluid.exec();
                    fluid_density_relaxation.exec(dt);
                    cylinder_wetting.exec(dt);

                    wave_making.exec(dt);
                    wall_boundary.updateCellLinkedList();
                    water_block_contact.updateConfiguration();

                    integ.stepBy(dt);

                    if (physical_time < release_time)
                    {
                        SimTK::State &state = integ.updAdvancedState();
                        state.updQ()[0] = displacement0[0];
                        state.updQ()[1] = displacement0[1];
                        state.updQ()[2] = 0.0;
                        state.updU()[0] = 0.0;
                        state.updU()[1] = 0.0;
                        state.updU()[2] = 0.0;
                        MBsystem.realize(state, SimTK::Stage::Velocity);
                    }

                    SimTK::State &state_for_update = integ.updAdvancedState();
                    force_on_bodies.clearAllBodyForces(state_for_update);
                    if (physical_time >= release_time)
                    {
                        force_on_bodies.setOneBodyForce(state_for_update, tethered_spot, force_on_tethered_spot.exec());
                    }
                    constraint_tethered_spot.exec();
                    //interpolation_front_center_position.exec();

                    relaxation_time += dt;
                    integration_time += dt;
                    physical_time += dt;
                }
                interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

                // 屏幕输出（仅时间信息，无数据文件）
                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << physical_time << "	Dt = " << Dt << "	dt = " << dt << "\n";
                }
                number_of_iterations++;

                // 更新链表和配置
                time_instance = TickCount::now();
                if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                {
                    particle_sorting.exec();
                }
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
            // 每个外层步结束后输出一次 VTP（释放前）
            body_states_recording.writeToFile();
        }
        // ========== 第二阶段：释放后（同时输出 VTP 和数据文件） ==========

        Real next_force_output = physical_time + output_interval_force;
        Real next_vtp_output = physical_time + output_interval_vtp;
        while (physical_time < end_time)
        {
            Real next_output_time = std::min(next_force_output, next_vtp_output);
            Real time_to_output = next_output_time - physical_time;
            if (time_to_output < 0)
                time_to_output = 0.0;
            Real target_time = std::min(time_to_output, end_time - physical_time);
            Real integration_time = 0.0;
            while (integration_time < target_time)
            {

                time_instance = TickCount::now();
                Real Dt = fluid_advection_time_step.exec();

                fluid_density_by_summation.exec();
                viscous_force.exec();
                transport_velocity_correction.exec();
                interval_computing_time_step += TickCount::now() - time_instance;

                time_instance = TickCount::now();
                Real relaxation_time = 0.0;
                Real dt = 0.0;
                viscous_force_from_fluid.exec();

                while (relaxation_time < Dt)
                {
                    dt = SMIN(SMIN(dt_thermal, fluid_acoustic_time_step.exec()), Dt);
                    fluid_pressure_relaxation.exec(dt);
                    pressure_force_from_fluid.exec();
                    fluid_density_relaxation.exec(dt);
                    cylinder_wetting.exec(dt);
                    wave_making.exec(dt);
                    wall_boundary.updateCellLinkedList();
                    water_block_contact.updateConfiguration();

                    integ.stepBy(dt);
                    SimTK::State &state_for_update = integ.updAdvancedState();
                    // 释放瞬间设置速度（仅一次）
                    if (!released && physical_time >= release_time)
                    {
                        state_for_update.updU()[0] = initial_angular_velocity;
                        state_for_update.updU()[1] = initial_speed * cos(initial_angle);
                        state_for_update.updU()[2] = initial_speed * sin(initial_angle);
                        MBsystem.realize(state_for_update, SimTK::Stage::Velocity);
                        released = true;
                    }

                    force_on_bodies.clearAllBodyForces(state_for_update);
                    if (physical_time >= release_time)
                    {
                        force_on_bodies.setOneBodyForce(state_for_update, tethered_spot, force_on_tethered_spot.exec());
                    }
                    constraint_tethered_spot.exec();

                    // interpolation_front_center_position.exec();

                    relaxation_time += dt;
                    integration_time += dt;
                    physical_time += dt;
                }
                interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;


                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << physical_time << "	Dt = " << Dt << "	dt = " << dt << "\n";
                    if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                    {
                        write_cylinder_wetting.writeToFile(number_of_iterations);
                        write_front_center_position.writeToFile(number_of_iterations);
                    }
                    if (number_of_iterations % restart_output_interval == 0)
                        restart_io.writeToFile(number_of_iterations);
                }
                number_of_iterations++;

                // 更新链表和配置（与第一阶段相同）
                time_instance = TickCount::now();
                if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                {
                    particle_sorting.exec();
                }
                water_block.updateCellLinkedList();
                cylinder.updateCellLinkedList();
                water_block_inner.updateConfiguration();
                cylinder_inner.updateConfiguration();
                cylinder_contact.updateConfiguration();
                water_block_complex.updateConfiguration();
                front_center_observer_contact.updateConfiguration();
                free_stream_surface_indicator.exec();
                interval_updating_configuration += TickCount::now() - time_instance;
                damping_wave.exec(Dt);
            }

            if (physical_time >= next_force_output)
            {
                TickCount t2 = TickCount::now();

                viscous_force_from_fluid.exec();
                pressure_force_from_fluid.exec();
                SimTK::State &output_state = integ.updAdvancedState();

                write_total_viscous_force_global.writeToFile(number_of_iterations); // 记录结果到文件
                write_total_pressure_force_global.writeToFile(number_of_iterations);

                Vec2d total_viscous_force_g = calculate_cylinder_total_viscous_force.exec();
                Vec2d total_pressure_force_g = calculate_cylinder_total_pressure_force.exec();
                Real current_rotation_from_matrix = getCylinderRotationAngle(tethered_spot, integ.getAdvancedState());
                Real total_rotation_angle = initial_rotation_angle + current_rotation_from_matrix;
                Vec2d viscous_local = transformGlobalForceToLocal(total_viscous_force_g, total_rotation_angle);
                Vec2d pressure_local = transformGlobalForceToLocal(total_pressure_force_g, total_rotation_angle);

                //  获取前段中心观测点位置
                auto &front_center_particles = front_center_observer.getBaseParticles();
                Vecd front_center_pos = front_center_particles.getVariableDataByName<Vecd>("Position")[0];

                summary_output.writeData(physical_time,
                                         total_viscous_force_g,  // 全局粘性力
                                         total_pressure_force_g, // 全局压力力
                                         viscous_local,          // 局部压力力X
                                         pressure_local,         // 新增局部力Y
                                         front_center_pos);
                TickCount t3 = TickCount::now();
                interval += t3 - t2;
                next_force_output += output_interval_force;
            }
            if (physical_time >= next_vtp_output)
            {
                TickCount t4 = TickCount::now();

                body_states_recording.writeToFile();

                TickCount t5 = TickCount::now();
                interval += t5 - t4;
                next_vtp_output += output_interval_vtp;
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



