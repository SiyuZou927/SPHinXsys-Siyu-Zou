/**
 * @file    water_entry_exit_3d.cpp
 * @brief   3D water entry and exit example with surface wetting considered.
 * @details This is a 3D FSI test case using M1.stl for spatial temporal identification approach.
 * @author  Siyu Zou Shuoguo Zhang and Xiangyu Hu
 */
#include "SimTKsimbody.h"
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

//----------------------------------------------------------------------
//  Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.2; /**< Water tank length. */
Real DW = 0.4; /**< Water tank width. */
Real DH = 1.0; /**< Water tank height. */
Real LH = 0.5; /**< Water column height. */

Real particle_spacing_ref = 0.0025;   /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
Vecd object_center(0.46, LH + 0.028, DW / 2); /**< Location of the object center. */
// 初始速度参数
Real initial_speed = 70.0;                     /**< Initial velocity magnitude (m/s). */
Real initial_angle = -20.0 * Pi / 180.0;       /**< Initial velocity angle (radians, negative = downward). */
Real initial_pitch_angle = -20.0 * Pi / 180.0; /**< Initial pitch angle about Z-axis(radians). */
std::string full_path_to_file = "./input/M1TR4Z02H05.stl"; // model file path
//----------------------------------------------------------------------
//  Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.00;                         /**< Fluid density (water). */
Real rho0_s = 2.063;                         /**< Object density. */
Real gravity_g = 9.81;                        /**< Gravity. */
Real U_max =10.95; /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                      /**< Reference sound speed. */
Real mu_f = 8.9e-7;                           /**< Water dynamics viscosity (Pa·s). */
//----------------------------------------------------------------------
//  Wetting parameters
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 75.0 * pow(particle_spacing_ref, 2); /**< Wetting coefficient. */
Real fluid_moisture = 1.0;                                   /**< fluid moisture. */
Real object_moisture = 0.0;                                  /**< object moisture. */
Real wall_moisture = 1.0;                                    /**< wall moisture. */
//----------------------------------------------------------------------
//  全局力到局部坐标系的转换
//----------------------------------------------------------------------
Vec3d transformGlobalForceToLocal(const Vec3d &global_force, Real rotation_angle_z)
{
    Real cos_theta = cos(rotation_angle_z);
    Real sin_theta = sin(rotation_angle_z);
    Vec3d local_force;
    local_force[0] = global_force[0] * cos_theta + global_force[1] * sin_theta;  // X方向
    local_force[1] = -global_force[0] * sin_theta + global_force[1] * cos_theta; // Y方向
    local_force[2] = global_force[2];
    return local_force;
}
//----------------------------------------------------------------------
//   观测点定义
//----------------------------------------------------------------------
StdVec<Vecd> leading_edge_location = {object_center}; /**< Leading edge observation point. */
//----------------------------------------------------------------------
//  汇总输出类（3D版本）
//----------------------------------------------------------------------
class SummaryOutput3D
{
  public:
    SummaryOutput3D(const std::string &filename)
        : output_file_(filename)
    {
        output_file_ << "Time[s]   "
                     << "LeadingEdge_X LeadingEdge_Y LeadingEdge_Z   "
                     << "ViscousForceGlobal_X ViscousForceGlobal_Y ViscousForceGlobal_Z   "
                     << "PressureForceGlobal_X PressureForceGlobal_Y PressureForceGlobal_Z   "
                     << "ViscousForceLocal_X ViscousForceLocal_Y ViscousForceLocal_Z   "
                     << "PressureForceLocal_X PressureForceLocal_Y PressureForceLocal_Z   "
                     << "TotalForceLocal_X TotalForceLocal_Y TotalForceLocal_Z\n"; 
    }

    void writeData(Real time,
                   const Vecd &leading_edge_position,
                   const Vecd &viscous_force_global,
                   const Vecd &pressure_force_global,
                   const Vecd &viscous_force_local,
                   const Vecd &pressure_force_local,
                   const Vecd &total_force_local) 
    {
        output_file_ << std::scientific << std::setprecision(9)
                     << time << " "
                     << leading_edge_position[0] << " " << leading_edge_position[1] << " " << leading_edge_position[2] << " "
                     << viscous_force_global[0] << " " << viscous_force_global[1] << " " << viscous_force_global[2] << " "
                     << pressure_force_global[0] << " " << pressure_force_global[1] << " " << pressure_force_global[2] << " "
                     << viscous_force_local[0] << " " << viscous_force_local[1] << " " << viscous_force_local[2] << " "
                     << pressure_force_local[0] << " " << pressure_force_local[1] << " " << pressure_force_local[2] << " "
                     << total_force_local[0] << " " << total_force_local[1] << " " << total_force_local[2] << "\n"; // 新增
        output_file_.flush();
    }
  private:
    std::ofstream output_file_;
};
//----------------------------------------------------------------------
//  Definition for water body 
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(Vecd(DL / 2, LH / 2, DW / 2)), Vecd(DL / 2, LH / 2, DW / 2));// 创建水体长方体
    }
};

//----------------------------------------------------------------------
//  Definition for wall body
//----------------------------------------------------------------------
 class WallBoundary : public ComplexShape {
public:
explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
{
    add<GeometricShapeBox>(Transform(Vecd(DL / 2, -BW / 2, DW / 2)), Vecd(DL / 2 + BW, BW / 2, DW / 2 + BW));       // bottom
    //add<GeometricShapeBox>(Transform(Vecd(-BW / 2, DH / 2, DW / 2)), Vecd(BW / 2, DH / 2 + BW, DW / 2 + BW));       // left
    //add<GeometricShapeBox>(Transform(Vecd(DL + BW / 2, DH / 2, DW / 2)), Vecd(BW / 2, DH / 2 + BW, DW / 2 + BW));   // right
    //add<GeometricShapeBox>(Transform(Vecd(DL / 2, DH / 2, -BW / 2)), Vecd(DL / 2 + BW, DH / 2 + BW, BW / 2));       // front
    //add<GeometricShapeBox>(Transform(Vecd(DL / 2, DH / 2, DW + BW / 2)), Vecd(DL / 2 + BW, DH / 2 + BW, BW / 2));   // back
    // add<GeometricShapeBox>(Transform(Vecd(DL/2, DH + BW/2, DW/2)),Vecd(DL/2 + BW, BW/2, DW/2 + BW));//top
}};
 //----------------------------------------------------------------------
 //  Definition for object body using STL
 //----------------------------------------------------------------------
 class FloatingObject : public ComplexShape
 {
   public:
     explicit FloatingObject(const std::string &shape_name) : ComplexShape(shape_name)
     {
         Vecd translation(0.0, 0.0, 0.0);
         add<TriangleMeshShapeSTL>(full_path_to_file, translation, 1.0);
     }
 };
 //----------------------------------------------------------------------
 // 创建Simbody约束区域的形状
 //----------------------------------------------------------------------
 std::shared_ptr<ComplexShape> createSimbodyConstrainShape3D()
 {
     auto constrain_shape = makeShared<ComplexShape>("SimbodyConstrainShape");
     Vecd translation(0.0, 0.0, 0.0);
     constrain_shape->add<TriangleMeshShapeSTL>(full_path_to_file, translation, 1.0);
     return constrain_shape;
 }
 //----------------------------------------------------------------------
 //  wetting body initial condition
 //----------------------------------------------------------------------
 template <class BodyType>
 class WettingBodyInitialCondition : public LocalDynamics
 {
   public:
     explicit WettingBodyInitialCondition(BodyType &sph_body)
         : LocalDynamics(sph_body),
           pos_(particles_->getVariableDataByName<Vecd>("Position")),
           phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)){};

     void update(size_t index_i, Real dt)
     {
         phi_[index_i] = fluid_moisture;
     }

   protected:
     Real *phi_;
     Vecd *pos_;
 };
 class WettingObjectBodyInitialCondition : public LocalDynamics
 {
   public:
     explicit WettingObjectBodyInitialCondition(SPHBody &sph_body)
         : LocalDynamics(sph_body),
           pos_(particles_->getVariableDataByName<Vecd>("Position")),
           phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {}
     void update(size_t index_i, Real dt)
     {
         phi_[index_i] = object_moisture; // 0.0
     }

   protected:
     Real *phi_;
     Vecd *pos_;
 };
 //// 为不同类型创建别名
using WettingFluidBodyInitialCondition = WettingBodyInitialCondition<FluidBody>;
using WettingWallBodyInitialCondition = WettingBodyInitialCondition<SolidBody>;
 //----------------------------------------------------------------------
 //	The diffusion model of wetting
 //----------------------------------------------------------------------
 using ObjectFluidDiffusionDirichlet =
     DiffusionRelaxationRK2<DiffusionRelaxation<Dirichlet<KernelGradientContact>, IsotropicDiffusion>>;

//----------------------------------------------------------------------
// 获取物体的旋转角度
//----------------------------------------------------------------------
Real getObjectRotationAngle(SimTK::MobilizedBody &mob, const SimTK::State &state)
{
    SimTK::Rotation rot = mob.getBodyRotation(state);
    SimTK::Vec3 angles = rot.convertRotationToBodyFixedXYZ();
    return angles[2]; // Z_axis rotation angle
}
//----------------------------------------------------------------------
// Set object initial velocity.
//----------------------------------------------------------------------
class ObjectInitialVelocity : public LocalDynamics
{
  public:
        ObjectInitialVelocity(SPHBody &sph_body) : LocalDynamics(sph_body),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")) {}
    void update(size_t index_i, Real dt)
    {
        vel_[index_i][0] = initial_speed * cos(initial_angle);
        vel_[index_i][1] = initial_speed * sin(initial_angle);
        vel_[index_i][2] = 0.0; // No initial velocity in Z direction
    }
  protected:
    Vecd *vel_;
};

//----------------------------------------------------------------------
//  Main program
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //  Build up an SPHSystem
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);

    // 运行配置
    //sph_system.setRunParticleRelaxation(false);
    //sph_system.setReloadParticles(true);
    sph_system.setRunParticleRelaxation(true);
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);

    //----------------------------------------------------------------------
    //  Creating bodies with corresponding materials and particles
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(
        ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody floating_object(sph_system, makeShared<FloatingObject>("FloatingObject"));
    floating_object.defineAdaptationRatios(1.15, 2.0);//Multi-resolution
    floating_object.defineBodyLevelSetShape();

    // 定义材料，包含湿表面扩散
    floating_object.defineClosure<Solid, IsotropicDiffusion>(
        rho0_s, ConstructArgs(diffusion_species_name, diffusion_coeff));
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? floating_object.generateParticles<BaseParticles, Reload>(floating_object.getName())
        : floating_object.generateParticles<BaseParticles, Lattice>();

    // 注册额外状态变量
    BaseParticles &object_particles = floating_object.getBaseParticles();

    object_particles.registerStateVariableData<Vecd>("Velocity");// 


    // 观测体
    ObserverBody leading_edge_observer(sph_system, "LeadingEdgeObserver");
    leading_edge_observer.generateParticles<ObserverParticles>(leading_edge_location);

    //----------------------------------------------------------------------
    //  Define body relation map
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation object_inner(floating_object);
    ContactRelation water_block_contact(water_block, {&wall_boundary, &floating_object});
    ContactRelation object_contact(floating_object, {&water_block});
    ContactRelation object_observer_contact(leading_edge_observer, {&floating_object});
    ContactRelation wetting_observer_contact(leading_edge_observer, {&floating_object});
    ContactRelation leading_edge_observer_contact(leading_edge_observer, {&floating_object});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);

     //----------------------------------------------------------------------
    // Run particle relaxation
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        std::cout << "Starting 3D particle relaxation..." << std::endl;
        InnerRelation object_inner(floating_object); //
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(floating_object);
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(floating_object); /** Write the body state to Vtp file. */
        ReloadParticleIO write_particle_reload_files(floating_object);        /** Write the particle reload files. */
        RelaxationStepInner relaxation_step_inner(object_inner);              /** A  Physics relaxation step. */
        //----------------------------------------------------------------------
        // Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0); // output inicial state
        //----------------------------------------------------------------------
        // Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        int max_iterations = 2000;

        while (ite_p < max_iterations)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "Relaxation steps for the 3D object N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "3D physics relaxation process finished!" << std::endl;

        write_particle_reload_files.writeToFile(0); // Output results
        return 0;
    }
    //----------------------------------------------------------------------
    //  Define the fluid dynamics
    //----------------------------------------------------------------------
    // 湿表面扩散模型
    GetDiffusionTimeStepSize get_diffusion_time_step(floating_object);
    ObjectFluidDiffusionDirichlet object_wetting(object_contact);

    // 初始化湿表面条件
    SimpleDynamics<WettingBodyInitialCondition<FluidBody>> wetting_water_initial_condition(water_block);
    SimpleDynamics<WettingBodyInitialCondition<SolidBody>>wetting_wall_initial_condition(wall_boundary);
    SimpleDynamics<WettingObjectBodyInitialCondition> wetting_object_initial_condition(floating_object);

    SimpleDynamics<ObjectInitialVelocity> object_set_initial_velocity(floating_object); // set object initial velocity

    Gravity gravity(Vecd(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);

    // 自由表面识别（使用带湿表面耦合的版本）
    InteractionWithUpdate<WettingCoupledSpatialTemporalFreeSurfaceIndicationComplex>
        free_stream_surface_indicator(water_block_inner, water_block_contact);


    // 法向计算
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> object_normal_direction(floating_object);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex>
        kernel_correction_complex(DynamicsArgs(water_block_inner, 0.5), water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);
    // 流体动力学
    //Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);

    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall>viscous_force(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>transport_velocity_correction(water_block_inner, water_block_contact);
    // periodic condition
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicAlongAxis periodic_along_z(water_block.getSPHBodyBounds(), zAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_block, periodic_along_x);
    PeriodicConditionUsingCellLinkedList periodic_condition_z(water_block, periodic_along_z);
    // 时间步控制
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep>fluid_acoustic_time_step(water_block);
    //----------------------------------------------------------------------
    //  Algorithms of FSI
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(object_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(fluid_density_relaxation)>>pressure_force_from_fluid(object_contact);

    // 力统计（全局坐标系）
    ReducedQuantityRecording<QuantitySummation<Vecd>>write_total_viscous_force(floating_object, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>>write_total_pressure_force(floating_object, "PressureForceFromFluid");

    ReduceDynamics<QuantitySummation<Vecd>> calculate_cylinder_total_pressure_force(floating_object, "PressureForceFromFluid");
    ReduceDynamics<QuantitySummation<Vecd>> calculate_cylinder_total_viscous_force(floating_object, "ViscousForceFromFluid");
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //  Building Simbody (3D)
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;                /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem); /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);  /** Mass properties of the fixed spot. */

    //SimTK::Body::Rigid fixed_spot_info(*structure_system.body_part_mass_properties_);
    SimTK::Body::Rigid fixed_spot_info(SimTK::MassProperties(1.0, SimTKVec3(0), SimTK::UnitInertia(1)));
    // 创建约束区域（使用与物体相同的形状）
    SolidBodyPartForSimbody structure_system(floating_object, createSimbodyConstrainShape3D());
    /** Mass properties of the constrained spot. */
    // SimTK::Body::Rigid structure_info(*structure_system.body_part_mass_properties_); // 质量属性
    Vecd tethering_point = structure_system.initial_mass_center_;       // 系留点使用实际质心
    SimTK::MassProperties object_mass_props(                            // set tethered spot mass properties
        structure_system.body_part_mass_properties_->getMass(),         // 保留原质量
        SimTKVec3(tethering_point[0], tethering_point[1], 0.0),         // 强制质心为tethering_point
        structure_system.body_part_mass_properties_->getUnitInertia()); // 保留原转动惯量
    SimTK::Body::Rigid tethered_spot_info(object_mass_props);    
        /** Mobility of the fixed spot. */
    SimTK::MobilizedBody::Weld fixed_spot(matter.Ground(), SimTK::Transform(SimTKVec3(tethering_point[0], tethering_point[1], 0.0)),
                                          fixed_spot_info, SimTK::Transform(SimTKVec3(0)));

    Vecd displacement0 = structure_system.initial_mass_center_ - tethering_point;
    SimTK::MobilizedBody::Planar structure_mob(matter.Ground(),
                                              SimTK::Transform(SimTKVec3(displacement0[0], displacement0[1], displacement0[2])),
                                              tethered_spot_info, SimTK::Transform(SimTKVec3(0)));
   //SimTK::MobilizedBody::Planar structure_mob(fixed_spot, // <--- 父级改为 fixed_spot
   //                                           SimTK::Transform(SimTKVec3(displacement0[0], displacement0[1], displacement0[2])),
   //                                           tethered_spot_info, SimTK::Transform(SimTKVec3(0)));

    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTKVec3(0.0, -gravity_g, 0.0), 0.0); // 重力
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);                                   // 离散力
    fixed_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.02));               // 可视化固定点
    tethered_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.4));             // 可视化约束区域
    SimTK::State state = MBsystem.realizeTopology(); // 初始化状态

    SimTK::Vec3 mobilizer_vel(0.0, initial_speed * cos(initial_angle), initial_speed * sin(initial_angle));
    structure_mob.setU(state, mobilizer_vel);

    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);

    // 输出Simbody状态
    std::cout << "Debug: Joint U = ("
              << structure_mob.getU(state)[0] << ", "
              << structure_mob.getU(state)[1] << ", "
              << structure_mob.getU(state)[2] << ")" << std::endl;
    std::cout << "Mass: " << structure_system.body_part_mass_properties_->getMass() << std::endl;
    std::cout << "Center of mass: " << structure_system.initial_mass_center_.transpose() << std::endl;
    //----------------------------------------------------------------------
    //  Coupling between SimBody and SPH
    //----------------------------------------------------------------------
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody> force_on_structure(structure_system, MBsystem, structure_mob, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_on_structure(structure_system, MBsystem, structure_mob, integ);
    //----------------------------------------------------------------------
    //  I/O operations and observations
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_states_recording.addToWrite<Real>(floating_object, diffusion_species_name);
    RestartIO restart_io(sph_system);
    ObservedQuantityRecording<Real> write_object_wetting("Phi", wetting_observer_contact);
    ObservedQuantityRecording<Vecd> write_leading_edge_position("Position", leading_edge_observer_contact);

    SummaryOutput3D summary_output("./output/SummaryOutput3D.dat"); // 创建汇总输出
    //----------------------------------------------------------------------
    //  Prepare the simulation
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_x.update_cell_linked_list_.exec();
    periodic_condition_z.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    object_normal_direction.exec();
    wetting_water_initial_condition.exec();
    wetting_wall_initial_condition.exec();
    wetting_object_initial_condition.exec();

    object_set_initial_velocity.exec();
    Real dt_thermal = get_diffusion_time_step.exec();
    free_stream_surface_indicator.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //  Time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1;
    int observation_sample_interval = screen_output_interval * 1;
    Real end_time = 0.009; 
    Real output_interval_vtp = end_time / 12.0;
    Real output_interval_force = end_time / 500.0;
    Real next_force_output = output_interval_force;
    Real next_vtp_output = output_interval_vtp;
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
    //  First output
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_object_wetting.writeToFile(number_of_iterations);
    write_leading_edge_position.writeToFile(number_of_iterations);
    std::cout << "Initialization complete. Starting simulation..." << std::endl;
    std::cout << "Object center: " << object_center.transpose() << std::endl;
    //----------------------------------------------------------------------
    //  Main loop
    //----------------------------------------------------------------------
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
            kernel_correction_complex.exec(); // with KGC correction
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real dt = 0.0;
            viscous_force_from_fluid.exec();
            while (relaxation_time < Dt)
            {
                /** inner loop for dual-time criteria time-stepping.*/
                dt = SMIN(SMIN(dt_thermal, fluid_acoustic_time_step.exec()), Dt);
                // 防止dt为0或负值
                if (dt < 1e-12)
                {
                    std::cerr << "Warning: dt too small: " << dt << std::endl;
                    dt = 1e-9;
                }

                fluid_pressure_relaxation.exec(dt);
                pressure_force_from_fluid.exec();
                fluid_density_relaxation.exec(dt);
                object_wetting.exec(dt);// 湿表面扩散

                 // 耦合多体动力学
                integ.stepBy(dt);
                SimTK::State &state_for_update = integ.updAdvancedState();
                force_on_bodies.clearAllBodyForces(state_for_update);
                force_on_bodies.setOneBodyForce(state_for_update, structure_mob,force_on_structure.exec());
                constraint_on_structure.exec();

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;
            // 输出和记录
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(6)
                          << "N=" << number_of_iterations
                          << " Time = " << physical_time
                          << " Dt = " << Dt << " dt = " << dt << "\n"
                          << " integration_time = " << integration_time<< "\n";

                if (number_of_iterations % observation_sample_interval == 0)
                {
                    write_object_wetting.writeToFile(number_of_iterations);
                    write_leading_edge_position.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            periodic_condition_x.bounding_.exec();
            periodic_condition_z.bounding_.exec();

            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            // 更新邻居列表和配置
                water_block.updateCellLinkedList();
                periodic_condition_x.update_cell_linked_list_.exec(); // 周期性修正
                periodic_condition_z.update_cell_linked_list_.exec();
                wall_boundary.updateCellLinkedList();
                floating_object.updateCellLinkedList();
                water_block_inner.updateConfiguration();
                object_inner.updateConfiguration();
                object_contact.updateConfiguration();
                water_block_complex.updateConfiguration();
                leading_edge_observer_contact.updateConfiguration();
                free_stream_surface_indicator.exec();
                interval_updating_configuration += TickCount::now() - time_instance;
        }

        if (physical_time >= next_force_output)
        {
            TickCount t2 = TickCount::now();
            viscous_force_from_fluid.exec();
            pressure_force_from_fluid.exec();
            auto &leading_edge_particles = leading_edge_observer.getBaseParticles(); // 获取前缘观测点数据
            Vecd leading_edge_pos = leading_edge_particles.getVariableDataByName<Vecd>("Position")[0];
             SimTK::State &output_state = integ.updAdvancedState();
            //  获取全局力
            write_total_viscous_force.writeToFile(number_of_iterations);
            write_total_pressure_force.writeToFile(number_of_iterations);

            // 计算粒子上的总力（全局坐标系）
            BaseParticles &object_particles = floating_object.getBaseParticles();
            Vec3d total_viscous_force_global = calculate_cylinder_total_viscous_force.exec();
            Vec3d total_pressure_force_global = calculate_cylinder_total_pressure_force.exec();

            SimTK::State current_state = integ.getAdvancedState(); // 获取当前物体的旋转角度（从Simbody状态）
            MBsystem.realize(current_state, SimTK::Stage::Velocity); // 确保状态已实现
            Real rotation_angle_z = initial_pitch_angle + getObjectRotationAngle(structure_mob, current_state);

            // 将全局力转换到局部坐标系
            Vecd viscous_force_local = transformGlobalForceToLocal(total_viscous_force_global, rotation_angle_z);
            Vecd pressure_force_local = transformGlobalForceToLocal(total_pressure_force_global, rotation_angle_z);

            Vecd total_force_local = viscous_force_local + pressure_force_local;

            summary_output.writeData(physical_time,
                                     leading_edge_pos,
                                     total_viscous_force_global,
                                     total_pressure_force_global,
                                     viscous_force_local,
                                     pressure_force_local,
                                     total_force_local); 
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
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
}


