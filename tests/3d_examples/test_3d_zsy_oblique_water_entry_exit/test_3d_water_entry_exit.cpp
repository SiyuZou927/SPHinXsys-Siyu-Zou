/**
 * @file    test_3d_water_entry_exit.cpp
 * @brief   3D water entry and exit example with surface wetting considered.
 * @details 3D fluid-structure interaction with STL geometry and wetting-coupled
 *          spatial-temporal free-surface identification.
 * @author  Siyu Zou Shuoguo Zhang and Xiangyu Hu
 */
#include "SimTKsimbody.h"
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//  Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.2; /**< Water tank length. */
Real DW = 0.4; /**< Water tank width. */
Real DH = 1.0; /**< Water tank height. */
Real LH = 0.5; /**< Water column height. */

Real particle_spacing_ref = 0.005;           /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;          /**< Thickness of the bottom wall. */
Vecd object_center(0.46, LH + 0.028, DW / 2); /**< Reference location used to initialize the observer. */

Real initial_speed = 50.0;                   /**< Initial velocity magnitude (m/s). */
Real initial_angle = -20.0 * Pi / 180.0;     /**< Initial velocity angle (radians, negative = downward). */
Real initial_pitch_angle = -20.0 * Pi / 180.0; /**< Initial pitch offset for reporting the local frame (radians). */
std::string full_path_to_file = "./input/M1TR4Z02H05.stl"; // STL geometry includes the initial position and orientation.
//----------------------------------------------------------------------
//  Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.00;         /**< Water density. */
Real rho0_s = 2.063;        /**< Object density. */
Real gravity_g = 9.81;      /**< Gravitational acceleration. */
Real U_max = initial_speed; /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;    /**< Reference sound speed. */
Real mu_f = 8.9e-7;         /**< Water dynamic viscosity, consistent with rho0_f. */
//----------------------------------------------------------------------
//  Wetting parameters
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 0.0 * pow(particle_spacing_ref, 2); /**< Wetting coefficient. */
Real fluid_moisture = 1.0;                               /**< Initial fluid moisture. */
Real object_moisture = 0.0;                              /**< Initial object moisture. */
Real wall_moisture = 1.0;                                /**< Wall moisture parameter. */
//----------------------------------------------------------------------
//  Transform global forces to the local frame rotated about Z.
//----------------------------------------------------------------------
Vec3d transformGlobalForceToLocal(const Vec3d &global_force, Real rotation_angle_z)
{
    Real cos_theta = cos(rotation_angle_z);
    Real sin_theta = sin(rotation_angle_z);
    Vec3d local_force;
    local_force[0] = global_force[0] * cos_theta + global_force[1] * sin_theta;  // Local X component.
    local_force[1] = -global_force[0] * sin_theta + global_force[1] * cos_theta; // Local Y component.
    local_force[2] = global_force[2];
    return local_force;
}

Real calculateHydrodynamicMomentZ(BaseParticles &particles, const Vec3d &mass_center)
{
    Vecd *position = particles.getVariableDataByName<Vecd>("Position");
    Vecd *viscous_force = particles.getVariableDataByName<Vecd>("ViscousForceFromFluid");
    Vecd *pressure_force = particles.getVariableDataByName<Vecd>("PressureForceFromFluid");
    Real moment_z = 0.0;
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
    {
        Vec3d displacement = position[i] - mass_center;
        Vec3d hydrodynamic_force = viscous_force[i] + pressure_force[i];
        moment_z += displacement[0] * hydrodynamic_force[1] -
                    displacement[1] * hydrodynamic_force[0];
    }
    return moment_z;
}

bool timeReached(Real current_time, Real event_time)
{
    Real tolerance = 100.0 * Eps * std::max(1.0, std::abs(event_time));
    return current_time >= event_time || event_time - current_time <= tolerance;
}
//----------------------------------------------------------------------
//  Observer reference location.
//----------------------------------------------------------------------
StdVec<Vecd> leading_edge_location = {object_center}; /**< Observer initialized at the reference location above. */
//----------------------------------------------------------------------
//  Summary output for the 3D case.
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
                     << "TotalForceLocal_X TotalForceLocal_Y TotalForceLocal_Z   "
                     << "MassCenter_X MassCenter_Y MassCenter_Z   "
                     << "Velocity_X Velocity_Y Velocity_Z   "
                     << "RotationAngle_Z AngularVelocity_Z HydrodynamicMoment_Z\n";
    }

    void writeData(Real time,
                   const Vecd &leading_edge_position,
                   const Vecd &viscous_force_global,
                   const Vecd &pressure_force_global,
                   const Vecd &viscous_force_local,
                   const Vecd &pressure_force_local,
                   const Vecd &total_force_local,
                   const Vec3d &mass_center,
                   const Vec3d &velocity,
                   Real rotation_angle_z,
                   Real angular_velocity_z,
                   Real hydrodynamic_moment_z)
    {
        output_file_ << std::scientific << std::setprecision(9)
                     << time << " "
                     << leading_edge_position[0] << " " << leading_edge_position[1] << " " << leading_edge_position[2] << " "
                     << viscous_force_global[0] << " " << viscous_force_global[1] << " " << viscous_force_global[2] << " "
                     << pressure_force_global[0] << " " << pressure_force_global[1] << " " << pressure_force_global[2] << " "
                     << viscous_force_local[0] << " " << viscous_force_local[1] << " " << viscous_force_local[2] << " "
                     << pressure_force_local[0] << " " << pressure_force_local[1] << " " << pressure_force_local[2] << " "
                     << total_force_local[0] << " " << total_force_local[1] << " " << total_force_local[2] << " "
                     << mass_center[0] << " " << mass_center[1] << " " << mass_center[2] << " "
                     << velocity[0] << " " << velocity[1] << " " << velocity[2] << " "
                     << rotation_angle_z << " " << angular_velocity_z << " "
                     << hydrodynamic_moment_z << "\n";
        output_file_.flush();
    }

  private:
    std::ofstream output_file_;
};
//----------------------------------------------------------------------
//  Water body geometry.
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(Vecd(DL / 2, LH / 2, DW / 2)), Vecd(DL / 2, LH / 2, DW / 2));
    }
};

//----------------------------------------------------------------------
//  Bottom wall geometry.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(Vecd(DL / 2, -BW / 2, DW / 2)), Vecd(DL / 2 + BW, BW / 2, DW / 2 + BW));
    }
};
//----------------------------------------------------------------------
//  Object geometry loaded from STL without an additional transform.
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
//  Simbody-constrained region using the same STL as the object.
//----------------------------------------------------------------------
std::shared_ptr<ComplexShape> createSimbodyConstrainShape3D()
{
    auto constrain_shape = makeShared<ComplexShape>("SimbodyConstrainShape");
    Vecd translation(0.0, 0.0, 0.0);
    constrain_shape->add<TriangleMeshShapeSTL>(full_path_to_file, translation, 1.0);
    return constrain_shape;
}
//----------------------------------------------------------------------
//  Initial wetting conditions.
//----------------------------------------------------------------------
template <class BodyType>
class WettingBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingBodyInitialCondition(BodyType &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name))
    {
    };

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
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name))
    {
    }

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = object_moisture;
    }

  protected:
    Real *phi_;
    Vecd *pos_;
};

using WettingFluidBodyInitialCondition = WettingBodyInitialCondition<FluidBody>;
using WettingWallBodyInitialCondition = WettingBodyInitialCondition<SolidBody>;
//----------------------------------------------------------------------
//  Wetting diffusion model.
//----------------------------------------------------------------------
using ObjectFluidDiffusionDirichlet =
    DiffusionRelaxationRK2<DiffusionRelaxation<Dirichlet<KernelGradientContact>, IsotropicDiffusion>>;

//----------------------------------------------------------------------
//  Z rotation from Simbody's body-fixed XYZ angles.
//----------------------------------------------------------------------
Real getObjectRotationAngle(SimTK::MobilizedBody &mob, const SimTK::State &state)
{
    SimTK::Rotation rot = mob.getBodyRotation(state);
    SimTK::Vec3 angles = rot.convertRotationToBodyFixedXYZ();
    return angles[2];
}
//----------------------------------------------------------------------
//  Set the initial object-particle velocity.
//----------------------------------------------------------------------
class ObjectInitialVelocity : public LocalDynamics
{
  public:
    ObjectInitialVelocity(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity"))
    {
    }

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][0] = initial_speed * cos(initial_angle);
        vel_[index_i][1] = initial_speed * sin(initial_angle);
        vel_[index_i][2] = 0.0; // No initial velocity in the Z direction.
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
    //  Set up the SPH system.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);

    // Run configuration.
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    // sph_system.setRunParticleRelaxation(true);
    // sph_system.setReloadParticles(false);
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
    floating_object.defineAdaptationRatios(1.15, 2.0); // Multi-resolution adaptation.
    floating_object.defineBodyLevelSetShape();

    // Solid material with wetting diffusion.
    floating_object.defineClosure<Solid, IsotropicDiffusion>(
        rho0_s, ConstructArgs(diffusion_species_name, diffusion_coeff));
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? floating_object.generateParticles<BaseParticles, Reload>(floating_object.getName())
        : floating_object.generateParticles<BaseParticles, Lattice>();

    // Register the additional state variable.
    BaseParticles &object_particles = floating_object.getBaseParticles();
    object_particles.registerStateVariableData<Vecd>("Velocity");

    // Observer body at the prescribed reference location.
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
    //  Combined relation used to update the water configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);

    //----------------------------------------------------------------------
    //  Run particle relaxation.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        std::cout << "Starting 3D particle relaxation..." << std::endl;
        InnerRelation object_inner(floating_object);
        //----------------------------------------------------------------------
        //  Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(floating_object);
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(floating_object); /** Write body states to VTP. */
        ReloadParticleIO write_particle_reload_files(floating_object);       /** Write particle reload files. */
        RelaxationStepInner relaxation_step_inner(object_inner);             /** Particle relaxation step. */
        //----------------------------------------------------------------------
        //  Initialize particle relaxation.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0); // Output the initial state.
        //----------------------------------------------------------------------
        //  Relax the object particles.
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

        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //  Define the fluid dynamics
    //----------------------------------------------------------------------
    // Wetting diffusion model.
    GetDiffusionTimeStepSize get_diffusion_time_step(floating_object);
    ObjectFluidDiffusionDirichlet object_wetting(object_contact);

    // Initial wetting conditions.
    SimpleDynamics<WettingBodyInitialCondition<FluidBody>> wetting_water_initial_condition(water_block);
    SimpleDynamics<WettingBodyInitialCondition<SolidBody>> wetting_wall_initial_condition(wall_boundary);
    SimpleDynamics<WettingObjectBodyInitialCondition> wetting_object_initial_condition(floating_object);

    SimpleDynamics<ObjectInitialVelocity> object_set_initial_velocity(floating_object);

    Gravity gravity(Vecd(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);

    // Wetting-coupled spatial-temporal free-surface identification.
    InteractionWithUpdate<WettingCoupledSpatialTemporalFreeSurfaceIndicationComplex>
        free_stream_surface_indicator(water_block_inner, water_block_contact);

    // Surface normal calculation.
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> object_normal_direction(floating_object);
    // InteractionWithUpdate<LinearGradientCorrectionMatrixComplex>
    //    kernel_correction_complex(DynamicsArgs(water_block_inner, 0.5), water_block_contact);
    // Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);

    // Fluid dynamics.
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);

    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);

    // Periodic boundaries along X and Z.
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicAlongAxis periodic_along_z(water_block.getSPHBodyBounds(), zAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_block, periodic_along_x);
    PeriodicConditionUsingCellLinkedList periodic_condition_z(water_block, periodic_along_z);
    // Time-step control.
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_time_step(water_block);
    //----------------------------------------------------------------------
    //  Algorithms of FSI
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(object_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(fluid_density_relaxation)>> pressure_force_from_fluid(object_contact);

    // Force totals in the global frame.
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_viscous_force(floating_object, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force(floating_object, "PressureForceFromFluid");

    ReduceDynamics<QuantitySummation<Vecd>> calculate_cylinder_total_pressure_force(floating_object, "PressureForceFromFluid");
    ReduceDynamics<QuantitySummation<Vecd>> calculate_cylinder_total_viscous_force(floating_object, "ViscousForceFromFluid");
    //----------------------------------------------------------------------
    //  Particle sorting for configuration updates.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //  Set up Simbody for planar rigid-body motion in the 3D domain.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;                /**< Multibody system. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem); /**< Bodies and mobilizers. */
    SimTK::GeneralForceSubsystem forces(MBsystem);  /**< Applied forces. */

    SimTK::Body::Rigid fixed_spot_info(SimTK::MassProperties(1.0, SimTKVec3(0), SimTK::UnitInertia(1)));
    // Constrain the SPH body part selected by the object STL.
    SolidBodyPartForSimbody structure_system(floating_object, createSimbodyConstrainShape3D());
    // Place the fixed reference frame at the initial body-part mass center.
    Vecd tethering_point = structure_system.initial_mass_center_;
    SimTK::Body::Rigid tethered_spot_info(*structure_system.body_part_mass_properties_);

    // Weld the reference body to ground.
    SimTK::MobilizedBody::Weld fixed_spot(matter.Ground(), SimTK::Transform(SimTKVec3(tethering_point[0], tethering_point[1], tethering_point[2])),
                                         fixed_spot_info, SimTK::Transform(SimTKVec3(0)));

    Vecd displacement0 = structure_system.initial_mass_center_ - tethering_point;

    // Allow X-Y translation and rotation about Z relative to the reference body.
    SimTK::MobilizedBody::Planar structure_mob(fixed_spot,
                                             SimTK::Transform(SimTKVec3(displacement0[0], displacement0[1], displacement0[2])),
                                             tethered_spot_info, SimTK::Transform(SimTKVec3(0)));

    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTKVec3(0.0, -gravity_g, 0.0), 0.0);
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter); // Discrete SPH coupling forces.
    fixed_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.02));    // Reference-body decoration.
    tethered_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.4)); // Object-body decoration.
    SimTK::State state = MBsystem.realizeTopology(); // Initialize the system state.

    SimTK::Vec3 mobilizer_vel(0.0, initial_speed * cos(initial_angle), initial_speed * sin(initial_angle));
    structure_mob.setU(state, mobilizer_vel);

    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);

    // Report the initial Simbody state.
    std::cout << "Debug: Joint U = ("
              << structure_mob.getU(state)[0] << ", "
              << structure_mob.getU(state)[1] << ", "
              << structure_mob.getU(state)[2] << ")" << std::endl;
    std::cout << "Mass: " << structure_system.body_part_mass_properties_->getMass() << std::endl;
    std::cout << "Center of mass: " << structure_system.initial_mass_center_.transpose() << std::endl;
    //----------------------------------------------------------------------
    //  Coupling between Simbody and SPH.
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

    SummaryOutput3D summary_output("./output/SummaryOutput3D.dat");
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
    Real end_time = 0.012;
    Real output_interval_vtp = end_time / 12.0;
    Real output_interval_force = end_time / 300.0;
    Real next_force_output = output_interval_force;
    Real next_vtp_output = output_interval_vtp;
    //----------------------------------------------------------------------
    //  Wall-clock timing statistics.
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
    while (!timeReached(physical_time, end_time))
    {
        Real next_event_time = std::min(end_time, std::min(next_force_output, next_vtp_output));
        Real target_time = timeReached(physical_time, next_event_time)
                               ? 0.0
                               : next_event_time - physical_time;
        Real integration_time = 0.0;

        while (integration_time < target_time && !timeReached(physical_time, end_time))
        {
            time_instance = TickCount::now();
            Real Dt = fluid_advection_time_step.exec();
            fluid_density_by_summation.exec();
            viscous_force.exec();
            // kernel_correction_complex.exec(); // With kernel gradient correction.
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real dt = 0.0;
            viscous_force_from_fluid.exec();
            while (relaxation_time < Dt && integration_time < target_time &&
                   !timeReached(physical_time, end_time))
            {
                // Inner loop for dual-criteria time stepping.
                Real remaining_relaxation_time = Dt - relaxation_time;
                Real remaining_integration_time = target_time - integration_time;
                Real remaining_physical_time = end_time - physical_time;
                if (remaining_relaxation_time <= TinyReal || remaining_integration_time <= TinyReal ||
                    remaining_physical_time <= TinyReal)
                    break;
                dt = SMIN(dt_thermal, fluid_acoustic_time_step.exec(), remaining_relaxation_time,
                          remaining_integration_time, remaining_physical_time);

                fluid_pressure_relaxation.exec(dt);
                pressure_force_from_fluid.exec();
                fluid_density_relaxation.exec(dt);
                object_wetting.exec(dt); // Wetting diffusion.

                // Couple SPH and multibody dynamics.
                integ.stepBy(dt);
                SimTK::State &state_for_update = integ.updAdvancedState();
                force_on_bodies.clearAllBodyForces(state_for_update);
                force_on_bodies.setOneBodyForce(state_for_update, structure_mob, force_on_structure.exec());
                constraint_on_structure.exec();

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;
            // Screen output and observer recording.
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(6)
                          << "N=" << number_of_iterations
                          << " Time = " << physical_time
                          << " Dt = " << Dt << " dt = " << dt << "\n"
                          << " integration_time = " << integration_time << "\n";

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
            // Update cell-linked lists, periodic entries, and neighbor configurations.
            water_block.updateCellLinkedList();
            periodic_condition_x.update_cell_linked_list_.exec();
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

        if (timeReached(physical_time, next_force_output))
        {
            TickCount t2 = TickCount::now();
            viscous_force_from_fluid.exec();
            pressure_force_from_fluid.exec();
            auto &leading_edge_particles = leading_edge_observer.getBaseParticles(); // Read observer-particle data.
            Vecd leading_edge_pos = leading_edge_particles.getVariableDataByName<Vecd>("Position")[0];
            SimTK::State &output_state = integ.updAdvancedState();
            // Record force totals in the global frame.
            write_total_viscous_force.writeToFile(number_of_iterations);
            write_total_pressure_force.writeToFile(number_of_iterations);

            // Sum hydrodynamic forces on the object particles in the global frame.
            BaseParticles &object_particles = floating_object.getBaseParticles();
            Vec3d total_viscous_force_global = calculate_cylinder_total_viscous_force.exec();
            Vec3d total_pressure_force_global = calculate_cylinder_total_pressure_force.exec();

            SimTK::State current_state = integ.getAdvancedState(); // Current rigid-body state.
            MBsystem.realize(current_state, SimTK::Stage::Velocity); // Realize positions and velocities.
            Real rotation_angle_z = initial_pitch_angle + getObjectRotationAngle(structure_mob, current_state);

            SimTKVec3 origin = structure_mob.getBodyOriginLocation(current_state);
            SimTKVec3 origin_velocity = structure_mob.getBodyOriginVelocity(current_state);
            Vec3d mass_center(origin[0], origin[1], origin[2]);
            Vec3d velocity(origin_velocity[0], origin_velocity[1], origin_velocity[2]);
            Real angular_velocity_z = structure_mob.getBodyAngularVelocity(current_state)[2];
            Real hydrodynamic_moment_z =
                calculateHydrodynamicMomentZ(object_particles, mass_center);

            // Transform global forces to the local frame.
            Vecd viscous_force_local = transformGlobalForceToLocal(total_viscous_force_global, rotation_angle_z);
            Vecd pressure_force_local = transformGlobalForceToLocal(total_pressure_force_global, rotation_angle_z);

            Vecd total_force_local = viscous_force_local + pressure_force_local;

            summary_output.writeData(physical_time,
                                     leading_edge_pos,
                                     total_viscous_force_global,
                                     total_pressure_force_global,
                                     viscous_force_local,
                                     pressure_force_local,
                                     total_force_local,
                                     mass_center,
                                     velocity,
                                     rotation_angle_z,
                                     angular_velocity_z,
                                     hydrodynamic_moment_z);
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
            next_force_output += output_interval_force;
        }
        if (timeReached(physical_time, next_vtp_output))
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
