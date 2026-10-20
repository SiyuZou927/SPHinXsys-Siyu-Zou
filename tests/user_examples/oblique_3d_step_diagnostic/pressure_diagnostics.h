#ifndef OBLIQUE_PRESSURE_DIAGNOSTICS_H
#define OBLIQUE_PRESSURE_DIAGNOSTICS_H

#include <array>
#include <cstdlib>
#include <fstream>
#include <iomanip>

// The base interaction remains authoritative. Diagnostics use the same state
// and cached contact geometry, and run only at scheduled recording events.
template <class FluidIntegration2ndHalfType>
class DiagnosticPressureForce
    : public solid_dynamics::PressureForceFromFluid<FluidIntegration2ndHalfType>
{
    using Base = solid_dynamics::PressureForceFromFluid<FluidIntegration2ndHalfType>;
    static constexpr size_t bins = 12;
    struct Terms
    {
        Vecd central = Vecd::Zero(), acceleration = Vecd::Zero();
        Vecd riemann = Vecd::Zero(), negative = Vecd::Zero();
        size_t pairs = 0, saturated = 0;
        Real max_distance_change = 0.0;
    };
    StdVec<Terms> terms_;
    StdVec<int> bin_;
    StdVec<bool> shoulder_;
    Vecd *position_, *viscous_;
    Real *volume_;
    SimTK::SpatialVec held_wrench_{SimTK::Vec3(0), SimTK::Vec3(0)};
    StdVec<Vecd *> fluid_position_;
    bool collecting_ = false;
    std::ofstream output_, wrench_, geometry_;

  public:
    explicit DiagnosticPressureForce(BaseContactRelation &relation)
        : Base(relation), terms_(this->particles_->TotalRealParticles()),
          bin_(terms_.size()), shoulder_(terms_.size()),
          position_(this->particles_->template getVariableDataByName<Vecd>("Position")),
          viscous_(this->particles_->template getVariableDataByName<Vecd>("ViscousForceFromFluid")),
          volume_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
          output_("./output/PressureDecomposition.csv"),
          wrench_("./output/AppliedWrench.csv"),
          geometry_("./output/DiagnosticGeometry.csv")
    {
        for (auto *p : this->contact_particles_)
            fluid_position_.push_back(p->template getVariableDataByName<Vecd>("Position"));
        output_ << "stage,time_s,simbody_time_s,region,particle_count,contact_particle_count,pair_count,saturated_pair_count,pressure_axial,central_axial,acceleration_axial,riemann_axial,viscous_axial,total_axial,negative_central_axial,pressure_moment_z,central_moment_z,acceleration_moment_z,riemann_moment_z,closure_norm,max_cached_distance_change_m\n";
        wrench_ << "time_s,simbody_time_s,dt_s,assigned_torque_x,assigned_torque_y,assigned_torque_z,assigned_force_x,assigned_force_y,assigned_force_z,used_torque_x,used_torque_y,used_torque_z,used_force_x,used_force_y,used_force_z\n";
        geometry_ << "particle_index,initial_x,initial_y,initial_z,station_m,radius_m,bin,annular_shoulder_mask,volume_m3,simbody_constrained\n";
        output_ << std::scientific << std::setprecision(17);
        wrench_ << std::scientific << std::setprecision(17);
        geometry_ << std::scientific << std::setprecision(17);
    }

    void initializeMaterialRegions(const IndexVector &constrained_indices)
    {
        StdVec<bool> constrained(terms_.size(), false);
        for (size_t i : constrained_indices) constrained[i] = true;
        Vecd axis(std::cos(initial_pitch_angle), std::sin(initial_pitch_angle), 0.0);
        for (size_t i = 0; i < terms_.size(); ++i)
        {
            Vecd d = object_center - position_[i];
            Real station = d.dot(axis);
            Real radius = (d - station * axis).norm();
            bin_[i] = std::max(0, std::min(int(bins) - 1, int(std::floor(station / 0.04))));
            // Supplementary overlapping mask; axial bins alone form a partition.
            shoulder_[i] = std::abs(station - 0.4029) <= 0.00375 &&
                           radius >= 0.016 && std::abs(this->n_[i].dot(axis)) >= 0.7;
            geometry_ << i << ',' << position_[i][0] << ',' << position_[i][1] << ','
                      << position_[i][2] << ',' << station << ',' << radius << ','
                      << bin_[i] << ',' << shoulder_[i] << ',' << volume_[i]
                      << ',' << constrained[i] << '\n';
        }
        geometry_.flush();
    }

    void collect(bool enabled) { collecting_ = enabled; }

    void interaction(size_t i, Real dt = 0.0)
    {
        Base::interaction(i, dt);
        if (!collecting_) return;
        Terms t;
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            const Neighborhood &neighbors = (*this->contact_configuration_[k])[i];
            for (size_t n = 0; n < neighbors.current_size_; ++n)
            {
                size_t j = neighbors.j_[n];
                const Vecd &e = neighbors.e_ij_[n];
                Real r = neighbors.r_ij_[n];
                Real external_acceleration =
                    (this->contact_force_prior_[k][j] / this->contact_mass_[k][j] - this->acc_ave_[i]).dot(e);
                Vecd normal = -SGN(e.dot(this->n_[i])) * this->n_[i];
                Real jump = 2.0 * (this->contact_vel_[k][j] - this->vel_ave_[i]).dot(normal);
                Real weight = -neighbors.dW_ij_[n] * this->contact_Vol_[k][j] * this->Vol_[i];
                Vecd central = (2.0 * this->contact_p_[k][j] * weight) * e;
                t.central += central;
                t.acceleration += (this->contact_rho_[k][j] * r *
                                   SMAX(Real(0), external_acceleration) * weight) * e;
                t.riemann += (this->riemann_solvers_[k].DissipativePJump(jump) * weight) * normal;
                if (this->contact_p_[k][j] < 0.0) t.negative += central;
                ++t.pairs;
                if (3.0 * SMAX(Real(0), jump) >= c_f) ++t.saturated;
                Vecd current_delta = position_[i] - fluid_position_[k][j];
                current_delta[0] -= DL * std::round(current_delta[0] / DL);
                current_delta[2] -= DW * std::round(current_delta[2] / DW);
                Real current_r = current_delta.norm();
                t.max_distance_change = std::max(t.max_distance_change, std::abs(current_r - r));
            }
        }
        terms_[i] = t;
    }

    void write(const char *stage, Real time, Real simbody_time, Real angle, const Vecd &center)
    {
        struct Sum
        {
            size_t particles = 0, contacted = 0, pairs = 0, saturated = 0;
            std::array<Vecd, 6> forces;
            std::array<Real, 4> moments{};
            Real max_change = 0.0;
            Sum() { for (auto &f : forces) f.setZero(); }
        };
        std::array<Sum, bins + 2> sums;
        for (size_t i = 0; i < terms_.size(); ++i)
        {
            const auto &t = terms_[i];
            std::array<Vecd, 6> f{this->force_from_fluid_[i], t.central, t.acceleration,
                                 t.riemann, viscous_[i], t.negative};
            Vecd arm = position_[i] - center;
            for (size_t region : {size_t(bin_[i]), bins}) add(sums[region], t, f, arm);
            if (shoulder_[i]) add(sums[bins + 1], t, f, arm);
        }
        Vecd axis(std::cos(angle), std::sin(angle), 0.0);
        for (size_t region = 0; region < sums.size(); ++region)
        {
            const auto &s = sums[region];
            output_ << stage << ',' << time << ',' << simbody_time << ',';
            if (region < bins) output_ << "axial_bin_" << region;
            else output_ << (region == bins ? "all" : "annular_shoulder_mask");
            output_ << ',' << s.particles << ',' << s.contacted << ',' << s.pairs << ',' << s.saturated;
            for (size_t j = 0; j < 5; ++j) output_ << ',' << s.forces[j].dot(axis);
            output_ << ',' << (s.forces[0] + s.forces[4]).dot(axis) << ',' << s.forces[5].dot(axis);
            for (Real m : s.moments) output_ << ',' << m;
            output_ << ',' << (s.forces[0] - s.forces[1] - s.forces[2] - s.forces[3]).norm()
                    << ',' << s.max_change << '\n';
        }
        output_.flush();
    }

    void writeAppliedWrench(Real time, Real simbody_time, Real dt,
                           const SimTK::SpatialVec &w, bool record)
    {
        if (record)
        {
            wrench_ << time << ',' << simbody_time << ',' << dt;
            for (const auto &value : {w, held_wrench_})
                for (int part = 0; part < 2; ++part)
                    for (int j = 0; j < 3; ++j) wrench_ << ',' << value[part][j];
            wrench_ << '\n';
            wrench_.flush();
        }
        held_wrench_ = w;
    }

  private:
    template <class SumType>
    static void add(SumType &sum, const Terms &t, const std::array<Vecd, 6> &forces, const Vecd &arm)
    {
        ++sum.particles;
        sum.contacted += t.pairs > 0;
        sum.pairs += t.pairs;
        sum.saturated += t.saturated;
        for (size_t j = 0; j < forces.size(); ++j) sum.forces[j] += forces[j];
        for (size_t j = 0; j < 4; ++j) sum.moments[j] += arm[0] * forces[j][1] - arm[1] * forces[j][0];
        sum.max_change = std::max(sum.max_change, t.max_distance_change);
    }
};

#endif
