// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef EFFECTIVEPOTENTIAL_ONE_LOOP_POTENTIAL_HPP_
#define EFFECTIVEPOTENTIAL_ONE_LOOP_POTENTIAL_HPP_

#include <vector>

#include <Eigen/Core>

#include "property.hpp"
#include "potential.hpp"

namespace EffectivePotential {

/** x * log(x) that safely treats x = 0 */
double xlogx(double);

/** Method for daisy corrections */
enum class DaisyMethod { None,
                         ArnoldEspinosa,
                         Parwani };

class OneLoopPotential : public Potential {
public:
  virtual double V0(Eigen::VectorXd phi) const = 0;
  /** Functions for squared field dependent masses, depending on:
      a vector of fields and for scalars a xi gauge parameter.  Note the
      latter will be unused for potentials implemented in a fixed gauge. */
  virtual std::vector<double> get_scalar_masses_sq(Eigen::VectorXd phi, double xi) const;
  virtual std::vector<double> get_fermion_masses_sq(Eigen::VectorXd phi) const { return {}; }
  virtual std::vector<double> get_vector_masses_sq(Eigen::VectorXd phi) const { return {}; }
  virtual std::vector<double> get_ghost_masses_sq(Eigen::VectorXd phi, double xi) const { return {}; }
  virtual std::vector<double> get_scalar_debye_sq(Eigen::VectorXd phi, double xi, double T) const { return {}; }
  virtual std::vector<double> get_scalar_thermal_sq(double T) const { return {}; }
  virtual std::vector<double> get_vector_debye_sq(Eigen::VectorXd phi, double T) const { return {}; }
  virtual std::vector<double> get_scalar_dofs() const;
  virtual std::vector<double> get_fermion_dofs() const { return {}; }
  virtual std::vector<double> get_vector_dofs() const { return {}; }
  virtual std::vector<double> get_ghost_dofs() const { return {}; }

  /** The Hessian matrix of tree-level potential */
  Eigen::MatrixXd d2V0_dx2(Eigen::VectorXd phi) const;
  /** The derivative of the gradient of potential with respect to temperature */
  Eigen::VectorXd d2V_dxdt(Eigen::VectorXd phi, double T) const;
  /** Finite-temperature effective potential */
  double V(Eigen::VectorXd phi, double T) const;
  /** Zero-temperature one-loop correction */
  virtual double V1(std::vector<double> scalar_masses_sq,
                    std::vector<double> fermion_masses_sq,
                    std::vector<double> vector_masses_sq,
                    std::vector<double> ghost_masses_sq) const;
  virtual double V1(Eigen::VectorXd phi, double T = 0.) const;
  /** Finite-temperature one-loop correction */
  double V1T(std::vector<double> scalar_masses_sq,
             std::vector<double> fermion_masses_sq,
             std::vector<double> vector_masses_sq,
             std::vector<double> ghost_masses_sq, double T) const;
  double V1T(Eigen::VectorXd phi, double T) const;
  /** Daisy corrections to potential */
  double daisy(std::vector<double> scalar_masses_sq,
               std::vector<double> scalar_debye_sq,
               std::vector<double> vector_masses_sq,
               std::vector<double> vector_debye_sq, double T) const;
  double daisy(Eigen::VectorXd phi, double T) const;

  /** Counter-term to potential */
  virtual double counter_term(Eigen::VectorXd phi, double T) const { return 0; }

  /** High-temperature expansion of potential */
  double VHT(Eigen::VectorXd phi, double T) const;

  /** Tree-level scalar masses */
  std::vector<double> get_tree_scalar_masses_sq(Eigen::VectorXd phi) const;

  /** One-loop scalar masses */
  std::vector<double> get_1l_scalar_masses_sq(Eigen::VectorXd phi, double T) const;

  /** The renormalization scale for the effective potential */
  void set_renormalization_scale(double Q) { renormalization_scale = Q; }
  double get_renormalization_scale() const { return renormalization_scale; }

  /** The gauge parameter \xi */
  void set_xi(double xi_in) { xi = xi_in; }
  double get_xi() const { return xi; }

  /** Treatment of the thermal masses */
  void set_daisy_method(DaisyMethod dm) { daisy_method = dm; }
  DaisyMethod get_daisy_method() const { return daisy_method; }

private:
  /** The renormalization scale for the effective potential */
  double renormalization_scale = 246.;
  /** The gauge parameter \xi */
  double xi = 0.;
  /** Treatment of the thermal masses */
  DaisyMethod daisy_method = DaisyMethod::ArnoldEspinosa;
};

} // namespace EffectivePotential

#endif // EFFECTIVEPOTENTIAL_ONE_LOOP_POTENTIAL_HPP_
