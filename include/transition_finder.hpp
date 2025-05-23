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

#ifndef PHASETRACER_TRANSITION_FINDER_HPP_
#define PHASETRACER_TRANSITION_FINDER_HPP_

#include <algorithm>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <boost/cstdint.hpp>

#include "phase_finder.hpp"
#include "overload.hpp"
#include "potential.hpp"
#include "action_calculator.hpp"
#include "transition_graph_util.hpp"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace TransitionGraph {
struct Path;
}

namespace PhaseTracer {

/** Information about root-finding for a transition */
enum Message { SUCCESS,
               NON_OVERLAPPING_T,
               ERROR };

struct Transition {
  /** Data about a particular transition */
  Message message;
  double TC;
  Phase true_phase;
  Phase false_phase;
  Eigen::VectorXd true_vacuum;
  Eigen::VectorXd false_vacuum;
  double gamma;
  std::vector<bool> changed;
  double delta_potential;
  size_t key;
  size_t id;

  double TN;
  Eigen::VectorXd true_vacuum_TN;
  Eigen::VectorXd false_vacuum_TN;
  bool subcritical = false;

  Transition(Message message) : message(message) {};

  Transition(Message message,
             double TC,
             Phase true_phase,
             Phase false_phase,
             Eigen::VectorXd true_vacuum,
             Eigen::VectorXd false_vacuum,
             double gamma,
             std::vector<bool> changed,
             double delta_potential,
             size_t key,
             size_t id) : message(message),
                          TC(TC),
                          true_phase(true_phase),
                          false_phase(false_phase),
                          true_vacuum(true_vacuum),
                          false_vacuum(false_vacuum),
                          gamma(gamma),
                          changed(changed),
                          delta_potential(delta_potential),
                          key(key),
                          id(id) {}

  void set_subcritical(bool subcritical_) {
    subcritical = subcritical_;
  }

  void set_nucleation(double T, Eigen::VectorXd true_vacuum, Eigen::VectorXd false_vacuum) {
    calculate_action = true;
    TN = T;
    true_vacuum_TN = true_vacuum;
    false_vacuum_TN = false_vacuum;
  }

  /** Pretty-printer for single transition */
  friend std::ostream &operator<<(std::ostream &o, const Transition &a) {
    if (a.message == SUCCESS) {
      o << "=== transition from phase " << a.false_phase.key;
      if (a.key == 0) {
        o << " to phase " << a.true_phase.key << " ===" << std::endl;
      } else {
        o << " to symmetric partner " << a.key
          << " of phase " << a.true_phase.key << " ===" << std::endl;
      }
      o << "changed = " << a.changed << std::endl
        << "TC = " << a.TC << std::endl
        << "false vacuum (TC) = " << a.false_vacuum << std::endl
        << "true vacuum (TC) = " << a.true_vacuum << std::endl
        << "gamma (TC) = " << a.gamma << std::endl
        << "delta potential (TC) = " << a.delta_potential << std::endl;

      if (a.calculate_action) {
        o << "TN = " << a.TN << std::endl
          << "false vacuum (TN) = " << a.true_vacuum_TN << std::endl
          << "true vacuum (TN) = " << a.false_vacuum_TN << std::endl;
      } else {
        o << "did not calculate action or check nucleation" << std::endl;
      }

      o << "transition was "
        << (a.subcritical ? "subcritical" : "not subcritical") << std::endl;

    } else {
      o << "=== failure. message =  " << a.message << " ===" << std::endl;
    }
    return o;
  }

private:
  bool calculate_action = false;
};

class TransitionFinder {
public:
  explicit TransitionFinder(PhaseFinder &pf_) : pf(pf_), ac(pf_.P) {
    calculate_action = false;
  }
  explicit TransitionFinder(PhaseFinder &pf_, ActionCalculator ac_) : pf(pf_), ac(ac_) {
    calculate_action = true;
  }
  virtual ~TransitionFinder() = default;

  /** Find all transitions between all phases */
  void find_transitions();

  double get_v_T();

  json get_S_T(int nS);

  /** Called from find_transitions; adds subcritical transitions to the transitions vector */
  void append_subcritical_transitions();

  /** Called from append_subcritical_transitions; checks whether there is a subcritical transition between two phases. */
  bool checkSubcriticalTransition(const std::vector<PhaseTracer::Phase> &phases, int i, int j, double Tmax,
                                  double energyAtTmax, bool checkFromNewPhase, std::vector<bool> &isTransitionedTo);

  /** Find all transition paths  */
  void find_transition_paths(const EffectivePotential::Potential &model, bool knownHighTPhase);

  /** Retrieve all transitions between all phases */
  std::vector<Transition> get_transitions() const { return transitions; }

  std::vector<Eigen::VectorXd> get_vacua_at_T(const Phase &phase1, const Phase &phase2, double T, size_t i_unique = 0) const;

  double get_action(const Eigen::VectorXd &vacuum_1, const Eigen::VectorXd &vacuum_2, double T) const;

  double get_action(const Phase &phase1, const Phase &phase2, double T, size_t i_unique = 0) const;

  std::vector<double> get_action(const Phase &phase1, const Phase &phase2, std::vector<double> T_list, size_t i_unique = 0) const;

  void write_action_to_text(const Phase &phase1, const Phase &phase2, std::vector<double> T_list, const std::string &filename, size_t i_unique = 0) const;

  void write_action_to_text(const Transition &tran, double T_min, double T_max, size_t n_step, const std::string &filename, size_t i_unique = 0) const;

  void write_action_to_text(const Transition &tran, const std::string &filename, size_t n_step = 50, size_t i_unique = 0) const;

  double get_Tnuc(const Phase &phase1, const Phase &phase2, size_t i_unique, double T_begin, double T_end) const;

  /** Retrieve all transition paths */
  std::vector<TransitionGraph::Path> get_transition_paths() const { return transition_paths; }

  /** Retrieve all phases */
  std::vector<Phase> get_phases() const { return pf.get_phases(); }

  /** Pretty-printer for set of transitions in this object */
  friend std::ostream &operator<<(std::ostream &o, const TransitionFinder &a);

  /** Object with phases and potential */
  PhaseFinder &pf;

private:
  ActionCalculator ac;

  /** Whether already calculated all transitions */
  bool calculated_transitions = false;

  /** Container for all transitions between any two phases */
  std::vector<Transition> transitions;

  /** Container for all transition paths. */
  std::vector<TransitionGraph::Path> transition_paths;

  /** Find transitions between two phases between two temperatures */
  std::vector<Transition> find_transition(const Phase &phase1, const Phase &phase2, double T1, double T2, size_t currentID) const;

  std::vector<Transition> symmetric_partners(const Phase &phase1, const Phase &phase2, double TC, size_t currentID) const;

  /** Find critical temperature between two phases */
  double find_critical_temperature(const Phase &phase1, const Phase &phase2, double T1, double T2) const;

  /** Find many transitions between two phases at a particular resolution */
  std::vector<Transition> divide_and_find_transition(const Phase &phase1, const Phase &phase2, double T1, double T2, size_t currentID) const;

  /** Find un-overlapped temperature region between the two phases*/
  std::vector<double> get_un_overlapped_T_range(const Phase &phase1, const Phase &phase2, double T1, double T2) const;

  /** Strength of phase transition for first N fields */
  double gamma(const Eigen::VectorXd &true_vacuum, const Eigen::VectorXd &false_vacuum, const double TC) const;

  /** Note which VEVs changed */
  std::vector<bool> changed(const Eigen::VectorXd &true_vacuum, const Eigen::VectorXd &false_vacuum) const;

  /** Number of scalar fields that could break electroweak symmetry */
  PROPERTY(int, n_ew_scalars, -1)

  /** Minimum separation between critical temperatures */
  PROPERTY(double, separation, 1.)
  /** Assume at most one critical temperature between two phases */
  PROPERTY(bool, assume_only_one_transition, true)
  /** Relative precision in critical temperature */
  //PROPERTY(double, TC_tol_rel, 1.e-4)
  PROPERTY(double, TC_tol_rel, 1.e-4)
  /** Maximum number of iterations when finding a critical temperature */
  PROPERTY(boost::uintmax_t, max_iter, 100)
  /** Relative tolerance for judging whether a field is changed during a transition */
  PROPERTY(double, change_rel_tol, 1.e-3)
  /** Absolute tolerance for judging whether a field is changed during a transition */
  PROPERTY(double, change_abs_tol, 1.e-3)
  //Change made below by Ansh
  //PROPERTY(double, Tnuc_step, 1)
  PROPERTY(double, Tnuc_step, 0.1)
  /** Relative precision in nucleation temperature */
  //PROPERTY(double, Tnuc_tol_rel, 1.e-3)
  PROPERTY(double, Tnuc_tol_rel, 1.e-4)
  
  PROPERTY(bool, calculate_action, false)

  /**
   * Whether we should check for subcritical transitions, i.e. transitions between phases when one phase is strictly of
   * lower energy than the other. This can occur when a new phase appears near a phase that is not of the highest energy
   * at that temperature.
   */
  PROPERTY(bool, check_subcritical_transitions, false)
};

} // namespace PhaseTracer

#endif // PHASETRACER_TRANSITION_FINDER_HPP_
