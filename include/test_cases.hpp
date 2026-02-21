#ifndef TEST_CASES_HPP
#define TEST_CASES_HPP

#include <cmath>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

using namespace dealii;

namespace WaveEquationProject {

// =========================================================================
// SCENARIO 1: MANUFACTURED SOLUTION (Dim-Independent)
// =========================================================================
// Description:
// Uses a smooth analytical function to verify the order of convergence.
// The solution is defined as a standing wave product of sines:
// u(x,t) = sin(pi*t) * product(sin(pi*x_i)) for i=0..dim-1
//
// In 2D: u = sin(pi*t) * sin(pi*x) * sin(pi*y)
// In 3D: u = sin(pi*t) * sin(pi*x) * sin(pi*y) * sin(pi*z)
// =========================================================================

/**
 * @brief Exact analytical solution u_ex(x,t).
 */
template <int dim> class ExactSolution_Case1 : public Function<dim> {
public:
  ExactSolution_Case1() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int /*component*/ = 0) const override {
    double time = this->get_time();

    // Time part: sin(pi * t)
    double result = std::sin(numbers::PI * time);

    // Spatial part: product of sin(pi * x_i) for all dimensions
    for (unsigned int d = 0; d < dim; ++d) {
      result *= std::sin(numbers::PI * p[d]);
    }
    return result;
  }
};

/**
 * @brief Initial Displacement u(x,0).
 * Since u contains sin(pi*t), at t=0 the displacement is uniformly zero.
 */
template <int dim> class InitialValuesU_Case1 : public Function<dim> {
public:
  virtual double value(const Point<dim> & /*p*/,
                       const unsigned int = 0) const override {
    return 0.0;
  }
};

/**
 * @brief Initial Velocity v(x,0) = du/dt(x,0).
 * Derivative: du/dt = pi * cos(pi*t) * spatial_part.
 * At t=0:     v(0)  = pi * 1.0 * spatial_part.
 */
template <int dim> class InitialValuesV_Case1 : public Function<dim> {
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int = 0) const override {
    double result = numbers::PI; // pi * cos(0)

    // Multiply by spatial sines
    for (unsigned int d = 0; d < dim; ++d) {
      result *= std::sin(numbers::PI * p[d]);
    }
    return result;
  }
};

/**
 * @brief Source Term f(x,t).
 * Calculated by substituting the exact solution into the wave equation:
 * f = d^2u/dt^2 - Laplacian(u)
 *
 * Derivation:
 * 1. Time derivative: d^2u/dt^2 = -pi^2 * u
 * 2. Laplacian:       Lapl(u)   = -dim * pi^2 * u  (each dimension adds a -pi^2
 * factor)
 * 3. Result:          f = (-pi^2 * u) - (-dim * pi^2 * u)
 * = (dim - 1) * pi^2 * u
 *
 * In 2D: f = 1 * pi^2 * u
 * In 3D: f = 2 * pi^2 * u
 */
template <int dim> class RightHandSide_Case1 : public Function<dim> {
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int = 0) const override {
    double time = this->get_time();

    // Calculate the coefficient (dim - 1) * pi^2
    double coefficient = (dim - 1.0) * numbers::PI * numbers::PI;

    // Reconstruct u(x,t) parts
    double spatial_part = 1.0;
    for (unsigned int d = 0; d < dim; ++d) {
      spatial_part *= std::sin(numbers::PI * p[d]);
    }

    return coefficient * std::sin(numbers::PI * time) * spatial_part;
  }
};

// SCENARIO 2: GAUSSIAN PULSE
// Physical Demo: A single pulse starting from the center.
// Demonstrates propagation and reflection against Dirichlet boundaries.

/**
 * @brief Initial Displacement: Gaussian bell curve centered at (0.5, 0.5).
 */
template <int dim> class InitialValuesU_Case2 : public Function<dim> {
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int = 0) const override {
    const Point<dim> center(0.5, 0.5);
    const double sigma = 0.05; // Controls the width of the pulse1
                               // Gaussian formula: exp(-r^2 / sigma^2)
    return std::exp(-(p.distance_square(center)) / (sigma * sigma));
  }
};

/**
 * @brief Initial Velocity: Zero (starts from rest).
 */
template <int dim> class InitialValuesV_Case2 : public Function<dim> {
public:
  virtual double value(const Point<dim> & /*p*/,
                       const unsigned int = 0) const override {
    return 0.0;
  }
};

/**
 * @brief Forcing Term: Zero (free evolution).
 */
template <int dim> class RightHandSide_Case2 : public Function<dim> {
public:
  virtual double value(const Point<dim> & /*p*/,
                       const unsigned int = 0) const override {
    return 0.0;
  }
};

// SCENARIO 3: DOUBLE PULSE COLLISION
// Two pulses starting at opposite corners to demonstrate superposition.
template <int dim> class InitialValuesU_Case3 : public Function<dim> {
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int = 0) const override {
    const Point<dim> center1(0.3, 0.3);
    const Point<dim> center2(0.7, 0.7);
    const double sigma = 0.08;

    // Linear superposition of two Gaussians
    double pulse1 = std::exp(-(p.distance_square(center1)) / (sigma * sigma));
    double pulse2 = std::exp(-(p.distance_square(center2)) / (sigma * sigma));

    return pulse1 + pulse2;
  }
};

// SCENARIO 4: RAINDROP SYMPHONY
// Multiple random drops with different amplitudes (positive/negative).
// Creates complex interference patterns.
template <int dim> class InitialValuesU_Raindrops : public Function<dim> {
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int = 0) const override {

    // Drop definition: {x, y, amplitude, sigma}
    const double drops[5][4] = {{0.50, 0.50, 1.0, 0.05},
                                {0.20, 0.20, -0.8, 0.04},
                                {0.80, 0.30, 0.6, 0.03},
                                {0.30, 0.70, -0.5, 0.03},
                                {0.75, 0.75, 0.7, 0.04}};

    double total_displacement = 0.0;

    for (int i = 0; i < 5; ++i) {
      const Point<dim> center(drops[i][0], drops[i][1]);
      const double amplitude = drops[i][2];
      const double sigma = drops[i][3];

      // sunm of Gaussians
      total_displacement +=
          amplitude * std::exp(-(p.distance_square(center)) / (sigma * sigma));
    }

    return total_displacement;
  }
};

// SCENARIO 5: RING IMPLOSION
// A Gaussian ring that collapses towards the center, creating a singularity.
template <int dim> class InitialValuesU_Ring : public Function<dim> {
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int = 0) const override {
    const double r = p.distance(Point<dim>(0.5, 0.5)); // Distance from center
    const double radius_ring = 0.35;                   // Radius of the ring
    const double sigma = 0.05;                         // Thickness

    // Gaussian function of the radius: peak at r == radius_ring
    return std::exp(-std::pow(r - radius_ring, 2) / (sigma * sigma));
  }
};

// =========================================================================
// SCENARIO 6: TIME-DEPENDENT BCs (The Pumping Wall)
// Domain: [0,1]^2.
// Boundary Condition:
// - Left wall (x=0): Moves with sin(10*t) * sin(pi*y)
// - Other walls: Fixed at 0
// =========================================================================

template <int dim> class BoundaryFunction_MovingWall : public Function<dim> {
public:
  BoundaryFunction_MovingWall() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int /*component*/ = 0) const override {
    double time = this->get_time();

    // Apply movement on left wall ONLY (x very close to 0)
    if (std::abs(p[0]) < 1e-10) {
      // Sinusoidal motion in time, moduled by a spatial sine
      // to keep angles still (avoid discontinuities)
      return 0.5 * std::sin(5.0 * time) * std::sin(numbers::PI * p[1]);
    }

    // Every other wall stays still
    return 0.0;
  }
};

/**
 * @brief Initial Velocity for Scenario 6.
 * Matches the time derivative of the moving wall at t=0 to prevent shocks.
 * Wall Position: g(t) = 0.5 * sin(5.0 * t) * sin(pi * y)
 * Wall Velocity: g'(t) = 2.5 * cos(5.0 * t) * sin(pi * y)
 * At t=0:        g'(0) = 2.5 * 1.0 * sin(pi * y)
 */
template <int dim> class InitialValuesV_MovingWall : public Function<dim> {
public:
  virtual double value(const Point<dim> &p,
                       const unsigned int /*component*/ = 0) const override {

    if (std::abs(p[0]) < 1e-10) {

      return 2.5 * std::sin(numbers::PI * p[1]);
    }
    return 0.0;
  }
};

} // namespace WaveEquationProject

#endif
