#ifndef TORUS_KNOT_H
#define TORUS_KNOT_H

#include <parametrics/gmpcurve.h>
#include <cmath>

// TorusKnot class definition inheriting from GMlib::PCurve
class TorusKnot : public GMlib::PCurve<float,3> {
    GM_SCENEOBJECT(TorusKnot)

  public:
    // Default constructor (no parameters needed)
    TorusKnot() {}

  protected:
    /*!
     *  eval(t, d, left):
     *  - Evaluates the torus knot at parameter `t`.
     *  - Computes position, first derivative, and second derivative.
     *  - Uses exact mathematical derivatives (no numerical approximation).
     */
    void eval(float t, int d, bool /*left*/ = true) const override {

      // Ensure _p has room for up to d derivatives (0 => just position)
      this->_p.setDim(d + 1);

      // Parameters defining the Torus Knot (p twists, q loops)
      float R = 2.0f;  // Major radius offset
      int   p = 2;      // Twists around torus axis
      int   q = 3;      // Loops through torus hole

      // 1. Compute the Position
      float x = ( R + std::cos(q * t) ) * std::cos(p * t);
      float y = ( R + std::cos(q * t) ) * std::sin(p * t);
      float z = std::sin(q * t);

      this->_p[0] = { x, y, z };

      // 2. Compute the First Derivative, if requested
      if(d > 0) {
        // dx/dt
        float dx = -p * (R + std::cos(q * t)) * std::sin(p * t) - q * std::sin(q * t) * std::cos(p * t);

        // dy/dt
        float dy = p * (R + std::cos(q * t)) * std::cos(p * t) - q * std::sin(q * t) * std::sin(p * t);

        // dz/dt
        float dz = q * std::cos(q * t);

        this->_p[1] = {dx, dy, dz};
      }

      // 3. Compute the Second Derivative, if requested
      if(d > 1) {
        // The second derivative components are calculated analytically.
        // Instead of blindly differentiating, they are derived carefully in steps.

        // Initialize xpp to 0.0f
        float xpp = 0.0f;
        {
          // Calculate partA and partB for xpp
          float partA = -p * (p * (R + std::cos(q * t)) * std::cos(p * t) - q * std::sin(q * t) * std::sin(p * t));
          float partB = -q * (q * std::cos(q * t) * std::cos(p * t) - p * std::sin(q * t) * std::sin(p * t));
          // Sum partA and partB to get xpp
          xpp = partA + partB;
        }

        // Initialize ypp to 0.0f
        float ypp = 0.0f;
        {
          // Calculate partC and partD for ypp
          float partC = p * (-p * (R + std::cos(q * t)) * std::sin(p * t) - q * std::sin(q * t) * std::cos(p * t));
          float partD = -q * (q * std::cos(q * t) * std::sin(p * t) + p * std::sin(q * t) * std::cos(p * t));
          // Sum partC and partD to get ypp
          ypp = partC + partD;
        }

        // Calculate zpp directly
        float zpp = -q * q * std::sin(q * t);

        // Assign the calculated values to the _p array
        this->_p[2] = {xpp, ypp, zpp};
      }
    }

    /*!
     *  getStartP()
     *  - Defines the start of the parametric domain.
     *  - 0 corresponds to the beginning of the torus knot.
     */
    float getStartP() const override {
      return 0.0f;
    }

    /*!
     *  getEndP()
     *  - Defines the end of the parametric domain.
     *  - The torus knot completes one full cycle at `t = 6π` for (p=2, q=3).
     */
    float getEndP() const override {
      return 6.0f * float(M_PI);
    }

    /*!
     *  isClosed()
     *  - The torus knot is a closed curve that loops back onto itself.
     */
    bool isClosed() const override {
      return true;
    }
};

#endif // TORUS_KNOT_H
