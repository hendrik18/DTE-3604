#ifndef CLOSED_SUBDIVISION_CURVE_H
#define CLOSED_SUBDIVISION_CURVE_H

#include <parametrics/gmpcurve.h>
#include <core/containers/gmdvector.h>

// ClosedSubdivisionCurve class definition inheriting from GMlib::PCurve
class ClosedSubdivisionCurve : public GMlib::PCurve<float, 3> {
  GM_SCENEOBJECT(ClosedSubdivisionCurve)

public:
  // Constructor: Initialize the curve with control points and subdivision degree
  ClosedSubdivisionCurve(const GMlib::DVector<GMlib::Vector<float, 3>> &controlPts, int degree)
      : _controlPoints(controlPts), _degree(degree) {

    // Constrain the parametric domain to [0, 1]
    this->setDomain(0.0f, 1.0f);

    // Compute subdivided points
    laneRiesenfeldSubdivision();
  }

  // Destructor (default)
  ~ClosedSubdivisionCurve() override = default;

  // PCurve interface overrides
  void eval(float t, int d, bool left = true) const override;
  float getStartP() const override { return 0.0f; }
  float getEndP() const override { return 1.0f; }
  bool isClosed() const override { return true; } // Mark as a closed curve

private:
  GMlib::DVector<GMlib::Vector<float, 3>> _controlPoints; // Original control polygon
  mutable GMlib::DVector<GMlib::Vector<float, 3>> _subdividedPoints; // Subdivided points
  int _degree; // Number of subdivision iterations

  // Perform Lane-Riesenfeld subdivision to refine the curve
  void laneRiesenfeldSubdivision();
};

/*!
 *  eval(float t, int d, bool left) const
 *
 *  - Maps parameter t ∈ [0,1] to an index in _subdividedPoints.
 *  - Interpolates linearly between adjacent points for a smooth curve.
 *  - Approximates the first derivative using finite differences if requested.
 */
void ClosedSubdivisionCurve::eval(float t, int d, bool /*left*/) const {

  // Ensure _p has space for position and derivatives
  this->_p.setDim(d + 1);

  // Map t ∈ [0,1] to an index in _subdividedPoints
  float scaled_t = t * (_subdividedPoints.getDim() - 1);
  int index = static_cast<int>(std::floor(scaled_t)) % _subdividedPoints.getDim();
  float alpha = scaled_t - index; // Fractional part for interpolation

  // Fetch adjacent points for interpolation
  GMlib::Vector<float, 3> p1 = _subdividedPoints[index];
  GMlib::Vector<float, 3> p2 = _subdividedPoints[(index + 1) % _subdividedPoints.getDim()];

  // Linearly interpolate between p1 and p2
  this->_p[0] = (1.0f - alpha) * p1 + alpha * p2;

  // Approximate the first derivative if d > 0 (central difference)
  if (d > 0) {
    int next = (index + 1) % _subdividedPoints.getDim();
    int prev = (index - 1 + _subdividedPoints.getDim()) % _subdividedPoints.getDim();
    this->_p[1] = (_subdividedPoints[next] - _subdividedPoints[prev]) * 0.5f;
  }
}

/*!
 *  laneRiesenfeldSubdivision()
 *
 *  - Implements the Lane-Riesenfeld subdivision algorithm for **closed** curves.
 *  - Inserts **midpoints** and applies **averaging passes** to generate a smooth result.
 *  - Ensures closure by explicitly setting the last point equal to the first.
 */
void ClosedSubdivisionCurve::laneRiesenfeldSubdivision() {

  // Start with the original control points
  GMlib::DVector<GMlib::Vector<float, 3>> points = _controlPoints;

  // Perform _degree_ iterations of Lane-Riesenfeld subdivision
  for (int iter = 0; iter < _degree; ++iter) {

    int numPoints = points.getDim();
    GMlib::DVector<GMlib::Vector<float, 3>> newPoints(2 * numPoints, GMlib::Vector<float, 3>(0.0f, 0.0f, 0.0f));

    // Step 1: Insert midpoints
    for (int i = 0; i < numPoints; ++i) {
      newPoints[2 * i] = points[i]; // Keep original point
      int nxt = (i + 1) % numPoints; // Wrap around for closed curve
      newPoints[2 * i + 1] = (points[i] + points[nxt]) * 0.5f; // Compute midpoint
    }

    // Step 2: Perform smoothing passes
    for (int avg = 1; avg < _degree; ++avg) {
      GMlib::DVector<GMlib::Vector<float, 3>> smoothedPoints(newPoints.getDim(), GMlib::Vector<float, 3>(0.0f, 0.0f, 0.0f));
      
      for (int i = 0; i < newPoints.getDim(); ++i) {
        int prev = (i - 1 + newPoints.getDim()) % newPoints.getDim();
        smoothedPoints[i] = (newPoints[i] + newPoints[prev]) * 0.5f; // Average with previous
      }
      
      newPoints = smoothedPoints; // Update newPoints after each pass
    }

    points = newPoints; // Update points after each iteration
  }

  // Store final refined points
  _subdividedPoints = points;

  // Ensure closure: explicitly set the last point to match the first
  if (_subdividedPoints.getDim() > 1) {
    _subdividedPoints[_subdividedPoints.getDim() - 1] = _subdividedPoints[0];
  }
}

#endif // CLOSED_SUBDIVISION_CURVE_H
