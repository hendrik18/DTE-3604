#include <parametrics/gmpcurve.h>
#include <core/containers/gmdvector.h>

// MyB_spline class definition inheriting from GMlib::PCurve
class MyB_spline : public GMlib::PCurve<float,3> {
    GM_SCENEOBJECT(MyB_spline)

public:
    // Constructor 1: Initialize with given control points
    MyB_spline(const GMlib::DVector<GMlib::Vector<float,3>>& c);
    
    // Constructor 2: Approximate a set of points using least squares
    MyB_spline(const GMlib::DVector<GMlib::Vector<float,3>>& p, int n);

protected:
    // Evaluate the curve at parameter t with d derivatives
    void eval(float t, int d, bool left = true) const override;
    
    // Return the first valid parameter value (avoid repeated knots at start)
    float getStartP() const override {
        return _knotVector[2]; // First non-repeated knot
    }
    
    // Return the last valid parameter value (avoid repeated knots at end)
    float getEndP() const override {
        return _knotVector[_knotVector.getDim() - 3]; // Last non-repeated knot
    }
    
    // Check if the curve is closed (always false for this B-spline)
    bool isClosed() const override {
        return false;
    }

private:
    GMlib::DVector<GMlib::Vector<float,3>> _controlPoints; // Control points defining the curve
    GMlib::DVector<float> _knotVector; // Knot vector defining parameter spacing

    // Generate a uniform knot vector for a 2nd-degree B-spline
    void generateKnotVector();
    
    // Compute control points using least squares fitting
    void leastSquaresFit(const GMlib::DVector<GMlib::Vector<float,3>>& p, int n);
    
    // Evaluate the basis function at index i, degree k, and parameter t
    float evaluateBasis(int i, int k, float t) const;
};

// Constructor: Create a B-spline from predefined control points
MyB_spline::MyB_spline(const GMlib::DVector<GMlib::Vector<float,3>>& c)
    : _controlPoints(c) {
    generateKnotVector(); // Generate knot vector for this set of control points
}

// Constructor: Approximate a given set of points using least squares fitting
MyB_spline::MyB_spline(const GMlib::DVector<GMlib::Vector<float,3>>& p, int n) {
    leastSquaresFit(p, n); // Fit control points using least squares
    generateKnotVector(); // Generate a corresponding knot vector
}

// Generate a uniform knot vector for a 2nd-degree (quadratic) B-spline
void MyB_spline::generateKnotVector() {
    int n = _controlPoints.getDim(); // Number of control points
    int k = 2; // Degree of the B-spline (quadratic)
    int m = n + k + 1; // Number of knots

    _knotVector.setDim(m);
    
    // First k+1 knots are set to 0
    for (int i = 0; i <= k; ++i) {
        _knotVector[i] = 0.0f;
    }
    
    // Middle knots are uniformly spaced
    for (int i = k + 1; i < m - (k + 1); ++i) {
        _knotVector[i] = static_cast<float>(i - k);
    }
    
    // Last k+1 knots are set to the maximum value
    float maxValue = static_cast<float>(m - 2 * (k + 1) + 1);
    for (int i = m - (k + 1); i < m; ++i) {
        _knotVector[i] = maxValue;
    }
}

// Compute control points using least squares fitting
void MyB_spline::leastSquaresFit(const GMlib::DVector<GMlib::Vector<float,3>>& p, int n) {
    int m = p.getDim(); // Number of input points
    int k = 2; // B-spline degree

    _controlPoints.setDim(n); // Allocate space for control points

    GMlib::DMatrix<float> N(m, n, 0.0f); // Basis function matrix
    GMlib::DVector<float> weights(m, 1.0f); // Weights for least squares

    // Compute basis functions for each input point
    for (int i = 0; i < m; ++i) {
        float t = static_cast<float>(i) / (m - 1); // Normalized parameter value
        int span = -1; // Find knot span

        // Find the knot span for the current parameter value t
        // The knot span is the interval [u_i, u_{i+1}) where the parameter t lies
        // This determines which basis functions are non-zero at t
        for (int j = k; j < n; ++j) {
            if (t >= _knotVector[j] && t < _knotVector[j + 1]) {
                span = j;
                break;
            }
        }
        if (span == -1) span = n - 1; // Handle boundary case

        GMlib::DVector<float> N_i(k + 1, 0.0f); // Basis function values
        N_i[0] = 1.0f;
        
        // Compute basis functions using recursion
        for (int d = 1; d <= k; ++d) {
            float left = t - _knotVector[span + 1 - d];
            float right = _knotVector[span + d] - t;
            float saved = 0.0f;
            
            // Update basis functions for the current degree
            for (int r = 0; r < d; ++r) {
                float temp = N_i[r] / (_knotVector[span + r + 1] - _knotVector[span + 1 - d + r]);
                N_i[r] = saved + right * temp;
                saved = left * temp;
            }
            N_i[d] = saved;
        }

        // Store the computed basis functions in the matrix N
        for (int j = 0; j <= k; ++j) {
            N[i][span - k + j] = N_i[j];
        }
    }
}

// Evaluate basis function recursively using the Cox–de Boor formula
float MyB_spline::evaluateBasis(int i, int degree, float t) const {
    // Base case: if degree is 0, check if t is within the knot span [u_i, u_{i+1})
    if (degree == 0) {
        if( (_knotVector[i] <= t && t < _knotVector[i+1]) ||
            (t == _knotVector[_knotVector.getDim()-1] && i == _controlPoints.getDim()-1) )
            return 1.0f;
        else
            return 0.0f;    }
    
    // Calculate the denominator for the first term
    float denom1 = _knotVector[i+degree] - _knotVector[i];
    // Calculate the first term of the Cox–de Boor recursion formula
    float term1  = (denom1 != 0.0f) ? (t - _knotVector[i]) / denom1 * evaluateBasis(i, degree - 1, t) : 0.0f;

    // Calculate the denominator for the second term
    float denom2 = _knotVector[i+degree+1] - _knotVector[i+1];
    // Calculate the second term of the Cox–de Boor recursion formula
    float term2  = (denom2 != 0.0f) ? (_knotVector[i+degree+1] - t) / denom2 * evaluateBasis(i+1, degree - 1, t) : 0.0f;

    // Return the sum of the two terms
    return term1 + term2;
}

// Evaluate the curve at parameter t using correct basis function evaluation
void MyB_spline::eval(float t, int d, bool left) const {
    this->_p.setDim(d+1);
    this->_p[0] = GMlib::Vector<float,3>(0.0f, 0.0f, 0.0f);
    
    int n = _controlPoints.getDim();
    int degree = 2; // B-spline degree
    
    // Sum over each control point multiplied by its corresponding basis function value.
    for (int i = 0; i < n; ++i) {
        float basisVal = evaluateBasis(i, degree, t);
        this->_p[0] += basisVal * _controlPoints[i];
    }
}