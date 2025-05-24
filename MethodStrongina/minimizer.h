#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <cmath>
#include <vector>
#include <limits>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

// Structure for representing a point of minimization
struct Point {
    vector<double> x;
    double y;
};

// Interface for the target function
class FunctionInterface {
public:
    // Computes the function value at point x
    virtual double ComputeFunction(const vector<double>& x) const = 0;
    // Returns the known optimum point
    virtual vector<double> GetOptimumPoint() const = 0;
    virtual ~FunctionInterface() {}
};

//
// Peano curve function declarations
//
extern int n1, nexp, l, iq, iu[10], iv[10]; // Global variables used by mapd and node

void mapd( double x, int m, double* y, int n, int key );
void node ( int is ); // Helper for mapd, could be static in .cpp if not called elsewhere
void resetMappingGlobals(int n);
vector<double> peanoMapping(double x, int m, int n, int key);


//
// Minimizer class (template, so definition is in the header)
//
template <typename T>
class Minimizer {
private:
    vector<double> leftBound;
    vector<double> rightBound;
    int iterationCount;
    double epsilon;
    double r_param;
    const T& function;
    ofstream logFile;
    int exitMainCount;
    int exitTestCount;
    vector<Point> points;

    // Flag for using mapping (multidimensional task)
    bool useMapping;
    int mappingOrder;  // order of mapping (m)
    int dimension;     // dimension of the original task
    int mappingKey;    // key for mapd function

public:
    // Constructor for one-dimensional case (without mapping)
    Minimizer(vector<double> a, vector<double> b, double eps, double r_val, const T& func)
        : leftBound(a), rightBound(b), epsilon(eps), r_param(r_val), function(func),
          exitMainCount(0), exitTestCount(0), useMapping(false)
    {
        logFile.open("minimization_log.txt", ios::out);
        if (!logFile.is_open()) {
            cerr << "Ошибка открытия файла журнала!" << endl;
        }
        // Initialize points
        points.push_back({ leftBound, function.ComputeFunction(leftBound) });
        points.push_back({ rightBound, function.ComputeFunction(rightBound) });
        logFile << "Начальные точки:\n";
        logPoint(points[0]);
        logPoint(points[1]);
    }

    // Constructor for multidimensional case (with mapping)
    Minimizer(vector<double> a, vector<double> b, double eps, double r_val, const T& func,
              int mappingOrder_, int dimension_, int mappingKey_)
        : leftBound(a), rightBound(b), epsilon(eps), r_param(r_val), function(func),
          exitMainCount(0), exitTestCount(0),
          useMapping(true), mappingOrder(mappingOrder_), dimension(dimension_), mappingKey(mappingKey_)
    {
        logFile.open("minimization_log.txt", ios::out);
        if (!logFile.is_open()) {
            cerr << "Ошибка открытия файла журнала!" << endl;
        }
        // For mapping, bounds are set for parameter x ∈ [0,1]
        vector<double> initPoint_a = { leftBound[0] }; // Ensure it's a vector
        points.push_back({ initPoint_a, calculateY(leftBound[0]) });
        vector<double> initPoint_b = { rightBound[0] }; // Ensure it's a vector
        points.push_back({ initPoint_b, calculateY(rightBound[0]) });
        logFile << "Начальные точки (одномерные для отображения):\n";
        logPoint(points[0]);
        logPoint(points[1]);
    }

    ~Minimizer() {
        if (logFile.is_open()) {
            logFile.close();
        }
    }

    // Calculate function value at point x
    double calculateY(double x_val) const { // Renamed x to x_val to avoid conflict
        if (!useMapping) {
            vector<double> pt = { x_val };
            return function.ComputeFunction(pt);
        } else {
            // Transform x_val via Peano curve to a point in ℝⁿ
            vector<double> mappedPoint = peanoMapping(x_val, mappingOrder, dimension, mappingKey);
            return function.ComputeFunction(mappedPoint);
        }
    }

    // Method to find the minimum
    vector<double> findMinimum() {
        iterationCount = 0;
        const int maxIterations = 10000; // Consider making this configurable
        int iteration = 0;
        double M, mVal, yNew;
        vector<double> xNew_vec(1); // For 1D search space, xNew is scalar but stored in vector<double> for Point
        double intervalSize = abs(rightBound[0] - leftBound[0]);

        // For multidimensional case, Hölder condition: alpha = 1/dimension, else alpha = 1.
        double alpha = useMapping ? (1.0 / static_cast<double>(dimension)) : 1.0;

        // Get actual optimum point (from file or analytically)
        vector<double> actualMinPoint = function.GetOptimumPoint();

        while (iteration < maxIterations) {
            M = calculateMaxSlope();
            mVal = (M > 0.0) ? r_param * M : 1.0;

            vector<double> intervalCharacteristics(points.size() - 1);
            for (size_t i = 0; i < points.size() - 1; ++i) {
                intervalCharacteristics[i] = calculateCharacteristic(i, mVal);
            }

            size_t maxCharacteristicIndex =
                distance(intervalCharacteristics.begin(),
                         max_element(intervalCharacteristics.begin(), intervalCharacteristics.end()));

            // Calculate new point as a one-dimensional value
            double x_left = points[maxCharacteristicIndex].x[0];
            double x_right = points[maxCharacteristicIndex + 1].x[0];
            double new_x_scalar = 0.5 * (x_left + x_right) - (points[maxCharacteristicIndex + 1].y - points[maxCharacteristicIndex].y) / (2.0 * mVal);

            // Clamp new_x_scalar to bounds
            if (new_x_scalar < leftBound[0])
                new_x_scalar = leftBound[0];
            if (new_x_scalar > rightBound[0])
                new_x_scalar = rightBound[0];

            xNew_vec[0] = new_x_scalar;
            yNew = calculateY(new_x_scalar);
            points.insert(points.begin() + maxCharacteristicIndex + 1, { xNew_vec, yNew });

            logFile << "Итерация " << iteration + 1 << ": ";
            logPoint(points[maxCharacteristicIndex + 1]);

            double threshold = epsilon * intervalSize;
            
            if (!useMapping) {
                // For 1D, the comparison is direct
                if (actualMinPoint.empty()) { // Handle case where GetOptimumPoint might return empty
                     // Fallback or alternative stop condition if actualMinPoint is not available
                    double currentInterval = abs(points[maxCharacteristicIndex + 1].x[0] - points[maxCharacteristicIndex].x[0]);
                     if (currentInterval <= epsilon) { // A simple interval based stop
                        exitTestCount++; // Or some other counter
                        break;
                     }
                } else if (pow(abs(new_x_scalar - actualMinPoint[0]), alpha) <= threshold) {
                    exitTestCount++;
                    break;
                }
            } else {
                // Stopping condition with Hölder consideration
                double currentInterval = abs(points[maxCharacteristicIndex + 1].x[0] - points[maxCharacteristicIndex].x[0]);
                if (currentInterval <= threshold * sqrt(5.0)) { 
                    exitMainCount++;

                    break;
                }
            }
            ++iteration;
        }

        // Find the point with the minimum value
        size_t minIndex = 0;
        double minY = points[0].y;
        for (size_t i = 1; i < points.size(); ++i) {
            if (points[i].y < minY) {
                minY = points[i].y;
                minIndex = i;
            }
        }
        iterationCount = iteration;
        if (useMapping)
            return peanoMapping(points[minIndex].x[0], mappingOrder, dimension, mappingKey);
        else
            return points[minIndex].x;
    }

    int GetIterationCount() const { return iterationCount; }
    int GetExitMainCount() const { return exitMainCount; }
    int GetExitTestCount() const { return exitTestCount; }

private:
    double calculateMaxSlope() const {
        double maxSlope = 0.0;
        for (size_t i = 1; i < points.size(); ++i) {
            double diffX = abs(points[i].x[0] - points[i - 1].x[0]);
            if (diffX > 1e-9) { // Avoid division by zero
                double slope = abs((points[i].y - points[i - 1].y) / diffX);
                maxSlope = std::max(maxSlope, slope);
            }
        }
        return maxSlope;
    }

    double calculateCharacteristic(size_t index, double mVal) const {
        double interval = points[index + 1].x[0] - points[index].x[0];
        if (abs(interval) < 1e-9) { // Avoid division by zero if interval is too small
            return -numeric_limits<double>::infinity(); // Or some other very small number
        }
        return mVal * interval + pow(points[index + 1].y - points[index].y, 2) / (mVal * interval)
               - 2.0 * (points[index + 1].y + points[index].y);
    }

    void logPoint(const Point& point) {
        if (!logFile.is_open()) return;
        logFile << "x: ";
        for (double val : point.x) {
            logFile << val << " ";
        }
        logFile << ", y: " << point.y << endl;
    }
};

#endif // MINIMIZER_H