#include <cmath>
#include <vector>
#include <limits>
#include <fstream>
#include <chrono>
#include <iostream>
#include <random>
#include <algorithm>
#include "HillProblem.hpp"
#include "ShekelProblem.hpp"
#include "HillProblemFamily.hpp"
#include "ShekelProblemFamily.hpp"
#include <cstdlib>
#include <string>
#include <filesystem>

using namespace std;

// Структура для представления точки (x, y)
struct Point {
    vector<double> x;
    double y;
};

// Абстрактный базовый класс для функций
class FunctionBase {
public:
    virtual double calculate(double x) const = 0;
};

// Конкретная функция для минимизации
class MyFunction : public FunctionBase {
public:
    double calculate(double x) const override {
        return sin(2.0 * x) + (x - 1.0) * (x - 1.0);
    }
};

// Класс для реализации метода поиска минимума
template <typename T>
class Minimizer {
private:
    vector<double> leftBound;
    vector<double> rightBound;
    int iterationCount;
    double epsilon;
    double r;
    const T& function;
    ofstream logFile;
    int exitMainCount;
    int exitTestCount;

    vector<Point> points;

public:
    Minimizer(vector<double> a, vector<double> b, double eps, double r, const T& func)
        : leftBound(a), rightBound(b), epsilon(eps), r(r), function(func),
        exitMainCount(0), exitTestCount(0)
        //exitMainHillCount(0), exitMainShekelCount(0), exitTestHillCount(0), exitTestShekelCount(0)
    {
        logFile.open("minimization_log.txt");
        if (!logFile.is_open()) {
            cerr << "Ошибка открытия файла журнала!" << endl;
        }

        // Добавляем начальные точки
        points.push_back({ leftBound, function.ComputeFunction(leftBound) });
        points.push_back({ rightBound, function.ComputeFunction(rightBound) });

        logFile << "Начальные точки:\n";
        logPoint(points[0]);
        logPoint(points[1]);
    }

    ~Minimizer() {
        logFile.close();
    }

    double calculateY(double x) const {
        std::vector<double> y = { x };
        return function.ComputeFunction(y);
    }

    vector<double> findMinimum() { 
        iterationCount = 0;
        const int maxIterations = 1000;
        int iteration = 0;
        double M, m, yNew;
        vector<double> xNew(leftBound.size());
        double intervalSize = abs(rightBound[0] - leftBound[0]);
        vector<double> actualMinPoint = function.GetOptimumPoint();

        while (iteration < maxIterations) {
            M = calculateMaxSlope();
            m = (M > 0.0) ? r * M : 1.0;

            vector<double> intervalCharacteristics(points.size() - 1);
            for (size_t i = 0; i < points.size() - 1; ++i) {
                intervalCharacteristics[i] = calculateCharacteristic(i, m);
            }

            size_t maxCharacteristicIndex =
                distance(intervalCharacteristics.begin(),
                    max_element(intervalCharacteristics.begin(), intervalCharacteristics.end()));

            for (size_t j = 0; j < leftBound.size(); j++) {
                xNew[j] = 0.5 * (points[maxCharacteristicIndex + 1].x[j] + points[maxCharacteristicIndex].x[j]) -
                    (points[maxCharacteristicIndex + 1].y - points[maxCharacteristicIndex].y) / (2.0 * m);

                if (xNew[j] < leftBound[j])
                    xNew[j] = leftBound[j];

                if (xNew[j] > rightBound[j])
                    xNew[j] = rightBound[j];
            }


            yNew = function.ComputeFunction(xNew);

            points.insert(points.begin() + maxCharacteristicIndex + 1, { xNew, yNew });


            logFile << "Итерация " << iteration + 1 << ": ";
            logPoint(points[maxCharacteristicIndex + 1]);



            if (abs(points[maxCharacteristicIndex + 1].x[0] - points[maxCharacteristicIndex].x[0]) <= epsilon * intervalSize / 2) {
                exitMainCount++;
                break;
            }
            if (abs(xNew[0] - actualMinPoint[0]) <= epsilon * intervalSize / 2.0) {
                exitTestCount++;
                break;
            }

            ++iteration;
        }


        size_t minIndex = 0;
        double minY = points[0].y;
        for (size_t i = 1; i < points.size(); ++i) {
            if (points[i].y < minY) {
                minY = points[i].y;
                minIndex = i;
            }
        }
        iterationCount = iteration;
        return points[minIndex].x; // Возвращаем  точку  минимума
    }
    int GetIterationCount() const { return iterationCount; }
    int GetExitMainCount() const { return exitMainCount; } // Геттер для счетчика 1
    int GetExitTestCount() const { return exitTestCount; } // Геттер для счетчика 2

private:


    double calculateMaxSlope() const {
        double maxSlope = 0.0;
        for (size_t i = 1; i < points.size(); ++i) {
            if (abs(points[i].x[0] - points[i - 1].x[0]) > 1e-9) { //защита от деления на ноль
                double slope = abs((points[i].y - points[i - 1].y) /
                    (points[i].x[0] - points[i - 1].x[0]));
                maxSlope = max(maxSlope, slope);
            }
        }
        return maxSlope;
    }

    double calculateCharacteristic(size_t index, double m) const {
        return m * (points[index + 1].x[0] - points[index].x[0]) +
            pow(points[index + 1].y - points[index].y, 2) /
            (m * (points[index + 1].x[0] - points[index].x[0])) -
            2.0 * (points[index + 1].y + points[index].y);
    }

    void logPoint(const Point& point) {
        logFile << "x: ";
        for (double val : point.x) {
            logFile << val << " ";
        }
        logFile << ", y: " << point.y << endl;
    }
};

int main() {
    setlocale(LC_ALL, "Russian");
    
    int exitMainHillCount = 0;
    int exitMainShekelCount = 0;
    int exitTestHillCount = 0;
    int exitTestShekelCount = 0;
    double epsilon = 0.0, r = 0.0;
    std::cout << "Введите точность(> 0)" << endl;
    std::cin >> epsilon;
    std::cout << "Введите r(2 < r < 4)" << endl;
    std::cin >> r;
    
    // Количество случайных функций для тестирования
    int numTests = 1000;
    std::cout << "Введите количество случайных функций для тестирования" << endl;
    std::cin >> numTests;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> hillDist(0, NUM_HILL_PROBLEMS - 1);
    uniform_int_distribution<> shekelDist(0, NUM_SHEKEL_PROBLEMS - 1);
    uniform_real_distribution<> pointDist(0.0, 10.0); // Диапазон для аргумента функций
    MyFunction func;

    THillProblemFamily hillFamily;
    TShekelProblemFamily shekelFamily;
    ofstream dataFile("plot_data.txt");
    if (!dataFile.is_open()) {
        cerr << "Ошибка открытия файла данных!" << endl;
        return 1;
    }

    for (int i = 0; i < numTests; ++i) {
        int hillIndex = hillDist(gen);
        THillProblem* hill = dynamic_cast<THillProblem*>(hillFamily[hillIndex]);
        double actualMinHill = hill->GetOptimumValue(); //минимум из файла
        vector<double> actualMinPointHill = hill->GetOptimumPoint();
        cout << "\nHill Problem " << hillIndex << endl;
        cout << "Фактический минимум (из файла): " << actualMinHill << 
            " в точке x = " << actualMinPointHill[0] << endl;

        vector<double> hillLeft = { 0.0 };
        vector<double> hillRight = { 1.0 };
        Minimizer<THillProblem> hillMinimizer(hillLeft, hillRight, epsilon, r , *hill);
        
        auto startTime = chrono::high_resolution_clock::now();
        vector<double> calculatedMinPointHill = hillMinimizer.findMinimum();
        auto endTime = chrono::high_resolution_clock::now();
        chrono::duration<double> durationHill = endTime - startTime;

        // Выводим только результаты Minimizer
        cout << "Посчитанный минимум: " << hill->ComputeFunction(calculatedMinPointHill) << " в точке x = " << calculatedMinPointHill[0] << endl;
        cout << "Количество итераций: " << hillMinimizer.GetIterationCount() << endl;
        cout << "Время выполнения: " << durationHill.count() * 1000 << " мс" << endl;
        exitMainHillCount += hillMinimizer.GetExitMainCount();
        exitTestHillCount += hillMinimizer.GetExitTestCount();
        dataFile << i + 1 << " " << hillMinimizer.GetIterationCount() << endl;

        int shekelIndex = shekelDist(gen);
        TShekelProblem* shekel = dynamic_cast<TShekelProblem*>(shekelFamily[shekelIndex]);

        double actualMinShekel = shekel->GetOptimumValue(); //минимум из файла
        vector<double> actualMinPointShekel = shekel->GetOptimumPoint();

        cout << "\nShekel Problem " << shekelIndex << endl;
        cout << "Фактический минимум (из файла): " << actualMinShekel << 
            " в точке x = " << actualMinPointShekel[0] << endl;


        vector<double> shekelLeft = { 0.0 };
        vector<double> shekelRight = { 10.0 };
        Minimizer<TShekelProblem> shekelMinimizer(shekelLeft, shekelRight, epsilon, r, *shekel);

        startTime = chrono::high_resolution_clock::now();
        vector<double> calculatedMinPointShekel = shekelMinimizer.findMinimum();
        endTime = chrono::high_resolution_clock::now();
        chrono::duration<double> durationShekel = endTime - startTime;

        cout << "Посчитанный минимум: " << shekel->ComputeFunction(calculatedMinPointShekel) << " в точке x = " << calculatedMinPointShekel[0] << endl;
        cout << "Количество итераций: " << shekelMinimizer.GetIterationCount() << endl;
        cout << "Время выполнения: " << durationShekel.count() * 1000 << " мс" << endl;
        exitMainShekelCount += shekelMinimizer.GetExitMainCount();
        exitTestShekelCount += shekelMinimizer.GetExitTestCount();
        dataFile << i + 1 << " " << shekelMinimizer.GetIterationCount() << endl;
        cout << "--------------------" << endl;
    }
    ofstream statsFile("stats.txt");
    if (!statsFile.is_open()) {
        cerr << "Ошибка открытия stats.txt" << endl;
        return 1;
    }
    dataFile.close();
    statsFile << exitMainHillCount << " " << exitMainShekelCount << " " <<
    exitTestHillCount << " " << exitTestShekelCount << endl;
    statsFile << epsilon << " " << r << " " << numTests << endl;
    statsFile.close();
    std::cout << "Полный путь к файлу данных: " << std::filesystem::absolute("plot_data.txt").string() << std::endl;
    string command = "python ../../plot_graph.py";
    int result = system(command.c_str());
    if (result != 0) {
       cerr << "Ошибка при запуске Python скрипта! Код ошибки: " << result << endl;
       return 1;
    }
    return 0;
}