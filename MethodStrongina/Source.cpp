#include <iostream> // For cout, cin
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <cstdlib>     // For system()
#include <filesystem>  // If used directly in main, else can be removed if not.
                       // Not directly used in the provided main, but kept for potential future use.

#include "minimizer.h" // Our new header

// Problem specific headers
#include "HillProblem.hpp"
#include "ShekelProblem.hpp"
#include "HillProblemFamily.hpp"
#include "ShekelProblemFamily.hpp"
#include "grishagin_function.hpp" // Assuming this defines TGrishaginProblem or similar
#include "GrishaginProblemFamily.hpp"

// using namespace std; // Already in minimizer.h, so available here

int main() {
    setlocale(LC_ALL, "Rus");

    int taskType;
    int numTests;
    double epsilon_val, r_val; // Renamed to avoid conflict with Minimizer's internal 'r'

    cout << "Выберите тип задачи:" << endl;
    cout << " 1 - Одномерная (Хилл/Шекеля)" << endl;
    cout << " 2 - Многомерная (Гришагина)" << endl;
    cout << "Ваш выбор: ";
    cin >> taskType;

    cout << "Введите количество функций для тестирования: ";
    cin >> numTests;

    cout << "Введите точность (> 0): ";
    cin >> epsilon_val;
    cout << "Введите параметр r (например, 2.5): ";
    cin >> r_val;

    // File for plotting data
    ofstream dataFile("plot_data.txt", ios::out);
    if (!dataFile.is_open()) {
        cerr << "Ошибка открытия файла plot_data.txt!" << endl;
        return 1;
    }

    // Counters for statistics
    int exitMainTotal = 0;
    int exitTestTotal = 0;
    double totalIterations = 0.0;

    random_device rd;
    mt19937 gen(rd());

    if (taskType == 1) {
        int subChoice;
        cout << "Выберите задачу:" << endl;
        cout << " 1 - Функция Хилла" << endl;
        cout << " 2 - Функция Шекеля" << endl;
        cout << "Ваш выбор: ";
        cin >> subChoice;

        THillProblemFamily hillFamily;
        TShekelProblemFamily shekelFamily;
        uniform_int_distribution<> hillDist(0, hillFamily.GetFamilySize() - 1);
        uniform_int_distribution<> shekelDist(0, shekelFamily.GetFamilySize() - 1);

        for (int i = 0; i < numTests; ++i) {
            if (subChoice == 1) { // Hill
                int index = hillDist(gen);
                THillProblem* hill = dynamic_cast<THillProblem*>(hillFamily[index]);
                if (!hill) { cerr << "Ошибка: hill problem is null" << endl; continue; }
                double actualMin = hill->GetOptimumValue();
                vector<double> actualMinPoint = hill->GetOptimumPoint();
                cout << "\nHill Problem " << index << endl;
                cout << "Фактический минимум (из файла): " << actualMin
                     << " в точке x = " << (actualMinPoint.empty() ? NAN : actualMinPoint[0]) << endl;
                vector<double> a = {0.0};
                vector<double> b = {1.0};
                Minimizer<THillProblem> minimizer(a, b, epsilon_val, r_val, *hill);
                auto startTime = chrono::high_resolution_clock::now();
                vector<double> computedMinPoint = minimizer.findMinimum();
                auto endTime = chrono::high_resolution_clock::now();
                chrono::duration<double> duration = endTime - startTime;
                cout << "Посчитанный минимум: " << hill->ComputeFunction(computedMinPoint)
                     << " в точке x = " << (computedMinPoint.empty() ? NAN : computedMinPoint[0]) << endl;
                cout << "Итераций: " << minimizer.GetIterationCount() << ", время: "
                     << duration.count() * 1000 << " мс" << endl;
                totalIterations += minimizer.GetIterationCount();
                exitMainTotal += minimizer.GetExitMainCount();
                exitTestTotal += minimizer.GetExitTestCount();
                dataFile << i + 1 << " " << minimizer.GetIterationCount() << endl;
            } else { // Shekel
                int index = shekelDist(gen);
                TShekelProblem* shekel = dynamic_cast<TShekelProblem*>(shekelFamily[index]);
                 if (!shekel) { cerr << "Ошибка: shekel problem is null" << endl; continue; }
                double actualMin = shekel->GetOptimumValue();
                vector<double> actualMinPoint = shekel->GetOptimumPoint();
                cout << "\nShekel Problem " << index << endl;
                cout << "Фактический минимум (из файла): " << actualMin
                     << " в точке x = " << (actualMinPoint.empty() ? NAN : actualMinPoint[0]) << endl;
                vector<double> a = {0.0};
                vector<double> b = {10.0};
                Minimizer<TShekelProblem> minimizer(a, b, epsilon_val, r_val, *shekel);
                auto startTime = chrono::high_resolution_clock::now();
                vector<double> computedMinPoint = minimizer.findMinimum();
                auto endTime = chrono::high_resolution_clock::now();
                chrono::duration<double> duration = endTime - startTime;
                cout << "Посчитанный минимум: " << shekel->ComputeFunction(computedMinPoint)
                     << " в точке x = " << (computedMinPoint.empty() ? NAN : computedMinPoint[0]) << endl;
                cout << "Итераций: " << minimizer.GetIterationCount() << ", время: "
                     << duration.count() * 1000 << " мс" << endl;
                totalIterations += minimizer.GetIterationCount();
                exitMainTotal += minimizer.GetExitMainCount();
                exitTestTotal += minimizer.GetExitTestCount();
                dataFile << i + 1 << " " << minimizer.GetIterationCount() << endl;
            }
        }
    } else if (taskType == 2) {
        // Multidimensional task: e.g., Grishagin problem family.
        TGrishaginProblemFamily grishFamily;
        // uniform_int_distribution<> grishDist(0, grishFamily.GetFamilySize() - 1); // FamilySize might be 100, indices 1-100
        
        for (int i = 0; i < numTests; ++i) {
            // int index = grishDist(gen) + 1; // if family indices are 1-based
            int index = (i % grishFamily.GetFamilySize()) + 1; // Cycle through problems if numTests > family size
            if (index > 100) index = 100; // Cap at 100 for Grishagin if that's the max
            
            IOptProblem* problemBase = grishFamily[index];
            if (!problemBase) {
                 cerr << "Ошибка: grishagin problem " << index << " is null" << endl;
                 continue;
            }
            // Assuming TGrishaginProblem inherits from FunctionInterface or is adaptable
            // For simplicity, let's assume grishFamily[index] returns a compatible type
            // If TGrishaginProblem is the concrete type:
            TGrishaginProblem* grishProblem = dynamic_cast<TGrishaginProblem*>(problemBase);
            if (!grishProblem) {
                cerr << "Ошибка: Не удалось привести к TGrishaginProblem для индекса " << index << endl;
                continue;
            }

            double actualMin = grishProblem->GetOptimumValue();
            vector<double> actualMinPoint = grishProblem->GetOptimumPoint();
            cout << "\nGrishagin Problem " << index << endl;
            cout << "Фактический минимум (из файла): " << actualMin;
            if (actualMinPoint.size() >= 2) {
                 cout << " в точке (x,y) = (" << actualMinPoint[0] << ", " << actualMinPoint[1] << ")" << endl;
            } else {
                cout << " (точка оптимума не полностью определена)" << endl;
            }

            vector<double> a = {0.0}; // Search space for the 1D parameter of Peano curve
            vector<double> b = {1.0};
            // The type for Minimizer should be the concrete problem type TGrishaginProblem
            Minimizer<TGrishaginProblem> minimizer(a, b, epsilon_val, r_val, *grishProblem, 10, 2, 1); // order=10, dim=2, key=1
            
            auto startTime = chrono::high_resolution_clock::now();
            vector<double> computedMinPoint = minimizer.findMinimum();
            auto endTime = chrono::high_resolution_clock::now();
            chrono::duration<double> duration = endTime - startTime;
            
            cout << "Посчитанный минимум: " << grishProblem->ComputeFunction(computedMinPoint);
            if (computedMinPoint.size() >= 2) {
                cout << " в точке (x,y) = (" << computedMinPoint[0] << ", " << computedMinPoint[1] << ")" << endl;
            } else if (!computedMinPoint.empty()){
                 cout << " в точке x_param = " << computedMinPoint[0] << " (многомерная точка не полностью определена)" << endl;
            } else {
                cout << " (точка не найдена)" << endl;
            }
            cout << "Итераций: " << minimizer.GetIterationCount() << ", время: "
                 << duration.count() * 1000 << " мс" << endl;
            totalIterations += minimizer.GetIterationCount();
            exitMainTotal += minimizer.GetExitMainCount();
            exitTestTotal += minimizer.GetExitTestCount();
            dataFile << i + 1 << " " << minimizer.GetIterationCount() << endl;
        }
    } else {
        cout << "Некорректный выбор." << endl;
        dataFile.close(); // Close file even on error
        return 1;
    }
    dataFile.close();

    // Write general statistics to stats.txt
    ofstream statsFile("stats.txt", ios::out);
    if (!statsFile.is_open()) {
        cerr << "Ошибка открытия файла stats.txt!" << endl;
        return 1;
    }
    statsFile << exitMainTotal << " " << exitTestTotal << endl;
    statsFile << epsilon_val << " " << r_val << " " << numTests << endl;
    statsFile.close();

    cout << "\nОбщая статистика:" << endl;
    if (numTests > 0) {
        cout << "Среднее число итераций: " << totalIterations / numTests << endl;
    } else {
        cout << "Среднее число итераций: N/A (0 тестов)" << endl;
    }
    cout << "Счетчик основного условия: " << exitMainTotal << endl;
    cout << "Счетчик дополнительного условия: " << exitTestTotal << endl;

    // Run Python script for plotting
    string command = "python ../../plot_graph.py"; // Adjust path if necessary
    #ifdef _WIN32
        // No special handling needed for system() on Windows for python usually
    #else // Linux/macOS
        // If python is python3 on your system:
        // command = "python3 ../../plot_graph.py";
    #endif
    
    cout << "Запуск скрипта построения графика: " << command << endl;
    int result = system(command.c_str());
    if (result != 0) {
       cerr << "Ошибка при запуске Python скрипта! Код ошибки: " << result << endl;
       // Don't return 1 here, as minimization might have been successful.
       // This is a post-processing step.
    }

    cout << "Нажмите Enter для выхода...";
    cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Clear remaining input
    cin.get(); // Wait for Enter key

    return 0;
}