// main.cpp (или Source.cpp)

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <fstream>
#include <cstdlib>     // Для system()
#include <limits>      // Для numeric_limits
#include <cmath>       // Для NAN, если используется

// Наш главный заголовочный файл с классом Minimizer
#include "minimizer.h"

// Заголовочные файлы для конкретных задач
#include "HillProblem.hpp"
#include "ShekelProblem.hpp"
#include "HillProblemFamily.hpp"
#include "ShekelProblemFamily.hpp"
#include "grishagin_function.hpp" // Предполагается, что определяет TGrishaginProblem
#include "GrishaginProblemFamily.hpp"

// using namespace std; // Уже есть в minimizer.h и доступно здесь

int main() {
    setlocale(LC_ALL, "Russian");
    
    std::ofstream main_log_file("minimization_log.txt", std::ios::out); // std::ios::out перезапишет старый файл
    if (!main_log_file.is_open()) {
        std::cerr << "КРИТИЧЕСКАЯ ОШИБКА: Не удалось открыть основной лог-файл minimization_log.txt!" << std::endl;
        return 1; // Завершаем, если лог не открылся
    }
    main_log_file << "--- Начало сессии логирования ---" << std::endl;

    int taskType_choice;
    int numTests_val;
    double epsilon_param, r_param_strongin;

    std::cout << "Выберите тип задачи:" << std::endl;
    std::cout << " 1 - Одномерная (Хилл/Шекеля)" << std::endl;
    std::cout << " 2 - Многомерная (Гришагина)" << std::endl;
    std::cout << "Ваш выбор: ";
    std::cin >> taskType_choice;

    std::cout << "Введите количество функций для тестирования: ";
    std::cin >> numTests_val;

    std::cout << "Введите точность эпсилон (> 0): ";
    std::cin >> epsilon_param;
    std::cout << "Введите параметр r метода Стронгина (например, >1, обычно 2.0-4.0): ";
    std::cin >> r_param_strongin;

    main_log_file << "Выбран тип задачи: " << taskType_choice << std::endl;
    main_log_file << "Количество тестов: " << numTests_val << ", Эпсилон: " << epsilon_param << ", r: " << r_param_strongin << std::endl;
    // Файл для записи данных для построения графика
    std::ofstream dataFile("plot_data.txt", std::ios::out);
    if (!dataFile.is_open()) {
        std::cerr << "Ошибка открытия файла plot_data.txt!" << std::endl;
        return 1;
    }

    // Счетчики для статистики
    int totalExitMainCount = 0;
    int totalExitTestCount = 0;
    double cumulativeIterations = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());

    if (taskType_choice == 1) {
        // --- Одномерные задачи ---
        int subChoice;
        std::cout << "Выберите одномерную задачу:" << std::endl;
        std::cout << " 1 - Функция Хилла" << std::endl;
        std::cout << " 2 - Функция Шекеля" << std::endl;
        std::cout << "Ваш выбор: ";
        std::cin >> subChoice;

        THillProblemFamily hillFamily;
        TShekelProblemFamily shekelFamily;
        // Распределения для случайного выбора индекса задачи из семейства
        std::uniform_int_distribution<> hillDist(0, hillFamily.GetFamilySize() - 1);
        std::uniform_int_distribution<> shekelDist(0, shekelFamily.GetFamilySize() - 1);

        for (int i = 0; i < numTests_val; ++i) {
            if (subChoice == 1) { // Задача Хилла
                int problem_idx = hillDist(gen);
                THillProblem* hillProblemInstance = dynamic_cast<THillProblem*>(hillFamily[problem_idx]);
                if (!hillProblemInstance) {
                    std::cerr << "Ошибка: не удалось получить экземпляр задачи Хилла " << problem_idx << std::endl;
                    continue;
                }

                double actualOptValue = hillProblemInstance->GetOptimumValue();
                std::vector<double> actualOptPointCoords = hillProblemInstance->GetOptimumPoint();

                std::cout << "\nЗадача Хилла, индекс " << problem_idx << std::endl;
                std::cout << "Известный оптимум (из файла): " << actualOptValue
                          << " в точке x = " << (actualOptPointCoords.empty() ? NAN : actualOptPointCoords[0]) << std::endl;

                // Границы для одномерной задачи Хилла (обычно [0,1])
                std::vector<double> problem_domain_a = {0.0};
                std::vector<double> problem_domain_b = {1.0};

                main_log_file << "\n--- Запуск теста 1D: " << (subChoice == 1 ? "Хилл" : "Шекель") 
                          << ", Индекс задачи " << problem_idx << " ---" << std::endl;
                std::cout << "DEBUG: Перед созданием 1D Minimizer..." << std::endl;
                Minimizer<THillProblem> minimizer(problem_domain_a, problem_domain_b, epsilon_param, r_param_strongin, 
                                                *hillProblemInstance, main_log_file); // Передаем main_log_file
                std::cout << "DEBUG: 1D Minimizer создан." << std::endl;

                auto startTime = std::chrono::high_resolution_clock::now();
                std::vector<double> computedMinCoords = minimizer.findMinimum();
                auto endTime = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = endTime - startTime;

                double computedOptValue = computedMinCoords.empty() ? NAN : hillProblemInstance->ComputeFunction(computedMinCoords);
                std::cout << "Найденный минимум: " << computedOptValue
                          << " в точке x = " << (computedMinCoords.empty() ? NAN : computedMinCoords[0]) << std::endl;
                std::cout << "Итераций: " << minimizer.GetIterationCount() << ", время: "
                          << duration.count() * 1000 << " мс" << std::endl;

                cumulativeIterations += minimizer.GetIterationCount();
                totalExitMainCount += minimizer.GetExitMainCount();
                totalExitTestCount += minimizer.GetExitTestCount();
                dataFile << i + 1 << " " << minimizer.GetIterationCount() << std::endl;

            } else { // Задача Шекеля
                int problem_idx = shekelDist(gen);
                TShekelProblem* shekelProblemInstance = dynamic_cast<TShekelProblem*>(shekelFamily[problem_idx]);
                 if (!shekelProblemInstance) {
                    std::cerr << "Ошибка: не удалось получить экземпляр задачи Шекеля " << problem_idx << std::endl;
                    continue;
                 }
                double actualOptValue = shekelProblemInstance->GetOptimumValue();
                std::vector<double> actualOptPointCoords = shekelProblemInstance->GetOptimumPoint();

                std::cout << "\nЗадача Шекеля, индекс " << problem_idx << std::endl;
                std::cout << "Известный оптимум (из файла): " << actualOptValue
                          << " в точке x = " << (actualOptPointCoords.empty() ? NAN : actualOptPointCoords[0]) << std::endl;

                // Границы для одномерной задачи Шекеля (обычно [0,10])
                std::vector<double> problem_domain_a = {0.0};
                std::vector<double> problem_domain_b = {10.0};
                
                main_log_file << "\n--- Запуск теста 1D: " << (subChoice == 1 ? "Хилл" : "Шекель") 
                          << ", Индекс задачи " << problem_idx << " ---" << std::endl;
                std::cout << "DEBUG: Перед созданием 1D Minimizer..." << std::endl;
                Minimizer<TShekelProblem> minimizer(problem_domain_a, problem_domain_b, epsilon_param, r_param_strongin, 
                                                *shekelProblemInstance, main_log_file); // Передаем main_log_file
                std::cout << "DEBUG: 1D Minimizer создан." << std::endl;

                //Minimizer<TShekelProblem> minimizer(problem_domain_a, problem_domain_b, epsilon_param, r_param_strongin, *shekelProblemInstance);

                auto startTime = std::chrono::high_resolution_clock::now();
                std::vector<double> computedMinCoords = minimizer.findMinimum();
                auto endTime = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = endTime - startTime;

                double computedOptValue = computedMinCoords.empty() ? NAN : shekelProblemInstance->ComputeFunction(computedMinCoords);
                std::cout << "Найденный минимум: " << computedOptValue
                          << " в точке x = " << (computedMinCoords.empty() ? NAN : computedMinCoords[0]) << std::endl;
                std::cout << "Итераций: " << minimizer.GetIterationCount() << ", время: "
                          << duration.count() * 1000 << " мс" << std::endl;

                cumulativeIterations += minimizer.GetIterationCount();
                totalExitMainCount += minimizer.GetExitMainCount();
                totalExitTestCount += minimizer.GetExitTestCount();
                dataFile << i + 1 << " " << minimizer.GetIterationCount() << std::endl;
            }
        }
    } else if (taskType_choice == 2) {
        // --- Многомерная задача: Гришагина ---
        TGrishaginProblemFamily grishaginFamily;
        // Распределение для случайного выбора индекса задачи (если нужно, обычно задачи Гришагина нумеруются 1-100)
        // std::uniform_int_distribution<> grishaginDist(1, grishaginFamily.GetFamilySize());

        for (int i = 0; i < numTests_val; ++i) {
            // int problem_idx = grishaginDist(gen); // Случайный выбор
            int problem_idx = (i % grishaginFamily.GetFamilySize()) + 1; // Последовательный выбор с циклическим повторением
            if (problem_idx > 100 && grishaginFamily.GetFamilySize() >= 100) problem_idx = 100; // Ограничение для стандартных 100 задач Гришагина

            TGrishaginProblem* grishaginProblemInstance = dynamic_cast<TGrishaginProblem*>(grishaginFamily[problem_idx]);
            if (!grishaginProblemInstance) {
                 std::cerr << "Ошибка: не удалось получить экземпляр задачи Гришагина " << problem_idx << std::endl;
                 continue;
            }

            double actualOptValue = grishaginProblemInstance->GetOptimumValue();
            std::vector<double> actualOptPointCoords = grishaginProblemInstance->GetOptimumPoint();

            std::cout << "\nЗадача Гришагина, индекс " << problem_idx << std::endl;
            std::cout << "Известный оптимум (из файла): " << actualOptValue;
            if (actualOptPointCoords.size() >= 2) {
                 std::cout << " в точке (y1,y2) = (" << actualOptPointCoords[0] << ", " << actualOptPointCoords[1] << ")" << std::endl;
            } else {
                std::cout << " (координаты точки оптимума не полностью определены)" << std::endl;
            }

            // Параметры для поиска с кривой Пеано
            double peano_param_min = 0.0; // Границы для одномерного параметра x_param
            double peano_param_max = 1.0;
            int peano_order_m = 15;       // Порядок кривой Пеано (m)
            int problem_dim_n = 2;        // Размерность исходной задачи Гришагина (n)
            int peano_key = 1;            // Ключ для функции mapd кривой Пеано

            main_log_file << "\n--- Запуск теста nD: Гришагин, Индекс задачи " << problem_idx << " ---" << std::endl;
            std::cout << "DEBUG: Перед созданием nD Minimizer для задачи " << problem_idx << "..." << std::endl;
            Minimizer<TGrishaginProblem> minimizer(
                peano_param_min, peano_param_max,
                epsilon_param, r_param_strongin, *grishaginProblemInstance,
                peano_order_m, problem_dim_n, peano_key,
                main_log_file // Передаем main_log_file
            );
            std::cout << "DEBUG: nD Minimizer создан." << std::endl;

            auto startTime = std::chrono::high_resolution_clock::now();
            std::vector<double> computedMinCoords = minimizer.findMinimum();
            auto endTime = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = endTime - startTime;

            double computedOptValue = computedMinCoords.empty() ? NAN : grishaginProblemInstance->ComputeFunction(computedMinCoords);
            std::cout << "Найденный минимум: " << computedOptValue;
            if (computedMinCoords.size() >= 2) {
                std::cout << " в точке (y1,y2) = (" << computedMinCoords[0] << ", " << computedMinCoords[1] << ")" << std::endl;
            } else if (!computedMinCoords.empty()){
                 std::cout << " (многомерные координаты не полностью определены, возможно только x_param найден)" << std::endl;
            } else {
                std::cout << " (точка не найдена)" << std::endl;
            }
            std::cout << "Итераций: " << minimizer.GetIterationCount() << ", время: "
                 << duration.count() * 1000 << " мс" << std::endl;

            cumulativeIterations += minimizer.GetIterationCount();
            totalExitMainCount += minimizer.GetExitMainCount();
            totalExitTestCount += minimizer.GetExitTestCount();
            dataFile << i + 1 << " " << minimizer.GetIterationCount() << std::endl;
        }
    } else {
        std::cout << "Некорректный выбор типа задачи." << std::endl;
        dataFile.close();
        return 1;
    }
    dataFile.close();

    // Запись общей статистики в файл stats.txt
    std::ofstream statsFile("stats.txt", std::ios::out);
    if (!statsFile.is_open()) {
        std::cerr << "Ошибка открытия файла stats.txt!" << std::endl;
        return 1; // Можно и не завершать программу, если статистика не критична
    }
    statsFile << totalExitMainCount << " " << totalExitTestCount << std::endl;
    statsFile << epsilon_param << " " << r_param_strongin << " " << numTests_val << std::endl;
    statsFile.close();

    std::cout << "\nОбщая статистика:" << std::endl;
    if (numTests_val > 0) {
        std::cout << "Среднее число итераций: " << cumulativeIterations / numTests_val << std::endl;
    } else {
        std::cout << "Среднее число итераций: N/A (0 тестов выполнено)" << std::endl;
    }
    std::cout << "Счетчик основного условия останова: " << totalExitMainCount << std::endl;
    std::cout << "Счетчик дополнительного условия останова: " << totalExitTestCount << std::endl;

    // Запуск Python-скрипта для построения графика
    std::string python_command = "python ../../plot_graph.py"; // Путь к скрипту может потребовать корректировки
    // Для Linux/macOS может потребоваться "python3"
    // #ifdef __linux__ || __APPLE__
    //  python_command = "python3 ../../plot_graph.py";
    // #endif

    std::cout << "Запуск скрипта для построения графика: " << python_command << std::endl;
    int script_result = system(python_command.c_str());
    if (script_result != 0) {
       std::cerr << "Ошибка при запуске Python скрипта! Код ошибки: " << script_result << std::endl;
       // Не завершаем программу из-за этого, так как основная работа выполнена
    }

    std::cout << "\nНажмите Enter для выхода...";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Очистка буфера ввода
    std::cin.get();

    main_log_file << "--- Конец сессии логирования ---" << std::endl;
    main_log_file.close(); // Закрываем файл в конце main
    std::cout << "DEBUG: Программа main завершает работу. Лог-файл закрыт." << std::endl;
    return 0;
}