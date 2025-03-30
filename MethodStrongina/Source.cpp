// source.cpp
#include <cmath>
#include <vector>
#include <limits>
#include <fstream>
#include <chrono>
#include <iostream>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <filesystem>

// Подключите ваши заголовочные файлы с задачами:
#include "HillProblem.hpp"
#include "ShekelProblem.hpp"
#include "HillProblemFamily.hpp"
#include "ShekelProblemFamily.hpp"
#include "grishagin_function.hpp"
#include "GrishaginProblemFamily.hpp"

using namespace std;

// Структура для представления точки минимизации
struct Point {
    vector<double> x; // одномерное представление параметра (при использовании отображения, x ∈ [0,1])
    double y;         // значение целевой функции в данной точке
};

// Интерфейс для целевой функции (если вам он нужен для унификации)
class FunctionInterface {
public:
    // Вычисляет значение функции в точке x
    virtual double ComputeFunction(const vector<double>& x) const = 0;
    // Возвращает известную точку оптимума
    virtual vector<double> GetOptimumPoint() const = 0;
    virtual ~FunctionInterface() {}
};

//
// Реализация кривой Пеано (ваша версия)
// Здесь приведён ваш код функции mapd и вспомогательной функции node.
// (Обратите внимание – глобальные переменные и массивы используются, как в вашем коде.)
//
int n1, nexp, l, iq, iu[10], iv[10];
void mapd( double x, int m, double* y, int n, int key )
{
    /* mapping y(x) : 1 - center, 2 - line, 3 - node */
    double d, mne, dd, dr;
    float p, r;
    int iw[11];
    int it, is, i, j, k;
    void node ( int );

    p = 0.0;
    n1 = n - 1;
    for ( nexp = 1, i = 0; i < n; nexp *= 2, i++ );
    d = x;
    r = 0.5;
    it = 0;
    dr = nexp;
    for ( mne = 1, i = 0; i < m; mne *= dr, i++ );
    for ( i = 0; i < n; i++ ) {
        iw[i] = 1; y[i] = 0.0;
    }
    if ( key == 2 ) {
        d = d * (1.0 - 1.0 / mne); k = 0;
    } else if ( key > 2 ) {
        dr = mne / nexp;
        dr = dr - fmod(dr, 1.0);
        dd = mne - dr;
        dr = d * dd;
        dd = dr - fmod(dr, 1.0);
        dr = dd + (dd - 1) / (nexp - 1);
        dd = dr - fmod(dr, 1.0);
        d = dd * (1.0 / mne);
    }
    for ( j = 0; j < m; j++ ) {
        iq = 0;
        if ( x == 1.0 ) {
            is = nexp - 1; d = 0.0;
        } else {
            d = d * nexp;
            is = d;
            d = d - is;
        }
        i = is;
        node(i);
        i = iu[0];
        iu[0] = iu[it];
        iu[it] = i;
        i = iv[0];
        iv[0] = iv[it];
        iv[it] = i;
        if ( l == 0 )
            l = it;
        else if ( l == it )
            l = 0;
        if ( (iq > 0) || ((iq == 0) && (is == 0)) )
            k = l;
        else if ( iq < 0 )
            k = ( it == n1 ) ? 0 : n1;
        r = r * 0.5;
        it = l;
        for ( i = 0; i < n; i++ ) {
            iu[i] = iu[i] * iw[i];
            iw[i] = -iv[i] * iw[i];
            p = r * iu[i];
            p = p + y[i];
            y[i] = p;
        }
    }
    if ( key == 2 ) {
        if ( is == (nexp - 1) ) i = -1;
        else i = 1;
        p = 2 * i * iu[k] * r * d;
        p = y[k] - p;
        y[k] = p;
    } else if ( key == 3 ) {
        for ( i = 0; i < n; i++ ) {
            p = r * iu[i];
            p = p + y[i];
            y[i] = p;
        }
    }
}
void node ( int is )
{
    /* calculate iu, iv, l by is */
    int n, i, j, k1, k2, iff;
    n = n1 + 1;
    if ( is == 0 ) {
        l = n1;
        for ( i = 0; i < n; i++ ) {
            iu[i] = -1; iv[i] = -1;
        }
    } else if ( is == (nexp - 1) ) {
        l = n1;
        iu[0] = 1;
        iv[0] = 1;
        for ( i = 1; i < n; i++ ) {
            iu[i] = -1; iv[i] = -1;
        }
        iv[n1] = 1;
    } else {
        iff = nexp;
        k1 = -1;
        for ( i = 0; i < n; i++ ) {
            iff = iff / 2;
            if ( is >= iff ) {
                if ( (is == iff) && (is != 1) ) { l = i; iq = -1; }
                is = is - iff;
                k2 = 1;
            } else {
                k2 = -1;
                if ( (is == (iff - 1)) && (is != 0) ) { l = i; iq = 1; }
            }
            j = -k1 * k2;
            iv[i] = j;
            iu[i] = j;
            k1 = k2;
        }
        iv[l] = iv[l] * iq;
        iv[n1] = -iv[n1];
    }
}

// Обёртка вокруг кривой Пеано: принимает одномерное значение x и возвращает вектор размерности n.
vector<double> peanoMapping(double x, int m, int n, int key) {
    vector<double> y(n, 0.0);
    mapd(x, m, y.data(), n, key);
    return y;
}

//
// Класс минимизатора по методу Стронгина.
// Шаблонный параметр T – тип задачи (например, THillProblem, TShekelProblem, TGrishaginProblem).
// Предполагается, что у классов задач есть методы ComputeFunction(vector<double>),
// GetOptimumPoint(), GetOptimumValue(), GetDimension(), GetConstraintsNumber() (если есть ограничения) и т.д.
//
template <typename T>
class Minimizer {
private:
    vector<double> leftBound;   // Границы для одномерного параметра x ∈ [leftBound[0], rightBound[0]]
    vector<double> rightBound;
    int iterationCount;
    double epsilon;
    double r;
    const T& function;
    ofstream logFile;
    int exitMainCount;
    int exitTestCount;
    vector<Point> points;       // Точки, где вычислялась функция (одномерное представление)

    // Параметры для режима отображения (многомерная задача)
    bool useMapping;   // если true, то решаем многомерную задачу посредством отображения
    int mappingOrder;  // порядок (m) отображения
    int dimension;     // размерность исходной задачи (например, 2 для Гришагина)
    int mappingKey;    // ключ для функции mapd

public:
    // Конструктор для одномерного случая (без отображения)
    Minimizer(vector<double> a, vector<double> b, double eps, double r, const T& func)
        : leftBound(a), rightBound(b), epsilon(eps), r(r), function(func),
          exitMainCount(0), exitTestCount(0), useMapping(false)
    {
        logFile.open("minimization_log.txt", ios::out);
        if (!logFile.is_open()) {
            cerr << "Ошибка открытия файла журнала!" << endl;
        }
        // Инициализируем точки
        points.push_back({ leftBound, function.ComputeFunction(leftBound) });
        points.push_back({ rightBound, function.ComputeFunction(rightBound) });
        logFile << "Начальные точки:\n";
        logPoint(points[0]);
        logPoint(points[1]);
    }

    // Конструктор для многомерного случая (с отображением)
    Minimizer(vector<double> a, vector<double> b, double eps, double r, const T& func,
              int mappingOrder_, int dimension_, int mappingKey_)
        : leftBound(a), rightBound(b), epsilon(eps), r(r), function(func),
          exitMainCount(0), exitTestCount(0),
          useMapping(true), mappingOrder(mappingOrder_), dimension(dimension_), mappingKey(mappingKey_)
    {
        logFile.open("minimization_log.txt", ios::out);
        if (!logFile.is_open()) {
            cerr << "Ошибка открытия файла журнала!" << endl;
        }
        // Здесь leftBound и rightBound задаются для параметра x ∈ [0,1]
        vector<double> initPoint = { leftBound[0] };
        points.push_back({ initPoint, calculateY(leftBound[0]) });
        initPoint = { rightBound[0] };
        points.push_back({ initPoint, calculateY(rightBound[0]) });
        logFile << "Начальные точки (одномерные для отображения):\n";
        logPoint(points[0]);
        logPoint(points[1]);
    }

    ~Minimizer() {
        logFile.close();
    }

    // Вычисление значения функции в точке x (одномерное значение)
    double calculateY(double x) const {
        if (!useMapping) {
            vector<double> pt = { x };
            return function.ComputeFunction(pt);
        } else {
            // Преобразуем x через кривую Пеано в точку в ℝⁿ
            vector<double> mappedPoint = peanoMapping(x, mappingOrder, dimension, mappingKey);
            return function.ComputeFunction(mappedPoint);
        }
    }

    // Метод поиска минимума
    vector<double> findMinimum() { 
        iterationCount = 0;
        const int maxIterations = 10000;
        int iteration = 0;
        double M, mVal, yNew;
        vector<double> xNew(leftBound.size());
        double intervalSize = abs(rightBound[0] - leftBound[0]);
        
        // Для многомерного случая условие Гёльдера: alpha = 1/dimension, иначе alpha = 1.
        double alpha = useMapping ? (1.0 / double(dimension)) : 1.0;
        
        // Получаем фактическую точку оптимума (из файла или аналитически)
        vector<double> actualMinPoint = function.GetOptimumPoint();

        while (iteration < maxIterations) {
            M = calculateMaxSlope();
            mVal = (M > 0.0) ? r * M : 1.0;

            vector<double> intervalCharacteristics(points.size() - 1);
            for (size_t i = 0; i < points.size() - 1; ++i) {
                intervalCharacteristics[i] = calculateCharacteristic(i, mVal);
            }

            size_t maxCharacteristicIndex =
                distance(intervalCharacteristics.begin(),
                         max_element(intervalCharacteristics.begin(), intervalCharacteristics.end()));

            // Вычисляем новую точку как одномерное значение
            double x_left = points[maxCharacteristicIndex].x[0];
            double x_right = points[maxCharacteristicIndex + 1].x[0];
            xNew[0] = 0.5 * (x_left + x_right) - (points[maxCharacteristicIndex + 1].y - points[maxCharacteristicIndex].y) / (2.0 * mVal);

            // Ограничиваем xNew границами
            if (xNew[0] < leftBound[0])
                xNew[0] = leftBound[0];
            if (xNew[0] > rightBound[0])
                xNew[0] = rightBound[0];

            yNew = calculateY(xNew[0]);
            points.insert(points.begin() + maxCharacteristicIndex + 1, { xNew, yNew });

            logFile << "Итерация " << iteration + 1 << ": ";
            logPoint(points[maxCharacteristicIndex + 1]);

            // Условие останова с учетом Гёльдера
            double currentInterval = abs(points[maxCharacteristicIndex + 1].x[0] - points[maxCharacteristicIndex].x[0]);
            double threshold = pow(epsilon * intervalSize / 2.0, 1.0 / alpha);
            if (currentInterval <= threshold) {
                exitMainCount++;
                break;
            }
            if (abs(xNew[0] - actualMinPoint[0]) <= threshold) {
                exitTestCount++;
                break;
            }

            ++iteration;
        }

        // Находим точку с минимальным значением
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
            // Возвращаем найденное одномерное значение, преобразованное в многомерное через peanoMapping
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
            if (diffX > 1e-9) {
                double slope = abs((points[i].y - points[i - 1].y) / diffX);
                maxSlope = max(maxSlope, slope);
            }
        }
        return maxSlope;
    }

    double calculateCharacteristic(size_t index, double mVal) const {
        double interval = points[index + 1].x[0] - points[index].x[0];
        return mVal * interval + pow(points[index + 1].y - points[index].y, 2) / (mVal * interval)
               - 2.0 * (points[index + 1].y + points[index].y);
    }

    void logPoint(const Point& point) {
        logFile << "x: ";
        for (double val : point.x) {
            logFile << val << " ";
        }
        logFile << ", y: " << point.y << endl;
    }
};

//
// MAIN
//
int main() {
    setlocale(LC_ALL, "Russian");
    
    int taskType;
    int numTests;
    double epsilon, r;
    
    cout << "Выберите тип задачи:" << endl;
    cout << " 1 - Одномерная (Хилл/Шекеля)" << endl;
    cout << " 2 - Многомерная (Гришагина)" << endl;
    cout << "Ваш выбор: ";
    cin >> taskType;
    
    cout << "Введите количество случайных функций для тестирования: ";
    cin >> numTests;
    
    cout << "Введите точность (> 0): ";
    cin >> epsilon;
    cout << "Введите параметр r (например, 2.5): ";
    cin >> r;
    
    // Файл для записи данных для построения графика
    ofstream dataFile("plot_data.txt", ios::out);
    if (!dataFile.is_open()) {
        cerr << "Ошибка открытия файла plot_data.txt!" << endl;
        return 1;
    }
    
    // Счетчики для статистики
    int exitMainTotal = 0;
    int exitTestTotal = 0;
    double totalIterations = 0.0;
    
    random_device rd;
    mt19937 gen(rd());
    
    if (taskType == 1) {
        // Одномерные задачи: случайным образом выбираем либо задачу Хилла, либо Шекеля.
        // Предполагается, что у вас есть константы NUM_HILL_PROBLEMS и NUM_SHEKEL_PROBLEMS,
        // а также классы THillProblemFamily и TShekelProblemFamily.
        int subChoice;
        cout << "Выберите задачу:" << endl;
        cout << " 1 - Функция Хилла" << endl;
        cout << " 2 - Функция Шекеля" << endl;
        cout << "Ваш выбор: ";
        cin >> subChoice;
        
        // Создаём семейство задач (предполагается, что соответствующий класс имеет метод GetFamilySize() и оператор [])
        THillProblemFamily hillFamily;
        TShekelProblemFamily shekelFamily;
        uniform_int_distribution<> hillDist(0, hillFamily.GetFamilySize() - 1);
        uniform_int_distribution<> shekelDist(0, shekelFamily.GetFamilySize() - 1);
        
        for (int i = 0; i < numTests; ++i) {
            if (subChoice == 1) { // Хилл
                int index = hillDist(gen);
                THillProblem* hill = dynamic_cast<THillProblem*>(hillFamily[index]);
                double actualMin = hill->GetOptimumValue();
                vector<double> actualMinPoint = hill->GetOptimumPoint();
                cout << "\nHill Problem " << index << endl;
                cout << "Фактический минимум (из файла): " << actualMin 
                     << " в точке x = " << actualMinPoint[0] << endl;
                vector<double> a = {0.0};
                vector<double> b = {1.0};
                Minimizer<THillProblem> minimizer(a, b, epsilon, r, *hill);
                auto startTime = chrono::high_resolution_clock::now();
                vector<double> computedMin = minimizer.findMinimum();
                auto endTime = chrono::high_resolution_clock::now();
                chrono::duration<double> duration = endTime - startTime;
                cout << "Посчитанный минимум: " << hill->ComputeFunction(computedMin)
                     << " в точке x = " << computedMin[0] << endl;
                cout << "Итераций: " << minimizer.GetIterationCount() << ", время: " 
                     << duration.count() * 1000 << " мс" << endl;
                totalIterations += minimizer.GetIterationCount();
                exitMainTotal += minimizer.GetExitMainCount();
                exitTestTotal += minimizer.GetExitTestCount();
                dataFile << i + 1 << " " << minimizer.GetIterationCount() << endl;
            } else { // Шекеля
                int index = shekelDist(gen);
                TShekelProblem* shekel = dynamic_cast<TShekelProblem*>(shekelFamily[index]);
                double actualMin = shekel->GetOptimumValue();
                vector<double> actualMinPoint = shekel->GetOptimumPoint();
                cout << "\nShekel Problem " << index << endl;
                cout << "Фактический минимум (из файла): " << actualMin 
                     << " в точке x = " << actualMinPoint[0] << endl;
                vector<double> a = {0.0};
                vector<double> b = {10.0};
                Minimizer<TShekelProblem> minimizer(a, b, epsilon, r, *shekel);
                auto startTime = chrono::high_resolution_clock::now();
                vector<double> computedMin = minimizer.findMinimum();
                auto endTime = chrono::high_resolution_clock::now();
                chrono::duration<double> duration = endTime - startTime;
                cout << "Посчитанный минимум: " << shekel->ComputeFunction(computedMin)
                     << " в точке x = " << computedMin[0] << endl;
                cout << "Итераций: " << minimizer.GetIterationCount() << ", время: " 
                     << duration.count() * 1000 << " мс" << endl;
                totalIterations += minimizer.GetIterationCount();
                exitMainTotal += minimizer.GetExitMainCount();
                exitTestTotal += minimizer.GetExitTestCount();
                dataFile << i + 1 << " " << minimizer.GetIterationCount() << endl;
            }
        }
    } else if (taskType == 2) {
        // Многомерная задача: например, семейство задач Гришагина.
        TGrishaginProblemFamily grishFamily;
        uniform_int_distribution<> grishDist(0, grishFamily.GetFamilySize() - 1);
        for (int i = 0; i < numTests; ++i) {
            int index = grishDist(gen);
            // Предполагается, что класс GrishaginProblem имеет методы GetOptimumPoint, GetOptimumValue и ComputeFunction
            auto grishProblem = grishFamily[index];
            double actualMin = grishProblem->GetOptimumValue();
            vector<double> actualMinPoint = grishProblem->GetOptimumPoint();
            cout << "\nGrishagin Problem " << index << endl;
            cout << "Фактический минимум (из файла): " << actualMin 
                 << " в точке (x,y) = (" << actualMinPoint[0] << ", " << actualMinPoint[1] << ")" << endl;
            // Для многомерной задачи задаём, что x ∈ [0,1] будет отображаться в ℝ², параметры отображения: order=3, dimension=2, key=1.
            vector<double> a = {0.0};
            vector<double> b = {1.0};
            Minimizer<decltype(*grishProblem)> minimizer(a, b, epsilon, r, *grishProblem, 3, 2, 1);
            auto startTime = chrono::high_resolution_clock::now();
            vector<double> computedMin = minimizer.findMinimum();
            auto endTime = chrono::high_resolution_clock::now();
            chrono::duration<double> duration = endTime - startTime;
            cout << "Посчитанный минимум: " << grishProblem->ComputeFunction(computedMin)
                 << " в точке (x,y) = (" << computedMin[0] << ", " << computedMin[1] << ")" << endl;
            cout << "Итераций: " << minimizer.GetIterationCount() << ", время: " 
                 << duration.count() * 1000 << " мс" << endl;
            totalIterations += minimizer.GetIterationCount();
            exitMainTotal += minimizer.GetExitMainCount();
            exitTestTotal += minimizer.GetExitTestCount();
            dataFile << i + 1 << " " << minimizer.GetIterationCount() << endl;
        }
    } else {
        cout << "Некорректный выбор." << endl;
        return 1;
    }
    dataFile.close();
    
    // Запись общей статистики в файл stats.txt
    ofstream statsFile("stats.txt", ios::out);
    if (!statsFile.is_open()) {
        cerr << "Ошибка открытия файла stats.txt!" << endl;
        return 1;
    }
    statsFile << exitMainTotal << " " << exitTestTotal << endl;
    statsFile << epsilon << " " << r << " " << numTests << endl;
    statsFile.close();
    
    cout << "\nОбщая статистика:" << endl;
    cout << "Среднее число итераций: " << totalIterations / numTests << endl;
    cout << "Счетчик основного условия: " << exitMainTotal << endl;
    cout << "Счетчик дополнительного условия: " << exitTestTotal << endl;
    
    // Запуск Python-скрипта для построения графика
    string command = "python plot_graph.py";
    int result = system(command.c_str());
    if (result != 0) {
       cerr << "Ошибка при запуске Python скрипта! Код ошибки: " << result << endl;
       return 1;
    }
    
    system("pause");
    return 0;
}
