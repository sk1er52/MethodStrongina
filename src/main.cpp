#include <chrono>
#include <clocale>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept> 
#include <string>
#include <vector>

// Наши заголовочные файлы
#include "common_types.hpp"
#include "dual_lipschitz_optimizer.hpp"
#include "json_logger.hpp"
#include "benchmark_runner.hpp"

// Семейства задач
#include "GrishaginProblemFamily.hpp"
#include "HillProblem.hpp"
#include "HillProblemFamily.hpp"
#include "ShekelProblem.hpp"
#include "ShekelProblemFamily.hpp"
#include "GKLSProblemFamily.hpp"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [mode] [params...]" << std::endl;
        return 1;
    }

    std::string mode = argv[1];

    if (mode == "--benchmark") {
        if (argc < 7) {
            std::cerr << "Usage: --benchmark [family: 1-Grishagin, 2-GKLS_Simple, 3-GKLS_Hard] [eps] [r_glob] [r_loc] [max_iter] [peano_m] [dim]" << std::endl;
            return 1;
        }
        try {
            int family_type = std::stoi(argv[2]);
            double eps = std::stod(argv[3]);
            double r_glob = std::stod(argv[4]);
            double r_loc = std::stod(argv[5]);
            int max_iter = std::stoi(argv[6]);
            int peano_m = std::stoi(argv[7]);
            int dim = (argc > 8) ? std::stoi(argv[8]) : 2;

            BenchmarkRunner::Config cfg_dl = {eps, r_loc, r_glob, max_iter, peano_m};
            BenchmarkRunner::Config cfg_classic = {eps, r_glob, r_glob, max_iter, peano_m}; 
            
            std::vector<int> res_agp_dl, res_agp, res_direct, res_crs2, res_isres, res_esch;

            if (family_type == 1) {
                TGrishaginProblemFamily family;
                res_agp_dl = BenchmarkRunner::RunMyAlgorithm(family, cfg_dl);
                res_agp    = BenchmarkRunner::RunMyAlgorithm(family, cfg_classic);
                res_direct = BenchmarkRunner::RunNLopt(family, nlopt::GN_DIRECT_L, cfg_dl);
                res_crs2   = BenchmarkRunner::RunNLopt(family, nlopt::GN_CRS2_LM, cfg_dl);
                res_isres  = BenchmarkRunner::RunNLopt(family, nlopt::GN_ISRES, cfg_dl);
                res_esch   = BenchmarkRunner::RunNLopt(family, nlopt::GN_ESCH, cfg_dl);
            } else if (family_type == 2) {
                TGKLSProblemFamily family(dim, Simple, TD);
                res_agp_dl = BenchmarkRunner::RunMyAlgorithm(family, cfg_dl);
                res_agp    = BenchmarkRunner::RunMyAlgorithm(family, cfg_classic);
                res_direct = BenchmarkRunner::RunNLopt(family, nlopt::GN_DIRECT_L, cfg_dl);
                res_crs2   = BenchmarkRunner::RunNLopt(family, nlopt::GN_CRS2_LM, cfg_dl);
                res_isres  = BenchmarkRunner::RunNLopt(family, nlopt::GN_ISRES, cfg_dl);
                res_esch   = BenchmarkRunner::RunNLopt(family, nlopt::GN_ESCH, cfg_dl);
            } else if (family_type == 3) {
                TGKLSProblemFamily family(dim, Hard, TD);
                res_agp_dl = BenchmarkRunner::RunMyAlgorithm(family, cfg_dl);
                res_agp    = BenchmarkRunner::RunMyAlgorithm(family, cfg_classic);
                res_direct = BenchmarkRunner::RunNLopt(family, nlopt::GN_DIRECT_L, cfg_dl);
                res_crs2   = BenchmarkRunner::RunNLopt(family, nlopt::GN_CRS2_LM, cfg_dl);
                res_isres  = BenchmarkRunner::RunNLopt(family, nlopt::GN_ISRES, cfg_dl);
                res_esch   = BenchmarkRunner::RunNLopt(family, nlopt::GN_ESCH, cfg_dl);
            }

            std::string family_str;
            int total_tasks = 0;
            if (family_type == 1) { family_str = "Grishagin"; total_tasks = 100; }
            else if (family_type == 2) { family_str = "GKLS Simple"; total_tasks = 100; }
            else if (family_type == 3) { family_str = "GKLS Hard"; total_tasks = 100; }

            // Вывод в обновленном формате JSON
            std::cout << "{" << std::endl;
            std::cout << "  \"info\": {" << std::endl;
            std::cout << "    \"family\": \"" << family_str << "\"," << std::endl;
            std::cout << "    \"total_tasks\": " << total_tasks << std::endl;
            std::cout << "  }," << std::endl;
            std::cout << "  \"parameters\": {" << std::endl;
            std::cout << "    \"r_glob\": " << r_glob << "," << std::endl;
            std::cout << "    \"r_loc\": " << r_loc << std::endl;
            std::cout << "  }," << std::endl;
            std::cout << "  \"results\": {" << std::endl;
            
            auto print_array = [](const std::string& name, const std::vector<int>& arr, bool last) {
                std::cout << "    \"" << name << "\": [";
                for(size_t i=0; i<arr.size(); ++i) std::cout << arr[i] << (i==arr.size()-1?"":",");
                std::cout << "]" << (last ? "" : ",") << std::endl;
            };

            print_array("agp_dl", res_agp_dl, false);
            print_array("agp_classic", res_agp, false);
            print_array("direct", res_direct, false);
            print_array("crs2", res_crs2, false);
            print_array("isres", res_isres, false);
            print_array("esch", res_esch, true);

            std::cout << "  }" << std::endl;
            std::cout << "}" << std::endl;

            return 0;
        } catch (const std::exception &e) {
            std::cerr << "Benchmark error: " << e.what() << std::endl;
            return 1;
        }
    }

    // --- РЕЖИМ 2: ЛИНИИ УРОВНЯ (Для визуализации в GUI) ---
    if (mode == "--get-level-lines") {
        if (argc < 4) return 1;
        int problem_idx = std::stoi(argv[2]);
        int grid_res = std::stoi(argv[3]);

        TGrishaginProblemFamily family;
        TGrishaginProblem *instance = dynamic_cast<TGrishaginProblem *>(family[problem_idx]);
        if (!instance) return 1;

        std::vector<double> x_grid(grid_res), y_grid(grid_res);
        std::vector<std::vector<double>> z_grid(grid_res, std::vector<double>(grid_res));

        for (int i = 0; i < grid_res; ++i) {
            y_grid[i] = 0.0 + i * (1.0 / (grid_res - 1));
            for (int j = 0; j < grid_res; ++j) {
                x_grid[j] = 0.0 + j * (1.0 / (grid_res - 1));
                z_grid[i][j] = instance->ComputeFunction({x_grid[j], y_grid[i]});
            }
        }
        JsonLogger::PrintLevelLines(x_grid, y_grid, z_grid);
        return 0;
    }

    // --- РЕЖИМ 3: ОДИНОЧНЫЙ ЗАПУСК (Стандартная минимизация) ---
    std::ofstream cpp_internal_log("minimization_log_cpp.txt", std::ios::out);
    
    int task_type = 0;
    int sub_choice = 0;
    int problem_idx = 0;
    double eps = 0.001;
    double r_glob = 3.0;
    double r_loc = 1.5;
    int peano_m = 10;
    int max_iter = 10000;

    try {
        task_type = std::stoi(argv[1]);
        eps = std::stod(argv[3]);
        r_glob = std::stod(argv[4]);
        max_iter = std::stoi(argv[6]);
        if (argc >= 8) r_loc = std::stod(argv[7]);

        if (task_type == 1) {
            sub_choice = std::stoi(argv[2]);
            problem_idx = std::stoi(argv[5]);
        } else if (task_type == 2) {
            problem_idx = std::stoi(argv[2]);
            peano_m = std::stoi(argv[5]);
        } else if (task_type == 3) {
            problem_idx = std::stoi(argv[2]);
            peano_m = 10;
        }
    } catch (...) { return 1; }

    std::string family_name;
    double found_z = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> found_y;
    int iterations = 0;
    double time_ms = 0;
    double known_opt_z = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> known_opt_y;
    std::vector<TrialPoint> history;

    if (task_type == 1) {
        if (sub_choice == 1) {
            family_name = "Hill";
            THillProblemFamily family;
            THillProblem *instance = dynamic_cast<THillProblem *>(family[problem_idx]);
            known_opt_z = instance->GetOptimumValue();
            known_opt_y = instance->GetOptimumPoint();
            DualLipschitzOptimizer<THillProblem> opt({0.0}, {1.0}, eps, r_loc, r_glob, *instance, cpp_internal_log, max_iter);
            auto start = std::chrono::high_resolution_clock::now();
            found_y = opt.Optimize();
            auto end = std::chrono::high_resolution_clock::now();
            time_ms = std::chrono::duration<double, std::milli>(end - start).count();
            iterations = opt.GetIterationCount();
            history = opt.GetTrialPointsHistory();
            if (!found_y.empty()) found_z = instance->ComputeFunction(found_y);
        } else {
            family_name = "Shekel";
            TShekelProblemFamily family;
            TShekelProblem *instance = dynamic_cast<TShekelProblem *>(family[problem_idx]);
            known_opt_z = instance->GetOptimumValue();
            known_opt_y = instance->GetOptimumPoint();
            DualLipschitzOptimizer<TShekelProblem> opt({0.0}, {10.0}, eps, r_loc, r_glob, *instance, cpp_internal_log, max_iter);
            auto start = std::chrono::high_resolution_clock::now();
            found_y = opt.Optimize();
            auto end = std::chrono::high_resolution_clock::now();
            time_ms = std::chrono::duration<double, std::milli>(end - start).count();
            iterations = opt.GetIterationCount();
            history = opt.GetTrialPointsHistory();
            if (!found_y.empty()) found_z = instance->ComputeFunction(found_y);
        }
    } else if (task_type == 2) {
        family_name = "Grishagin";
        TGrishaginProblemFamily family;
        TGrishaginProblem *instance = dynamic_cast<TGrishaginProblem *>(family[problem_idx]);
        known_opt_z = instance->GetOptimumValue();
        known_opt_y = instance->GetOptimumPoint();
        DualLipschitzOptimizer<TGrishaginProblem> opt(0.0, 1.0, eps, r_loc, r_glob, *instance, peano_m, 2, 1, cpp_internal_log, max_iter);
        auto start = std::chrono::high_resolution_clock::now();
        found_y = opt.Optimize();
        auto end = std::chrono::high_resolution_clock::now();
        time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        iterations = opt.GetIterationCount();
        history = opt.GetTrialPointsHistory();
        if (!found_y.empty()) found_z = instance->ComputeFunction(found_y);
    } else if (task_type == 3) {
        family_name = "GKLS";
        int dim = std::stoi(argv[5]);
        TGKLSProblemFamily family(dim, Simple, TD);
        TGKLSProblem *instance = dynamic_cast<TGKLSProblem *>(family[problem_idx]);
        known_opt_z = instance->GetOptimumValue();
        known_opt_y = instance->GetOptimumPoint();
        DualLipschitzOptimizer<TGKLSProblem> opt(0.0, 1.0, eps, r_loc, r_glob, *instance, peano_m, dim, 1, cpp_internal_log, max_iter);
        auto start = std::chrono::high_resolution_clock::now();
        found_y = opt.Optimize();
        auto end = std::chrono::high_resolution_clock::now();
        time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        iterations = opt.GetIterationCount();
        history = opt.GetTrialPointsHistory();
        if (!found_y.empty()) found_z = instance->ComputeFunction(found_y);
    }

    JsonLogger::PrintMinimizationResult(family_name, problem_idx, found_z, found_y, iterations, time_ms, 0, 0, known_opt_z, known_opt_y, history);
    cpp_internal_log.close();
    return 0;
}