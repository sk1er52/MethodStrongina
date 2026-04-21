#ifndef BENCHMARK_RUNNER_HPP
#define BENCHMARK_RUNNER_HPP

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <nlopt.hpp>
#include "dual_lipschitz_optimizer.hpp"
#include "IOptProblemFamily.hpp"

// Структура для отслеживания истории внутри NLopt
struct NLoptContext {
    IOptProblem* problem;
    double target_z;
    double eps;
    int current_eval;
    int first_hit_eval; // Здесь запомним итерацию успеха
};

// Обертка для NLopt
inline double nlopt_objective_wrapper(const std::vector<double> &x, std::vector<double> &grad, void *f_data) {
    NLoptContext* ctx = static_cast<NLoptContext*>(f_data);
    ctx->current_eval++; 
    
    double val = ctx->problem->ComputeFunction(x);
    
    // Если мы еще не находили минимум, проверяем текущую точку
    if (ctx->first_hit_eval == -1) {
        // Критерий: попали в окрестность известного глобального минимума
        if (std::abs(val - ctx->target_z) <= ctx->eps * 10.0) {
            ctx->first_hit_eval = ctx->current_eval; // ЗАПОМНИЛИ!
        }
    }
    
    return val;
}

class BenchmarkRunner {
public:
    struct Config {
        double eps;
        double r_loc;
        double r_glob;
        int max_iter;
        int peano_m;
    };

    static std::vector<int> RunMyAlgorithm(IOptProblemFamily& family, const Config& cfg) {
        std::vector<int> results;
        std::ostream null_stream(nullptr); 

        for (int i = 0; i < family.GetFamilySize(); ++i) {
            IOptProblem* problem = family[i];
            int dim = problem->GetDimension();

            DualLipschitzOptimizer<IOptProblem> opt(0.0, 1.0, cfg.eps, cfg.r_loc, cfg.r_glob, 
                                                   *problem, cfg.peano_m, dim, 1, null_stream, cfg.max_iter);
            
            // Даем алгоритму отработать по его правилам
            opt.Optimize();
            
            // АНАЛИЗИРУЕМ ИСТОРИЮ ЭКСПЕРИМЕНТОВ
            const auto& history = opt.GetTrialPointsHistory();
            int hit_iter = cfg.max_iter + 1; // По умолчанию - не решено
            
            for (size_t k = 0; k < history.size(); ++k) {
                if (!std::isnan(history[k].z_value)) {
                    if (IsSolved(history[k].z_value, problem->GetOptimumValue(), cfg.eps)) {
                        hit_iter = k + 1; // Нашли! Запоминаем номер итерации (1-based)
                        break;            // Дальше историю можно не смотреть
                    }
                }
            }
            results.push_back(hit_iter);
        }
        return results;
    }

    static std::vector<int> RunNLopt(IOptProblemFamily& family, nlopt::algorithm alg, const Config& cfg) {
        std::vector<int> results;
        for (int i = 0; i < family.GetFamilySize(); ++i) {
            IOptProblem* problem = family[i];
            int dim = problem->GetDimension();

            // Инициализируем контекст: first_hit_eval = -1 (еще не найдено)
            NLoptContext ctx = {problem, problem->GetOptimumValue(), cfg.eps, 0, -1};

            nlopt::opt opt(alg, dim);
            opt.set_min_objective(nlopt_objective_wrapper, &ctx);
            
            std::vector<double> lb, ub;
            problem->GetBounds(lb, ub);
            opt.set_lower_bounds(lb);
            opt.set_upper_bounds(ub);
            
            opt.set_xtol_abs(cfg.eps);
            opt.set_maxeval(cfg.max_iter);

            std::vector<double> x(dim);
            for (int j = 0; j < dim; ++j) x[j] = (lb[j] + ub[j]) / 2.0;
            
            double min_f;
            try {
                opt.optimize(x, min_f);
            } catch (...) {
                // Игнорируем внутренние ошибки NLopt
            }

            // АНАЛИЗИРУЕМ РЕЗУЛЬТАТ
            if (ctx.first_hit_eval != -1) {
                // Алгоритм наткнулся на минимум в процессе работы
                results.push_back(ctx.first_hit_eval);
            } else {
                // Так и не нашел за max_iter
                results.push_back(cfg.max_iter + 1);
            }
        }
        return results;
    }

private:
    static bool IsSolved(double found_z, double known_opt_z, double eps) {
        return std::abs(found_z - known_opt_z) <= eps * 10.0; 
    }
};

#endif