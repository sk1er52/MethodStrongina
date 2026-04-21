import json
import matplotlib.pyplot as plt
import numpy as np
import sys

def get_performance_curve(data, max_iter, total_tasks):
    solved_tasks = np.sort([x for x in data if x <= max_iter])
    
    if len(solved_tasks) == 0:
        return [0, max_iter], [0, 0], None
    
    x = [0] + list(solved_tasks) + [max_iter]
    y_steps = np.arange(1, len(solved_tasks) + 1) / total_tasks * 100
    final_percent = (len(solved_tasks) / total_tasks) * 100
    y = [0] + list(y_steps) + [final_percent]
    
    iter_100_percent = solved_tasks[-1] if len(solved_tasks) == total_tasks else None
    return x, y, iter_100_percent

try:
    with open('benchmark_results.json', 'r') as f:
        data = json.load(f)
        info = data['info']
        params = data['parameters']
        results = data['results']
except Exception as e:
    print(f"Ошибка чтения JSON: {e}. Пересобери C++ код!")
    sys.exit()

max_iter = 5000 # Убедись, что совпадает с лимитом из C++
total_tasks = info['total_tasks']
family_name = info['family']

plt.figure(figsize=(13, 8))

label_agp_dl = f"AGP-DL ($r_{{glob}}={params['r_glob']}, r_{{loc}}={params['r_loc']}$)"
label_agp = f"AGP Classic ($r={params['r_glob']}$)"

styles = {
    'agp_dl': {'label': label_agp_dl, 'color': 'red', 'linewidth': 2.5, 'zorder': 10},
    'agp_classic': {'label': label_agp, 'color': 'orange', 'linewidth': 2, 'linestyle': '--', 'zorder': 9},
    'direct': {'label': 'NLopt DIRECT', 'color': 'blue', 'linewidth': 2, 'linestyle': ':'},
    'crs2': {'label': 'NLopt CRS2 (Stochastic)', 'color': 'green', 'linewidth': 2, 'linestyle': '-.'},
    'isres': {'label': 'NLopt ISRES (Evolutionary)', 'color': 'purple', 'linewidth': 1.5, 'linestyle': '--'},
    'esch': {'label': 'NLopt ESCH (Evolutionary)', 'color': 'brown', 'linewidth': 1.5, 'linestyle': '-.'}
}

# Для раздвигания подписей 100%
annotation_offsets = [(0, 25), (0, 45), (0, 65), (0, 85)]
offset_idx = 0

for algo_name, iter_counts in results.items():
    x, y, iter_100 = get_performance_curve(iter_counts, max_iter, total_tasks)
    style = styles.get(algo_name, {'label': algo_name})
    
    plt.step(x, y, where='post', **style)
    
    if iter_100 is not None:
        plt.plot(iter_100, 100, marker='o', markersize=6, color=style['color'], zorder=11)
        
        # Рисуем выноску со стрелочкой, чтобы текст не слипался
        offset = annotation_offsets[offset_idx % len(annotation_offsets)]
        plt.annotate(f'{iter_100} ит.', 
                     xy=(iter_100, 100), 
                     xytext=offset, 
                     textcoords='offset points', 
                     ha='center', va='bottom',
                     color=style['color'], 
                     fontsize=10, fontweight='bold',
                     arrowprops=dict(arrowstyle='->', color=style['color'], lw=1.5),
                     bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=style['color'], alpha=0.9))
        offset_idx += 1

title_str = f"Операционные характеристики алгоритмов\nКласс задач: {family_name} (Всего задач: {total_tasks})"
plt.title(title_str, fontsize=15, pad=15)
plt.xlabel("Количество испытаний (итераций)", fontsize=12)
plt.ylabel("Процент решенных задач (%)", fontsize=12)
plt.xlim(0, max_iter)
plt.ylim(0, 115) # Увеличили Y, чтобы влезли стрелочки
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend(fontsize=11, loc='lower right')

plt.tight_layout()
plt.savefig("benchmark_chart.png", dpi=300)
plt.show()