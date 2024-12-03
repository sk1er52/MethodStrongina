import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import interp1d

# Читаем данные из файла
data_file_path = os.path.join("plot_data.txt")
stats_file_path = os.path.join("stats.txt")
try:
    data = np.loadtxt(data_file_path)
except OSError as e:
    print(f"Ошибка при чтении файла: {e}")
    exit()
    
# Разделяем данные для Хилла и Шекеля
problem_numbers = data[::2, 0].astype(int)
hill_iterations = data[::2, 1]
shekel_iterations = data[1::2, 1]

num_problems = len(problem_numbers) # Общее число задач 

# --- Сортировка данных ---
# Для Хилла
hill_data = np.column_stack((hill_iterations, problem_numbers))  # Объединяем итерации и номера задач
hill_data = hill_data[hill_data[:, 0].argsort()] # Сортируем по количеству итераций
hill_iterations_sorted = hill_data[:, 0]
hill_problem_numbers_sorted = hill_data[:, 0]

# Для Шекеля (аналогично)
shekel_data = np.column_stack((shekel_iterations, problem_numbers))
shekel_data = shekel_data[shekel_data[:, 0].argsort()]
shekel_iterations_sorted = shekel_data[:, 0]
shekel_problem_numbers_sorted = shekel_data[:, 0]
# --- Конец сортировки ---

# Создаем массивы для хранения процента/количества решенных задач
hill_solved = []
shekel_solved = []

for i in range(num_problems):
    hill_solved_count = i + 1  # Теперь просто i + 1, так как данные отсортированы
    hill_solved.append(hill_solved_count * 100 / num_problems)

    shekel_solved_count = i + 1
    shekel_solved.append(shekel_solved_count * 100 / num_problems)
    
# --- Вычисление среднего количества итераций ---
avg_hill_iterations = np.mean(hill_iterations)
avg_shekel_iterations = np.mean(shekel_iterations)
    
# Интерполяция для Хилла
f_hill = interp1d(hill_iterations_sorted, hill_solved, kind='linear') #  kind='linear' для линейной интерполяции
x_hill_new = np.linspace(hill_iterations_sorted.min(), hill_iterations_sorted.max(), 500) # 500 новых точек
y_hill_new = f_hill(x_hill_new)



# Интерполяция для Шекеля (аналогично)
f_shekel = interp1d(shekel_iterations_sorted, shekel_solved, kind='linear')
x_shekel_new = np.linspace(shekel_iterations_sorted.min(), shekel_iterations_sorted.max(), 500)
y_shekel_new = f_shekel(x_shekel_new)

try:
    with open(stats_file_path, 'r') as f:
        exit_counts = list(map(int, f.readline().split()))
        epsilon, r, numTests = map(float, f.readline().split())
        hill_exit1_count, shekel_exit1_count, hill_exit2_count, shekel_exit2_count = exit_counts
except OSError as e:
    print(f"Ошибка при чтении файла статистики: {e}")
    exit()
    
# --- Выводим информацию на экран ---
print("Информация о выходе из цикла:")
print(f"  Hill (Основное условие): {hill_exit1_count}")
print(f"  Hill (Доп. условие): {hill_exit2_count}")
print(f"  Shekel (Основное условие): {shekel_exit1_count}")
print(f"  Shekel (Доп. условие): {shekel_exit2_count}")
print("\nПараметры:")
print(f"  Epsilon: {epsilon}")
print(f"  r: {r}")
print(f"  Количество тестов: {int(numTests)}")

# Создаем два отдельных графика
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# --- Вывод r и epsilon НАД графиками ---
fig.suptitle(f"Параметры: r = {r}, epsilon = {epsilon}", fontsize=14)

# График для Хилла
ax1.plot(hill_iterations_sorted, hill_solved, marker='o', markersize=3, label='Реальные данные') # Маркеры  для  реальных  данных
#  Интерполяция
f_hill = interp1d(hill_iterations_sorted, hill_solved, kind='linear')
x_hill_new = np.linspace(hill_iterations_sorted.min(), hill_iterations_sorted.max(), 500)
y_hill_new = f_hill(x_hill_new)
ax1.plot(x_hill_new, y_hill_new, label='Интерполяция')

ax1.set_xlabel("Количество итераций")
ax1.set_ylabel("Процент решенных задач")
ax1.set_title("Функция Хилла")
ax1.grid(True)
ax1.legend()

# График для Шекеля
ax2.plot(shekel_iterations_sorted, shekel_solved, marker='o', markersize=3, label='Реальные данные') # Маркеры  для  реальных  данных
#  Интерполяция
f_shekel = interp1d(shekel_iterations_sorted, shekel_solved, kind='linear')
x_shekel_new = np.linspace(shekel_iterations_sorted.min(), shekel_iterations_sorted.max(), 500)
y_shekel_new = f_shekel(x_shekel_new)
ax2.plot(x_shekel_new, y_shekel_new, label='Интерполяция')

ax2.set_xlabel("Количество итераций")
ax2.set_ylabel("Процент решенных задач")
ax2.set_title("Функция Шекеля")
ax2.grid(True)
ax2.legend()

# --- Добавление текста на графики ---
ax1.text(0.05, 0.95, f"Среднее число итераций: {avg_hill_iterations:.2f}\nОсновное условие 1: {hill_exit1_count}, Тест 2: {hill_exit2_count}",
        transform=ax1.transAxes, verticalalignment='top')

ax2.text(0.05, 0.95, f"Среднее число итераций: {avg_shekel_iterations:.2f}\nОсновное условие 1: {shekel_exit1_count}, Тест 2: {shekel_exit2_count}",
        transform=ax2.transAxes, verticalalignment='top')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()