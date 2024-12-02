import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import interp1d

# Читаем данные из файла
data_file_path = os.path.join("plot_data.txt")
try:
    data = np.loadtxt(data_file_path)
except OSError as e:
    print(f"Ошибка при чтении файла: {e}")
    exit()


# Разделяем данные для Хилла и Шекеля
problem_numbers = data[::2, 0]
hill_iterations = data[::2, 1]
shekel_iterations = data[1::2, 1]

num_problems = len(problem_numbers) # Общее число задач (в вашем примере 5)

# --- Сортировка данных ---
# Для Хилла
hill_data = np.column_stack((hill_iterations, problem_numbers))  # Объединяем итерации и номера задач
hill_data = hill_data[hill_data[:, 0].argsort()] # Сортируем по количеству итераций
hill_iterations_sorted = hill_data[:, 0]
hill_problem_numbers_sorted = hill_data[:, 1]

# Для Шекеля (аналогично)
shekel_data = np.column_stack((shekel_iterations, problem_numbers))
shekel_data = shekel_data[shekel_data[:, 0].argsort()]
shekel_iterations_sorted = shekel_data[:, 0]
shekel_problem_numbers_sorted = shekel_data[:, 1]
# --- Конец сортировки ---

# Создаем массивы для хранения процента/количества решенных задач
hill_solved = []
shekel_solved = []

for i in range(num_problems):
    hill_solved_count = i + 1  # Теперь просто i + 1, так как данные отсортированы
    hill_solved.append(hill_solved_count * 100 / num_problems)

    shekel_solved_count = i + 1
    shekel_solved.append(shekel_solved_count * 100 / num_problems)
    
# Интерполяция для Хилла
f_hill = interp1d(hill_iterations_sorted, hill_solved, kind='linear') #  kind='linear' для линейной интерполяции
x_hill_new = np.linspace(hill_iterations_sorted.min(), hill_iterations_sorted.max(), 500) # 500 новых точек
y_hill_new = f_hill(x_hill_new)



# Интерполяция для Шекеля (аналогично)
f_shekel = interp1d(shekel_iterations_sorted, shekel_solved, kind='linear')
x_shekel_new = np.linspace(shekel_iterations_sorted.min(), shekel_iterations_sorted.max(), 500)
y_shekel_new = f_shekel(x_shekel_new)

# Создаем два отдельных графика
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

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

plt.tight_layout()
plt.show()