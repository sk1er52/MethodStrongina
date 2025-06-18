import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.interpolate import interp1d
import os # Для проверки существования файлов

# --- Константы и пути к файлам ---
DATA_FILE_PATH = "plot_data.txt"
STATS_FILE_PATH = "stats.txt"

# --- Проверка существования файлов ---
if not os.path.exists(DATA_FILE_PATH):
    print(f"Ошибка: Файл данных '{DATA_FILE_PATH}' не найден.")
    exit()
if not os.path.exists(STATS_FILE_PATH):
    print(f"Ошибка: Файл статистики '{STATS_FILE_PATH}' не найден.")
    exit()

# --- Чтение данных из файла plot_data.txt ---
try:
    # Читаем только второй столбец (количество итераций)
    data = np.loadtxt(DATA_FILE_PATH, usecols=(1,))
    # Если был только один тест, loadtxt вернет скаляр, преобразуем в массив
    if data.ndim == 0:
        iterations = np.array([data])
    else:
        iterations = data
    num_tests_from_data = len(iterations)
    if num_tests_from_data == 0:
        print(f"Ошибка: Файл данных '{DATA_FILE_PATH}' пуст.")
        exit()
except OSError as e:
    print(f"Ошибка при чтении файла данных '{DATA_FILE_PATH}': {e}")
    exit()
except ValueError as e:
     print(f"Ошибка формата данных в '{DATA_FILE_PATH}': {e}")
     print("Возможно, файл пуст, содержит некорректные значения или неверное количество столбцов.")
     exit()
except IndexError:
     print(f"Ошибка: Не удалось прочитать второй столбец из '{DATA_FILE_PATH}'. Убедитесь, что файл содержит как минимум два столбца.")
     exit()


# --- Чтение статистики из файла stats.txt ---
try:
    with open(STATS_FILE_PATH, 'r') as f:
        # Читаем первую строку: общие счетчики выхода
        line1 = f.readline().split()
        if len(line1) != 2:
             raise ValueError(f"Ожидалось 2 значения счетчика выхода в первой строке {STATS_FILE_PATH}, найдено {len(line1)}")
        exit_main_count, exit_test_count = map(int, line1)

        # Читаем вторую строку: параметры epsilon, r, numTests
        line2 = f.readline().split()
        if len(line2) != 3:
             raise ValueError(f"Ожидалось 3 параметра (epsilon, r, numTests) во второй строке {STATS_FILE_PATH}, найдено {len(line2)}")
        epsilon, r, num_tests_from_stats = map(float, line2)
        num_tests_from_stats = int(num_tests_from_stats) # Количество тестов должно быть целым

        # (Опционально) Проверим, совпадает ли количество тестов
        if num_tests_from_data != num_tests_from_stats:
            print(f"Предупреждение: Количество тестов в '{DATA_FILE_PATH}' ({num_tests_from_data}) "
                  f"не совпадает с количеством в '{STATS_FILE_PATH}' ({num_tests_from_stats}).")
            print(f"График будет построен на основе {num_tests_from_data} тестов из '{DATA_FILE_PATH}'.")
        num_tests = num_tests_from_data # Используем фактическое количество результатов

except OSError as e:
    print(f"Ошибка при чтении файла статистики '{STATS_FILE_PATH}': {e}")
    exit()
except (ValueError, IndexError) as e:
     print(f"Ошибка формата или содержимого файла статистики '{STATS_FILE_PATH}': {e}")
     exit()

# --- Определение названия задачи (Заглушка!) ---
# Рекомендуется модифицировать C++ код, чтобы он записывал имя задачи
# (например, "Hill", "Shekel", "Grishagin") в stats.txt (например, на 3 строку),
# а затем раскомментировать и адаптировать блок чтения имени задачи здесь.
problem_name = "График сходимости метода" # Общее название
# try:
#     with open(STATS_FILE_PATH, 'r') as f:
#         lines = f.readlines()
#         if len(lines) >= 3:
#             name_line = lines[2].strip()
#             if name_line: # Если строка не пустая
#                 problem_name = name_line
# except Exception as e:
#      print(f"Не удалось прочитать имя задачи из {STATS_FILE_PATH}: {e}. Используется имя по умолчанию.")
#      pass # Оставить имя по умолчанию


# --- Функция для получения точек верхней границы ---
def get_cumulative_points(x, y):
    """
    Находит уникальные значения x и соответствующие им максимальные значения y.
    Это создает точки для верхней границы данных.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    # Сначала получим уникальные значения x (уже отсортированные)
    unique_x = np.unique(x)
    max_y_for_unique_x = []

    # Для каждого уникального x найдем максимальный y
    for x_val in unique_x:
        # Находим все y, соответствующие текущему x_val в исходных данных
        corresponding_y = y[x == x_val]
        if corresponding_y.size > 0: # Проверка, что массив не пустой
             max_y_for_unique_x.append(np.max(corresponding_y))
        # Иначе пропускаем это значение x (не должно происходить при норм. данных)

    return unique_x, np.array(max_y_for_unique_x)

# --- Подготовка данных для графика ---

# Сортируем итерации по возрастанию
iterations_sorted = np.sort(iterations)

# Создаем ось Y: процент решенных задач (от 0 до 100)
percent_solved = np.linspace(0, 100, num=num_tests, endpoint=True)

# Находим точки для верхней границы (макс. % для каждой итерации)
plot_iterations, plot_percentages = get_cumulative_points(iterations_sorted, percent_solved)

# --- ДОБАВЛЕНО: Добавляем точку (0, 0) в начало, если нужно ---
if len(plot_iterations) == 0:
    # Если данных нет, создаем точку (0,0) для пустого графика
    plot_iterations = np.array([0])
    plot_percentages = np.array([0])
elif plot_iterations[0] != 0:
    # Если первая точка не при x=0, добавляем (0,0) в начало
    plot_iterations = np.insert(plot_iterations, 0, 0)
    plot_percentages = np.insert(plot_percentages, 0, 0)
elif plot_percentages[0] != 0:
    # Если первая точка при x=0, но y > 0 (не должно быть с linspace),
    # принудительно ставим y=0
    plot_percentages[0] = 0
    
# --- Интерполяция (если достаточно точек) ---
can_interpolate = False
x_interp = np.array([]) # Инициализация пустыми массивами
y_interp = np.array([])
# Используем plot_iterations и plot_percentages для интерполяции
if len(plot_iterations) >= 2:
    try:
        # Убедимся, что точки уникальны по X для interp1d
        # (get_cumulative_points уже должна это гарантировать)
        interp_x = plot_iterations
        interp_y = plot_percentages

        # Проверяем, что после отбора уникальных осталось хотя бы 2 точки
        if len(interp_x) >= 2:
             f_interp = interp1d(interp_x, interp_y, kind='linear', fill_value="extrapolate")
             # Генерируем точки для плавной кривой интерполяции
             x_min = interp_x.min()
             x_max = interp_x.max()
             # Проверяем, что min и max не совпадают
             if x_max > x_min:
                 x_interp = np.linspace(x_min, x_max, 500)
                 y_interp = f_interp(x_interp)
                 # Ограничим интерполированные значения Y диапазоном [0, 100]
                 y_interp = np.clip(y_interp, 0, 100)
                 can_interpolate = True
             else: # Если все точки имеют одинаковое значение x
                 x_interp = interp_x
                 y_interp = interp_y
                 can_interpolate = True # Отобразим линию, соединяющую точки
                 print("Все тесты завершились за одинаковое число итераций. Интерполяция как линия.")
        else:
             print("Недостаточно уникальных точек по оси X для построения интерполяции (нужно минимум 2).")
             x_interp = interp_x # Оставляем исходные точки
             y_interp = interp_y

    except Exception as e:
        print(f"Ошибка во время интерполяции: {e}")
        can_interpolate = False
        x_interp = plot_iterations # На случай ошибки оставляем исходные точки
        y_interp = plot_percentages
else:
    print("Недостаточно точек для построения интерполяции (нужно минимум 2).")
    x_interp = plot_iterations # Оставляем исходные точки
    y_interp = plot_percentages


# --- Вывод информации на экран ---
print("\n" + "="*30)
print("      Статистика Запуска")
print("="*30)
print(f"Параметры:")
print(f"  Точность (epsilon): {epsilon}")
print(f"  Параметр надежности (r): {r}")
print(f"  Количество тестов: {num_tests}")
print(f"\nИнформация о выходе из цикла:")
print(f"  Выход по основному условию (exitMainCount): {exit_main_count}")
print(f"  Выход по доп. условию (exitTestCount):   {exit_test_count}")
mean_iterations_val = iterations.mean() if num_tests > 0 else 0
print(f"\nСреднее число итераций: {mean_iterations_val:.2f}")
print("="*30 + "\n")

# --- Построение графика ---
fig, ax = plt.subplots(figsize=(10, 7)) # Один график

# Общий заголовок с параметрами
fig.suptitle(f"Параметры: r = {r:.2f}, epsilon = {epsilon:.4f}", fontsize=14)

# 1. Точки реальных данных (ВЕРХНЯЯ ГРАНИЦА)
# Используем plot_iterations и plot_percentages
ax.plot(plot_iterations, plot_percentages,
        marker='o', linestyle='', markersize=5, # Немного увеличим маркер
        color='blue', zorder=3, # zorder=3 помещает точки поверх линии
        label=f'Результаты {num_tests} тестов (точки)')

# 2. Интерполяционная кривая (если была построена или есть точки)
if can_interpolate and len(x_interp) > 0:
    ax.plot(x_interp, y_interp,
            linestyle='-', color='red', linewidth=2, # Немного утолщим линию
            zorder=2, # zorder=2 помещает линию под точками
            label='Интерполяция')
elif len(plot_iterations) > 0: # Если интерполяции нет, но есть точки
     print("Отрисовка без интерполяционной кривой.")
     # Можем просто не рисовать красную линию


ax.set_xlabel("Количество итераций")
ax.set_ylabel("Процент решенных задач (%)")
ax.set_title(problem_name)
ax.grid(True, linestyle='--', alpha=0.6)
# Проверяем, есть ли что добавить в легенду, чтобы не вызывать ее для пустого графика
handles, labels = ax.get_legend_handles_labels()
if handles:
    ax.legend(loc='lower right') # Разместим легенду внизу справа

# Установим пределы осей для лучшего вида
ax.set_ylim(-5, 105) # Начнем чуть ниже 0% и закончим чуть выше 100%
max_iter = plot_iterations.max() if len(plot_iterations) > 0 else 100
min_x_limit = -max_iter*0.02 if max_iter > 0 else -1 # Начнем чуть левее 0
max_x_limit = max_iter * 1.05 if max_iter > 0 else 10 # Закончить немного правее максимума
ax.set_xlim(left=min_x_limit,
            right=max_x_limit)

# Добавление текста со статистикой на график
stats_text = (
    f"Среднее итер.: {mean_iterations_val:.2f}\n"
    f"Выход (осн.): {exit_main_count}\n"
    f"Выход (доп.): {exit_test_count}"
)
ax.text(0.03, 0.97, stats_text,
        transform=ax.transAxes, verticalalignment='top', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.4', fc='white', alpha=0.8))

plt.tight_layout(rect=[0, 0.03, 1, 0.93]) # Отрегулируем отступы, учитывая suptitle
plt.show()

print("График построен и отображен.")