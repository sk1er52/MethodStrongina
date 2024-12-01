import matplotlib.pyplot as plt
import numpy as np
import os

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


# Создаем два отдельных графика
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6)) #  figsize  -  размер  графиков


# График для Хилла
ax1.plot(problem_numbers, hill_iterations)
ax1.set_xlabel("Количество решенных задач")
ax1.set_ylabel("Количество итераций")
ax1.set_title("Функция Хилла")
ax1.grid(True)


# График для  Шекеля
ax2.plot(problem_numbers, shekel_iterations)
ax2.set_xlabel("Количество решенных задач")
ax2.set_ylabel("Количество итераций")
ax2.set_title("Функция Шекеля")
ax2.grid(True)



#  Отображаем  графики
plt.tight_layout()  #  Чтобы  графики  не  перекрывались
plt.show()