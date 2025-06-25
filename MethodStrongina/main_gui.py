import sys
import subprocess
import json
import os
import random 
import numpy as np 
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
                             QPushButton, QComboBox, QTextEdit, QMainWindow, QGridLayout, QGroupBox,
                             QSizePolicy, QSpinBox, QMessageBox) 
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtCore import Qt

import matplotlib.pyplot as plt

# Предполагаемые размеры семейств задач (УТОЧНИТЕ ЭТИ ЗНАЧЕНИЯ!)
FAMILY_SIZES = {
    "Hill": 10,      # Пример, сколько задач в THillProblemFamily
    "Shekel": 10,    # Пример, сколько задач в TShekelProblemFamily
    "Grishagin": 100 # Обычно 100 задач Гришагина
}


class OptimizerUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Метод Стронгина UI (PyQt)")
        self.setGeometry(100, 100, 900, 800)

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QHBoxLayout(self.central_widget) 

        self.left_panel_widget = QWidget()
        self.left_panel_layout = QVBoxLayout(self.left_panel_widget)
        
        self.right_panel_widget = QWidget() 
        self.right_panel_layout = QVBoxLayout(self.right_panel_widget)
        self.right_panel_widget.setFixedWidth(280) 

        self.main_layout.addWidget(self.left_panel_widget, 2) 
        self.main_layout.addWidget(self.right_panel_widget, 1) 

        script_dir = os.path.dirname(os.path.abspath(__file__))
        # --- ВАЖНО: Настройте этот блок для поиска вашего .exe ---
        potential_paths = [
            os.path.join(script_dir, "build", "bin", "Release", "OptimizationProject.exe"),
            os.path.join(script_dir, "build", "bin", "Debug", "OptimizationProject.exe"),
            os.path.join(script_dir, "build", "Release", "OptimizationProject.exe"),
            os.path.join(script_dir, "build", "Debug", "OptimizationProject.exe"),
            os.path.join(script_dir, "Release", "OptimizationProject.exe"), 
            os.path.join(script_dir, "Debug", "OptimizationProject.exe"),   
            os.path.join(script_dir, "OptimizationProject.exe"), 
            "OptimizationProject.exe" 
        ]
        base_project_dir = os.path.normpath(os.path.join(script_dir, "..")) # Если скрипт в подкаталоге проекта
        potential_paths.insert(0, os.path.join(base_project_dir, "build", "Release", "OptimizationProject.exe"))
        potential_paths.insert(1, os.path.join(base_project_dir, "build", "Debug", "OptimizationProject.exe"))
        potential_paths.insert(2, os.path.join(base_project_dir, "build", "bin", "Release", "OptimizationProject.exe"))
        potential_paths.insert(3, os.path.join(base_project_dir, "build", "bin", "Debug", "OptimizationProject.exe"))
        potential_paths.insert(4, os.path.join(base_project_dir, "build", "bin", "OptimizationProject.exe"))


        self.cpp_executable_path = ""
        for path_check in potential_paths:
            normalized_path = os.path.normpath(path_check)
            if os.path.isfile(normalized_path):
                self.cpp_executable_path = normalized_path
                print(f"INFO: Найден C++ exe: {self.cpp_executable_path}")
                break
        
        if not self.cpp_executable_path:
            print(f"КРИТИЧЕСКАЯ ОШИБКА: Исполняемый файл C++ НЕ НАЙДЕН. Проверены варианты, например: {potential_paths[0]}")
            # QMessageBox будет показан в __main__ после создания окна.

        self.current_plot_data_file = None 
        self.last_successful_run_data = None 
        self.all_runs_data_for_stats = [] 

        self._create_input_widgets_definitions()
        self._create_plot_buttons() 
        self._create_input_widgets_layout()      
        self._create_output_widgets()            
        
        self._on_parameter_changed()

    def _create_input_widgets_definitions(self):
        self.task_type_combo = QComboBox()
        self.task_type_combo.addItems(["1: Одномерная", "2: Многомерная (Гришагина)"])
        self.task_type_combo.setCurrentIndex(1) 
        self.task_type_combo.currentIndexChanged.connect(self._on_parameter_changed)

        self.sub_choice_label = QLabel("Подтип (1D):")
        self.sub_choice_combo = QComboBox()
        self.sub_choice_combo.addItems(["1: Хилл", "2: Шекель"])
        self.sub_choice_combo.currentIndexChanged.connect(self._on_parameter_changed)

        self.num_runs_spinbox = QSpinBox()
        self.num_runs_spinbox.setMinimum(1); self.num_runs_spinbox.setMaximum(1000); self.num_runs_spinbox.setValue(1)
        self.num_runs_spinbox.valueChanged.connect(self._on_parameter_changed)

        self.problem_idx_label = QLabel("Индекс задачи (1-based):")
        self.problem_idx_entry = QLineEdit("1") 
        self.problem_idx_entry.setValidator(QIntValidator(1, max(FAMILY_SIZES.values(), default=100)))
        
        self.epsilon_entry = QLineEdit("0.01") 
        self.epsilon_entry.setValidator(QDoubleValidator(0.000000001, 1.0, 9, notation=QDoubleValidator.StandardNotation))
        
        self.r_param_entry = QLineEdit("3.0")
        self.r_param_entry.setValidator(QDoubleValidator(0.1, 100.0, 2, notation=QDoubleValidator.StandardNotation))
        
        self.peano_order_label = QLabel("Порядок Пеано (m):")
        self.peano_order_entry = QLineEdit("12") 
        self.peano_order_entry.setValidator(QIntValidator(1, 30)) 
        
        self.max_iter_label = QLabel("Макс. итераций:")
        self.max_iter_entry = QLineEdit("1000") 
        self.max_iter_entry.setValidator(QIntValidator(1, 200000))
        
        self.calculate_button = QPushButton("Посчитать")
        self.calculate_button.clicked.connect(self.run_calculation_handler)

    def _create_input_widgets_layout(self):
        input_group = QGroupBox("Параметры задачи")
        self.input_layout = QGridLayout() 
        current_row = 0
        self.input_layout.addWidget(QLabel("Тип задачи:"), current_row, 0, Qt.AlignTop)
        self.input_layout.addWidget(self.task_type_combo, current_row, 1); current_row += 1
        self.input_layout.addWidget(self.sub_choice_label, current_row, 0, Qt.AlignTop)
        self.input_layout.addWidget(self.sub_choice_combo, current_row, 1); current_row += 1
        self.input_layout.addWidget(QLabel("Количество запусков:"), current_row, 0, Qt.AlignTop)
        self.input_layout.addWidget(self.num_runs_spinbox, current_row, 1); current_row += 1
        self.input_layout.addWidget(self.problem_idx_label, current_row, 0, Qt.AlignTop)
        self.input_layout.addWidget(self.problem_idx_entry, current_row, 1); current_row += 1
        self.input_layout.addWidget(QLabel("Эпсилон:"), current_row, 0, Qt.AlignTop)
        self.input_layout.addWidget(self.epsilon_entry, current_row, 1); current_row += 1
        self.input_layout.addWidget(QLabel("Параметр r:"), current_row, 0, Qt.AlignTop)
        self.input_layout.addWidget(self.r_param_entry, current_row, 1); current_row += 1
        self.input_layout.addWidget(self.peano_order_label, current_row, 0, Qt.AlignTop)
        self.input_layout.addWidget(self.peano_order_entry, current_row, 1); current_row += 1
        self.input_layout.addWidget(self.max_iter_label, current_row, 0, Qt.AlignTop)
        self.input_layout.addWidget(self.max_iter_entry, current_row, 1); current_row += 1
        self.input_layout.addWidget(self.calculate_button, current_row, 0, 1, 2)
        input_group.setLayout(self.input_layout)
        self.left_panel_layout.addWidget(input_group)
        self.left_panel_layout.addStretch(1)

    def _create_output_widgets(self):
        output_group = QGroupBox("Результаты вычислений")
        output_layout = QVBoxLayout()
        self.result_display = QTextEdit(); self.result_display.setReadOnly(True)
        self.result_display.setPlaceholderText("Результаты расчета появятся здесь...")
        self.result_display.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        output_layout.addWidget(self.result_display); output_group.setLayout(output_layout)
        self.left_panel_layout.addWidget(output_group)
        self.left_panel_layout.setStretchFactor(output_group, 1)

    def _create_plot_buttons(self): 
        plot_group = QGroupBox("Графики")
        plot_layout = QVBoxLayout()
        
        self.single_run_plot_button = QPushButton("График сходимости (plot_data.txt)")
        self.single_run_plot_button.clicked.connect(self.show_single_run_convergence_plot)
        plot_layout.addWidget(self.single_run_plot_button)

        self.plot_level_lines_button = QPushButton("Линии уровня и точки (2D)")
        self.plot_level_lines_button.clicked.connect(self.show_level_lines_plot)
        self.plot_level_lines_button.setEnabled(False) 
        plot_layout.addWidget(self.plot_level_lines_button)
        
        self.op_char_plot_button = QPushButton("График рабочих характеристик (серия)")
        self.op_char_plot_button.clicked.connect(self.show_operating_characteristics_plot)
        self.op_char_plot_button.setEnabled(False)
        plot_layout.addWidget(self.op_char_plot_button)

        plot_group.setLayout(plot_layout)
        self.right_panel_layout.addWidget(plot_group)
        self.right_panel_layout.addStretch(1)

    def _on_parameter_changed(self):
        is_1d_type = (self.task_type_combo.currentText().startswith("1:"))
        self.sub_choice_label.setVisible(is_1d_type)
        self.sub_choice_combo.setVisible(is_1d_type)
        self.peano_order_label.setVisible(not is_1d_type)
        self.peano_order_entry.setVisible(not is_1d_type)
        
        num_runs = self.num_runs_spinbox.value()
        is_single_run = (num_runs == 1)
        self.problem_idx_label.setEnabled(is_single_run)
        self.problem_idx_entry.setEnabled(is_single_run)
        
        can_plot_levels = (not is_1d_type and 
                           is_single_run and 
                           self.last_successful_run_data is not None and
                           self.last_successful_run_data.get("problem_family") == "Grishagin")
        self.plot_level_lines_button.setEnabled(can_plot_levels)
        self.op_char_plot_button.setEnabled(num_runs > 1 and bool(self.all_runs_data_for_stats))
        self.single_run_plot_button.setEnabled(is_single_run)

        max_idx_val = 1 
        if is_1d_type:
            family_key = self.sub_choice_combo.currentText().split(':')[1].strip()
            max_idx_val = FAMILY_SIZES.get(family_key, 10) 
            self.problem_idx_label.setText(f"Индекс задачи (1D, 1-{max_idx_val}):")
        else:
            family_key = "Grishagin"
            max_idx_val = FAMILY_SIZES.get(family_key, 100)
            self.problem_idx_label.setText(f"Индекс задачи ({family_key}, 1-{max_idx_val}):")
        
        current_idx_validator = self.problem_idx_entry.validator()
        if isinstance(current_idx_validator, QIntValidator):
            current_idx_validator.setTop(max_idx_val if max_idx_val > 0 else 1)

    def show_single_run_convergence_plot(self):
        if not self.current_plot_data_file or not os.path.exists(self.current_plot_data_file):
            self.result_display.append(f"Файл данных для графика ({self.current_plot_data_file}) не найден. Сначала выполните расчет.")
            return
        try:
            plot_x, plot_y = [], []
            with open(self.current_plot_data_file, "r") as f:
                for line_num, line_content in enumerate(f,1):
                    parts = line_content.strip().split()
                    if len(parts) == 2:
                        try: plot_x.append(int(parts[0])); plot_y.append(float(parts[1])) 
                        except ValueError: self.result_display.append(f"Предупреждение: Не удалось прочитать строку {line_num} данных для графика: {line_content.strip()}")
            if not plot_x or not plot_y: self.result_display.append(f"Нет данных для построения в файле {self.current_plot_data_file}"); return
            plt.figure(f"График из plot_data.txt ({os.path.basename(self.current_plot_data_file)})")
            plt.plot(plot_x, plot_y, marker='o', linestyle='-')
            plt.xlabel("Номер точки в файле (из plot_data.txt)") 
            plt.ylabel("Количество итераций (из plot_data.txt)") 
            plt.title(f"График данных из {os.path.basename(self.current_plot_data_file)}")
            plt.grid(True); plt.tight_layout(); plt.show() 
            self.result_display.append(f"График из {self.current_plot_data_file} построен.")
        except Exception as e: self.result_display.append(f"Ошибка при построении графика из plot_data.txt: {type(e).__name__} - {e}")
    
    def show_operating_characteristics_plot(self):
        if not self.all_runs_data_for_stats:
            self.result_display.append("Нет данных для построения графика характеристик. Выполните серию расчетов (Количество запусков > 1).")
            return
        iterations_list = [rd.get("iterations",0) for rd in self.all_runs_data_for_stats if isinstance(rd,dict) and rd.get("iterations") is not None]
        if not iterations_list:
            self.result_display.append("Нет данных по итерациям для построения графика характеристик."); return

        sorted_iterations = np.sort(iterations_list)
        y_values = np.arange(1, len(sorted_iterations) + 1) / len(sorted_iterations) * 100 
        
        plt.figure("Рабочие характеристики (Operating Characteristics)")
        plt.plot(sorted_iterations, y_values, marker='.', linestyle='-', drawstyle='steps-post') 
        plt.xlabel("Количество итераций (NFE)"); plt.ylabel("Доля решенных задач P(NFE), %") 
        plt.title(f"Рабочие характеристики для {len(sorted_iterations)} запусков")
        plt.xlim(left=0); plt.ylim(bottom=0, top=105) 
        plt.grid(True, which='both', linestyle='--', linewidth='0.5', color='gray')
        plt.minorticks_on()
        try:
            epsilon_ui = self.epsilon_entry.text(); r_param_ui = self.r_param_entry.text()
            peano_m_ui = self.peano_order_entry.text() if not self.task_type_combo.currentText().startswith("1:") else "N/A"
            task_type_for_plot_full = self.task_type_combo.currentText()
            if task_type_for_plot_full.startswith("1:"): task_type_for_plot = self.sub_choice_combo.currentText().split(':')[1].strip()
            else: task_type_for_plot = task_type_for_plot_full.split(':')[1].strip()
            
            params_text = (f"Eps={str(epsilon_ui)}, R={str(r_param_ui)}, Peano M={str(peano_m_ui)}\n"
                           f"Задач: {len(sorted_iterations)}, Тип: {task_type_for_plot}")
            plt.text(0.98, 0.02, params_text, transform=plt.gca().transAxes, fontsize=8, 
                     verticalalignment='bottom', horizontalalignment='right',
                     bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.7))
        except Exception as e_text: print(f"Ошибка при добавлении текста на график: {e_text}")
        plt.tight_layout(); plt.show()
        self.result_display.append("График рабочих характеристик построен.")

    def show_level_lines_plot(self):
        self.result_display.append("\nЗапрос данных для линий уровня..."); QApplication.processEvents()
        
        if self.num_runs_spinbox.value() != 1 or not self.last_successful_run_data: 
            self.result_display.append("Ошибка: Линии уровня строятся только для одиночного успешного расчета задачи Гришагина.")
            return
        
        problem_idx_cpp = self.last_successful_run_data.get("problem_index_cpp")
        problem_family = self.last_successful_run_data.get("problem_family")

        if problem_family != "Grishagin" or problem_idx_cpp is None: 
            self.result_display.append("Ошибка: Линии уровня доступны только для успешно решенной задачи Гришагина.")
            return

        grid_resolution = 100 
        command_level_lines = [self.cpp_executable_path, "--get-level-lines", str(problem_idx_cpp), str(grid_resolution)]
        
        self.result_display.append(f"Команда для C++ (линии уровня): {' '.join(command_level_lines)}"); QApplication.processEvents()
        try:
            exe_dir = os.path.dirname(os.path.abspath(self.cpp_executable_path))
            if not exe_dir: exe_dir = os.getcwd()
            
            creation_flags_val = subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0
            level_lines_process = subprocess.run(command_level_lines, capture_output=True, text=True, 
                                                 encoding='utf-8', errors='replace', cwd=exe_dir, 
                                                 timeout=60, creationflags=creation_flags_val)

            if level_lines_process.returncode == 0 and level_lines_process.stdout:
                try:
                    level_lines_data = json.loads(level_lines_process.stdout)
                    if level_lines_data.get("type") != "level_lines_data": 
                        self.result_display.append(f"Ошибка: C++ вернул JSON неизвестного типа для линий уровня: {level_lines_data.get('type')}")
                        return
                    
                    x_grid_coords = np.array(level_lines_data.get("x_coords_grid",[]))
                    y_grid_coords = np.array(level_lines_data.get("y_coords_grid",[]))
                    z_grid_values = np.array(level_lines_data.get("z_values_grid",[]))

                    if not (x_grid_coords.size > 0 and y_grid_coords.size > 0 and z_grid_values.size > 0 and \
                            z_grid_values.shape == (len(y_grid_coords), len(x_grid_coords))):
                        self.result_display.append(f"Ошибка: Некорректные размеры массивов данных для линий уровня. XG:{x_grid_coords.shape}, YG:{y_grid_coords.shape}, ZG:{z_grid_values.shape}")
                        return
                    
                    plt.figure(f"Линии уровня: {problem_family} Задача {problem_idx_cpp+1}", figsize=(8.5,7))
                    X_mesh, Y_mesh = np.meshgrid(x_grid_coords, y_grid_coords)
                    
                    contour_fill = plt.contourf(X_mesh, Y_mesh, z_grid_values, levels=30, cmap='viridis') 
                    plt.colorbar(contour_fill, label='Значение Z')
                    plt.contour(X_mesh, Y_mesh, z_grid_values, levels=contour_fill.levels, colors='black', linewidths=0.5, alpha=0.7)
                    
                    trial_history_y = self.last_successful_run_data.get("trial_history_y_coords", []) 
                    if trial_history_y and all(isinstance(p, list) and len(p) == 2 for p in trial_history_y):
                        y1_trials = [p[0] for p in trial_history_y]
                        y2_trials = [p[1] for p in trial_history_y]
                        plt.scatter(y1_trials, y2_trials, s=15, color='black', alpha=0.5, label=f'Исслед.т. ({len(trial_history_y)})', zorder=3)

                    found_y = self.last_successful_run_data.get("found_minimum_y", [])
                    if len(found_y) == 2:
                        plt.scatter(found_y[0], found_y[1], s=150, color='red', marker='*', edgecolor='black', label='Найден min', zorder=5)
                    
                    known_y = self.last_successful_run_data.get("known_optimum_y", [])
                    if len(known_y) == 2:
                        plt.scatter(known_y[0], known_y[1], s=150, color='lime', marker='P', edgecolor='black', label='Известен opt', zorder=5)
                    
                    plt.xlabel("Y1"); plt.ylabel("Y2")
                    plt.title(f"Линии уровня: {problem_family} Задача {problem_idx_cpp+1}")
                    plt.legend(loc='upper right'); plt.axis('scaled'); 
                    if x_grid_coords.size > 0 : plt.xlim(x_grid_coords.min(), x_grid_coords.max())
                    if y_grid_coords.size > 0 : plt.ylim(y_grid_coords.min(), y_grid_coords.max())
                    plt.tight_layout(); plt.show()
                    self.result_display.append("График линий уровня построен.")

                except json.JSONDecodeError as e_json: 
                    self.result_display.append(f"Ошибка парсинга JSON (линии уровня): {e_json}\nStdout C++:\n{level_lines_process.stdout}")
                except Exception as e_plot: 
                    self.result_display.append(f"Ошибка при построении графика линий уровня: {type(e_plot).__name__} - {e_plot}")
            else:
                self.result_display.append(f"Ошибка C++ при генерации данных для линий уровня (код {level_lines_process.returncode}):\nStderr:\n{level_lines_process.stderr}\nStdout:\n{level_lines_process.stdout}")
        except Exception as e_general: 
            self.result_display.append(f"Общая ошибка при запросе линий уровня: {type(e_general).__name__} - {e_general}")

    def run_calculation_handler(self):
        self.result_display.clear(); self.all_runs_data_for_stats.clear(); 
        self.last_successful_run_data=None 
        self._on_parameter_changed() 

        num_runs = self.num_runs_spinbox.value()
        self.result_display.append(f"Запускается серия из {num_runs} расчетов..."); QApplication.processEvents()
        
        for i in range(num_runs):
            self.result_display.append(f"\n--- Расчет {i+1} из {num_runs} ---"); QApplication.processEvents()
            task_type_str = self.task_type_combo.currentText().split(":")[0]
            sub_choice_str, current_problem_family_key = "", ""
            
            if task_type_str == "1": 
                sub_choice_str = self.sub_choice_combo.currentText().split(":")[0]
                current_problem_family_key = self.sub_choice_combo.currentText().split(':')[1].strip()
            elif task_type_str == "2": 
                current_problem_family_key = "Grishagin"
            else: self.result_display.append("  Неизвестный тип задачи."); continue
            
            max_idx_for_family = FAMILY_SIZES.get(current_problem_family_key, 1)
            if max_idx_for_family <= 0: max_idx_for_family = 1 

            current_problem_idx_ui = int(self.problem_idx_entry.text()) if num_runs == 1 else random.randint(1, max_idx_for_family)
            
            if num_runs > 1: self.result_display.append(f"  Случайно выбран: {current_problem_family_key}, индекс (1-based): {current_problem_idx_ui}")
            
            stdout, stderr, ret_code = self._execute_single_cpp_run(
                task_type_str, sub_choice_str, current_problem_idx_ui, 
                self.epsilon_entry.text(), self.r_param_entry.text(),
                self.peano_order_entry.text(), self.max_iter_entry.text())
            
            if stdout is None and stderr is None: 
                self.result_display.append("  Критическая ошибка подготовки C++ вызова (stdout/stderr is None)."); continue

            if ret_code == 0 and stdout:
                try:
                    cpp_results = json.loads(stdout)
                    cpp_results['problem_index_ui_selected'] = current_problem_idx_ui 
                    self.all_runs_data_for_stats.append(cpp_results)
                    if num_runs == 1: 
                        self.last_successful_run_data = cpp_results 
                    
                    problem_index_from_cpp = cpp_results.get('problem_index_cpp', -1)
                    # index_to_display_ui = problem_index_from_cpp + 1 if problem_index_from_cpp != -1 else "N/A"

                    result_str = (f"  Семейство: {cpp_results.get('problem_family','N/A')}\n"
                                  f"  Индекс (1-based UI): {current_problem_idx_ui} (C++ 0-based: {problem_index_from_cpp})\n"
                                  f"  Найденный Z: {cpp_results.get('found_minimum_z','N/A')}\n"
                                  f"  Найденный Y: {cpp_results.get('found_minimum_y','N/A')}\n"
                                  f"  Известный Z*: {cpp_results.get('known_optimum_z','N/A')}\n"
                                  f"  Известный Y*: {cpp_results.get('known_optimum_y','N/A')}\n"
                                  f"  Итераций: {cpp_results.get('iterations','N/A')}\n"
                                  f"  Время: {cpp_results.get('time_ms','N/A')} мс\n"
                                  f"  Выход осн.усл.: {cpp_results.get('exit_main_count','N/A')}\n"
                                  f"  Выход доп.усл.: {cpp_results.get('exit_test_count','N/A')}\n")
                    # Если C++ передает историю точек для линий уровня (добавить ключ "trial_history_y_coords" в JSON C++)
                    # if "trial_history_y_coords" in cpp_results:
                    #     result_str += f"  Получено точек истории: {len(cpp_results['trial_history_y_coords'])}\n"

                    self.result_display.append("  Расчет успешно завершен.\n" + result_str)
                    
                    if os.path.exists(self.current_plot_data_file): 
                        self.result_display.append(f"  Файл {self.current_plot_data_file} обновлен.")
                    else: 
                        self.result_display.append(f"  ВНИМАНИЕ: Файл {self.current_plot_data_file} не был создан C++ программой.")
                
                except json.JSONDecodeError as e: 
                    self.result_display.append(f"  Ошибка парсинга JSON: {e}\nStdout:\n{stdout}\nStderr:\n{stderr}")
                except Exception as e: 
                    self.result_display.append(f"  Ошибка обработки вывода C++: {e}\nStdout:\n{stdout}\nStderr:{stderr}")
            else: 
                self.result_display.append(f"  Ошибка C++ (код {ret_code}):\nStderr:{stderr}\nStdout:{stdout}")
            
            QApplication.processEvents()
            if num_runs > 1 and i < num_runs-1: self.result_display.append("-" * 30)
        
        if num_runs > 1 and self.all_runs_data_for_stats: 
            self._display_aggregated_stats()
        
        self._on_parameter_changed() 
        self.result_display.append("\nВсе расчеты завершены.")

    def _execute_single_cpp_run(self, task_type_str, sub_choice_str, problem_idx_ui, 
                                epsilon_str_ui, r_param_str_ui, peano_order_m_str_ui, max_iter_str_ui):
        try:
            current_exe_path = self.cpp_executable_path
            if not os.path.isfile(current_exe_path):
                script_dir = os.path.dirname(os.path.abspath(__file__)); basename_exe = os.path.basename(self.cpp_executable_path) if self.cpp_executable_path else "OptimizationProject.exe"
                potential_paths = [os.path.join(script_dir, basename_exe), os.path.normpath(os.path.join(script_dir,"..","build","Release",basename_exe)), os.path.normpath(os.path.join(script_dir,"..","build","Debug",basename_exe)), os.path.normpath(os.path.join(script_dir,"Release",basename_exe)), os.path.normpath(os.path.join(script_dir,"Debug",basename_exe)), os.path.normpath(os.path.join(script_dir,"build","bin","Release",basename_exe)), os.path.normpath(os.path.join(script_dir,"build","bin","Debug",basename_exe)), os.path.normpath(os.path.join(script_dir,"build","bin",basename_exe))]
                base_project_dir = os.path.normpath(os.path.join(script_dir, ".."))
                potential_paths.insert(0, os.path.join(base_project_dir, "build", "Release", basename_exe))
                potential_paths.insert(1, os.path.join(base_project_dir, "build", "Debug", basename_exe))
                potential_paths.insert(2, os.path.join(base_project_dir, "build", "bin", "Release", basename_exe))
                potential_paths.insert(3, os.path.join(base_project_dir, "build", "bin", "Debug", basename_exe))
                found_exe = False
                for p_path in potential_paths:
                    if os.path.isfile(p_path): current_exe_path = os.path.abspath(p_path); found_exe = True; break
                if not found_exe: self.result_display.append(f"  КРИТ.ОШИБКА: exe C++ НЕ НАЙДЕН при попытке запуска. Проверенный путь: {self.cpp_executable_path}"); return None, "C++ exe not found during execution", -1
            
            # self.result_display.append(f"  Используется C++ exe: {current_exe_path}")

            problem_idx_cpp = problem_idx_ui - 1 
            epsilon_val_str, r_param_val_str = epsilon_str_ui.replace(',','.'), r_param_str_ui.replace(',','.')
            command = [current_exe_path, task_type_str]
            if task_type_str=="1": 
                command.extend([sub_choice_str,epsilon_val_str,r_param_val_str,str(problem_idx_cpp),max_iter_str_ui])
            elif task_type_str=="2": 
                command.extend([str(problem_idx_cpp),epsilon_val_str,r_param_val_str,peano_order_m_str_ui,max_iter_str_ui])
            else: 
                return None, "Invalid task_type_str in _execute_single_cpp_run", -5 
            
            # self.result_display.append(f"  DEBUG Команда: {' '.join(command)}"); QApplication.processEvents()
            exe_dir = os.path.dirname(os.path.abspath(current_exe_path))
            if not exe_dir: exe_dir=os.getcwd() 
            self.current_plot_data_file = os.path.join(exe_dir, "plot_data.txt") 
            
            cf = subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0
            cp = subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace',cwd=exe_dir,timeout=600,creationflags=cf)
            return cp.stdout, cp.stderr, cp.returncode
        except subprocess.TimeoutExpired: self.result_display.append("  Ошибка: Расчет C++ таймаут."); return None, "Timeout", -2 
        except FileNotFoundError: self.result_display.append(f"  КРИТ.ОШИБКА: Не удалось запустить C++. Файл не найден: {current_exe_path if 'current_exe_path' in locals() else self.cpp_executable_path}"); return None, f"C++ exe {current_exe_path if 'current_exe_path' in locals() else self.cpp_executable_path} not found by OS", -4
        except Exception as e: self.result_display.append(f"  Общая ошибка вызова C++: {type(e).__name__} - {e}"); return None, str(e), -3

    def _display_aggregated_stats(self): # Восстановленный метод
        if not self.all_runs_data_for_stats:
            self.result_display.append("\nНет данных для агрегированной статистики.")
            return

        self.result_display.append("\n--- Агрегированная статистика ---")
        iterations_list, time_ms_list, z_diffs_list = [], [], []
        total_exit_main, total_exit_test = 0, 0
        num_valid_runs_for_stats = 0 # Используем отдельный счетчик для статистики

        for run_data in self.all_runs_data_for_stats:
            if isinstance(run_data, dict) and run_data.get("iterations") is not None:
                num_valid_runs_for_stats += 1
                iterations_list.append(run_data.get("iterations", 0))
                time_ms_list.append(run_data.get("time_ms", 0.0))
                total_exit_main += run_data.get("exit_main_count", 0)
                total_exit_test += run_data.get("exit_test_count", 0)
                
                found_z = run_data.get("found_minimum_z")
                known_z = run_data.get("known_optimum_z")

                if found_z is not None and not (isinstance(found_z, str) and found_z.lower() == "null") and \
                   known_z is not None and not (isinstance(known_z, str) and known_z.lower() == "null"):
                    try:
                        z_diffs_list.append(abs(float(found_z) - float(known_z)))
                    except ValueError:
                        self.result_display.append(f"  Предупреждение: не удалось вычислить разницу Z для запуска с Z_found={found_z}, Z_known={known_z}")
        
        if num_valid_runs_for_stats > 0:
            self.result_display.append(f"Всего корректно обработанных запусков для статистики: {num_valid_runs_for_stats}")
            if iterations_list:
                self.result_display.append(
                    f"Итерации: Ср={np.mean(iterations_list):.2f}, Мед={np.median(iterations_list):.0f}, "
                    f"Мин={np.min(iterations_list)}, Макс={np.max(iterations_list)}, СтОткл={np.std(iterations_list):.2f}"
                )
            if time_ms_list:
                self.result_display.append(
                    f"Время(мс): Ср={np.mean(time_ms_list):.2f}, Мед={np.median(time_ms_list):.2f}, "
                    f"Мин={np.min(time_ms_list):.2f}, Макс={np.max(time_ms_list):.2f}, Сумм={np.sum(time_ms_list):.2f}"
                )
            if z_diffs_list:
                 self.result_display.append(
                    f"|Z_found-Z*|: Ср={np.mean(z_diffs_list):.4f}, Мед={np.median(z_diffs_list):.4f}, "
                    f"Мин={np.min(z_diffs_list):.4f}, Макс={np.max(z_diffs_list):.4f}, СтОткл={np.std(z_diffs_list):.4f}"
                )
            else:
                self.result_display.append("Нет данных для статистики по разнице Z.")

            self.result_display.append(f"Всего выходов по осн. условию: {total_exit_main} ({total_exit_main/num_valid_runs_for_stats*100:.1f}%)")
            self.result_display.append(f"Всего выходов по доп. условию: {total_exit_test} ({total_exit_test/num_valid_runs_for_stats*100:.1f}%)")
            
            exits_by_limit = num_valid_runs_for_stats - (total_exit_main + total_exit_test)
            if exits_by_limit > 0 : # Показывать только если есть такие выходы
                 self.result_display.append(f"Выходов по лимиту итераций: {exits_by_limit} ({exits_by_limit/num_valid_runs_for_stats*100:.1f}%)")
        else:
            self.result_display.append("Не было успешных запусков с данными для статистики.")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ui = OptimizerUI()
    if not ui.cpp_executable_path or not os.path.isfile(ui.cpp_executable_path):
        QMessageBox.critical(ui, "Критическая ошибка",
                             f"Исполняемый файл C++ 'OptimizationProject.exe' не найден при старте.\n"
                             f"Пожалуйста, проверьте его наличие и правильность пути в переменной "
                             f"self.cpp_executable_path в скрипте Python.\n"
                             f"Попытка найти по пути: {ui.cpp_executable_path if ui.cpp_executable_path else 'не определен'}")
        # sys.exit(1) 
    ui.show()
    sys.exit(app.exec_())