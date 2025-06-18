import sys
import subprocess
import json
import os
import random # Для случайного выбора
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
                             QPushButton, QComboBox, QTextEdit, QMainWindow, QGridLayout, QGroupBox,
                             QSizePolicy, QSpinBox, QMessageBox) 
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtCore import Qt

# Используем Matplotlib для построения графика
import matplotlib.pyplot as plt
import numpy as np 


# Предполагаемые размеры семейств задач (для случайного выбора)
FAMILY_SIZES = {
    "Hill": 1000,    # Пример, замените на реальный размер вашего THillProblemFamily
    "Shekel": 1000,  # Пример, замените на реальный размер вашего TShekelProblemFamily
    "Grishagin": 100 
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
        self.right_panel_widget.setFixedWidth(250) 

        self.main_layout.addWidget(self.left_panel_widget, 2) 
        self.main_layout.addWidget(self.right_panel_widget, 1) 

        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.cpp_executable_path = os.path.normpath(os.path.join(script_dir, "..", "build", "Release", "OptimizationProject.exe"))
        if not os.path.isfile(self.cpp_executable_path):
            debug_path = os.path.normpath(os.path.join(script_dir, "..", "build", "Debug", "OptimizationProject.exe"))
            if os.path.isfile(debug_path): self.cpp_executable_path = debug_path
            else:
                build_exe_path = os.path.normpath(os.path.join(script_dir, "Release", "OptimizationProject.exe")) 
                if os.path.isfile(build_exe_path): self.cpp_executable_path = build_exe_path
                else:
                    build_exe_path_debug = os.path.normpath(os.path.join(script_dir, "Debug", "OptimizationProject.exe"))
                    if os.path.isfile(build_exe_path_debug): self.cpp_executable_path = build_exe_path_debug
                    elif '__file__' in locals(): # Проверка, чтобы QMessageBox не вызывался при импорте как модуля
                        QMessageBox.critical(self, "Критическая ошибка",
                                             f"Исполняемый файл C++ 'OptimizationProject.exe' не найден.\n"
                                             f"Проверьте его наличие и правильность пути в скрипте Python.\n"
                                             f"Ожидаемый путь (пример): {self.cpp_executable_path}")

        # --- ИНИЦИАЛИЗИРУЕМ АТРИБУТЫ ДАННЫХ ДО ВЫЗОВА _on_parameter_changed ---
        self.current_plot_data_file = None 
        self.last_successful_run_data = None 
        self.all_runs_data_for_stats = [] 
        # ------------------------------------------------------------------

        self._create_input_widgets_definitions()
        self._create_plot_buttons() # Кнопки создаются до того, как к ним обратится _on_parameter_changed
        self._create_input_widgets_layout() # Размещаем виджеты ввода
        self._create_output_widgets()
        
        self._on_parameter_changed() # Теперь этот вызов безопасен

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
        self.max_iter_entry.setValidator(QIntValidator(10, 200000))
        
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

    def _create_plot_buttons(self): 
        plot_group = QGroupBox("Графики")
        plot_layout = QVBoxLayout()
        self.plot_convergence_button = QPushButton("График сходимости")
        self.plot_convergence_button.clicked.connect(self.show_matplotlib_plot_from_file)
        plot_layout.addWidget(self.plot_convergence_button)
        self.plot_level_lines_button = QPushButton("Линии уровня и точки (2D)")
        self.plot_level_lines_button.clicked.connect(self.show_level_lines_plot)
        self.plot_level_lines_button.setEnabled(False) 
        plot_layout.addWidget(self.plot_level_lines_button)
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
        self.problem_idx_label.setEnabled(num_runs == 1)
        self.problem_idx_entry.setEnabled(num_runs == 1)
        
        # Теперь self.last_successful_run_data существует, даже если None
        self.plot_level_lines_button.setEnabled(not is_1d_type and num_runs == 1 and self.last_successful_run_data is not None)

        max_idx_val = 1 
        if is_1d_type:
            family_key = self.sub_choice_combo.currentText().split(':')[1].strip()
            max_idx_val = FAMILY_SIZES.get(family_key, 1000) 
            self.problem_idx_label.setText(f"Индекс задачи (1D, 1-{max_idx_val}):")
        else:
            max_idx_val = FAMILY_SIZES.get("Grishagin", 100)
            self.problem_idx_label.setText(f"Индекс задачи (Гришагина, 1-{max_idx_val}):")
        current_idx_validator = self.problem_idx_entry.validator()
        if isinstance(current_idx_validator, QIntValidator):
            current_idx_validator.setTop(max_idx_val if max_idx_val > 0 else 1)

    def _create_output_widgets(self):
        output_group = QGroupBox("Результаты вычислений")
        output_layout = QVBoxLayout()
        self.result_display = QTextEdit(); self.result_display.setReadOnly(True)
        self.result_display.setPlaceholderText("Результаты расчета появятся здесь...")
        self.result_display.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        output_layout.addWidget(self.result_display); output_group.setLayout(output_layout)
        self.left_panel_layout.addWidget(output_group)
        self.left_panel_layout.setStretchFactor(output_group, 1)

    def show_matplotlib_plot_from_file(self):
        num_runs = self.num_runs_spinbox.value()
        if num_runs > 1 and self.all_runs_data_for_stats:
            iterations_list = [rd.get("iterations",0) for rd in self.all_runs_data_for_stats if isinstance(rd,dict) and rd.get("iterations") is not None]
            if iterations_list:
                plt.figure("Распределение итераций"); plt.hist(iterations_list, bins='auto', color='skyblue', rwidth=0.85)
                plt.xlabel("Кол-во итераций"); plt.ylabel("Частота"); plt.title(f"Распределение итераций ({len(iterations_list)} зап.)")
                plt.grid(axis='y', alpha=0.75); plt.tight_layout(); plt.show()
                self.result_display.append("Гистограмма итераций построена.")
            else: self.result_display.append("Нет данных для гистограммы итераций.")
        elif self.current_plot_data_file and os.path.exists(self.current_plot_data_file):
            try:
                px,py=[],[]
                with open(self.current_plot_data_file,"r") as f:
                    for i,l in enumerate(f,1):
                        pts=l.strip().split()
                        if len(pts)==2:
                            try:px.append(int(pts[0]));py.append(float(pts[1]))
                            except ValueError: self.result_display.append(f"Предупр.: плохая строка {i} в plot_data.txt: {l.strip()}")
                if not px or not py: self.result_display.append(f"Нет данных в {self.current_plot_data_file}");return
                plt.figure(f"График ({os.path.basename(self.current_plot_data_file)})");plt.plot(px,py,marker='o');plt.xlabel("Точка");plt.ylabel("Итерации");plt.grid(True);plt.show()
                self.result_display.append(f"График из {self.current_plot_data_file} построен.")
            except Exception as e:self.result_display.append(f"Ошибка графика из файла: {e}")
        else:self.result_display.append(f"Файл данных ({self.current_plot_data_file}) не найден. Расчет?")

    def show_level_lines_plot(self):
        self.result_display.append("\nЗапрос данных для линий уровня..."); QApplication.processEvents()
        if not self.last_successful_run_data: self.result_display.append("Ошибка: Сначала выполните успешный расчет."); return
        
        problem_idx_cpp = self.last_successful_run_data.get("problem_index_cpp")
        problem_family = self.last_successful_run_data.get("problem_family")
        if problem_family != "Grishagin" or problem_idx_cpp is None: self.result_display.append("Ошибка: Линии уровня только для Гришагина."); return

        grid_res = 50
        cmd_ll = [self.cpp_executable_path, "--get-level-lines", str(problem_idx_cpp), str(grid_res)]
        self.result_display.append(f"Команда (линии уровня): {' '.join(cmd_ll)}"); QApplication.processEvents()
        try:
            exe_dir = os.path.dirname(os.path.abspath(self.cpp_executable_path));
            if not exe_dir: exe_dir = os.getcwd()
            cf = subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0
            ll_proc = subprocess.run(cmd_ll, capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=exe_dir, timeout=60, creationflags=cf)
            if ll_proc.returncode == 0:
                try:
                    ll_data = json.loads(ll_proc.stdout)
                    if ll_data.get("type") != "level_lines_data": self.result_display.append("Ошибка: C++ вернул неверный JSON для линий уровня."); return
                    xg,yg,zg = np.array(ll_data.get("x_coords_grid",[])), np.array(ll_data.get("y_coords_grid",[])), np.array(ll_data.get("z_values_grid",[]))
                    if not (xg.size>0 and yg.size>0 and zg.size>0 and zg.shape==(len(yg),len(xg))):
                        self.result_display.append(f"Ошибка: неверные размеры массивов для линий уровня. XG:{xg.shape},YG:{yg.shape},ZG:{zg.shape}"); return
                    
                    trials, found, known = self.last_successful_run_data.get("trial_history_y_coords",[]), self.last_successful_run_data.get("found_minimum_y",[]), self.last_successful_run_data.get("known_optimum_y",[])
                    plt.figure(f"Линии уровня: {problem_family} {problem_idx_cpp+1}", figsize=(8.5,7)); X,Y=np.meshgrid(xg,yg)
                    cp = plt.contourf(X,Y,zg,levels=30,cmap='viridis'); plt.colorbar(cp,label='Z'); plt.contour(X,Y,zg,levels=cp.levels,colors='k',linewidths=0.5)
                    if trials: tx=[p[0] for p in trials if len(p)==2]; ty=[p[1] for p in trials if len(p)==2]; plt.scatter(tx,ty,s=10,c='k',alpha=0.6,label=f'Исслед.т.({len(tx)})')
                    if len(found)==2: plt.scatter(found[0],found[1],s=120,c='r',marker='*',edgecolor='k',label='Найден min',zorder=5)
                    if len(known)==2: plt.scatter(known[0],known[1],s=120,c='lime',marker='P',edgecolor='k',label='Известен opt',zorder=5)
                    plt.xlabel("Y1");plt.ylabel("Y2");plt.title(f"Линии уровня: {problem_family} Задача {problem_idx_cpp+1}");plt.legend();plt.axis('scaled');plt.xlim(xg.min(),xg.max());plt.ylim(yg.min(),yg.max());plt.tight_layout();plt.show()
                    self.result_display.append("График линий уровня построен.")
                except json.JSONDecodeError as e: self.result_display.append(f"Ошибка JSON (линии уровня): {e}\n{ll_proc.stdout}")
                except Exception as ep: self.result_display.append(f"Ошибка построения линий уровня: {type(ep).__name__} - {ep}")
            else: self.result_display.append(f"Ошибка C++ (линии уровня, код {ll_proc.returncode}):\n{ll_proc.stderr}\n{ll_proc.stdout}")
        except Exception as e: self.result_display.append(f"Общая ошибка (линии уровня): {type(e).__name__} - {e}")

    def run_calculation_handler(self):
        self.result_display.clear(); self.all_runs_data_for_stats.clear(); self.last_successful_run_data=None; self.plot_level_lines_button.setEnabled(False)
        num_runs = self.num_runs_spinbox.value()
        self.result_display.append(f"Запускается серия из {num_runs} расчетов..."); QApplication.processEvents()
        for i in range(num_runs):
            self.result_display.append(f"\n--- Расчет {i+1} из {num_runs} ---"); QApplication.processEvents()
            task_type, sub_choice, family_key = self.task_type_combo.currentText().split(":")[0], "", ""
            if task_type == "1": sub_choice, family_key = self.sub_choice_combo.currentText().split(":")[0], self.sub_choice_combo.currentText().split(':')[1].strip()
            elif task_type == "2": family_key = "Grishagin"
            
            problem_idx_ui = int(self.problem_idx_entry.text()) if num_runs == 1 else random.randint(1, FAMILY_SIZES.get(family_key,1) if FAMILY_SIZES.get(family_key,1)>0 else 1)
            if num_runs > 1: self.result_display.append(f"  Случ.выбран тип: {family_key}, индекс(1-based): {problem_idx_ui}")
            
            stdout, stderr, ret_code = self._execute_single_cpp_run(task_type, sub_choice, problem_idx_ui, self.epsilon_entry.text(), self.r_param_entry.text(), self.peano_order_entry.text(), self.max_iter_entry.text())
            if stdout is None and stderr is None: self.result_display.append("  Крит.ошибка подготовки C++ вызова."); continue
            if ret_code == 0 and stdout:
                try:
                    res = json.loads(stdout); res['problem_index_ui_used'] = problem_idx_ui
                    self.all_runs_data_for_stats.append(res); self.last_successful_run_data = res
                    idx_disp = res.get('problem_index_cpp',-1)+1
                    self.result_display.append("  Расчет успешно завершен.\n"+f"  Семейство: {res.get('problem_family','N/A')}\n  Индекс (1-based UI): {problem_idx_ui} (C++ 0-based: {res.get('problem_index_cpp','N/A')})\n  Найденный Z: {res.get('found_minimum_z','N/A')}\n  Найденный Y: {res.get('found_minimum_y','N/A')}\n  Известный Z*: {res.get('known_optimum_z','N/A')}\n  Известный Y*: {res.get('known_optimum_y','N/A')}\n  Итераций: {res.get('iterations','N/A')}\n  Время: {res.get('time_ms','N/A')} мс\n  Выход осн.усл.: {res.get('exit_main_count','N/A')}\n  Выход доп.усл.: {res.get('exit_test_count','N/A')}\n")
                    if os.path.exists(self.current_plot_data_file): self.result_display.append(f"  Файл {self.current_plot_data_file} обновлен.")
                    else: self.result_display.append(f"  ВНИМАНИЕ: Файл {self.current_plot_data_file} не создан C++.")
                except json.JSONDecodeError as e: self.result_display.append(f"  Ошибка парсинга JSON: {e}\nStdout:{stdout}\nStderr:{stderr}")
                except Exception as e: self.result_display.append(f"  Ошибка обработки вывода C++: {e}\nStdout:{stdout}\nStderr:{stderr}")
            else: self.result_display.append(f"  Ошибка C++ (код {ret_code}):\nStderr:{stderr}\nStdout:{stdout}")
            QApplication.processEvents()
            if num_runs > 1 and i < num_runs-1: self.result_display.append("-" * 30)
        if num_runs > 1: self._display_aggregated_stats()
        if self.last_successful_run_data and self.last_successful_run_data.get("problem_family")=="Grishagin": self.plot_level_lines_button.setEnabled(True)
        else: self.plot_level_lines_button.setEnabled(False)
        self.result_display.append("\nВсе расчеты завершены.")

    def _execute_single_cpp_run(self, task_type_str, sub_choice_str, problem_idx_ui, epsilon_str_ui, r_param_str_ui, peano_order_m_str_ui, max_iter_str_ui):
        try:
            current_exe_path = self.cpp_executable_path
            if not os.path.isfile(current_exe_path):
                script_dir = os.path.dirname(os.path.abspath(__file__)); basename_exe = os.path.basename(self.cpp_executable_path) if self.cpp_executable_path else "OptimizationProject.exe"
                potential_paths = [os.path.join(script_dir, basename_exe), os.path.normpath(os.path.join(script_dir,"..","build","Release",basename_exe)), os.path.normpath(os.path.join(script_dir,"..","build","Debug",basename_exe)), os.path.normpath(os.path.join(script_dir,"Release",basename_exe)), os.path.normpath(os.path.join(script_dir,"Debug",basename_exe)), os.path.normpath(os.path.join(script_dir,"build","bin","Release",basename_exe)), os.path.normpath(os.path.join(script_dir,"build","bin","Debug",basename_exe)), os.path.normpath(os.path.join(script_dir,"build","bin",basename_exe))]
                found_exe = False
                for p_path in potential_paths:
                    if os.path.isfile(p_path): current_exe_path = p_path; found_exe = True; self.result_display.append(f"  Найден C++ exe: {current_exe_path}"); break
                if not found_exe: self.result_display.append(f"КРИТ.ОШИБКА: exe C++ НЕ НАЙДЕН. Путь: {self.cpp_executable_path}"); return None, "C++ exe not found", -1

            problem_idx_cpp = problem_idx_ui - 1 
            epsilon_val_str, r_param_val_str = epsilon_str_ui.replace(',','.'), r_param_str_ui.replace(',','.')
            command = [current_exe_path, task_type_str]
            if task_type_str=="1": command.extend([sub_choice_str,epsilon_val_str,r_param_val_str,str(problem_idx_cpp),max_iter_str_ui])
            elif task_type_str=="2": command.extend([str(problem_idx_cpp),epsilon_val_str,r_param_val_str,peano_order_m_str_ui,max_iter_str_ui])
            else: return None, "Invalid task_type_str", -5 
            
            self.result_display.append(f"  Команда для C++: {' '.join(command)}"); QApplication.processEvents()
            exe_dir = os.path.dirname(os.path.abspath(current_exe_path));
            if not exe_dir: exe_dir=os.getcwd() 
            self.current_plot_data_file = os.path.join(exe_dir, "plot_data.txt") 
            
            cf = subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0
            cp = subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace',cwd=exe_dir,timeout=600,creationflags=cf) # Увеличен таймаут
            return cp.stdout, cp.stderr, cp.returncode
        except subprocess.TimeoutExpired: self.result_display.append("  Ошибка: Расчет C++ таймаут."); return None, "Timeout", -2 
        except FileNotFoundError: self.result_display.append(f"  КРИТ.ОШИБКА: Не удалось запустить C++. Файл не найден: {self.cpp_executable_path}"); return None, f"C++ exe {self.cpp_executable_path} not found by OS", -4
        except Exception as e: self.result_display.append(f"  Общая ошибка вызова C++: {type(e).__name__} - {e}"); return None, str(e), -3
            
    def _display_aggregated_stats(self):
        if not self.all_runs_data_for_stats: self.result_display.append("\nНет данных для агрег.статистики."); return
        self.result_display.append("\n--- Агрегированная статистика ---")
        iters,times,z_diffs = [],[],[]; main_exits,test_exits,valid_runs = 0,0,0
        for d in self.all_runs_data_for_stats:
            if isinstance(d,dict) and d.get("iterations") is not None:
                valid_runs+=1; iters.append(d.get("iterations",0)); times.append(d.get("time_ms",0))
                main_exits+=d.get("exit_main_count",0); test_exits+=d.get("exit_test_count",0)
                fz,kz=d.get("found_minimum_z"),d.get("known_optimum_z")
                if fz is not None and kz is not None and not(isinstance(fz,str)and fz=="null") and not(isinstance(kz,str)and kz=="null"):
                    try: z_diffs.append(abs(float(fz)-float(kz)))
                    except ValueError: pass
        if valid_runs>0:
            self.result_display.append(f"Всего корректно обработанных запусков: {valid_runs}")
            if iters: self.result_display.append(f"Итерации: Avg={np.mean(iters):.2f}, Std={np.std(iters):.2f}, Min={np.min(iters)}, Max={np.max(iters)}")
            if times: self.result_display.append(f"Время(мс): Avg={np.mean(times):.2f}, Std={np.std(times):.2f}, Min={np.min(times):.2f}, Max={np.max(times):.2f}")
            if z_diffs: self.result_display.append(f"|Z_found-Z*|: Avg={np.mean(z_diffs):.4f}, Std={np.std(z_diffs):.4f}, Min={np.min(z_diffs):.4f}, Max={np.max(z_diffs):.4f}")
            self.result_display.append(f"Всего выходов осн.усл.: {main_exits} ({main_exits/valid_runs*100:.1f}%)")
            self.result_display.append(f"Всего выходов доп.усл.: {test_exits} ({test_exits/valid_runs*100:.1f}%)")
            if(main_exits+test_exits)<valid_runs: self.result_display.append(f"Выходов по лимиту итераций: {valid_runs-(main_exits+test_exits)}")
        else: self.result_display.append("Не было успешных запусков для статистики.")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ui = OptimizerUI()
    ui.show()
    sys.exit(app.exec_())