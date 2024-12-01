cmake_minimum_required(VERSION 3.10)
project(OptimizationProject)

# Найдите математическую библиотеку
find_library(M_LIB m)

# Добавьте флаг компиляции для C++11 (или более поздней версии)
set(CMAKE_CXX_STANDARD 17)

# Укажите исходные файлы вашего проекта
set(SOURCE_FILES
    MethodStrongina/Source.cpp # Замените на имя вашего исходного файла
    Hill/HillProblem.cpp       # Предполагается, что вы создадите этот файл
    Shekel/ShekelProblem.cpp    # Предполагается, что вы создадите этот файл
    src/IOptProblem.cpp
    src/IOptProblemFamily.cpp
    src/IGeneralOptProblem.cpp
    src/IGeneralOptProblemFamily.cpp
    src/IConstrainedOptProblem.cpp
    src/IConstrainedOptProblemFamily.cpp
)

# Создайте исполняемый файл
add_executable(OptimizationProject ${SOURCE_FILES})

# Подключите математическую библиотеку, если она найдена
if(M_LIB)
    target_link_libraries(OptimizationProject ${M_LIB})
endif()


# Добавьте каталоги include (если необходимо) - если ваши .hpp файлы находятся
# в отдельном каталоге, раскомментируйте следующую строку и укажите путь
target_include_directories(OptimizationProject PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Hill" "${CMAKE_CURRENT_SOURCE_DIR}/Shekel" "${CMAKE_CURRENT_SOURCE_DIR}/src")

# (Опционально) создайте библиотеку, если вы хотите использовать функции Хилла и Шекеля
# в других проектах
# add_library(OptimizationFunctions STATIC HillProblem.cpp ShekelProblem.cpp)
# target_link_libraries(OptimizationProject OptimizationFunctions)