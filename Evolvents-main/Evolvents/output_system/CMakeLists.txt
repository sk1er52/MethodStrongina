cmake_minimum_required(VERSION 3.10)
project(EvolventsLib LANGUAGES CXX)

# --- Исходные файлы библиотеки ---
set(EVOLVENTS_SOURCES
    evolvent.cpp
    extended.cpp
    exception.cpp
    output_system.cpp
    messages.cpp
)

add_library(evolvents_lib STATIC ${EVOLVENTS_SOURCES})
target_include_directories(evolvents_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
if(WIN32)
    message(STATUS "Configuring EvolventsLib dependencies for Windows (MPFR/MPIR)")

    # Определяем архитектуру (x64 или x86)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(ARCH_DIR "x64")
        message(STATUS "Detected 64-bit build")
    else()
        set(ARCH_DIR "x86")
        message(STATUS "Detected 32-bit build")
    endif()

    # --- Настройка MPIR (GMP) ---
    set(MPIR_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mpir/include)
    set(MPIR_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/mpir/lib/win/vc14/${ARCH_DIR}/mpir.lib) # Путь к .lib

    if(EXISTS ${MPIR_INCLUDE_DIR} AND EXISTS ${MPIR_LIB_PATH})
        add_library(mpir STATIC IMPORTED GLOBAL) # GLOBAL чтобы было видно в родительском проекте если надо
        set_target_properties(mpir PROPERTIES
            IMPORTED_LOCATION "${MPIR_LIB_PATH}"
            INTERFACE_INCLUDE_DIRECTORIES "${MPIR_INCLUDE_DIR}"
        )
        message(STATUS "Found MPIR: ${MPIR_LIB_PATH}")
        # Линкуем библиотеку evolvents_lib с mpir
        target_link_libraries(evolvents_lib PRIVATE mpir)
    else()
        message(FATAL_ERROR "MPIR include or library not found at expected paths:\nInclude: ${MPIR_INCLUDE_DIR}\nLibrary: ${MPIR_LIB_PATH}")
    endif()

    # --- Настройка MPFR ---
    set(MPFR_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mpfr/include)
    set(MPFR_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/mpfr/lib/win/vc14/${ARCH_DIR}/mpfr.lib) # Путь к .lib

    if(EXISTS ${MPFR_INCLUDE_DIR} AND EXISTS ${MPFR_LIB_PATH})
        add_library(mpfr STATIC IMPORTED GLOBAL)
        set_target_properties(mpfr PROPERTIES
            IMPORTED_LOCATION "${MPFR_LIB_PATH}"
            INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}"
        )
        message(STATUS "Found MPFR: ${MPFR_LIB_PATH}")
        # Линкуем библиотеку evolvents_lib с mpfr
        target_link_libraries(evolvents_lib PRIVATE mpfr)
    else()
        message(FATAL_ERROR "MPFR include or library not found at expected paths:\nInclude: ${MPFR_INCLUDE_DIR}\nLibrary: ${MPFR_LIB_PATH}")
    endif()

    # Добавляем include директории для MPIR и MPFR, чтобы evolvents_lib их видела при компиляции
    # Хотя INTERFACE_INCLUDE_DIRECTORIES должно быть достаточно, но для надежности
    target_include_directories(evolvents_lib PRIVATE
         ${MPIR_INCLUDE_DIR}
         ${MPFR_INCLUDE_DIR}
    )

else()
    message(WARNING "Pre-built MPFR/MPIR are provided only for Windows. Build might fail on other platforms unless system libraries are found.")
    message(FATAL_ERROR "Automatic configuration for non-Windows MPFR/MPIR is not implemented in this example. Please install them and adjust CMakeLists.txt or provide pre-built libraries.")
endif()