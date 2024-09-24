cmake_minimum_required(VERSION 3.10)

# A projekt neve
project(FFTProject)

# C szabvány beállítása
set(CMAKE_C_STANDARD 11)

# Forrásfájlok
set(SRC_FILES src/fft.c)
set(TEST_FILES tests/test_fft.c)

# Főprogram létrehozása
add_executable(fft_program ${SRC_FILES})

# Tesztprogram létrehozása
add_executable(test_fft ${SRC_FILES} ${TEST_FILES})

# Linkelés a Check keretrendszerhez
find_package(Check REQUIRED)
target_link_libraries(test_fft Check::Check)