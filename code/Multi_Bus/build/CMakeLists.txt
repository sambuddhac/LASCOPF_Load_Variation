##############################################KNITRO CMake##################################################################################################
# make a directory build in KTR, go there and run “cmake ..”.
# may need to specify generator in Windows:
# e.g.: for Visual Studio 2013 on x64, run cmake -G "Visual Studio 12 Win64" ..
# need separate build folders for Win32 and Win64
# The makefile will be created there. It uses KNITRODIR environment variable.
# Specify the minimum version for CMake

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)

PROJECT(LASCOPFDV)

# Set the output folder where your program will be created
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-std=c++0x")
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
#set(CMAKE_BINARY_DIR /usr/local/lib)
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
# Project's name
INCLUDE_DIRECTORIES($ENV{GUROBI_HOME}/include)
#INCLUDE_DIRECTORIES($ENV{GLPKDIR}/include)
SET(GUR_LIB libgurobi80.so)
SET(GUR_LIBR libgurobi_c++.a)
SET(GUR_LIBRA libgurobi_g++5.2.a)
#SET(GLPK_LIB libglpk.so)
# The following folder will be included
#include_directories("/usr/local/include")
# set the path to the library folderibrary folder
#link_directories(/usr/local/lib)

# link the libraries to the executable
ADD_EXECUTABLE(
LASCOPFDV
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/mainserialLASCOPFLV.cpp  
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/network.cpp
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/node.cpp
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/generator.cpp
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/gensolverFirst.cpp
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/gensolverIntermediate.cpp
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/gensolverLast.cpp
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/transl.cpp 
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/tlinecont.cpp
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/load.cpp
/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/code/Multi_Bus/loadcont.cpp
)
TARGET_LINK_LIBRARIES(LASCOPFDV /home/samie/opt/gurobi801/linux64/lib/${GUR_LIBR};/home/samie/opt/gurobi801/linux64/lib/${GUR_LIBRA};/home/samie/opt/gurobi801/linux64/lib/${GUR_LIB}   ${CMAKE_DL_LIBS})
#TARGET_LINK_LIBRARIES(HorInvestCoord /home/samie/opt/gurobi801/linux64/lib/${GUR_LIBR}  ${CMAKE_DL_LIBS})
#TARGET_LINK_LIBRARIES(HorInvestCoord /usr/local/lib/${GLPK_LIB}  ${CMAKE_DL_LIBS})
