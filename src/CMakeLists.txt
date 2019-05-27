project(multitopo)
cmake_minimum_required(VERSION 3.0)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 ")

SET(EIGEN_INCLUDE_DIRS "/opt/local/include/eigen3")
SET(EIGEN_INCLUDE_DIRS2 "/opt/local/include/")

SET(GUROBI_INCLUDE_DIRS "/Library/gurobi811/mac64/include/")
SET(GUROBI_LIB_DIRS "/Library/gurobi811/mac64/lib/")
SET(GUROBI_LIB gurobi_c++ gurobi81)


SET(TETGEN_INCLUDE_DIRS "../external/tetgen/")
SET(TETGEN_LIB_DIRS "../external/tetgen/")
SET(TETGEN_LIB tet)

include_directories(${EIGEN_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS2} ${GUROBI_INCLUDE_DIRS} ${TETGEN_INCLUDE_DIRS})


aux_source_directory(. MAIN)
aux_source_directory(./MainProcess MP)
aux_source_directory(./MMCore MM)
aux_source_directory(./Utility UTL)
aux_source_directory(./DMCore/Arrangement ARR)
aux_source_directory(./DMCore/GroupingManager GM)
aux_source_directory(./DMCore/Utility DMU)

aux_source_directory(./DMCore/LuCtr2Suf LU)
aux_source_directory(./DMCore/LuCtr2Suf/ContourHandler LUC)
aux_source_directory(./DMCore/LuCtr2Suf/JuFair LUJ)
aux_source_directory(./DMCore/LuCtr2Suf/Math LUM)
aux_source_directory(./DMCore/LuCtr2Suf/SpacePartitioner LUS)
aux_source_directory(./DMCore/LuCtr2Suf/Util LUU)

include_directories(. ./MainProcess ./MMCore ./Utility ./DMCore ./DMCore/Arrangement ./DMCore/GroupingManager ./DMCore/Utility)
include_directories(./DMCore/LuCtr2Suf ./DMCore/LuCtr2Suf/ContourHandler ./DMCore/LuCtr2Suf/JuFair ./DMCore/LuCtr2Suf/Math ./DMCore/LuCtr2Suf/SpacePartitioner ./DMCore/LuCtr2Suf/Util)

LINK_DIRECTORIES(${GUROBI_LIB_DIRS} ${TETGEN_LIB_DIRS})
add_executable(${PROJECT_NAME} ${MAIN} ${MP} ${MM} ${UTL} ${ARR} ${GM} ${DMU} ${LU} ${LUC} ${LUJ} ${LUM} ${LUS} ${LUU})

target_link_libraries(${PROJECT_NAME} ${GUROBI_LIB} ${TETGEN_LIB})
