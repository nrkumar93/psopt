cmake_minimum_required (VERSION 3.10)
project(cartpole LANGUAGES CXX)

set(MANIPULATOR_MDL_SHIP_DIR "/home/gaussian/cmu_ri_phd/phd_research/csbpl_mex/nabil/insat/rbt_matlab_test/simscape_model/models/manipulator_ip_ship")

MACRO(HEADER_DIRECTORIES return_list)
    FILE(GLOB_RECURSE new_list "${MANIPULATOR_MDL_SHIP_DIR}/*.h")
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

HEADER_DIRECTORIES(header_dir_list)

include_directories(
        ${CARTPOLE_INCLUDE_DIR}
        ${header_dir_list}) # Recursive

add_executable(${PROJECT_NAME} $<$<BOOL:$<TARGET_EXISTS:PSOPT_SNOPT_interface>>:$<TARGET_OBJECTS:PSOPT_SNOPT_interface>>
        dynamics.cpp
        cartpole.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt adolc Eigen3::Eigen)

if(${WITH_SNOPT_INTERFACE})
    enable_language(Fortran)
    target_link_libraries(${PROJECT_NAME} PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)
endif()
