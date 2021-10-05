cmake_minimum_required (VERSION 3.10)
project(mj_pendulum LANGUAGES CXX)

#set(MUJOCO_LIB_DIR "/home/gaussian/cmu_ri_phd/phd_misc/mujoco200_linux/bin")
#link_directories(${MUJOCO_LIB_DIR})

add_executable(${PROJECT_NAME} $<$<BOOL:$<TARGET_EXISTS:PSOPT_SNOPT_interface>>:$<TARGET_OBJECTS:PSOPT_SNOPT_interface>>
        mj_pendulum.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt adolc Eigen3::Eigen)

if(${WITH_SNOPT_INTERFACE})
    enable_language(Fortran)
    target_link_libraries(${PROJECT_NAME} PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)
endif()
