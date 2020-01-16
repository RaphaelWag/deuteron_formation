cmake_minimum_required( VERSION 3.10 )
project( mySimulation )

set( CMAKE_CXX_STANDARD 17 )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info")

#add_executable( mySimulation ALICE_data_filter.cpp myParticle.cpp)
#add_executable( mySimulation ALICE_np_weighting.cpp mySimulation.cpp myParticle.cpp histogram.cpp )
add_executable( mySimulation ALICE_LT_weighting.cpp mySimulation.cpp myParticle.cpp histogram.cpp )
#add_executable( mySimulation ALICE_coalescence.cpp mySimulation.cpp myParticle.cpp histogram.cpp )
#add_executable( mySimulation ALICE_coalescence_multiplicity.cpp mySimulation.cpp myParticle.cpp histogram.cpp )
#add_executable( mySimulation ALICE_coalescence_plot.cpp mySimulation.cpp myParticle.cpp histogram.cpp )
#add_executable( mySimulation p_x_correlation.cpp mySimulation.cpp myParticle.cpp histogram2d.cpp )
#add_executable( mySimulation ALICE_y_kinematics.cpp mySimulation.cpp myParticle.cpp histogram2d.cpp )