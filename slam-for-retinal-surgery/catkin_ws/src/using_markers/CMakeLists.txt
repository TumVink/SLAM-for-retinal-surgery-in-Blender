cmake_minimum_required(VERSION 2.8.3)
project(using_markers)


find_package(catkin REQUIRED COMPONENTS
  roscpp
  visualization_msgs
)









catkin_package(

)


include_directories(
# include
  ${catkin_INCLUDE_DIRS}
)





## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
add_executable(basic_shapes src/basic_shapes.cpp)
target_link_libraries(basic_shapes ${catkin_LIBRARIES})
add_executable(points_and_lines src/points_and_lines.cpp)
target_link_libraries(points_and_lines ${catkin_LIBRARIES})
