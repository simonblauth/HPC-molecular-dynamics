# Defining the executable target
add_executable(08 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(08 PUBLIC my_md_lib)

add_input_file(08 cluster_923.xyz)
add_input_file(08 cluster_3871.xyz)
