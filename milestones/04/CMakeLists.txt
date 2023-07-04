# Defining the executable target
add_executable(04 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(04 PUBLIC my_md_lib)

add_input_file(04 lj54.xyz)
