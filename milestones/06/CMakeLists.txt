# Defining the executable target
add_executable(06 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(06 PUBLIC my_md_lib)
