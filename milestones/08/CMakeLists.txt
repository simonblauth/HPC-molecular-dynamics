# Defining the executable target
add_executable(08 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(08 PUBLIC my_md_lib)
