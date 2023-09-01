# Defining the executable target
add_executable(07 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(07 PUBLIC my_md_lib)
