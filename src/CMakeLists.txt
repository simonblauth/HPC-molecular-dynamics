# List of headers
set(MY_MD_HEADERS
  hello.h
)

# List of implementation files
set(MY_MD_CPP
  hello.cpp
)

# Create a static library
add_library(my_md_lib STATIC ${MY_MD_CPP} ${MY_MD_HEADERS})

# Include the current directory as INTERFACE
# This will propagate the include to any target linking against the lib
target_include_directories(my_md_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Link against Eigen (automatically sets up the includes)
target_link_libraries(my_md_lib PUBLIC Eigen3::Eigen)

# Set up MPI includes and library linking
# This also propagates to further targets
if (MPI_FOUND)
  target_compile_options(my_md_lib PUBLIC -Wall -Wextra -Wpedantic)
  target_compile_definitions(my_md_lib INTERFACE USE_MPI)
  target_include_directories(my_md_lib PUBLIC ${MPI_CXX_INCLUDE_DIRS} SYSTEM)
  target_link_libraries(my_md_lib PUBLIC MPI::MPI_CXX)
endif()
