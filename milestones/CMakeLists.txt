# Function that helps copy input files for milestones
function(add_input_file milestone input_file)
  add_custom_command(
    TARGET ${milestone} POST_BUILD
    DEPENDS ${input_file}
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/${input_file}
    $<TARGET_FILE_DIR:${milestone}>
  )
endfunction()

add_subdirectory(01)
add_subdirectory(04)
add_subdirectory(05)
add_subdirectory(06)

# Future milestones
# add_subdirectory(07)
# add_subdirectory(08)
# add_subdirectory(09)

