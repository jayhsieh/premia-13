project(check_help)
add_executable(check_help check_help.cpp)
target_link_libraries(check_help ${LIBS} ${GEN_LIBS} ${PYTHON_LIBRARIES})

add_custom_command(TARGET check_help POST_BUILD       
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        $<TARGET_FILE:check_help>    
        "${OUTPUT_DIR}/$<TARGET_FILE_NAME:check_help>")
