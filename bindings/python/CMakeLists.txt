project(python)

set(OUTPUT_DIR ${OUTPUT_DIR}/python)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${OUTPUT_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/interop.py.cmake 
    ${OUTPUT_DIR}/premia/interop.py)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/common.py
    ${OUTPUT_DIR}/premia/common.py)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/time_all_methods.py
    ${OUTPUT_DIR}/time_all_methods.py)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/get_all_method_results.py
    ${OUTPUT_DIR}/get_all_method_results.py)


add_library(pypremia SHARED pypremia.cpp)
target_link_libraries(pypremia ${LIBS} ${GEN_LIBS} ${PYTHON_LIBRARIES})

    

if (MSVC)
	set_target_properties(pypremia PROPERTIES SUFFIX ".pyd")
	set(pypremia_suffix ".pyd")
elseif (UNIX)
	set_target_properties(pypremia PROPERTIES SUFFIX ".so")
endif()
set_target_properties(pypremia PROPERTIES PREFIX "lib")
add_custom_command(TARGET pypremia POST_BUILD       
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:pypremia>    
        "${OUTPUT_DIR}/premia/$<TARGET_FILE_NAME:pypremia>")

add_executable(pygen pygen.cpp)
target_link_libraries(pygen ${LIBS} ${GEN_LIBS} ${PYTHON_LIBRARIES})

add_custom_target(install_python       
    COMMAND $<TARGET_FILE:pygen> --output-dir="${OUTPUT_DIR}/premia" --dll-dir="${OUTPUT_DIR}/premia" --data-dir="${PREMIA_DATA_DIR}" --template-dir="${CMAKE_CURRENT_SOURCE_DIR}"
    DEPENDS pygen pypremia)

