project(karrigell)

add_executable(ksgen ksgen.cpp)
target_link_libraries(ksgen ${LIBS} ${GEN_LIBS} ${PYTHON_LIBRARIES})

if (KARRIGELL_INSTALLATION_PATH)
   add_custom_target(install_karrigell        
       COMMAND $<TARGET_FILE:ksgen> --output-dir="${KARRIGELL_INSTALLATION_PATH}" --pypremia-dir="${OUTPUT_DIR}/python" --data-dir="${PREMIA_DATA_DIR}" --template-dir="${CMAKE_CURRENT_SOURCE_DIR}/templates"
       DEPENDS ksgen install_python)

   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.cmake
       ${KARRIGELL_INSTALLATION_PATH}/data/www/conf.py)

endif (KARRIGELL_INSTALLATION_PATH)