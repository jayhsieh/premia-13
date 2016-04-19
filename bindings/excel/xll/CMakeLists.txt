set(xll_sources 
	ExcelFunc.cpp
	OptionPrice.cpp
	xlwoptionprice.cpp
)

ADD_MSVC_PRECOMPILED_HEADER("stdafx.h" "stdafx.cpp" xll_sources)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/import_xlpremia.h.cmake 
    ${CMAKE_CURRENT_BINARY_DIR}/import_xlpremia.h)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(xll SHARED ${xll_sources})
set_target_properties(xll PROPERTIES FOLDER "Excel")

set_target_properties(xll PROPERTIES SUFFIX ".xll")
set_target_properties(xll PROPERTIES OUTPUT_NAME "premia")

target_link_libraries(xll pnl premia xlw xlpremia)
add_dependencies(xll xllgen)

add_custom_command(TARGET xll POST_BUILD       
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        $<TARGET_FILE:xll> "${OUTPUT_DIR}")
