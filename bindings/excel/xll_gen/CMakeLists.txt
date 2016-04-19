set(xllgen_sources 
	ExcelFuncGenerator.cpp
)

ADD_MSVC_PRECOMPILED_HEADER("stdafx.h" "stdafx.cpp" xllgen_sources)

add_executable(xllgen ${xllgen_sources})
set_target_properties(xllgen PROPERTIES FOLDER "Excel")


target_link_libraries(xllgen pnl premia)

add_custom_command(TARGET xllgen POST_BUILD       
    COMMAND $<TARGET_FILE:xllgen> ${XLL_GENERATED_SOURCE})
