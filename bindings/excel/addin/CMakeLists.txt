set(xlpremia_sources 
	Excel_2.rc
	Excel_2.cpp
	OptionPricer.cpp
	PremiaAddin.cpp
)

ADD_MSVC_PRECOMPILED_HEADER("stdafx.h" "stdafx.cpp" xlpremia_sources)

add_library(xlpremia SHARED ${xlpremia_sources})

set_target_properties(xlpremia PROPERTIES FOLDER "Excel")

target_link_libraries(xlpremia pnl premia)

add_custom_command(TARGET xlpremia POST_BUILD       
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        $<TARGET_FILE:xlpremia> "${OUTPUT_DIR}/")
