include_directories(
    ${PROJECT_SOURCE_DIR}/Src
    )
add_subdirectory (math)
add_library (libcommon OBJECT
	chk.c
	ftools.c
	method.c
	numfunc.c
	error_msg.c
	test.c
	timeinfo.c
	tools.c
	var.c
	model.c
	option.c
	enums.c
)

install(FILES 
    var.h tools.h optype.h error_msg.h enums.h
    DESTINATION include)
