macro(oledgf_local_options)
    option(oledgf_WARNINGS_AS_ERRORS "Treat Warnings As Errors" OFF)

    if(PROJECT_IS_TOP_LEVEL)
     include(cmake/StandardProjectSettings.cmake)
    endif()

    add_library(oledgf_warnings INTERFACE)
    add_library(oledgf_options INTERFACE)

    include(cmake/CompilerWarnings.cmake)
    oledgf_set_project_warnings(
        oledgf_warnings
        ${oledgf_WARNINGS_AS_ERRORS}
        ""
        ""
        ""
        "")

endmacro()
