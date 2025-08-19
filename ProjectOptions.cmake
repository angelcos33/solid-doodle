macro(oledgf_local_options)
    option(oledgf_WARNINGS_AS_ERRORS "Treat Warnings As Errors" OFF)
    option(oledgf_ENABLE_CPPCHECK "Enable cpp-check analysis" OFF)

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
    
    include(cmake/StaticAnalyzers.cmake)
    if (oledgf_ENABLE_CPPCHECK)
      oledgf_enable_cppcheck(${oledgf_WARNINGS_AS_ERRORS}
      "") 
    endif()

endmacro()
