include(cmake/CPM.cmake)

function(oledgf_setup_dependencies)

if(NOT TARGET fmtlib::fmtlib)
    CPMAddPackage("gh:fmtlib/fmt#9.1.0")
endif()
    
endfunction()
