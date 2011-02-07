if(MINIMAL_BUILD)
  set(tools
  )
else(MINIMAL_BUILD)
  set(tools
  )
endif(MINIMAL_BUILD)

INCLUDE(CheckFunctionExists)
SET(GETOPT getopt)
CHECK_FUNCTION_EXISTS(${GETOPT} GETOPT_FOUND)

if(BUILD_SHARED)
  foreach(tool ${tools})
    add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} ${BABEL_LIBRARY})
    if(MSVC)
      set_target_properties(${tool} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    endif(MSVC)
    install(TARGETS ${tool}
                    RUNTIME DESTINATION bin
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib
    )
  endforeach(tool)

  foreach(tool confab calcrmsd)
    if (NOT GETOPT_FOUND)
      add_executable(${tool} ${tool}.cpp getopt.c)
    else(NOT GETOPT_FOUND)
      add_executable(${tool} ${tool}.cpp)
    endif(NOT GETOPT_FOUND)
    target_link_libraries(${tool} ${BABEL_LIBRARY})
    if(MSVC)
      set_target_properties(${tool} PROPERTIES
          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    endif(MSVC)
    install(TARGETS ${tool}
                    RUNTIME DESTINATION bin
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib
    )
  endforeach(tool)

else(BUILD_SHARED)
  if (NOT MSVC)
    if (BUILD_MIXED)
      #
      # Relocatable binaries on linux using -static-libgcc
      #
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
      # make sure the linker will find the static libstdc++
      file(REMOVE "${CMAKE_BINARY_DIR}/libstdc++.a")
      execute_process(COMMAND "g++" "-print-file-name=libstdc++.a"
                      OUTPUT_VARIABLE libstdcpp
                      OUTPUT_STRIP_TRAILING_WHITESPACE)
      execute_process(COMMAND "ln" "-s" "${libstdcpp}"
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
      link_directories(${CMAKE_BINARY_DIR})
    else(BUILD_MIXED)
      #
      # Default static build on linux.
      #
      set(CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
    endif(BUILD_MIXED)
  endif(NOT MSVC)

  # build a static babel only
  add_executable(babel babel.cpp static_plugins.cpp)
  target_link_libraries(babel ${BABEL_LIBRARY} ${ZLIB_LIBRARY})
    if(MSVC)
      set_target_properties(babel PROPERTIES
          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    endif(MSVC)
    install(TARGETS babel
                    RUNTIME DESTINATION bin
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib
    )

endif(BUILD_SHARED)
