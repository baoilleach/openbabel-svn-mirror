
file(GLOB openbabel_headers openbabel/*.h)
file(GLOB math_headers      openbabel/math/*.h)
file(GLOB stereo_headers      openbabel/stereo/*.h)
file(GLOB depict_headers      openbabel/stereo/*.h)

install(FILES inchi_api.h          DESTINATION include/inchi)
install(FILES chemdrawcdx.h        DESTINATION ${OB_INCLUDE_DIRS}/openbabel)
install(FILES ${openbabel_headers} DESTINATION ${OB_INCLUDE_DIRS}/openbabel)
install(FILES ${math_headers}      DESTINATION ${OB_INCLUDE_DIRS}/openbabel/math)
install(FILES ${stereo_headers}    DESTINATION ${OB_INCLUDE_DIRS}/openbabel/stereo)
