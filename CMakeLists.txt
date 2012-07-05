CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (nuclear_network)
ENABLE_LANGUAGE(CXX)

SET (CMAKE_VERBOSE_MAKEFILE true)

SET (CMAKE_MODULE_PATH ${nuclear_network_SOURCE_DIR}/cmake/modules
${CMAKE_MODULE_PATH})
INCLUDE(MacroOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out-of-source
build.")

ADD_SUBDIRECTORY (src)
