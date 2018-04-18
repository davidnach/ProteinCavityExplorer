cmake_minimum_required(VERSION 2.6)

project(XMLPARSE)

set ( XMLPARSESRC xmlParser.cpp )
set ( XMLPARSEINC xmlParser.h )
set ( XMLPARSE_README README )

# make the library
add_library(XmlParse ${XMLPARSESRC})

# install the library into the current xmlParser directory
INSTALL(TARGETS XmlParse
  ARCHIVE DESTINATION ${XMLPARSER_DIR}
)

add_custom_target(dist-level1-external-xmlParser
  COMMENT "Adding xlmParser files into:"
  COMMENT "${RIGDYNCURATION_DIST_LEVEL1_DIR}/external/xmlParser")

add_dependencies(level1dist dist-level1-external-xmlParser)

add_custom_command(
  TARGET dist-level1-external-xmlParser POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  ${XMLPARSER_DIR} 
  ${RIGDYNCURATION_DIST_LEVEL1_DIR}/external/xmlParser
)

