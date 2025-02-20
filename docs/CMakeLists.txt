cmake_minimum_required(VERSION 3.18)
project(Ikarus-Docs)

find_program(
  MKDOCS_EXECUTABLE
  NAMES mkdocs
  DOC "MkDocs documentation generation tool (http://www.mkdocs.org)" REQUIRED
)

add_subdirectory(website/doxygen)

add_custom_target(
  site
  COMMAND mkdocs build --config-file mkdocs.insiders.yml
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_dependencies(site doxygen_ikarus)

file(GLOB_RECURSE markdownFiles *.md)

add_custom_target(
  localSite
  COMMAND ${MKDOCS_EXECUTABLE} serve
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  DEPENDS mkdocs-macros.py mkdocs.yml ${markdownFiles}
)

add_dependencies(localSite doxygen_ikarus)

file(COPY mkdocs.yml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY mkdocs.insiders.yml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY overrides DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY website DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY literature.bib DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY mkdocs-macros.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
