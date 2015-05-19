################################################################################
 #   MindTheGap: Integrated detection and assembly of insertion variants
 #   A tool from the GATB (Genome Assembly Tool Box)
 #   Copyright (C) 2014  INRIA
 #   Authors: C.Lemaitre, G. Rizk
 #
 #  This program is free software: you can redistribute it and/or modify
 #  it under the terms of the GNU Affero General Public License as
 #  published by the Free Software Foundation, either version 3 of the
 #  License, or (at your option) any later version.
 #
 #  This program is distributed in the hope that it will be useful,
 #  but WITHOUT ANY WARRANTY; without even the implied warranty of
 #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #  GNU Affero General Public License for more details.
 #
 #  You should have received a copy of the GNU Affero General Public License
 #  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ################################################################################

project(MindTheGap)

cmake_minimum_required (VERSION 2.6)

################################################################################
# Define cmake modules directory
################################################################################
FOREACH (path "cmake" "../../cmake")
    IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
        SET (CMAKE_MODULE_PATH  "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
    ENDIF()
ENDFOREACH(path)

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
 SET (GATB_CORE_EXCLUDE_TESTS     1)
 SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

#line commented below : we want to keep the tools dbgh5, dbginfo, etc.
#SET (GATB_CORE_EXCLUDE_TOOLS     1)


# GATB CORE
include (GatbCore)

################################################################################
# TOOL
################################################################################

# we get compilation definitions from the gatb-core part
add_definitions (${gatb-core-flags})

# SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra" )
SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# we give the headers directories from : 
#       - from project source 
#       - from GATB-CORE source
#       - from dsk source 
include_directories (include ${gatb-core-includes})

# we define the files to be compiled
file (GLOB_RECURSE  ProjectFiles  src/*)

# we define the artifact to be built: the project binary
add_executable        (${PROJECT_NAME}  src/main.cpp ${ProjectFiles})

# we define which libraries to be linked with project binary
target_link_libraries (${PROJECT_NAME}  ${gatb-core-libraries})


################################################################################
#  PACKAGING 
################################################################################

# We chose the kind of archive we want to generate
SET (CPACK_GENERATOR            "TGZ")
SET (CPACK_SOURCE_GENERATOR     "TGZ")

# We copy the project binary to the 'bin' directory
INSTALL (TARGETS ${PROJECT_NAME} DESTINATION bin)

# We copy other stuff 
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data          DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc           DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/README        DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE       DESTINATION .)

# For the source archive, we exclude some unwanted directories.
SET (CPACK_SOURCE_IGNORE_FILES  ${CPACK_SOURCE_IGNORE_FILES} "/build*" "/distrib*"   "/test*" "../../cmake/Delivery.cmake" "genere_distrib.sh"  ".gitignore$")

# To be done at the end. 
INCLUDE (Delivery)

