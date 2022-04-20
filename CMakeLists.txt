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

cmake_minimum_required(VERSION 3.1)

################################################################################
# The version number.
################################################################################
# The default version number is the latest official build
SET (gatb-tool_VERSION_MAJOR 2)
SET (gatb-tool_VERSION_MINOR 3)
SET (gatb-tool_VERSION_PATCH 0)

# But, it is possible to define another release number during a local build
IF (DEFINED MAJOR)
    SET (gatb-tool_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-tool_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-tool_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-tool-version ${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH})

# However, continuous integration has priority over local compilation
IF (DEFINED JENKINS_TAG)
    SET (gatb-tool-version ${JENKINS_TAG})
ENDIF()

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# SUPPORTED KMER SIZES
################################################################################

# One can uncomment this line and set the wanted values
#set (KSIZE_LIST "32   64   96  128  160  192  224  256")

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
#SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
# TOOL
################################################################################


# We also set a flag for TR1 management
if (use_new_cxx)
    set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DUSE_NEW_CXX ")
endif()
 
message("-- cxx: ${CMAKE_CXX_FLAGS}")



# we get compilation definitions from the gatb-core part
add_definitions (${gatb-core-flags})

# SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra" )
SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# we give the headers directories from : 
#       - from project source 
#       - from GATB-CORE source
#       - from dsk source 
set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories (${PROGRAM_SOURCE_DIR} ${gatb-core-includes})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# we define the files to be compiled
file (GLOB  ProjectFiles  src/*)


# we define the artifact to be built: the project binary
add_executable        (${PROJECT_NAME}  src/main.cpp ${ProjectFiles})
add_executable(nwalign src/nwAlign/nwalign.cpp)

# we define which libraries to be linked with project binary
target_link_libraries (${PROJECT_NAME}  ${gatb-core-libraries})
target_link_libraries (nwalign  ${gatb-core-libraries})



################################################################################
#  PACKAGING 
################################################################################

# We set the version number
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
SET (CPACK_PACKAGE_VENDOR               "Genscale team (INRIA)")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-tool_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-tool_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-tool_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-tool-version}")  

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# We ignore unwated files for the source archive
SET (CPACK_SOURCE_IGNORE_FILES          
    "^${PROJECT_SOURCE_DIR}/\\.git/"     ;
    "^${PROJECT_SOURCE_DIR}/\\.gitmodules"  ;
    "^${PROJECT_SOURCE_DIR}/\\.gitignore";
    "^${PROJECT_SOURCE_DIR}/build/"  ;
    "^${GATB_CORE_HOME}/\\.cproject" ;
    "^${GATB_CORE_HOME}/\\.git/"     ;
    "^${GATB_CORE_HOME}/\\.project"  ;
    "^${GATB_CORE_HOME}/\\.gitignore";
    "^${GATB_CORE_HOME}/doc/design"  ;
    "^${GATB_CORE_HOME}/DELIVERY.md"
)

# We copy the project binary to the 'bin' directory
INSTALL (TARGETS   ${PROJECT_NAME}                           DESTINATION bin)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test          DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data          DESTINATION .)
#INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc           DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/README.md     DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/CHANGELOG.md     DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE       DESTINATION .)

# We include the "bin" tag into binary archive file name
set (CPACK_PACKAGE_FILE_NAME  ${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CMAKE_SYSTEM_NAME})

include (CPack)
