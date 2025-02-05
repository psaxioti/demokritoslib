cmake_minimum_required(VERSION 3.16)

project(demokritos VERSION 1.0)

option(STATIC "Build static executable" OFF)
if(STATIC)
   set(LINKING_TYPE STATIC)
else()
   set(LINKING_TYPE SHARED)
endif(STATIC)


#find_package(GSL REQUIRED)

set(RC_DEPENDS "")
# If you want to make the symbol name look prettier, just use relative path as the input
function(add_resource input)
  string(MAKE_C_IDENTIFIER ${input} input_identifier)
  set(res_intermediate_dir ${CMAKE_CURRENT_BINARY_DIR}/resources)
  set(res_with_null_output "${res_intermediate_dir}/${input}")
  set(output "${res_intermediate_dir}/${input_identifier}.o")
  # Add null-terminated character to the file
  add_custom_command(
      DEPENDS ${input}
      OUTPUT ${res_with_null_output}
      COMMAND ${CMAKE_COMMAND} -E copy ${input} ${res_with_null_output};
#      COMMAND echo -n '\\0' >> ${res_with_null_output}
      WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
  )
  add_custom_command(
      DEPENDS ${res_with_null_output}
      OUTPUT ${output}
      COMMAND ${CMAKE_LINKER} -z noexecstack --relocatable --format binary --output ${output} ${input}
      WORKING_DIRECTORY ${res_intermediate_dir}
  )
  set(RC_DEPENDS ${RC_DEPENDS} ${output} PARENT_SCOPE)
endfunction()
# Resource file list
add_custom_target( rc ALL DEPENDS ${RC_DEPENDS} )



# Download AME2020 masses file
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/data/masses.dat)
   file(DOWNLOAD https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt ${CMAKE_SOURCE_DIR}/data/masses.dat)
endif()
#add_resource("data/masses.dat" "demmasses")

# Download RIPL3 adundances file
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/data/abundance.dat)
   file(DOWNLOAD https://www-nds.iaea.org/RIPL/masses/abundance.dat ${CMAKE_SOURCE_DIR}/data/abundance.dat)
endif()



include_directories(include)


set(RC_DEPENDS "")
add_resource("data/masses.dat")
add_resource("data/abundance.dat")

add_library(demmasses ${LINKING_TYPE}
   src/MassFunctions.cc
   src/Element.cc
   src/Isotope.cc
   ${RC_DEPENDS}
)

target_sources(demmasses PUBLIC FILE_SET HEADERS BASE_DIRS include FILES include/MassFunctions.hh include/Element.hh include/Isotope.hh)
install(TARGETS demmasses
    COMPONENT demmasses
    LIBRARY
    FILE_SET HEADERS
)

set(RC_DEPENDS "")

add_library(demreaction ${LINKING_TYPE}
   src/Reaction.cc
   src/ReactionList.cc
   ${RC_DEPENDS}
)

target_sources(demreaction PUBLIC FILE_SET HEADERS BASE_DIRS include FILES include/Reaction.hh include/ReactionList.hh)
install(TARGETS demreaction
    COMPONENT demreaction
    LIBRARY
    FILE_SET HEADERS
)

set(RC_DEPENDS "")
add_resource("data/scoef.95a")
add_resource("data/scoef.95b")
add_resource("data/srim.txt")
add_resource("data/chu.dat")

add_library(demstopping ${LINKING_TYPE}
   src/Stopping.cc
   ${RC_DEPENDS}
)

target_sources(demstopping PUBLIC FILE_SET HEADERS BASE_DIRS include FILES include/Stopping.hh)
install(TARGETS demstopping
    COMPONENT demstopping
    LIBRARY
    FILE_SET HEADERS
)

set(RC_DEPENDS "")

add_library(demtarget ${LINKING_TYPE}
   src/Target.cc
   src/TargetElement.cc
   src/TargetLayer.cc
   ${RC_DEPENDS}
)

target_sources(demtarget PUBLIC FILE_SET HEADERS BASE_DIRS include FILES include/Target.hh include/FitStateEnum.hh)
install(TARGETS demtarget
    COMPONENT demtarget
    LIBRARY
    FILE_SET HEADERS
)

# set(THREADS_PREFER_PTHREAD_FLAG ON)
# find_package(Threads REQUIRED)
# set (CMAKE_CXX_STANDARD 11)

# set(RC_DEPENDS "")

# add_library(demthreadpool ${LINKING_TYPE}
#    src/ThreadPool.cc
#    ${RC_DEPENDS}
# )

# target_sources(demthreadpool PUBLIC FILE_SET HEADERS BASE_DIRS include FILES include/ThreadPool.hh)
# target_link_libraries(demthreadpool ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_PTHREAD_LIBS_INIT})
# install(TARGETS demthreadpool
#     COMPONENT demthreadpool
#     LIBRARY
#     FILE_SET HEADERS
# )
