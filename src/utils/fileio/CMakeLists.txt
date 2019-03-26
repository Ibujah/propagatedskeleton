set(LIBRARY_NAME ${FILEIO_LIB})

include_directories(${CMAKE_SOURCE_DIR}/src/lib
					${Boost_INCLUDE_DIR})
set(SOURCE_FILES   SkeletonFile.cpp
				   BoundaryFile2D.cpp)

# make the library
add_library(
    ${LIBRARY_NAME}
    ${INCLUDE_FILES}
    ${SOURCE_FILES}
	)

target_link_libraries(${LIBRARY_NAME} ${OpenCV_LIBS})

SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES LINKER_LANGUAGE CXX)

