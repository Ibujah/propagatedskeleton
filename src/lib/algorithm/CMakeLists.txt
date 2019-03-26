set(LIBRARY_NAME ${ALGORITHM_LIB})

include_directories(${CMAKE_SOURCE_DIR}/src/lib
					${VOROPLUSPLUS_INCLUDE_DIR}
					${OpenCV_INCLUDE_DIRS}
					${Boost_INCLUDE_DIR}
					${NLOPT_INCLUDE_DIR})
					
set(SOURCE_FILES    extractboundary/NaiveBoundary.cpp
					skeletonization/propagation/SpherePropagation2D.cpp
					skeletonization/propagation/DistanceField.cpp
					skeletonization/propagation/MovingCenter.cpp
					skeletonization/voronoi/VoronoiSkeleton2D.cpp
					graphoperation/ConnectedComponents.cpp
					skinning/Filling.cpp
					evaluation/ReprojError.cpp
					pruning/ScaleAxisTransform.cpp
					pruning/LambdaMedialAxis.cpp
					pruning/ThetaMedialAxis.cpp)
# make the library
add_library(
    ${LIBRARY_NAME}
    ${INCLUDE_FILES}
    ${SOURCE_FILES}
	)

SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(${LIBRARY_NAME} ${MATHTOOLS_LIB}
									  ${SHAPE_LIB}
									  ${FIELD_LIB}
									  ${SKELETON_LIB}
									  ${BOUNDARY_LIB}
									  ${CAMERA_LIB}
									  ${VOROPLUSPLUS_LIBRARY}
									  ${NLOPT_LIBRARY})

