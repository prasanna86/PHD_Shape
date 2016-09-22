## util apps
Add_Executable(MapsEllipsoidWithSource apps/MapsEllipsoidWithSource.cxx)

Target_Link_Libraries(MapsEllipsoidWithSource
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
)

ADD_EXECUTABLE(extractlabel apps/extractlabel.cxx)

TARGET_LINK_LIBRARIES(extractlabel
 ${ITK_LIBRARIES}
)

ADD_EXECUTABLE(surftransform apps/surftransform.cxx apps/surfio.cxx)

TARGET_LINK_LIBRARIES(surftransform
 ${ITK_LIBRARIES}
 ${VTK_LIBRARIES})

ADD_EXECUTABLE(surfvolume apps/surfvolume.cxx apps/surfio.cxx)

TARGET_LINK_LIBRARIES(surfvolume
 ${VTK_LIBRARIES})

ADD_EXECUTABLE(surfdecimate apps/surfdecimate.cxx apps/surfio.cxx)

TARGET_LINK_LIBRARIES(surfdecimate
 ${VTK_LIBRARIES})
