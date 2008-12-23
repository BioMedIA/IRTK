PROJECT(RVIEW_FLTK)

LINK_LIBRARIES(rview++ registration++ transformation++ image++ geometry++
common++ recipes)

SET(RVIEW_FLTK_SRCS
Fl_Simple_Group.cc
Fl_Simple_Browser.cc
Fl_RView.cc
Fl_RViewUI.cc
Fl_RView_Landmark.cc
Fl_RView_Registration.cc
Fl_RView_Histogram.cc
Fl_RView_Transformation.cc
Fl_RView_Visualisation.cc
Fl_RView_Segmentation.cc
main.cc)

ADD_EXECUTABLE(rview ${RVIEW_FLTK_SRCS})

INSTALL_TARGETS(/bin rview)

IF(APPLE)
  FIND_PROGRAM(APPLE_REZ Rez /Developer/Tools)
  ADD_CUSTOM_COMMAND(TARGET rview
    POST_BUILD
    COMMAND ${APPLE_REZ}
    ARGS ${FLTK_INCLUDE_DIR}/mac.r -o ${EXECUTABLE_OUTPUT_PATH}/rview)
ENDIF(APPLE)
