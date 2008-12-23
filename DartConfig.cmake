# Dart server to submit results (used by client)
SET (DROP_METHOD "scp")
SET (DROP_SITE "itk.doc.ic.ac.uk")
SET (DROP_LOCATION "/vol/vipcvs/Dart/projectincoming")
#SET (DROP_SITE_USER "rc3")
#SET (DROP_SITE_PASSWORD "insight-tester@somewhere.com")
SET (TRIGGER_SITE "http://www.doc.ic.ac.uk/~rc3/Dashboard/TestingResults/Dart.pl")

# Dart server configuration 
SET (CVS_WEB_URL "http://www.doc.ic.ac.uk/~dr/cvsweb/cvsweb.cgi/project/")
SET (CVS_WEB_CVSROOT "/vol/vipcvs/CVS/project/&f=h")

#OPTION(BUILD_DOXYGEN "Build source documentation using doxygen" "On")
#SET (DOXYGEN_URL "http://www.doc.ic.ac.uk/~rc3/Dashboard/Documentation")
#SET (DOXYGEN_CONFIG "${PROJECT_BINARY_DIR}/documentation/doxyfile" )

#SET (USE_GNATS "On")
#SET (GNATS_WEB_URL "http://${DROP_SITE}/cgi-bin/gnatsweb.pl/Insight/")
