# Gives access to SET_DEFAULT_AND_FROM_ENV function
INCLUDE(SetDefaultAndFromEnv)

#Settings for CTest/CDash
IF(NOT DEFINED CTEST_DROP_METHOD)
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_METHOD "https")
ENDIF()

IF(CTEST_DROP_METHOD STREQUAL "https")
  #The normal default is ${HOST_TYPE}-${COMPILER_VERSION}-${BUILD_DIR}
  #we are over-riding it by using set.
  SET(CTEST_BUILD_NAME "Linux-HW2-<user_not_set>")

  #Not sure if these does anything?
  SET_DEFAULT_AND_FROM_ENV(CTEST_PROJECT_NAME "HW2")
  SET_DEFAULT_AND_FROM_ENV(CTEST_TRIGGER_SITE "")

  #Can you submit to projects that you are not a member of?
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE "cdash-ners590.aura.arc-ts.umich.edu")
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_LOCATION "/submit.php?project=HW2")
  SET_DEFAULT_AND_FROM_ENV(CTEST_DROP_SITE_CDASH TRUE)
ENDIF()
