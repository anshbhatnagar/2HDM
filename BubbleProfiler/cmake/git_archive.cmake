if(NOT GIT_ARCHIVE_REPO)
  set(GIT_ARCHIVE_REPO ..)
endif()

if(NOT GIT_FOUND)
  find_package(Git QUIET)
endif()

if(NOT GIT_FOUND)
  message(FATAL_ERROR "Git executable not found")
endif()

if(NOT GIT_ARCHIVE_REF)
  set(GIT_ARCHIVE_REF "HEAD")
endif()

execute_process(COMMAND
  "${GIT_EXECUTABLE}"
  describe
  --tags
  "${GIT_ARCHIVE_REF}"
  WORKING_DIRECTORY "${GIT_ARCHIVE_REPO}"
  RESULT_VARIABLE GIT_DESCRIBE_ERROR
  OUTPUT_VARIABLE GIT_DESCRIBE_OUTPUT
  ERROR_VARIABLE GIT_DESCRIBE_ERROR_MSG
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

if(GIT_DESCRIBE_ERROR)
  message(FATAL_ERROR "${GIT_DESCRIBE_ERROR_MSG}")
endif()

if(GIT_ARCHIVE_PREFIX_ROOT)
  set(GIT_ARCHIVE_PREFIX ${GIT_ARCHIVE_PREFIX_ROOT}-${GIT_DESCRIBE_OUTPUT})
endif()

if(NOT GIT_ARCHIVE_OUTPUT_ROOT)
  message(FATAL_ERROR "Output file must be given")
endif()

if(NOT GIT_ARCHIVE_EXTENSION)
  set(GIT_ARCHIVE_EXTENSION "tar.gz")
endif()

set(GIT_ARCHIVE_OUTPUT ${GIT_ARCHIVE_OUTPUT_ROOT}-${GIT_DESCRIBE_OUTPUT}.${GIT_ARCHIVE_EXTENSION})

execute_process(COMMAND
  "${GIT_EXECUTABLE}"
  archive
  --worktree-attributes
  --prefix=${GIT_ARCHIVE_PREFIX}/
  --output=${GIT_ARCHIVE_OUTPUT}
  "${GIT_ARCHIVE_REF}"
  WORKING_DIRECTORY "${GIT_ARCHIVE_REPO}"
  RESULT_VARIABLE GIT_ARCHIVE_ERROR
  ERROR_VARIABLE GIT_ARCHIVE_ERROR_MSG)

if(GIT_ARCHIVE_ERROR)
  message(FATAL_ERROR "${GIT_ARCHIVE_ERROR_MSG}")
else()
  file(MD5 "${GIT_ARCHIVE_OUTPUT}" GIT_ARCHIVE_OUTPUT_MD5)
  file(WRITE ${GIT_ARCHIVE_OUTPUT}.md5 ${GIT_ARCHIVE_OUTPUT_MD5})
endif()
