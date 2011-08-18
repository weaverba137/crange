dnl
dnl This m4 code is designed to extract information from svn:keywords properties.
dnl

# AC_SVN_REVISION
# ---------------
AC_DEFUN([AC_SVN_REVISION],[$Revision$])

# AC_SVN_HEADURL
# --------------
AC_DEFUN([AC_SVN_HEADURL], [$HeadURL$])

# AC_SVN_AFTER_COLON
# ------------------
AC_DEFUN([AC_SVN_AFTER_COLON],[m4_substr($1,m4_eval(m4_index($1,[:])+2))])

# AC_SVN_BEFORE_SPACE
# -------------------
AC_DEFUN([AC_SVN_BEFORE_SPACE],[m4_substr($1,0,m4_index($1,[ ]))])

# AC_SVN_BEFORE_SLASH
# -------------------
AC_DEFUN([AC_SVN_BEFORE_SLASH],[m4_substr($1,0,m4_index($1,[/]))])

# AC_SVN_GET_INFO
# ---------------
AC_DEFUN([AC_SVN_GET_INFO],[AC_SVN_BEFORE_SPACE(AC_SVN_AFTER_COLON($1))])

# AC_SVN_GET_TAG
# --------------
AC_DEFUN([AC_SVN_GET_TAG],[m4_index(GET_INFO($1),[tags])])

# AC_SVN_GET_BRANCH
# -----------------
AC_DEFUN([AC_SVN_GET_BRANCH],[m4_index(GET_INFO($1),[branches])])

# AC_SVN_GET_TRUNK
# ----------------
AC_DEFUN([AC_SVN_GET_TRUNK],[m4_index(GET_INFO($1),[trunk])])

# AC_SVN_GET_VERSION
# ------------------
# Return a suitable version number.
AC_DEFUN([AC_SVN_GET_VERSION],
[BEFORE_SLASH(m4_substr(AC_SVN_GET_INFO($1),
              m4_eval(AC_SVN_GET_BRANCH($1)+9)))dnl
])# AC_SVN_GET_VERSION

