dnl
dnl This m4 code is designed to extract information from svn:keywords properties.
dnl
m4_divert(-1)
AC_DEFUN([REVISION],[$Revision$])
AC_DEFUN([HEADURL], [$HeadURL$])
AC_DEFUN([AC_SVN_AFTER_COLON],[m4_substr($1,m4_eval(m4_index($1,[:])+2))])
AC_DEFUN([AC_SVN_BEFORE_SPACE],[m4_substr($1,0,m4_index($1,[ ]))])
AC_DEFUN([AC_SVN_BEFORE_SLASH],[m4_substr($1,0,m4_index($1,[/]))])
AC_DEFUN([AC_SVN_GET_INFO],[AC_SVN_BEFORE_SPACE(AC_SVN_AFTER_COLON($1))])
AC_DEFUN([AC_SVN_GET_TAG],[m4_index(GET_INFO($1),[tags])])
AC_DEFUN([AC_SVN_GET_BRANCH],[m4_index(GET_INFO($1),[branches])])
AC_DEFUN([AC_SVN_GET_TRUNK],[m4_index(GET_INFO($1),[trunk])])
AC_DEFUN([AC_SVN_GET_VERSION],[BEFORE_SLASH(m4_substr(AC_SVN_GET_INFO($1),m4_eval(AC_SVN_GET_BRANCH($1)+9)))])
m4_divert(0)dnl
dnl AC_SVN_GET_INFO(REVISION)
dnl AC_SVN_GET_VERSION(HEADURL)
