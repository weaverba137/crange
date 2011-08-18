dnl
dnl This m4 code is designed to extract information from svn:keywords properties.
dnl
divert(-1)
define(`REVISION',`$Revision$')
define(`HEADURL', `$HeadURL$')
define(`AFTER_COLON',`substr($1,eval(index($1,`:')+2))')
define(`BEFORE_SPACE',`substr($1,0,index($1,` '))')
define(`BEFORE_SLASH',`substr($1,0,index($1,`/'))')
define(`GET_INFO',`BEFORE_SPACE(AFTER_COLON($1))')
define(`GET_TAG',`index(GET_INFO($1),`tags')')
define(`GET_BRANCH',`index(GET_INFO($1),`branches')')
define(`GET_TRUNK',`index(GET_INFO($1),`trunk')')
define(`GET_VERSION',`BEFORE_SLASH(substr(GET_INFO($1),eval(GET_BRANCH($1)+9)))')
divert(0)dnl
dnl GET_INFO(REVISION)
dnl GET_VERSION(HEADURL)
