dnl this is a test
divert(-1)
define(`REVISION',`$Revision: 1234 $')
define(`HEADURL', `$HeadURL: http://svn.example.com/svn/crange/tags/1.2.3/m4/svn.m4 $')
define(`AFTER_COLON',`substr($1,eval(index($1,`:')+2))')
define(`BEFORE_SPACE',`substr($1,0,index($1,` '))')
define(`BEFORE_SLASH',`substr($1,0,index($1,`/'))')
define(`GET_REVISION',`BEFORE_SPACE(AFTER_COLON($1))')
define(`GET_TAG',`BEFORE_SLASH(substr(GET_REVISION($1),eval(index(GET_REVISION($1),`tags')+5)))')
divert(0)dnl
GET_REVISION(REVISION)
GET_TAG(HEADURL)
