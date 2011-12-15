AC_PREREQ([2.63])
# TODO: change name (test), version (0.0.1) and contact email address
# (me@me.com)
AC_INIT([test], [0.0.1], [me@me.com])
AC_CONFIG_AUX_DIR([build-aux])
AM_SILENT_RULES([yes])
# TODO: test.hpp may need to change to actual file
AC_CONFIG_SRCDIR([src/test.hpp]) 
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
# Uncomment for creating libraries
# AC_PROG_LIBTOOL


# Checks for programs.
AC_PROG_CXX

# Check for yaggo
AC_ARG_VAR([YAGGO], [Yaggo switch parser generator])
AS_IF([test "x$YAGGO" = "x"], [AC_PATH_PROG([YAGGO], [yaggo], [false])])


# Checks for libraries.

# Checks for header files.
#AC_CHECK_HEADERS([stdlib.h string.h])

# Uncomment to get libraries via pkg-config
# PKG_CHECK_MODULES([JELLYFISH], [jellyfish-1.1])


# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_TYPE_SIZE_T
#AC_TYPE_SSIZE_T
#AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
#AC_CHECK_FUNCS([])

AC_OUTPUT([Makefile])
