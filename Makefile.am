# Don't require all the GNU mandated files
#
AUTOMAKE_OPTIONS = 1.11 subdir-objects parallel-tests color-tests dist-zip dist-bzip2 dist-xz foreign

topincludedir = @includedir@/Multipole

AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)

ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = multipole
multipole_SOURCES = multipole.c

check_PROGRAMS = multipole

check-local: multipole
	./multipole 65536 4

# Headers which are installed
pkginclude_HEADERS = multipole.h

# Additional files to distribute
EXTRA_DIST = ChangeLog version.sh
