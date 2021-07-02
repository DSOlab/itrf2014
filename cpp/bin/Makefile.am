bin_PROGRAMS = itrftool

itrftool_CXXFLAGS = \
	-std=c++17 \
	-g \
	-pg \
	-Wall \
	-Wextra \
	-Werror \
	-pedantic \
	-W \
	-Wshadow \
	-Winline \
	-Wdisabled-optimization \
	-DDEBUG \
	-Wno-error=inline \
	-I$(top_srcdir)/src

itrftool_SOURCES = \
	itrf_tool.cpp

AM_LDFLAGS = -lsinex -lggdatetime -lggeodesy

itrftool_LDADD = $(top_srcdir)/src/libitrfcpp.la $(AM_LIBS)