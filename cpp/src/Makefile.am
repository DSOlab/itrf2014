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
	itrf_tool.cpp \
	itrf_tools.hpp \
	itrf_tools.cpp

AM_LDFLAGS = -lggdatetime -lggeodesy
