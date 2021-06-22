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
	-I$(top_srcdir)/src \
	-L$(top_srcdir)/src

itrftool_SOURCES = \
	itrf_tool.cpp

AM_LDFLAGS = -lggdatetime -lggeodesy -litrfcpp
