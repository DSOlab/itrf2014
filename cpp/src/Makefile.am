lib_LTLIBRARIES = libitrfcpp.la

libitrfcpp_la_CXXFLAGS = \
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
	-Wno-error=inline

dist_libitrfcpp_la_SOURCES = \
	itrf_tools.cpp \
        psd_parser.cpp \
        ssc_parser.cpp

dist_include_HEADERS = \
	itrf_details.hpp \
	itrf_tools.hpp

#AM_LDFLAGS = -lggdatetime -lggeodesy