ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := libfsa.a
SOURCES  := argument_parser.cpp logger.cpp overlap.cpp read_store.cpp sequence.cpp\
             utility.cpp file_io.cpp seq_io.cpp overlap_store.cpp \
   			 ./simple_align.cpp overlap_filter.cpp overlap_stat.cpp overlap_trim.cpp

TGT_CXXFLAGS := ${TGT_CXXFLAGS} -U_GLIBCXX_PARALLEL -std=c++11 -Wall -O3 -D_FILE_OFFSET_BITS=64 ${CXXFLAGS}
SRC_INCDIRS  := . -I${PREFIX}/include

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lontcns -lz
TGT_PREREQS := libontcns.a

SUBMAKEFILES :=
