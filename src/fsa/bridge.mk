ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := fsa_ctg_bridge
SOURCES  := ./prog/fsa_ctg_bridge.cpp contig_bridge.cpp contig_graph.cpp contig_link.cpp contig_link_store.cpp

SRC_INCDIRS  := . 
TGT_CXXFLAGS := -U_GLIBCXX_PARALLEL -std=c++11 -Wall -O3 -D_FILE_OFFSET_BITS=64 ${CXXFLAGS}

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lfsa -lz -lpthread
TGT_PREREQS := libfsa.a

SUBMAKEFILES :=
