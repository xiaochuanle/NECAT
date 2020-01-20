ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := fsa_ol_filter
SOURCES  := ./prog/fsa_ol_filter.cpp

SRC_INCDIRS  := . 
TGT_CXXFLAGS := -U_GLIBCXX_PARALLEL -std=c++11 -Wall -O3 -D_FILE_OFFSET_BITS=64

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lfsa
TGT_PREREQS := libfsa.a

SUBMAKEFILES :=
