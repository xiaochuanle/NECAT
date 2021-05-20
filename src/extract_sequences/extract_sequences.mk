ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := oc2elr
SOURCES  := extract_sequences.c 

SRC_INCDIRS  := . 

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lontcns -lz
TGT_PREREQS := libontcns.a

SUBMAKEFILES :=
