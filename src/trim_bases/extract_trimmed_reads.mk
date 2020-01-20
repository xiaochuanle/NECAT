ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := oc2etr
SOURCES  := extract_trimmed_reads.c largest_cover_range.c range_list.c detect_chimeric_reads.c

SRC_INCDIRS  := . 

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lontcns
TGT_PREREQS := libontcns.a

SUBMAKEFILES :=
