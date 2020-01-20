ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := ctgcns
SOURCES  := main.c \
	pm4.c \
	load_ctg_read_ids.c \
	cns_one_ctg.c \
	chain_dp.c \
	mem_finder.c \
	cns_ctg_subseq.c \
	fc_correct_one_read.c \
	small_object_alloc.c \

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lontcns
TGT_PREREQS := libontcns.a

SUBMAKEFILES :=
