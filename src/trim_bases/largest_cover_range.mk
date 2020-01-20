ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := oc2lcr
SOURCES  := pm4_aux.c \
	detect_chimeric_reads.c \
	range_list.c \
	largest_cover_range.c \
	largest_cover_range_main.c \
	../consensus/cns_options.c \
	../consensus/consensus_aux.c \
	../consensus/consensus_one_partition.c \
	../consensus/consensus_one_read.c \
	../consensus/error_estimate.c \
	../consensus/overlaps_pool.c \
	../consensus/read_id_pool.c

SRC_INCDIRS  := . 

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lontcns
TGT_PREREQS := libontcns.a

SUBMAKEFILES :=
