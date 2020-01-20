ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libontcns.a

SOURCES      := \
	./common/check_nonrepeat_suffix.cpp \
	./common/cns_seq.c \
	./common/gapped_candidate.c \
	./common/ontcns_aux.c \
	./common/ontcns_defs.c \
	./common/m4_record.c \
	./common/makedb_aux.c \
	./common/map_aux.c \
	./common/map_options.c \
	./common/nst_nt4_table.c \
	./common/oc_assert.c \
	./common/packed_db.c \
	./common/record_reader.c \
	./common/record_writer.c \
	./common/soa.c \
	./common/symdust.cpp \
	./edlib/edlib.c \
	./edlib/edlib_wrapper.c \
	./gapped_align/blockwise_align.c \
	./gapped_align/edlib_ex.c \
	./gapped_align/edlib_ex_aux.c \
	./gapped_align/oc_aligner.c \
	./gapped_align/align.c \
	./gapped_align/oc_daligner.c \
	./klib/kalloc.c \
	./klib/kstring.c \
	./lookup_table/hash_list_bucket_sort.c \
	./lookup_table/lookup_table.c \
	./partition_candidates/pcan_aux.c \
	./tasc/align_tags.c \
	./tasc/cbcns.c \
	./tasc/cns_aux.c \
	./tasc/align_tags.c \
	./word_finder/chain_dp.c \
	./word_finder/word_finder_aux.c \
	./word_finder/word_finder.c

SRC_INCDIRS  := common \

SUBMAKEFILES := ./test/main.mk \
	./makedb/main.mk \
	./pm_one_volume/main.mk \
	./pairwise_mapping/main.mk \
	./partition_candidates/main.mk \
	./consensus/main.mk \
	./sequence_length_stats/main.mk \
	./split_long_reads/main.mk \
	./asm_pm/asmpm.mk \
	./partition_m4/main.mk \
	./reference_mapping/main.mk \
	./trim_bases/pm4.mk \
	./trim_bases/largest_cover_range.mk \
	./trim_bases/trim_bases.mk \
	./trim_bases/order_results.mk \
	./trim_bases/extract_trimmed_reads.mk \
	./trim_bases_accurate0/pm4.mk \
	./trim_bases_accurate0/largest_cover_range.mk \
	./trim_bases_accurate0/trim_bases.mk \
	./trim_bases_accurate/largest_cover_range.mk \
	./trim_bases_accurate/pm4.mk \
	./trim_bases_accurate/trim_bases.mk \
	./extract_sequences/extract_sequences.mk \
	./test/main.mk \
	./renumber_sequences/main.mk \
	./reorder_cns_reads/main.mk \
	./ctgpm/split_ctgs.mk \
	./ctgpm/fix_can_info.mk \
	./ctgpm/ctgpm.mk \
	./preprocess_raw_reads/main.mk \
	./fsa/fsa.mk \
	./fsa/filter.mk \
	./fsa/assemble.mk \
	./fsa/bridge.mk \
	./fsa/rd_stat.mk \
	./fsa/rd_extract.mk \
	./fsa/rd_tools.mk \
	./pipeline/main.mk \
	./reference_mapping/rm_one_vol_main.mk \
	./ctg_cns/filter_m4.mk \
	./ctg_cns/main.mk \
	./ctg_cns/pm4_main.mk \
	../tool/main.mk 

