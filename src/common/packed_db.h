#ifndef PACKED_DB_H
#define PACKED_DB_H

#include "ontcns_aux.h"
#include "ontcns_defs.h"
#include "../klib/kseq.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"

KSEQ_DECLARE(gzFile)

typedef struct
{
    idx offset;
    idx size;
	idx hdr_offset;
    int platform;
} SequenceInfo;

void
make_packed_db_name(const char* prefix, kstring_t* packed_db_name);

typedef struct {
    u8*                 	m_pac;
    idx                 	m_db_size;
    idx                 	m_max_db_size;
    kvec_t(SequenceInfo)	m_seq_info;   
	kstring_t				m_hdr;
} PackedDB;

#define PDB_NUM_SEQS(pdb) 			(kv_size((pdb)->m_seq_info))
#define PDB_SEQ_OFFSET(pdb, i)		((kv_A((pdb)->m_seq_info, (i))).offset)
#define PDB_SEQ_SIZE(pdb, i)		(kv_A((pdb)->m_seq_info, (i)).size)
#define PDB_GET_CHAR(pdb, i)		(_get_pac((pdb)->m_pac, (i)))
#define PDB_SET_CHAR(pdb, i, c)		(_set_pac((pdb)->m_pac, i, c))
#define PDB_SIZE(pdb)				((pdb)->m_db_size)
#define PDB_SEQ_NAME(pdb, i)		(kstr_str((pdb)->m_hdr) + kv_A((pdb)->m_seq_info, i).hdr_offset)

PackedDB*
new_PackedDB();

PackedDB*
free_PackedDB(PackedDB* pdb);

void
init_packed_db(PackedDB* pdb);

void
destroy_packed_db(PackedDB* pdb);

void
pdb_clear(PackedDB* pdb);

idx
pdb_size(PackedDB* pdb);

idx
pdb_num_seqs(PackedDB* pdb);

int
pdb_seq_platform(PackedDB* pdb, const idx i);

idx 
pdb_seq_offset(PackedDB* pdb, const idx i);

idx
pdb_seq_size(PackedDB* pdb, const idx i);

idx 
pdb_offset_to_id(PackedDB* pdb, const idx offset);

void
pdb_enlarge_size(PackedDB* pdb, const idx added_size);

void
pdb_add_string_seq(PackedDB* pdb, kstring_t* seq, int platform);

void
pdb_add_one_seq(PackedDB* pdb, kseq_t* seq, int platform);

void
pdb_extract_subsequence(PackedDB* pdb, const idx i, const idx from, const idx to, const int strand, kstring_t* s);

void
pdb_extract_sequence(PackedDB* pdb, const idx i, const int strand, kstring_t* s);

void
pdb_print_info(PackedDB* pdb);

void
pdb_dump(PackedDB* pdb, const char* path);

void
pdb_load(PackedDB* pdb, const char* path, const int platform);

void
pdb_merge(PackedDB* from, PackedDB* to);

typedef enum 
{
    eFasta,
    eFastq,
    ePac,
    eEmptyFile,
    eUnknown
} DBType;

DBType
detect_db_type(const char* path);

#endif // PACKED_DB_H
