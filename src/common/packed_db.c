#include "packed_db.h"

#include "nst_nt4_table.h"
#include "oc_assert.h"
#include "ontcns_aux.h"

static const char* pac_header = "ontcns_pac_header_hofuwhogfuewo";

static const char*
extract_hdr_id(kseq_t* read, size_t* id_size)
{
	const char* name = kstr_str(read->name);
	size_t i = 0;
	for(; i < kstr_size(read->name); ++i) {
		if (isspace(name[i])) break;
	}
	*id_size = i;
	return name;
}

void
make_packed_db_name(const char* prefix, kstring_t* packed_db_name)
{
    kstr_clear(*packed_db_name);
    kputs(prefix, packed_db_name);
    kputs(".pac", packed_db_name);
}

static inline idx
size2bytes(const idx size)
{
	idx bytes = (size + 3) >> 2;
	return bytes;
}

void
write_pac_header(FILE* out)
{
    FWRITE(pac_header, sizeof(char), strlen(pac_header), out);
}

BOOL
validate_pac_header(FILE* in)
{
    char buf[256];
    size_t s = strlen(pac_header);
	FREAD(buf, sizeof(char), strlen(pac_header), in);
    buf[s] = '\0';
    BOOL r = (strcmp(buf, pac_header) == 0) ? TRUE : FALSE;
    return r;
}

BOOL
validate_pac_header_path(const char* path)
{
	if (FILE_SIZE(path) < strlen(pac_header)) return FALSE;
	
	DFOPEN(in, path, "rb");
	BOOL r = validate_pac_header(in);
	FCLOSE(in);
	return r;
}

DBType
detect_db_type(const char* path)
{
    if (access(path, F_OK) == -1) {
        OC_ERROR("File '%s' does not exist", path);
    }
    idx s = FILE_SIZE(path);
    if (s == 0) {
        return eEmptyFile;
    }
    if (validate_pac_header_path(path)) {
        return ePac;
    }

    DBType dbt = eUnknown;
    DGZ_OPEN(in, path, "r");
    char c;
    GZ_READ(in, &c, 1);
    GZ_CLOSE(in);
    if (c == '>') dbt = eFasta;
    if (c == '@') dbt = eFastq;
    return dbt;
}

PackedDB*
new_PackedDB()
{
	PackedDB* pdb = (PackedDB*)malloc(sizeof(PackedDB));
	pdb->m_pac = 0;
	pdb->m_db_size = 0;
	pdb->m_max_db_size = 0;
	kv_init(pdb->m_seq_info);
	kstr_init(pdb->m_hdr);
	return pdb;
}

PackedDB*
free_PackedDB(PackedDB* pdb)
{
	free(pdb->m_pac);
	kv_destroy(pdb->m_seq_info);
	free_kstring(pdb->m_hdr);
	free(pdb);
	return 0;
}

void
init_packed_db(PackedDB* pdb)
{
    pdb->m_pac = NULL;
    pdb->m_db_size = 0;
    pdb->m_max_db_size = 0;
    kv_init(pdb->m_seq_info);
	kstr_init(pdb->m_hdr);
}

void
destroy_packed_db(PackedDB* pdb)
{
    if (pdb->m_pac) free(pdb->m_pac);
    pdb->m_db_size = 0;
    pdb->m_max_db_size = 0;
    kv_destroy(pdb->m_seq_info);
	free_kstring(pdb->m_hdr);
}

void
pdb_clear(PackedDB* pdb)
{
	kv_clear(pdb->m_seq_info);
	kstr_clear(pdb->m_hdr);
	if (pdb->m_db_size) {
		size_t pac_bytes = size2bytes(pdb->m_db_size);
		memset(pdb->m_pac, 0, pac_bytes);
		pdb->m_db_size = 0;
	}
}

idx
pdb_size(PackedDB* pdb)
{
    return pdb->m_db_size;
}

idx
pdb_num_seqs(PackedDB* pdb)
{
    return kv_size(pdb->m_seq_info);
}

int
pdb_seq_platform(PackedDB* pdb, const idx i) 
{
    return (kv_A(pdb->m_seq_info, i)).platform;
}

idx 
pdb_seq_offset(PackedDB* pdb, const idx i)
{
    return (kv_A(pdb->m_seq_info, i)).offset;
}

inline idx
pdb_seq_size(PackedDB* pdb, const idx i)
{
    return (kv_A(pdb->m_seq_info, i)).size;
}

idx 
pdb_offset_to_id(PackedDB* pdb, const idx offset)
{
	idx ns = pdb_num_seqs(pdb);
	idx left = 0, mid = 0, right = ns;
	while (left < right) {
		mid = (left + right) >> 1;
		if (offset >= pdb_seq_offset(pdb, mid)) {
			if (mid == ns - 1) break;
			if (offset < pdb_seq_offset(pdb, mid + 1)) break;
			left = mid + 1;
		} else {
			right = mid;
		}
	}
	return mid;
}

void
pdb_enlarge_size(PackedDB* pdb, const idx added_size)
{
	if (pdb->m_db_size + added_size > pdb->m_max_db_size) {
		idx new_size = (pdb->m_max_db_size) ? pdb->m_max_db_size : 4096;
		while (pdb->m_db_size + added_size > new_size) new_size *= 2;
		idx bytes = size2bytes(new_size);
		u8* new_pac = (u8*)calloc(bytes, 1);
		if (pdb->m_pac) {
			memcpy(new_pac, pdb->m_pac, size2bytes(pdb->m_db_size));
			free(pdb->m_pac);
		}
		pdb->m_pac = new_pac;
		pdb->m_max_db_size = new_size;
	}
}

void
pdb_add_string_seq(PackedDB* pdb, kstring_t* seq, int platform)
{
	const char* s = kstr_str(*seq);
	const size_t ssize = kstr_size(*seq);
	SequenceInfo si;
	si.offset = pdb->m_db_size;
	si.size = ssize;
	si.hdr_offset = (idx)-1;
	si.platform = platform;
	kv_push(SequenceInfo, pdb->m_seq_info, si);
	
	pdb_enlarge_size(pdb, ssize);
	for (size_t i = 0; i != ssize; ++i) {
		int oc = s[i];
		u8 ec = nst_nt4_table[oc];
		_set_pac(pdb->m_pac, pdb->m_db_size, ec);
		++pdb->m_db_size;
	}
}

void
pdb_add_one_seq(PackedDB* pdb, kseq_t* seq, int platform)
{
	const char* s = kseq_cstr(*seq);
	const size_t ssize = kseq_size(*seq);
	SequenceInfo si;
	si.offset = pdb->m_db_size;
	si.size = ssize;
	si.hdr_offset = kstr_size(pdb->m_hdr);
	si.platform = platform;
	kv_push(SequenceInfo, pdb->m_seq_info, si);
	const char* name;
	size_t name_size;
	name = extract_hdr_id(seq, &name_size);
	kputsn(name, name_size, &pdb->m_hdr);
	kputc('\0', &pdb->m_hdr);
	
	pdb_enlarge_size(pdb, ssize);
	for (size_t i = 0; i != ssize; ++i) {
		int oc = s[i];
		u8 ec = nst_nt4_table[oc];
		_set_pac(pdb->m_pac, pdb->m_db_size, ec);
		++pdb->m_db_size;
	}
}

void
pdb_extract_subsequence(PackedDB* pdb, const idx i, const idx from, const idx to, const int strand, kstring_t* seq)
{
	const idx s = pdb_seq_offset(pdb, i) + from;
	const idx e = pdb_seq_offset(pdb, i) + to;
	ks_resize(seq, to - from);
	idx pos = 0;
	if (strand == FWD) {
		for (idx k = s; k < e; ++k) {
			u8 c = _get_pac(pdb->m_pac, k);
			kstr_A(*seq, pos) = c;
			++pos;
		}
	} else {
		for (idx k = e; k != s; --k) {
			u8 c = _get_pac(pdb->m_pac, k - 1);
			c = 3 - c;
			kstr_A(*seq, pos) = c;
			++pos;
		}
	}
}

void
pdb_extract_sequence(PackedDB* pdb, const idx i, const int strand, kstring_t* seq)
{
	pdb_extract_subsequence(pdb, i, 0, pdb_seq_size(pdb, i), strand, seq);
}

void
pdb_print_info(PackedDB* pdb)
{
	printf("Number of sequences: %ld\n", pdb_num_seqs(pdb));
	printf("Number of bps: %ld\n", pdb_size(pdb));
}

void
pdb_dump(PackedDB* pdb, const char* path)
{
	DFOPEN(out, path, "wb");
	write_pac_header(out);
	/// 1) number of sequence
	const idx ns = pdb_num_seqs(pdb);
	FWRITE(&ns, sizeof(idx), 1, out);
	/// 2) db size
	FWRITE(&pdb->m_db_size, sizeof(idx), 1, out);
	/// 3) seq info
	FWRITE(kv_data(pdb->m_seq_info), sizeof(SequenceInfo), ns, out);
	/// 4) hdr size
	size_t hdr_size = kstr_size(pdb->m_hdr);
	FWRITE(&hdr_size, sizeof(size_t), 1, out);
	/// 5) hdr
	const char* hdr = kstr_str(pdb->m_hdr);
	if (hdr_size) {
		FWRITE(hdr, 1, hdr_size, out);
	}
	/// 6) pac
	idx pac_bytes = size2bytes(pdb->m_db_size);
	FWRITE(pdb->m_pac, sizeof(u8), pac_bytes, out);
	FCLOSE(out);
}

static void
pdb_load_pac(PackedDB* pdb, const char* path)
{	
	DFOPEN(in, path, "rb");
	if (!validate_pac_header(in)) {
		OC_ERROR("Invalid pac format database: '%s'", path);
	}
	/// 1) number of sequence
	idx ns;
	FREAD(&ns, sizeof(idx), 1, in);
	/// 2) db size
	idx db_size;
	FREAD(&db_size, sizeof(idx), 1, in);
	/// 3) seq info
	kv_resize(SequenceInfo, pdb->m_seq_info, ns);
	FREAD(kv_data(pdb->m_seq_info), sizeof(SequenceInfo), ns, in);
	/// 4) hdr size
	size_t hdr_size;
	FREAD(&hdr_size, sizeof(size_t), 1, in);
	/// 5) hdr
	if (hdr_size) {
		ks_resize(&pdb->m_hdr, hdr_size);
		char* hdr = kstr_str(pdb->m_hdr);
		FREAD(hdr, 1, hdr_size, in);
	}
	/// 6) pac
	pdb_enlarge_size(pdb, db_size);
	idx pac_bytes = size2bytes(db_size);
	FREAD(pdb->m_pac, sizeof(u8), pac_bytes, in);
	pdb->m_db_size = db_size;
	
	FCLOSE(in);
}

static void
pdb_load_fasta(PackedDB* pdb, const char* path, const int platform)
{
	DBType dbt = detect_db_type(path);
	if (dbt == eEmptyFile) return;
	if (dbt == eUnknown) OC_ERROR("unknown file format '%s'", path);
	idx db_size = FILE_SIZE(path);
	if (dbt == eFastq) db_size /= 2;
	db_size += 1000000;
	pdb_enlarge_size(pdb, db_size);
	
	DGZ_OPEN(in, path, "r");
	kseq_t* read = kseq_init(in);
	while (kseq_read(read) >= 0) {
		pdb_add_one_seq(pdb, read, platform);
	}
	GZ_CLOSE(in);
	kseq_destroy(read);
}

static void
pdb_load_worker(PackedDB* pdb, const char* path, const int platform)
{
	DBType dbt = detect_db_type(path);
	if (dbt == eEmptyFile) return;
	if (dbt == eUnknown) {
		OC_ERROR("unknown file format '%s'", path);
	}
	if (dbt == ePac) {
		pdb_load_pac(pdb, path);
	} else {
		pdb_load_fasta(pdb, path, platform);
	}
}

void
pdb_load(PackedDB* pdb, const char* path, const int platform)
{
	pdb_load_worker(pdb, path, platform);
}

void
pdb_merge(PackedDB* from, PackedDB* to)
{
	const idx to_db_size = PDB_SIZE(to);
	const idx to_hdr_size = kstr_size(to->m_hdr);
	SequenceInfo si;
	for (size_t i = 0; i < PDB_NUM_SEQS(from); ++i) {
		si = kv_A(from->m_seq_info, i);
		si.offset += to_db_size;
		si.hdr_offset += to_hdr_size;
		kv_push(SequenceInfo, to->m_seq_info, si);
		
		const char* hdr = PDB_SEQ_NAME(from, i);
		const size_t n = strlen(hdr) + 1;
		kputsn(hdr, n, &to->m_hdr);
	}
	
	idx from_size = PDB_SIZE(from);
	pdb_enlarge_size(to, from_size);
	idx idx_to = to_db_size;
	for (idx i = 0; i < from_size; ++i) {
		u8 c = PDB_GET_CHAR(from, i);
		PDB_SET_CHAR(to, idx_to, c);
		++idx_to;
	}
	to->m_db_size = idx_to;
	
	oc_assert(to->m_db_size <= to->m_max_db_size);
}
