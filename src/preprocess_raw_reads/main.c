#include "../common/ontcns_aux.h"
#include "../common/nst_nt4_table.h"
#include "../common/packed_db.h"
#include "../common/check_nonrepeat_suffix.h"

#include "truncate_end_spaces.h"

static void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s file_list_in out_folder file_list_out min_seq_length\n", prog);
}

/*
static void
process_fasta_file(const char* path, int* vid, const char* out_folder, FILE* list_out)
{
	char out_path[2048];
	sprintf(out_path, "%s/raw_reads_%d.fasta", out_folder, *vid);
	fprintf(list_out, "%s\n", out_path);
	char link_cmd[2048];
	sprintf(link_cmd, "ln -s %s %s", path, out_path);
	SYSTEM(link_cmd);
}
*/

static void
print_kstring(FILE* out, kstring_t* s)
{
	for (size_t i = 0; i < kstr_size(*s); ++i) {
		fprintf(out, "%c", kstr_A(*s, i));
	}
	fprintf(out, "\n");
}

static BOOL
check_data_line(kstring_t* seq, const char* filename, const int linenum)
{
	const char* s = kstr_data(*seq);
	const size_t n = kstr_size(*seq);
	BOOL r = TRUE;
	for (size_t i = 0; i < n; ++i) {
		int c = s[i];
		c = nst_nt4_table[c];
		if (c > 3) {
			r = FALSE;
			break;
		}
	}
	if (!r) {
		fprintf(stderr, "detect a misline at (%s, %d):\n", filename, linenum);
		print_kstring(stderr, seq);
	}
	return r;
}

#define KGetLine(line) kstr_clear(line); r = kgetline_gz(&(line), in); ++linenum

static void
process_fasta_file(const char* path, 
				   int* vid, 
				   const char* out_folder, 
				   FILE* list_out,
				   FILE* garbage_out, 
				   const int min_seq_length)
{
	char out_path[2048];
	sprintf(out_path, "%s/raw_reads_%d.fasta", out_folder, *vid);
	++(*vid);
	fprintf(list_out, "%s\n", out_path);
	DFOPEN(out, out_path, "w");
	new_kstring(hdr);
	new_kstring(seq);
	DGZ_OPEN(in, path, "r");
	int linenum = 1;
	int r;
	
	while (1) {
		BOOL is_done = FALSE;
		while (1) {
			KGetLine(hdr);
			if (r == EOF) {
				is_done = TRUE;
				break;
			}
			if (kstr_size(hdr) > 0 && kstr_A(hdr, 0) == '>') break;
		}
		if (is_done) break;
		KGetLine(seq);
		if (r == EOF) break;
		if (!check_data_line(&seq, path, linenum)) continue;
		
		if (kstr_size(seq) >= min_seq_length) { //&& is_nonrepeat_sequence(kstr_str(seq), kstr_size(seq))) {
			print_kstring(out, &hdr);
			print_kstring(out, &seq);
		} else {
			print_kstring(garbage_out, &hdr);
			print_kstring(garbage_out, &seq);
		}
	}
	
	FCLOSE(out);
	GZ_CLOSE(in);
	free_kstring(hdr);
	free_kstring(seq);
}

static void
process_fastq_file(const char* path, int* vid, const char* out_folder, FILE* list_out, FILE* garbage_out, const int min_seq_length)
{
	OC_LOG("process fastq file %s", path);
	new_kstring(hdr);
	new_kstring(seq);
	new_kstring(plus);
	new_kstring(qual);
	DGZ_OPEN(in, path, "r");
	int linenum = 1;
	int r;
	char out_path[2048];
	sprintf(out_path, "%s/raw_reads_%d.fastq", out_folder, *vid);
	++(*vid);
	fprintf(list_out, "%s\n", out_path);
	DFOPEN(out, out_path, "w");
	
	while (1) {
		BOOL is_done = FALSE;
		while (1) {
			KGetLine(hdr);
			if (r == EOF) {
				is_done = TRUE;
				break;
			}
			if (kstr_size(hdr) > 0 && kstr_A(hdr, 0) == '@') break;
		}
		if (is_done) break;
		KGetLine(seq);
		if (r == EOF) break;
		if (!check_data_line(&seq, path, linenum)) continue;
		KGetLine(plus);
		if (r == EOF) break;
		if (kstr_size(plus) == 0 || kstr_A(plus, 0) != '+') {
			fprintf(stderr, "Invalid third line of a fastq sequence at (%s, %d):\n", path, linenum);
			print_kstring(stderr, &plus);
			continue;
		}
		KGetLine(qual);
		if (r == EOF) break;
		if (kstr_size(qual) != kstr_size(seq)) {
			fprintf(stderr, "sequence and quality line do not have the same length at (%s, %d)\n", path, linenum);
			print_kstring(stderr, &seq);
			print_kstring(stderr, &qual);
			continue;
		}
		if (kstr_size(seq) >= min_seq_length) {
			print_kstring(out, &hdr);
			print_kstring(out, &seq);
			fprintf(out, "+\n");
			print_kstring(out, &qual);
		} else {
			print_kstring(garbage_out, &hdr);
			print_kstring(garbage_out, &seq);
			fprintf(garbage_out, "+\n");
			print_kstring(garbage_out, &qual);
		}
	}
	FCLOSE(out);
	GZ_CLOSE(in);
	free_kstring(hdr);
	free_kstring(seq);
	free_kstring(plus);
	free_kstring(qual);
}

FILE* open_garbage_file(const char* out_foler)
{
	char path[2048];
	sprintf(path, "%s/garbage_reads.fastq", out_foler);
	DFOPEN(out, path, "w");
	return out;
}

int main(int argc, char* argv[])
{
	if (argc != 5) {
		print_usage(argv[0]);
		return 1;
	}
	
	const char* file_list_in_path = argv[1];
	const char* out_folder = argv[2];
	const char* file_list_out_path = argv[3];
	const int min_seq_length = atoi(argv[4]);
	char line[2048];
	DGZ_OPEN(in, file_list_in_path, "r");
	DFOPEN(out, file_list_out_path, "w");
	FILE* garbage_out = open_garbage_file(out_folder);
	int vid = 0;
	while (gzgets(in, line, 2048)) {
		truncate_end_spaces(line);
		size_t n = strlen(line);
		if (n == 0) continue;
		if (line[n - 1] == '\n') {
			if (n > 1 && line[n - 2] == '\r') line[n - 2] = '\0';
			else line[n - 1] = '\0';
		}
		
		DBType t = detect_db_type(line);
		if (t == eEmptyFile) continue;
		if (t == eUnknown) {
			OC_ERROR("unknown file format: %s", line);
		}
		if (t == eFasta) {
			process_fasta_file(line, &vid, out_folder, out, garbage_out, min_seq_length);
		} else if (t == eFastq) {
			process_fastq_file(line, &vid, out_folder, out, garbage_out, min_seq_length);
		}
	}
	FCLOSE(out);
	FCLOSE(garbage_out);
	GZ_CLOSE(in);
}
