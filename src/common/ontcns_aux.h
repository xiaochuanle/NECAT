#ifndef ONTCNS_AUX_H
#define ONTCNS_AUX_H

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "ontcns_defs.h"
#include "../klib/ksort.h"
#include "../klib/kstring.h"

/// log info

void
what_is_time_now(char now[]);

void
oc_dump_info(FILE* out,
			const char* level,
			const char* file,
			const int line,
			const char* fmt,
			...);

#define OC_LOG(fmt, args...) oc_dump_info(stdout, "INFO", __FILE__, __LINE__, fmt, ##args)
#define OC_LOG_FULL(file, line, fmt, args...) oc_dump_info(stdout, "INFO", file, line, fmt, ##args)
#define OC_ERROR(fmt, args...) do { oc_dump_info(stderr, "ERROR", __FILE__, __LINE__, fmt, ##args); exit(1); } while(0)
#define OC_ERROR_FULL(file, line, fmt, args...) oc_dump_info(stderr, "ERROR", file, line, fmt, ##args); abort();

/// gzfile
int 
safe_gzread(const char* src_file,
		    const int src_line,
		    gzFile file, 
		    void *ptr, 
		    unsigned int len);

int 
err_gzread(gzFile file, void* ptr, unsigned int len);

gzFile 
safe_gzopen(const char* src_file, const int src_line, const char *fn, const char *mode);

int
safe_gzclose(const char* src_file, const int src_line, gzFile file);

#define GZ_OPEN(file, path, mode) file = safe_gzopen(__FILE__, __LINE__, path, mode)
#define DGZ_OPEN(file, path, mode) gzFile file; GZ_OPEN(file, path, mode)
#define GZ_CLOSE(file) safe_gzclose(__FILE__, __LINE__, file)
#define GZ_READ(file, ptr, len) safe_gzread(__FILE__, __LINE__, file, ptr, len)

/// file

#define SAFE_SCANF(input_func, stream, nread, format, ...) \
	do { \
	int __nread = input_func(stream, format, __VA_ARGS__); \
		if (nread != __nread) { \
		OC_ERROR("read error: should read %d iterms, but only %d has been read", nread, __nread); \
	} \
} while(0)

FILE*
safe_fopen(const char* src_file,
		   const int src_line,
		   const char* path,
		   const char* mode);

size_t
safe_fwrite(const char* src_file,
			const int src_line,
			const void* ptr,
			size_t size,
			size_t nmemb,
			FILE* stream);

size_t
safe_fread(const char* src_file,
		   const int src_line,
		   void* ptr,
		   size_t size,
		   size_t nmemb,
		   FILE* stream);

int 
safe_fclose(const char* src_file,
			const int src_line,
			FILE* stream);

#define FOPEN(file, path, mode) file = safe_fopen(__FILE__, __LINE__, path, mode)
#define DFOPEN(file, path, mode); FILE* file; FOPEN(file, path, mode)
#define FWRITE(ptr, size, n, stream) safe_fwrite(__FILE__, __LINE__, ptr, size, n, stream)
#define FREAD(ptr, size, n, stream) safe_fread(__FILE__, __LINE__, ptr, size, n, stream)
#define FCLOSE(file) safe_fclose(__FILE__, __LINE__, file)

/// timing 

double 
time_diff(const struct timeval* start, const struct timeval* end);

#define TIMING_START(title) \
	struct timeval __oc_start; \
	gettimeofday(&__oc_start, NULL); \
	OC_LOG("'%s' BEGINS", title)
   
#define TIMING_END(title) \
	struct timeval __oc_end; \
	gettimeofday(&__oc_end, NULL); \
	double __oc_time = time_diff(&__oc_start, &__oc_end); \
	OC_LOG("'%s' takes %.2lf secs.", title, __oc_time);

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

idx
file_size(const char* src_path, const int src_line, const char* path);

#define FILE_SIZE(path) file_size(__FILE__, __LINE__, path)

char*
ontcns_getline(char* buf, int bufsize, FILE* stream);

typedef int (*KstringOutputFunc)(kstring_t *s, const char *fmt, ...);
typedef int (*BufferOutputFunc)(char* s, const char* fmt, ...);
typedef int (*StreamOutputFunc)(FILE* s, const char* fmt, ...);

#define SYSTEM(cmd) do { if (system(cmd) < 0) OC_ERROR("Failed to running '%s'", cmd); } while(0)

#define OC_MIN(a, b) ((a) < (b) ? (a) : (b))
#define OC_MAX(a, b) ((a) > (b) ? (a) : (b))

#define HBN_MAX_PATH_LEN 2048

#endif // ONTCNS_AUX_H
