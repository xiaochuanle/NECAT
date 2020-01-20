#include "ontcns_aux.h"

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/stat.h>

#include "../klib/kseq.h"
KSEQ_INIT2(, gzFile, err_gzread)

int fileno(FILE*);

void
what_is_time_now(char now[])
{
    time_t ltime;
    time(&ltime);
    char* ltime_str = ctime(&ltime);
    size_t n = strlen(ltime_str) - 1;
    strncpy(now, ltime_str, n);
    now[n] = '\0';
}

void
oc_dump_info(FILE* out,
        const char* level,
        const char* file,
        const int line,
        const char* fmt,
        ...)
{
    char now[256];
    what_is_time_now(now);
    fprintf(out, "[%s] %s: ", now, level);
	va_list ap;
	va_start(ap, fmt);
	vfprintf(out, fmt, ap);
	va_end(ap);
    if (file) fprintf(out, " (%s, %d)", file, line);
	fprintf(out, "\n");
}

void
oc_dump_info_simple(FILE* out, const char* fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vfprintf(out, fmt, ap);
	va_end(ap);
	fprintf(out, "\n");
}

int safe_gzread(const char* src_file,
			   const int src_line,
			   gzFile file, 
			   void *ptr, 
			   unsigned int len)
{
	int ret = gzread(file, ptr, len);

	if (ret < 0)
	{
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		const char* why = Z_ERRNO == errnum ? strerror(errno) : msg;
		OC_ERROR_FULL(src_file, src_line, "gzread error: %s", why);
	}

	return ret;
}

int 
err_gzread(gzFile file, void* ptr, unsigned int len)
{
	return safe_gzread(__FILE__, __LINE__, file, ptr, len);
}

gzFile 
safe_gzopen(const char* src_file,
				const int src_line,
				const char *fn, 
				const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0) {
		fp = gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
		/* According to zlib.h, this is the only reason gzdopen can fail */
		if (!fp) OC_ERROR_FULL(src_file, src_line, "failed to open file '%s' with mode '%s': Out of memory", fn, mode);
		return fp;
	}
	if ((fp = gzopen(fn, mode)) == 0) {
		const char* why = errno ? strerror(errno) : "Out of memory";
		OC_ERROR_FULL(src_file, src_line, "failed to open file '%s'with mode '%s': %s", fn, mode, why);
	}
	return fp;
}

int 
safe_gzclose(const char* src_file,
				const int src_line,
				gzFile file)
{
	int ret = gzclose(file);
	if (Z_OK != ret)
	{
		const char* why = Z_ERRNO == ret ? strerror(errno) : zError(ret);
		OC_ERROR_FULL(src_file, src_line, "gzclose error: %s", why);
	}

	return ret;
}

FILE*
safe_fopen(const char* src_file,
			 const int src_line,
			 const char* path,
			 const char* mode)
{
	FILE *fp = 0;
	if (strcmp(path, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(path, mode)) == 0) {
		const char* why = strerror(errno);
		OC_ERROR_FULL(src_file, src_line, "failed to open file '%s' with mode '%s': %s", path, mode, why);
	}
	return fp;
}

size_t
safe_fwrite(const char* src_file,
			const int src_line,
			const void* ptr,
			size_t size,
			size_t nmemb,
			FILE* stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb) {
		OC_ERROR_FULL(src_file, src_line, "fwrite error: %s", strerror(errno));
	}
	return ret;
}

size_t
safe_fread(const char* src_file,
		   const int src_line,
		   void* ptr,
		   size_t size,
		   size_t nmemb,
		   FILE* stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb) {
		const char* why = ferror(stream) ? strerror(errno) : "Unexpected end of file";
		OC_ERROR_FULL(src_file, src_line, "fread error: %s", why);
	}
	return ret;
}

int 
safe_fclose(const char* src_file,
			const int src_line,
			FILE* stream)
{
	int ret = fclose(stream);
	if (ret != 0) {
		OC_ERROR_FULL(src_file, src_line, "fclose error: %s", strerror(errno));
	}
	return ret;
}

idx
file_size(const char* src_path, const int src_line, const char* path)
{
	struct stat statbuf;
	if (stat(path, &statbuf) == -1) {
		OC_ERROR_FULL(src_path, src_line, "failed to stat file '%s'", path);
	}
	return statbuf.st_size;
}

double 
time_diff(const struct timeval* start, const struct timeval* end)
{
	double d = end->tv_sec - start->tv_sec;
	d += 1.0 * (end->tv_usec - start->tv_usec) / 1e6;
	return d;
}

char*
ontcns_getline(char* buf, int bufsize, FILE* stream)
{
	if (!fgets(buf, bufsize, stream)) return 0;
	size_t n = strlen(buf);
	if (buf[n - 1] == '\n') buf[n - 1] = '\0';
	return buf;
}
