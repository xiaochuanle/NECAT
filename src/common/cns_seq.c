#include "cns_seq.h"

#include "ontcns_defs.h"

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>

const char* ontsa_hdr_prefix = "ontsa_id";

BOOL 
is_ontsa_hdr(const char* hdr)
{
	const size_t N = strlen(ontsa_hdr_prefix);
	size_t n = strlen(hdr);
	if (n < N) return FALSE;
	for (size_t i = 0; i < N; ++i) {
		if (ontsa_hdr_prefix[i] != hdr[i]) return FALSE;
	}
	return TRUE;
}

int
extract_ontsa_id(const char* hdr)
{
	int id = 0;
	const char* p = hdr;
	while(p != NULL && (!isdigit(*p)))++p;
	assert(p != NULL);
	while(p != NULL && isdigit(*p)) {
		id = id * 10 + *p - '0';
		++p;
	}
	return id;
}
