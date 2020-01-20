#include "oc_assert.h"

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

void
assertion_fail_handler(const char* expr, const char* file, const int line, const char* fmt, ...)
{
	fprintf(stderr, "Assertion Fail At '%s: %d'\n", file, line);
	fprintf(stderr, "Expression: '%s'\n", expr);
	if (!fmt) return;
	fprintf(stderr, "Context Information: '");
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fprintf(stderr, "'\n");
	exit(1);
}
