#ifndef OC_ASSERT_H
#define OC_ASSERT_H

	  
void
assertion_fail_handler(const char* expr, const char* file, const int line, const char* fmt, ...);

#define __oc_assert(expr, ...) \
	do { \
		if (!(expr)) { \
			assertion_fail_handler(#expr, __VA_ARGS__, NULL); \
			exit(1); \
		} \
	} while(0)
   
#define oc_assert(expr, args...) __oc_assert(expr, __FILE__, __LINE__, ##args);

#endif // OC_ASSERT_H
