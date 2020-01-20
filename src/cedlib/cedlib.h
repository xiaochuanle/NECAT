#ifndef CEDLIB_H
#define CEDLIB_H

#include "cedlib_defs.h"

void
cEdlibAlign(const char* const query, 
			const int queryLength,
			const char* const target, 
			const int targetLength,
			const EdlibAlignConfig config,
			EdlibData* edlibData1,
			EdlibData* edlibData2,
			void* km,
			unsigned char** alignment,
			int* alignmentLength); 

#endif // CEDLIB_H
