#ifndef PCAN_OPTIONS_H
#define PCAN_OPTIONS_H

typedef struct {
	int batch_size;
	int num_output_files;
	int num_threads;
} PcanOptions;

int
parse_PcanOptions(int argc, char* argv[], PcanOptions* options);

void
print_PcanOptions(const PcanOptions* options);

void
describe_PcanOptions();

#endif // PCAN_OPTIONS_H
