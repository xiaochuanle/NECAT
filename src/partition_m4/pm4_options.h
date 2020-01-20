#ifndef PM4_OPTIONS_H
#define PM4_OPTIONS_H

typedef struct {
	int batch_size;
	int num_output_files;
} PM4Options;

int
parse_PM4Options(int argc, char* argv[], PM4Options* options);

void
print_PM4Options(const PM4Options* options);

void
describe_PM4Options();

#endif // PM4_OPTIONS_H
