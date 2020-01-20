#include "../common/map_options.h"
#include "../common/makedb_aux.h"
#include "../klib/kstring.h"
#include "../common/ontcns_aux.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s [options] wrk-dir output\n", prog);
	fprintf(stderr, "\n");
	fprintf(stderr, "OPTIONS AND DESCRIPTIONS:\n");
	describe_MapOptions(&sDefaultPairwiseMapingOptions);
}

void
make_volume_result_name(const char* wrk_dir, const int vid, kstring_t* name)
{
	kstr_clear(*name);
	size_t n = strlen(wrk_dir);
	kputs(wrk_dir, name);
	if (wrk_dir[n - 1] != '/') kputc('/', name);
	ksprintf(name, "pm_result_%d", vid);
	kputc('\0', name);
}

void
make_pm_command(MapOptions* options, const char* wrk_dir, const int vid, const char* output, kstring_t* cmd)
{
	kstr_clear(*cmd);
	kputs("oc2pmov ", cmd);
	MapOptions2String(options, cmd);
	kputc(' ', cmd);
	ksprintf(cmd, "%s ", wrk_dir);
	ksprintf(cmd, "%d ", vid);
	ksprintf(cmd, "%s", output);
	kputc('\0', cmd);
}

void
make_merge_results_command(const int vid, const char* volume_results_name, const char* output, kstring_t* cmd)
{
	kstr_clear(*cmd);
	if (vid == 0) {
		ksprintf(cmd, "cat %s > %s", volume_results_name, output);
	} else {
		ksprintf(cmd, "cat %s >> %s", volume_results_name, output);
	}
	kputc('\0', cmd);
}

static BOOL
job_is_finised(const char* wrk_dir, const int vid)
{
	char name[2048];
	sprintf(name, "%s/pm%d.finished", wrk_dir, vid);
	return !access(name, F_OK);
}

static void
set_job_finished(const char* wrk_dir, const int vid)
{
	char name[2048];
	sprintf(name, "%s/pm%d.finished", wrk_dir, vid);
	FILE* file = fopen(name, "w");
	fclose(file);
}

static void
rm_volume_result(const char* result_path)
{
	char name[2048];
	sprintf(name, "rm -f %s", result_path);
	SYSTEM(name);
}

int main(int argc, char* argv[])
{
	TIMING_START(__func__);
	
	MapOptions options = sDefaultPairwiseMapingOptions;
	if (argc < 3 || (parse_MapOptions(argc - 2, argv, &options) != ARG_PARSE_SUCCESS)) {
		print_usage(argv[0]);
		return 1;
	}
	
	const char* wrk_dir = argv[argc - 2];
	const char* output = argv[argc - 1];
	int num_volumes = load_num_volumes(wrk_dir);
	
	new_kstring(cmd);
	new_kstring(result_name);
	for (int i = 0; i < num_volumes; ++i) {
		if (job_is_finised(wrk_dir, i)) continue;
		make_volume_result_name(wrk_dir, i, &result_name);
		make_pm_command(&options, wrk_dir, i, kstr_str(result_name), &cmd);
		const char* cmd_str = kstr_str(cmd);
		fprintf(stdout, "Running command '%s'\n", cmd_str);
		SYSTEM(cmd_str);
		set_job_finished(wrk_dir, i);
	}
	
	for (int i = 0; i < num_volumes; ++i) {
		make_volume_result_name(wrk_dir, i, &result_name);
		make_merge_results_command(i, kstr_str(result_name), output, &cmd);
		const char* cmd_str = kstr_str(cmd);
		SYSTEM(cmd_str);
		rm_volume_result(kstr_str(result_name));
	}
	
	free_kstring(cmd);
	free_kstring(result_name);
	TIMING_END(__func__);
}
