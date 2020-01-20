#include "asm_pm_common.c"

static void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s -Pwork_folder -Tnum_threads -Sstart_vid -Eend_vid\n", prog);
}

static int
parse_params(int argc, char* argv[], char** work_folder, int *nthreads, int *svid, int *evid)
{
	int set_wf = 0, set_nt = 0, set_svid = 0, set_evid = 0;
	for (int i = 1; i < argc; ++i) {
		int c = argv[i][0];
		size_t n = strlen(argv[i]);
		if (c != '-' || n < 3) {
			OC_LOG("Invalid option: %s", argv[i]);
			return ARG_PARSE_FAIL;
		}
		c = argv[i][1];
		switch (c) {
			case 'P':
				*work_folder = argv[i] + 2;
				set_wf = 1;
				break;
			case 'T':
				*nthreads = atoi(argv[i] + 2);
				set_nt = 1;
				break;
			case 'S':
				*svid = atoi(argv[i] + 2);
				set_svid = 1;
				break;
			case 'E':
				*evid = atoi(argv[i] + 2);
				set_evid = 1;
				break;
			default:
				OC_LOG("Invalid option: %s", argv[i]);
				return ARG_PARSE_FAIL;
				break;
		}
	}
	if (!set_wf) {
		OC_LOG("Work Folder (-P) is not set");
		return ARG_PARSE_FAIL;
	}
	if (!set_nt) {
		OC_LOG("CPU Threads (-T) is not set");
		return ARG_PARSE_FAIL;
	}
	if (!set_svid) {
		OC_LOG("Start Vid (-S) is not set");
		return ARG_PARSE_FAIL;
	}
	if (!set_evid) {
		OC_LOG("End Vid (-E) is not set");
		return ARG_PARSE_FAIL;
	}
	return ARG_PARSE_SUCCESS;
}

void fileidconvert(char *fileid,int noid){
	int i,j;
	char str1[50],str2[50];
	sprintf(str1,"%d",noid);
	j=strlen(str1);
	for(i=0;i<6-j;i++)str2[i]='0';
	str2[i]='\0';
	strcat(str2,str1);
	strcpy(fileid,str2);
}

int main(int argc, char* argv[])
{
	char* work_folder = 0;
	int nthreads = 0, svid = 0, evid = 0;
	if (parse_params(argc, argv, &work_folder, &nthreads, &svid, &evid) != ARG_PARSE_SUCCESS) {
		print_usage(argv[0]);
		return 1;
	}
	
	char tempstr[2048], tempstr1[1000];
	FILE* in;
	int* vrsid = (int*)malloc(sizeof(int) * evid);
	int* vreid = (int*)malloc(sizeof(int) * evid);
	sprintf(tempstr, "%s/ovlprep", work_folder);
	FOPEN(in, tempstr, "r");
	for (int i = 0; i < evid; ++i) {
		SAFE_SCANF(fscanf, in, 6, " %s %s %s %d %s %d\n", tempstr, tempstr, tempstr, vrsid + i, tempstr, vreid + i);
	}
	FCLOSE(in);
	
	fileidconvert(tempstr1, svid);
	sprintf(tempstr, "%s/%s.fasta", work_folder, tempstr1);
	reference = new_PackedDB();
OC_LOG("reference_path: %s", tempstr);
	pdb_load(reference, tempstr, TECH_NANOPORE);
OC_LOG("number of subjects: %d", PDB_NUM_SEQS(reference));
OC_LOG("number of bps: %d", PDB_SIZE(reference));
	lktbl = build_lookup_table(reference, seed_len, 1000, nthreads);
	reference_start_id = vrsid[svid - 1] - 1;
	
	sprintf(tempstr,"%s/%d_%d.r", work_folder, svid, evid);
	FOPEN(out, tempstr, "w");
	
	char job_name[1024];
	pthread_t jobs[nthreads];
	alloc_m4_sink();
	for (int i = svid; i <= evid; ++i) {
		sprintf(job_name, "pairwise mapping v%d vs v%d", svid, i);
		TIMING_START(job_name);
		fileidconvert(tempstr1, i);
		sprintf(tempstr, "%s/%s.fasta", work_folder, tempstr1);
		reads = new_PackedDB();
		pdb_load(reads, tempstr, TECH_NANOPORE);
		read_start_id = vrsid[i - 1] - 1;
		runnumber=0;
		readcount = PDB_NUM_SEQS(reads);
		terminalnum = (readcount + PLL - 1) / PLL;
		
		for (int k = 0; k < nthreads; ++k) {
			pthread_create(jobs + k, NULL, pairwise_mapping, NULL);
		}
		for (int k = 0; k < nthreads; ++k) {
			pthread_join(jobs[k], NULL);
		}
		dump_m4_sink();
		reads = free_PackedDB(reads);
		TIMING_END(job_name);
	}
	
	free_m4_sink();
	FCLOSE(out);
	lktbl = destroy_lookup_table(lktbl);
	reference = free_PackedDB(reference);
}
