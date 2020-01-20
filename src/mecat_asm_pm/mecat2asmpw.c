/* 
 * File:   pairwise_pthread.c
 * Author: Chuan-Le Xiao 
 * Email: xiaochuanle@126.com
 *
 * Created on September 1, 2015, 11:28 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
//#include <Windows.h>
#define RM 100000
#define ZV 1000
#define DN 500
#define BC 10
#define SM 60
#define SI 61
#define MAXC 100
#define MAXSTR 1000000000
#define SVM 200000
#define PLL 500
#define ErrorRate 0.10
typedef struct{
    int readno,readlen;
    char *seqloc;
} ReadFasta;
typedef struct{
	int readno,length;
    char *mem;
}readmemory;

typedef struct{
	int curnum,maxnum;
    int *pre;
}seedindex;

typedef struct {    
    int aln_str_size,dist,aln_q_s,aln_q_e,aln_t_s,aln_t_e;
    char q_aln_str[2500];
    char t_aln_str[2500];
} alignment;

typedef struct {
    int pre_k,x1,y1,x2,y2;
} d_path_data;

typedef struct {
    int d,k,pre_k,x1,y1,x2,y2;
} d_path_data2;
typedef struct {
    int x,y;
} path_point;

typedef struct{
	int loc1,loc2,left1,left2,right1,right2,score,num1,num2,readno,readstart;
	char chain;
}canidate_save;

typedef struct{
 char left_store1[RM],left_store2[RM],right_store1[RM],right_store2[RM],out_store1[RM],out_store2[RM];
// char tempstr1[RM],tempstr2[RM];
}output_store;

struct Back_List{
	short int score,loczhi[SM],seedno[SM],seednum;
	int index;
};

//whole var.
pthread_t *thread; //???????????????
int threadnum;
FILE **outfile;
pthread_mutex_t mutilock; //????????
int runnumber=0,runthreadnum=0, readcount,terminalnum;
int *countin,**databaseindex,*allloc,sumcount;
int seed_len,*llocation,seqcount,curreadcount;
char *STRMEM;
readmemory *indexread;
ReadFasta *readinfo;
char *savework,workpath[300];

int compare_d_path(const void * a, const void * b)
{
    const d_path_data2 * arg1 = (d_path_data2 *)a;
    const d_path_data2 * arg2 = (d_path_data2 *)b;
    if (arg1->d - arg2->d == 0) {
        return  arg1->k - arg2->k;
    } else {
        return arg1->d - arg2->d;
    }
}


d_path_data2 * get_dpath_idx( int d, int k, unsigned long max_idx, d_path_data2 * base) {
    d_path_data2 d_tmp;
    d_path_data2 *rtn;
    d_tmp.d = d;
    d_tmp.k = k;
    rtn = (d_path_data2 *)  bsearch( &d_tmp, base, max_idx, sizeof(d_path_data2), compare_d_path);
    //printf("dp %ld %ld %ld %ld %ld %ld %ld\n", (rtn)->d, (rtn)->k, (rtn)->x1, (rtn)->y1, (rtn)->x2, (rtn)->y2, (rtn)->pre_k); 
    return rtn;
}

int align(char * query_seq, char * target_seq,int band_tolerance,int get_aln_str,alignment * align_rtn,int * V,int * U,d_path_data2 *d_path,path_point * aln_path) {
    int k_offset,d,k, k2,best_m,min_k, new_min_k,max_k, new_max_k,pre_k,x, y;
    int ck,cd,cx, cy, nx, ny,max_d,band_size,q_len,t_len;
    unsigned long d_path_idx = 0,max_idx = 0;
    int aln_path_idx,aln_pos,i,aligned=0;
	d_path_data2 * d_path_aux;
	q_len=strlen(query_seq);t_len=strlen(target_seq);
    max_d = (int) (ErrorRate*(q_len + t_len));
    band_size = band_tolerance * 2;
    k_offset = max_d;
    align_rtn->aln_str_size = 0;align_rtn->aln_q_s = 0;align_rtn->aln_q_e = 0;align_rtn->aln_t_s = 0;align_rtn->aln_t_e = 0;
    best_m = -1;min_k = 0; max_k = 0;d_path_idx = 0; max_idx = 0;
    for (d = 0; d < max_d; d ++ ) {
        if(max_k - min_k > band_size)break;
        for (k = min_k; k <= max_k;  k += 2) {
            if ( k == min_k || k != max_k && V[ k - 1 + k_offset ] < V[ k + 1 + k_offset] ) {pre_k = k + 1;x = V[ k + 1 + k_offset];} 
			else {pre_k = k - 1;x = V[ k - 1 + k_offset] + 1;}
            y = x - k;
            d_path[d_path_idx].d = d;d_path[d_path_idx].k = k;d_path[d_path_idx].x1 = x;d_path[d_path_idx].y1 = y;

            while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){x++;y++;}
            d_path[d_path_idx].x2 = x;d_path[d_path_idx].y2 = y;d_path[d_path_idx].pre_k = pre_k;d_path_idx ++;
            V[ k + k_offset ] = x;U[ k + k_offset ] = x + y;
            if( x + y > best_m) {best_m = x + y;}
            if ( x >= q_len || y >= t_len) {
                aligned = 1;
                max_idx = d_path_idx;
                break;
            }
        }
        
        // For banding
        new_min_k = max_k;
        new_max_k = min_k;
        for (k2 = min_k; k2 <= max_k;  k2 += 2){
            if (U[ k2 + k_offset] >= best_m - band_tolerance ) {
                if( k2 < new_min_k )new_min_k = k2;
                if( k2 > new_max_k )new_max_k = k2;
            }
        }
        max_k = new_max_k + 1;
        min_k = new_min_k - 1;
         
        if (aligned == 1) {
            align_rtn->aln_q_e = x;align_rtn->aln_t_e = y;
            align_rtn->dist = d;align_rtn->aln_str_size = (x + y + d) / 2;
            align_rtn->aln_q_s = 0;align_rtn->aln_t_s = 0;
			//qsort(d_path, max_idx, sizeof(d_path_data2), compare_d_path);
            if (get_aln_str > 0) {
                cd = d;ck = k;aln_path_idx = 0;
                while (cd >= 0 && aln_path_idx < q_len + t_len + 1) {    
                    d_path_aux = (d_path_data2 *) get_dpath_idx( cd, ck, max_idx, d_path);
                    aln_path[aln_path_idx].x = d_path_aux -> x2;aln_path[aln_path_idx].y = d_path_aux -> y2;aln_path_idx ++;
                    aln_path[aln_path_idx].x = d_path_aux -> x1;aln_path[aln_path_idx].y = d_path_aux -> y1;aln_path_idx ++;
                    ck = d_path_aux -> pre_k;cd -= 1;
                }
                aln_path_idx --;
                cx = aln_path[aln_path_idx].x;cy = aln_path[aln_path_idx].y;align_rtn->aln_q_s = cx;align_rtn->aln_t_s = cy;
                aln_pos = 0;
                while ( aln_path_idx > 0 ) {
                    aln_path_idx --;
                    nx = aln_path[aln_path_idx].x;ny = aln_path[aln_path_idx].y;
                    if (cx == nx && cy == ny){
                        continue;
                    }
                    if (nx == cx && ny != cy){ //advance in y
                        for (i = 0; i <  ny - cy; i++)align_rtn->q_aln_str[aln_pos + i] = '-';
                       for (i = 0; i <  ny - cy; i++)align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                        aln_pos += ny - cy;} 
					else if (nx != cx && ny == cy){ //advance in x
                        for (i = 0; i <  nx - cx; i++)align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                        for (i = 0; i <  nx - cx; i++)align_rtn->t_aln_str[aln_pos + i] = '-';
                        aln_pos += nx - cx;
                    } 
					else {
                        for (i = 0; i <  nx - cx; i++)align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                        for (i = 0; i <  ny - cy; i++)align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                        aln_pos += ny - cy;
                    }
                    cx = nx;cy = ny;
                }
                align_rtn->aln_str_size = aln_pos;
            }
            break;
        }
    }
if(align_rtn->aln_q_e==q_len||align_rtn->aln_t_e==t_len)return(1);
else return(0);
}

void string_check(char *seq1,char *seq2,char *str1,char *str2){
	//char seq1[500]="GACCGCCGGACAGCCCACAAACACAACAGCATTTGGCGTATTTCCCGTCAAAGGACTGCGAGTGGGACCGCGCACCGATTTATAGAGTAACGGTGGGACTTACCCCCGACGACTAGAGG";
	//char seq2[500]="GACCGCCGGACAGGCCACAAACACAAATCCGCGAGGCGTATTCCTGTCAAAGGGACTACGGCCCAGTGGGACGCCGCACGACTATATAGTAGTAAGGTGTGCTTTACCCGACGCCCTAGAGG";
	//char str1[700]="GACCGCCGGACAG-CCCACAAACACAACA---GC-ATTTGGCGTATTTCCC-GTCAAAGG-ACTG-CG----AGTGGGACCGC-GCACCGA-T-TTATAG-AGTAACGGTG-GGACTT-ACCCCCGACGAC--TAGAGG";
	//char str2[700]="GACCGCCGGACAGGCC-ACAAACACAA-ATCCGCGA---GGCGTATT-CC-TGTCAAAGGGACT-ACGGCCCAGTGGGAC-GCCGCAC-GACTAT-ATAGTAGTAA-GGTGTG--CTTTACCC--GACG-CCCTAGAGG";
	int len1,len2,slen1,slen2,loc1,loc2,eit,k,s,j,ii,ii2,slen11;
	char one1[RM];
	len1=strlen(seq1)-1;len2=strlen(seq2)-1;
	slen11=strlen(str1)-1;
	for(slen1=strlen(str1)-1,loc1=0,loc2=0,eit=0;slen1>-1;slen1--){
					if(str1[slen1]!='-')loc1++; 
					else if(loc1<=len1&&loc2<=len2&&seq1[len1-loc1]==seq2[len2-loc2]){
						k=1;
						while(loc1+k<=len1&&loc2+k<=len2&&seq1[len1-loc1-k]==seq2[len2-loc2-k])k++;
			            s=0;j=slen1;
						while(s<k){
							if(str1[j]!='-'){str1[j]='-';s++;}
							j--;
						}
						s=0;j=slen1;
						while(s<k){
							if(str2[j]!='-'){str2[j]='-';s++;}
							j--;
						}
						for(s=0,j=slen1;s<k;j--){str1[j]=seq1[len1-loc1-s];str2[j]=seq2[len2-loc2-s];s++;}
						if(str1[slen1]!='-')loc1++;
					}
					/*for(ii=0,ii2=0;ii<=slen11;ii++)if(str1[ii]!='-')one1[ii2++]=str1[ii];
					one1[ii2]='\0';
					if(strcmp(one1,seq1)!=0){
						printf("xiao");
					}
					if(slen1==2246){
						printf("xiao");
					}*/
					if(str2[slen1]!='-')loc2++;
					else if(str1[slen1]!='-'&&(loc1-1<=len1&&loc2<=len2&&seq1[len1-loc1+1]==seq2[len2-loc2])){
						k=1;
						while(loc1+k-1<=len1&&loc2+k<=len2&&seq1[len1-loc1+1-k]==seq2[len2-loc2-k])k++;
			            s=0;j=slen1;
						while(s<k){
							if(str1[j]!='-'){str1[j]='-';s++;}
							j--;
						}
						s=0;j=slen1;
						while(s<k){
							if(str2[j]!='-'){str2[j]='-';s++;}
							j--;
						}
						for(s=0,j=slen1;s<k;j--){
							str1[j]=seq1[len1-loc1+1-s];
							str2[j]=seq2[len2-loc2-s];s++;}
						if(str2[slen1]!='-')loc2++;
					}
					else if(str1[slen1]=='-'&&(loc1-1<=len1&&loc2<=len2&&seq1[len1-loc1]==seq2[len2-loc2])){
						k=1;
						while(loc1+k<=len1&&loc2+k<=len2&&seq1[len1-loc1-k]==seq2[len2-loc2-k])k++;
			            s=0;j=slen1;
						while(s<k){
							if(str1[j]!='-'){str1[j]='-';s++;}
							j--;
						}
						s=0;j=slen1;
						while(s<k){
							if(str2[j]!='-'){str2[j]='-';s++;}
							j--;
						}
						for(s=0,j=slen1;s<k;j--){
							str1[j]=seq1[len1-loc1-s];
							str2[j]=seq2[len2-loc2-s];s++;}
						   if(str2[slen1]!='-')loc2++;

					}
				/*for(ii=0,ii2=0;ii<=slen11;ii++)if(str1[ii]!='-')one1[ii2++]=str1[ii];
					one1[ii2]='\0';
					if(strcmp(one1,seq1)!=0){
						printf("xiao");
					}*/

				}
}

int binary( int *a, int key, int n )
{
int left = 0, right = n-1, mid;
mid = ( left + right ) / 2;
if(a[right]<key) return right;
while( left <= right && a[mid] != key )
{
if( a[mid] < key &&a[mid+1]>key) return mid;
if( a[mid] < key &&a[mid+1]==key) return mid+1;
else if( a[mid] < key &&a[mid+1]<key) left = mid+1;
else if( a[mid] > key && a[mid-1]<=key ) return mid-1;
else if( a[mid] > key && a[mid-1]>key ) right = mid-1;
mid = ( left + right ) / 2;
}
if( a[mid] == key )return mid;
}

unsigned short atcttrans(char c){
if(c=='A'||c=='a')return 0;
else if(c=='T'||c=='T')return 1;
else if(c=='C'||c=='c')return 2;
else if(c=='G'||c=='g')return 3;
else return 4;
}


int sumvalue_x(int *intarry,int count){
	int i,sumval=0;
	for(i=0;i<count;i++){
            if(intarry[i]>0&&intarry[i]<257)sumval=sumval+intarry[i];
            else if(intarry[i]>256)intarry[i]=0;
        }
	return(sumval);
}

int transnum_buchang(char *seqm,int *value,int *endn,int len_str,int readnum){
   int eit=0,temp;
	int i,j,start,num;
	 num=(len_str-readnum)/BC+1;
	 *endn=(len_str-readnum)%BC;
	 //if((len_str%readnum)>0)num=num+1;
	 for(i=0;i<num;i++){
		 eit=0;start=i*BC;
		 //if(i==num-1){eit=0;start=len_str-readnum;}
		 for(j=0;j<readnum;j++){
            temp=atcttrans(seqm[start+j]);
			if(temp==4){eit=-1;break;}
			 eit=eit<<2;
			 eit=eit+temp;
		 }
		 value[i]=eit;
	 }
	 return(num);
}



int str2num(char *str){
	int sum=0,i;
	char ch;
	i=0;ch=str[i];
	while(ch!='\0'){sum=sum*10+(ch-'0');i++;ch=str[i];}
	return(sum);
}

void insert_loc(struct Back_List *spr,int loc,int seedn,float len){
	int list_loc[SI],list_score[SI],list_seed[SI],i,j,minval,mini;
	for(i=0;i<SM;i++){list_loc[i]=spr->loczhi[i];list_seed[i]=spr->seedno[i];list_score[i]=0;}
	list_loc[SM]=loc;list_seed[SM]=seedn;list_score[SM]=0;
	mini=-1;
	minval=10000;
	for(i=0;i<SM;i++)for(j=i+1;j<SI;j++)if(list_seed[j]-list_seed[i]>0&&list_loc[j]-list_loc[i]>0&&fabs((list_loc[j]-list_loc[i])/((list_seed[j]-list_seed[i])*len)-1.0)<ErrorRate){list_score[i]++;list_score[j]++;}
	for(i=0;i<SI;i++)if(minval>list_score[i]){minval=list_score[i];mini=i;}
	if(minval==SM){spr->loczhi[SM-1]=loc;spr->seedno[SM-1]=seedn;}
	else if(minval<SM&&mini<SM){
		for(i=mini;i<SM;i++){spr->loczhi[i]=list_loc[i+1];spr->seedno[i]=list_seed[i+1];}
		spr->score--;
	}
}
int find_location(int *t_loc,int *t_seedn,int *t_score,int *loc,int k,int *rep_loc,float len,int read_len1){
	int i,j,maxval=0,maxi,rep=0,lasti,tempi;
	for(i=0;i<k;i++)t_score[i]=0;
	for(i=0;i<k-1;i++)for(j=i+1,tempi=t_seedn[i];j<k;j++)if(tempi!=t_seedn[j]&&t_seedn[j]-t_seedn[i]>0&&t_loc[j]-t_loc[i]>0&&t_loc[j]-t_loc[i]<read_len1&&fabs((t_loc[j]-t_loc[i])/((t_seedn[j]-t_seedn[i])*len)-1)<0.10){t_score[i]++;t_score[j]++;tempi=t_seedn[j];}
	
	for(i=0;i<k;i++){
		if(maxval<t_score[i]){maxval=t_score[i];maxi=i;rep=0;}
		else if(maxval==t_score[i]){rep++;lasti=i;}
	}
	for(i=0;i<4;i++)loc[i]=0;
	if(maxval>=5&&rep==maxval){loc[0]=t_loc[maxi],loc[1]=t_seedn[maxi];*rep_loc=maxi;loc[2]=t_loc[lasti],loc[3]=t_seedn[lasti];return(1);}
	else if(maxval>=5&&rep!=maxval){
		for(j=0;j<maxi;j++)if(t_seedn[maxi]-t_seedn[j]>0&&t_loc[maxi]-t_loc[j]>0&&t_loc[maxi]-t_loc[j]<read_len1&&fabs((t_loc[maxi]-t_loc[j])/((t_seedn[maxi]-t_seedn[j])*len)-1)<0.10){
				if(loc[0]==0){loc[0]=t_loc[j];loc[1]=t_seedn[j];*rep_loc=j;}
				else { loc[2]=t_loc[j];loc[3]=t_seedn[j];}
		}
		j=maxi;
		if(loc[0]==0){loc[0]=t_loc[j];loc[1]=t_seedn[j];*rep_loc=j;}
		else { loc[2]=t_loc[j];loc[3]=t_seedn[j];}
		for(j=maxi+1;j<k;j++)if(t_seedn[j]-t_seedn[maxi]>0&&t_loc[j]-t_loc[maxi]>0&&t_loc[j]-t_loc[maxi]<=read_len1&&fabs((t_loc[j]-t_loc[maxi])/((t_seedn[j]-t_seedn[maxi])*len)-1)<0.10){
				if(loc[0]==0){loc[0]=t_loc[j];loc[1]=t_seedn[j];*rep_loc=j;}
				else { loc[2]=t_loc[j];loc[3]=t_seedn[j];}
				}
		return(1);
		}
	else return(0);
}

int load_read(int readcount,char *STRMEM,char *filename,int *ll,int qstart){
	int readlen,readno,count=0,sum=0,i;
	char onedata[RM],*pre,tempstr[300];
	FILE *fq;
	fq=fopen(filename,"r");
	count=0;
	pre=STRMEM;
	while(fscanf(fq,">%[^\n]s",tempstr)!=EOF&&fscanf(fq,"%s\n",onedata)!=EOF){
		readno=qstart+count;
		readlen=strlen(onedata);
       for(i=0;i<readlen;i++)if(onedata[i]>='a')onedata[i]=toupper(onedata[i]);
		ll[count]=sum;
		strcpy(pre,onedata);
		indexread[count].mem=pre;
		indexread[count].readno=readno;
		indexread[count].length=readlen;
		sum=sum+readlen+1;
		pre=pre+readlen+1;count++;
	}
	fclose(fq);
	return(sum);
}
int compare_readindex(const void * a, const void * b)
{
    const readmemory * arg1 = (readmemory *)a;
    const readmemory * arg2 = (readmemory *)b;
    if (arg1->length - arg2->length == 0) {
        return  arg1->readno - arg2->readno;
    } else {
        return arg2->length-arg1->length;
    }
}


void creat_ref_index(char *seq,int seqcount){
     unsigned int eit,temp;
     int i, start, indexcount=0,leftnum=0;
     //
    if(seed_len==14)indexcount=268435456;
    else if(seed_len==13)indexcount=67108864;
    else if(seed_len==12)indexcount=16777216;
    else if(seed_len==11)indexcount=4194304;
    else if(seed_len==10)indexcount=1048576;
    else if(seed_len==9)indexcount=262144;
    else if(seed_len==8)indexcount=65536;
    else if(seed_len==7)indexcount=16384;
    else if(seed_len==6)indexcount=4096;
    leftnum=34-2*seed_len;  
 //printf("Constructing look-up table...\n");
    leftnum=34-2*seed_len;  
 //printf("Constructing look-up table...\n");
    countin=(int *)malloc((indexcount)*sizeof(int));
    for(i=0;i<indexcount;i++)countin[i]=0;
    
// Count the number
  eit=0;start=0;
 for(i=0;i<seqcount;i++){
	if(seq[i]=='N'||(temp=atcttrans(seq[i]))==4){
	     eit=0;
	     start=0;
	     continue;
	  }
   temp=atcttrans(seq[i]);
   if(start<seed_len-1){eit=eit<<2;
	           eit=eit+temp;
               start=start+1;
   }
   else if(start>=seed_len-1){
	 eit=eit<<2;
	 eit=eit+temp;
	 start=start+1;
	 countin[eit]=countin[eit]+1;
	 	 eit=eit<<leftnum;
	     eit=eit>>leftnum;
   }
 }


 //Max_index
sumcount=sumvalue_x(countin,indexcount);
allloc=(int *)malloc(sumcount*sizeof(int));
databaseindex=(int **)malloc((indexcount)*sizeof(int *));
 //allocate memory
sumcount=0;
 for(i=0;i<indexcount;i++){
	 if(countin[i]>0){
		databaseindex[i]=allloc+sumcount;
		 sumcount=sumcount+countin[i];
	     countin[i]=0;                           
         }
	 else databaseindex[i]=NULL;
 }
   
	// printf("xiao");//10834098
	 
//constructing the look-up table
  eit=0;start=0;
for(i=0;i<seqcount;i++){
	if(seq[i]=='N'||(temp=atcttrans(seq[i]))==4){
	eit=0;start=0;
		continue;
	}
   temp=atcttrans(seq[i]);
   if(start<seed_len-1){	 
         eit=eit<<2;
	             eit=eit+temp;
               start=start+1;
   }
   else if(start>=seed_len-1){
	 eit=eit<<2;
	 eit=eit+temp;
	 start=start+1;

     if(databaseindex[eit]!=NULL){
	 countin[eit]=countin[eit]+1;
	 databaseindex[eit][countin[eit]-1]=i+2-seed_len;
	 }
         eit=eit<<leftnum;
	 eit=eit>>leftnum;
   }
}

}



void pairwise_mapping(int threadint){
 int cleave_num,read_len,missreal,s_k;
 int mvalue[50000],*leadarray,loc_flag,flag_end,u_k;
 int count1=0,i,j,k,read_name;
 struct Back_List *database,*temp_spr,*temp_spr1;
 int location_loc[4],repeat_loc,*index_list,*index_spr;
 short int *index_score,*index_ss;
 int temp_list[150],temp_seedn[150],temp_score[150],start_loc;
 int endnum,loc_count,missall,ii,eit;
 FILE *out,*fid;
 char tempstr[300],onedata1[RM],onedata2[RM],*onedata,seq1[2500],seq2[2500],*seq_pr1,*seq_pr2,FR,*seq; 
 int sci=0,loc,templong,localnum,read_i,read_end;
 int left_loc,left_loc1,right_loc=0,right_loc1,cc1,readno,fileid,canidatenum,loc_seed,loc_list;
 int length1,gg,num1,num2;
 int low,high,mid,seedcount,lread_count;
 canidate_save canidate_loc[MAXC],canidate_temp;
 readmemory tempread;
 alignment strvalue;
 int longstr1[2000],longstr2[2000],readstart1,readend1,left_length1,right_length1,left_length2,right_length2,align_flag;
 d_path_data2 *d_path;
 path_point aln_path[5000];
 output_store *resultstore,*resultstore1;
 float jscore;
 // m5 various
 int numMatch,numMismatch,numIns,numDel;
 char matchPattern[RM],strand1,strand2;
 //get 
 seq=STRMEM;
 lread_count=curreadcount;   
   
   //sprintf(tempstr,"%s/%d.srel",workpath,threadint);
  // fid=fopen(tempstr, "a");
   j=seqcount/ZV+5;
   index_list=(int *)malloc(j*sizeof(int));
   index_score=(short int *)malloc(j*sizeof(short int));
   database=(struct Back_List *)malloc(j*sizeof(struct Back_List));
   for(i=0,temp_spr=database;i<j;temp_spr++,i++){temp_spr->score=0;temp_spr->index=-1;}
   d_path=(d_path_data2*)malloc(600*700*2*sizeof(d_path_data2));
   resultstore=(output_store *)malloc(1*sizeof(output_store));
   resultstore1=(output_store *)malloc(1*sizeof(output_store));
   //t1=(float)clock();
   fileid=1;
  while(fileid){
        pthread_mutex_lock(&mutilock);
        localnum=runnumber;
        runnumber++;
        //printf("localnum=%d\n",localnum);
        pthread_mutex_unlock(&mutilock);
        if(localnum>=terminalnum){
            fileid=0;
            break;
        }
       if(localnum==terminalnum-1)read_end=readcount;
	   else read_end = (localnum + 1)*PLL;
     for(read_i=localnum*PLL;read_i<read_end;read_i++){
          read_name=readinfo[read_i].readno;
          read_len=readinfo[read_i].readlen;
          strcpy(onedata1,readinfo[read_i].seqloc);
         // printf("%d\n",read_name);
		  //if(read_name!=1548)continue;
          canidatenum=0;
	   for(ii=1;ii<=2;ii++){
		  if(ii==1)onedata=onedata1;
		  else if(ii==2){
			  strcpy(onedata2,onedata1);
			  onedata=onedata2;
			  for(j=read_len-1,i=0;j>i;j--,i++){FR=onedata[i];onedata[i]=onedata[j];onedata[j]=FR;}
			  for(i=0;i<read_len;i++){
				FR=onedata[i];
				switch(FR){
				case 'A':{onedata[i]='T';break;}
				case 'T':{onedata[i]='A';break;}
				case 'C':{onedata[i]='G';break;}
				case 'G':{onedata[i]='C';break;}
				}
			}
		  }
	  loc_flag=0;
	  endnum=0;
	  loc_count=0;
	  missall=0;
      read_len=strlen(onedata);
	  cleave_num=transnum_buchang(onedata,mvalue,&endnum,read_len,seed_len);
	  //printf("x1\n");
	  j=0;missreal=0;
	  index_spr=index_list;
	  index_ss=index_score;
          endnum=0;
            for(k=0; k<cleave_num; k++)if(mvalue[k]>=0)
                {
                    count1=countin[mvalue[k]];
                    leadarray=databaseindex[mvalue[k]];
                    for(i=0; i<count1; i++,leadarray++)
                    {
                        templong=(*leadarray)/ZV;
                        u_k=(*leadarray)%ZV;
                        if(templong>=0)
                        {
                            temp_spr=database+templong;
                            if(temp_spr->score==0||temp_spr->seednum<k+1)
                            {
                                loc=++(temp_spr->score);
                                if(loc<=SM)
                                {
                                    temp_spr->loczhi[loc-1]=u_k;
                                    temp_spr->seedno[loc-1]=k+1;
                                }
                                //else insert_loc(temp_spr,u_k,k+1,BC);
                                if(templong>0)s_k=temp_spr->score+(temp_spr-1)->score;
                                else s_k=temp_spr->score;
                                if(endnum<s_k)endnum=s_k;
                                if(temp_spr->index==-1)
                                {
                                    *(index_spr++)=templong;
                                    *(index_ss++)=s_k;
                                    temp_spr->index=j;
                                    j++;
                                }
                                else index_score[temp_spr->index]=s_k;
                            }
                            temp_spr->seednum=k+1;
                        }
                    }
                }
//get most mapping canidate location
         cc1=j;
         for(i=0,index_spr=index_list,index_ss=index_score;i<cc1;i++,index_spr++,index_ss++)if(*index_ss>10){
		                     temp_spr=database+*index_spr;
				     if(temp_spr->score==0)continue;
		                     s_k=temp_spr->score;
				     start_loc=(*index_spr)*ZV;
			          if(*index_spr>0){loc=(temp_spr-1)->score;if(loc>0)start_loc=(*index_spr-1)*ZV;}
				  else loc=0;
		              if(loc==0)for(j=0,u_k=0;j<s_k&&j<SM;j++){temp_list[u_k]=temp_spr->loczhi[j];temp_seedn[u_k]=temp_spr->seedno[j];u_k++;}
		              else {
							k=loc;
							u_k=0;
							temp_spr1=temp_spr-1;
							for(j=0;j<k&&j<SM;j++){temp_list[u_k]=temp_spr1->loczhi[j];temp_seedn[u_k]=temp_spr1->seedno[j];u_k++;}
							for(j=0;j<s_k&&j<SM;j++){temp_list[u_k]=temp_spr->loczhi[j]+ZV;temp_seedn[u_k]=temp_spr->seedno[j];u_k++;}
							 }
						    flag_end=find_location(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len);
							if(flag_end==0)continue;
							if(temp_score[repeat_loc]<6)continue;
							canidate_temp.score=temp_score[repeat_loc];
							loc_seed=temp_seedn[repeat_loc];
							//loc_list=temp_list[repeat_loc];
							location_loc[0]=start_loc+location_loc[0];
							loc_list=location_loc[0];
							readno=binary(llocation,location_loc[0],lread_count);
							readstart1=llocation[readno]; readend1=llocation[readno+1];
							length1=readend1-readstart1;
							readno=readno;
							if(indexread[readno].readno>read_name)continue;
							if(indexread[readno].readno==read_name){
								u_k=readstart1/ZV;temp_spr=database+u_k;s_k=readstart1%ZV;
								for(j=0,k=0;j<temp_spr->score&&j<SM;j++)if(temp_spr->loczhi[j]<s_k){temp_spr->loczhi[k]=temp_spr->loczhi[j];k++;}
								temp_spr->score=k;
								for(temp_spr++,u_k++,k=readend1/ZV;u_k<k;u_k++,temp_spr++)temp_spr->score=0;
								for(j=0,k=0,s_k=readend1%ZV;j<temp_spr->score&&j<SM;j++)if(temp_spr->loczhi[j]>s_k){temp_spr->loczhi[k]=temp_spr->loczhi[j];k++;}
								temp_spr->score=k;
							}
							else {
								  canidate_temp.readno=readno;
							      canidate_temp.readstart=readstart1;
					              location_loc[1]=(location_loc[1]-1)*BC;
					              left_length1=location_loc[0]-readstart1+seed_len-1;right_length1=readend1-location_loc[0];
					              left_length2=location_loc[1]+seed_len-1;right_length2=read_len-location_loc[1];
					              if(left_length1>=left_length2)num1=left_length2;
					              else num1=left_length1;
					              if(right_length1>=right_length2)num2=right_length2;
					              else num2=right_length1;
					              if(num1+num2<400)continue;
								   seedcount=0;
								  //find all left seed 
								  canidate_temp.loc1=location_loc[0];canidate_temp.num1=num1;
								  canidate_temp.loc2=location_loc[1];canidate_temp.num2=num2;
								  canidate_temp.left1=left_length1;canidate_temp.left2=left_length2;
								   canidate_temp.right1=right_length1;canidate_temp.right2=right_length2;
								  for(u_k=*index_spr-2,k=num1/ZV,temp_spr1=temp_spr-2;u_k>=0&&k>=0;temp_spr1--,k--,u_k--)if(temp_spr1->score>0){
									  start_loc=u_k*ZV;
									  for(j=0,s_k=0;j<temp_spr1->score;j++)if(fabs((loc_list-start_loc-temp_spr1->loczhi[j])/((loc_seed-temp_spr1->seedno[j])*BC*1.0)-1.0)<0.10){
										  seedcount++;s_k++;
									  }
									  if(s_k*1.0/temp_spr1->score>0.4)temp_spr1->score=0;
								  }
								  //find all right seed
							     for(u_k=*index_spr+1,k=num2/ZV,temp_spr1=temp_spr+1;k>0;temp_spr1++,k--,u_k++)if(temp_spr1->score>0){
									  start_loc=u_k*ZV;
									  for(j=0,s_k=0;j<temp_spr1->score;j++)if(fabs((start_loc+temp_spr1->loczhi[j]-loc_list)/((temp_spr1->seedno[j]-loc_seed)*BC*1.0)-1.0)<0.10){
										  seedcount++;s_k++;
									  }
									  if(s_k*1.0/temp_spr1->score>0.4)temp_spr1->score=0;
								  }
								 canidate_temp.score=canidate_temp.score+seedcount;
								 if(ii==1)canidate_temp.chain='F';
								 else canidate_temp.chain='R';
								 //insert canidate position or delete this position
						          low=0;high=canidatenum-1;
						          while(low<=high){
							           mid=(low+high)/2;
							           if(mid>=canidatenum||canidate_loc[mid].score<canidate_temp.score)high=mid-1;
							           else low=mid+1;
						         }
                                 if(canidatenum<MAXC)for(u_k=canidatenum-1;u_k>high;u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
								 else for(u_k=canidatenum-2;u_k>high;u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
								 if(high+1<MAXC)canidate_loc[high+1]=canidate_temp;
							     if(canidatenum<MAXC)canidatenum++;
								 else canidatenum=MAXC;
							}
		 }
 		 for(i=0,index_spr=index_list;i<cc1;i++,index_spr++){database[*index_spr].score=0;database[*index_spr].index=-1;}
		}


 for(i=0;i<canidatenum;i++){  		                 
						//left alignment search
	                    location_loc[0]=canidate_loc[i].loc1;
			    location_loc[1]=canidate_loc[i].loc2;
			    readno=canidate_loc[i].readno;
			    num1=canidate_loc[i].num1;
			    num2=canidate_loc[i].num2;
				left_length1=canidate_loc[i].left1;left_length2=canidate_loc[i].left2;
				right_length1=canidate_loc[i].right1;right_length2=canidate_loc[i].right2;
			    if(canidate_loc[i].chain=='F')onedata=onedata1;
			    else if(canidate_loc[i].chain=='R')onedata=onedata2;
                            readstart1=canidate_loc[i].readstart;
                            seq_pr1=seq+location_loc[0]+seed_len-2;seq_pr2=onedata+location_loc[1]+seed_len-1;
			    left_loc1=0;left_loc=0;
			    resultstore->left_store1[0]='\0';resultstore->left_store2[0]='\0';
			    flag_end=1;
			    while(flag_end){
				if(num1>600){
				     for(s_k=0;s_k<DN;s_k++,seq_pr1--,seq_pr2--){seq1[s_k]=*seq_pr1;seq2[s_k]=*seq_pr2;}
					seq1[s_k]='\0';seq2[s_k]='\0';
                                }
				else {
                                    flag_end=0;
				    for(s_k=0;s_k<num1;s_k++,seq_pr1--,seq_pr2--){seq1[s_k]=*seq_pr1;seq2[s_k]=*seq_pr2;}
				    seq1[s_k]='\0';seq2[s_k]='\0';
                                }
								
				for(loc=0;loc<2000;loc++){longstr1[loc]=0;longstr2[loc]=0;}
                                align_flag=align(seq1,seq2,ErrorRate*s_k,400,&strvalue,longstr1,longstr2,d_path,aln_path);
				if(align_flag==1){
					strvalue.q_aln_str[strvalue.aln_str_size]='\0';strvalue.t_aln_str[strvalue.aln_str_size]='\0';
					for(k=strvalue.aln_str_size-1,loc=0,sci=0,eit=0;k>-1&&eit<4;k--){
						if(strvalue.q_aln_str[k]!='-')loc++;
						if(strvalue.t_aln_str[k]!='-')sci++;
						if(strvalue.q_aln_str[k]==strvalue.t_aln_str[k])eit++;
						else eit=0;
                                        }
								     
					if(flag_end==1){
						loc=DN-strvalue.aln_q_e+loc;sci=DN-strvalue.aln_t_e+sci;
						if(loc==DN)align_flag=0;
						seq_pr1=seq_pr1+loc;seq_pr2=seq_pr2+sci;
						strvalue.q_aln_str[k+1]='\0';strvalue.t_aln_str[k+1]='\0';
						left_loc1=left_loc1+DN-loc;left_loc=left_loc+DN-sci;			
                                        }
					else {
						loc=num1-strvalue.aln_q_e;sci=num1-strvalue.aln_t_e;
						if(loc==num1)align_flag=0;
						left_loc1=left_loc1+num1-loc;left_loc=left_loc+num1-sci;
						seq_pr1=seq_pr1+loc+1;seq_pr2=seq_pr2+sci+1;
									//strvalue.q_aln_str[k+7]='\0';strvalue.t_aln_str[k+7]='\0';			
                                        }
								
					if(align_flag==1){
					strcat(resultstore->left_store1,strvalue.q_aln_str);
					strcat(resultstore->left_store2,strvalue.t_aln_str);
					if(left_length1-left_loc1>=left_length2-left_loc)num1=left_length2-left_loc;
					else num1=left_length1-left_loc1;
					}
				}
								
			if(align_flag==0)break;
                            }
							
		//if(flag_end==1&&align_flag==0)continue;
		//right alignment search
		right_loc1=0;right_loc=0;seq_pr1=seq+location_loc[0]-1;seq_pr2=onedata+location_loc[1];
		resultstore->right_store1[0]='\0';resultstore->right_store2[0]='\0';
		flag_end=1;
		while(flag_end){
			if(num2>600){
				for(s_k=0;s_k<DN;s_k++,seq_pr1++,seq_pr2++){seq1[s_k]=*seq_pr1;seq2[s_k]=*seq_pr2;}
				seq1[s_k]='\0';seq2[s_k]='\0';
                        }
                        else {
                               flag_end=0;
			       for(s_k=0;s_k<num2;s_k++,seq_pr1++,seq_pr2++){seq1[s_k]=*seq_pr1;seq2[s_k]=*seq_pr2;}
				seq1[s_k]='\0';seq2[s_k]='\0';
                        }
		for(loc=0;loc<2000;loc++){longstr1[loc]=0;longstr2[loc]=0;}
                align_flag=align(seq1,seq2,ErrorRate*s_k,400,&strvalue,longstr1,longstr2,d_path,aln_path);
		if(align_flag==1){
			strvalue.q_aln_str[strvalue.aln_str_size]='\0';strvalue.t_aln_str[strvalue.aln_str_size]='\0';
			for(k=strvalue.aln_str_size-1,loc=0,sci=0,eit=0;k>-1&&eit<4;k--){
				if(strvalue.q_aln_str[k]!='-')loc++;
				if(strvalue.t_aln_str[k]!='-')sci++;
				if(strvalue.q_aln_str[k]==strvalue.t_aln_str[k])eit++;
				else eit=0;
                        }
                        if(flag_end==1){
				loc=DN-strvalue.aln_q_e+loc;sci=DN-strvalue.aln_t_e+sci;
				//right_loc1=right_loc1+DN-loc;right_loc=right_loc+DN-sci;
				if(loc==DN)align_flag=0;
				seq_pr1=seq_pr1-loc;seq_pr2=seq_pr2-sci;
				strvalue.q_aln_str[k+1]='\0';strvalue.t_aln_str[k+1]='\0';
				right_loc1=right_loc1+DN-loc;right_loc=right_loc+DN-sci;
                        }
                        else {
				loc=num2-strvalue.aln_q_e;sci=num2-strvalue.aln_t_e;
				if(loc==num2)align_flag=0;
				right_loc1=right_loc1+num2-loc;right_loc=right_loc+num2-sci;
				seq_pr1=seq_pr1-loc;seq_pr2=seq_pr2-sci;
				//strvalue.q_aln_str[k+7]='\0';strvalue.t_aln_str[k+7]='\0';
                        }
			if(align_flag==1){
			        strcat(resultstore->right_store1,strvalue.q_aln_str);
				strcat(resultstore->right_store2,strvalue.t_aln_str);
				if(right_length1-right_loc1>=right_length2-right_loc)num2=right_length2-right_loc;
				else num2=right_length1-right_loc1;
                        }
                }
		if(align_flag==0)break;
                }
		//if(flag_end==1&&align_flag==0)continue;
		//yushu=strlen(resultstore->left_store1)+strlen(resultstore->right_store1);
                //collect left
                
        
		
		s_k=strlen(resultstore->left_store1);
		for(j=0,loc=0,k=0;j<s_k;j++){
			if(resultstore->left_store1[j]!='-'){resultstore->out_store1[loc]=resultstore->left_store1[j];loc++;}
			if(resultstore->left_store2[j]!='-'){resultstore->out_store2[k]=resultstore->left_store2[j];k++;}
                }
		resultstore->out_store1[loc]='\0';
		resultstore->out_store2[k]='\0';
		//printf("%s\n%s\n%s\n%s\n",resultstore->out_store1,resultstore->out_store2,resultstore->left_store1,resultstore->left_store2);
		//if(read_name==1548&&readno==1341){
		//	printf("xiao");
	//	}
		string_check(resultstore->out_store1,resultstore->out_store2,resultstore->left_store1,resultstore->left_store2);
		//printf("%s\n%s\n",resultstore->left_store1,resultstore->left_store2);
		/*s_k=strlen(resultstore->left_store1);
		for(j=0,loc=0,k=0;j<s_k;j++){
			if(resultstore->left_store1[j]!='-'){resultstore1->out_store1[loc]=resultstore->left_store1[j];loc++;}
			if(resultstore->left_store2[j]!='-'){resultstore1->out_store2[k]=resultstore->left_store2[j];k++;}
                }
		resultstore1->out_store1[loc]='\0';
		resultstore1->out_store2[k]='\0';		
        if(strcmp(resultstore1->out_store1,resultstore->out_store1)!=0||strcmp(resultstore1->out_store2,resultstore->out_store2)!=0){
			k=strlen(resultstore1->out_store1)-1;        u_k=0;
			while(resultstore1->out_store1[k]==resultstore->out_store1[k]){k--;u_k++;}
			printf("xiao"); 
		}*/






		 /*printf("%d\t%d\n",indexread[readno].readno,read_name);
			if(readno==5&&read_name==306){
				printf("xiao");
			}*/


		//output result
		u_k=strlen(resultstore->left_store1);
		for(j=u_k-1,loc=0,k=0,eit=0;j>-1;j--,k++){
			FR=resultstore->left_store1[j];resultstore->out_store1[k]=FR;
			if(FR!='-')loc++;
			FR=resultstore->left_store2[j];resultstore->out_store2[k]=FR;
			if(FR!='-')eit++;
                }
		resultstore->out_store1[k]='\0';resultstore->out_store2[k]='\0';
		//resultstore->out_store1[k-seed_len]='\0';resultstore->out_store2[k-seed_len]='\0';
		if(u_k==seed_len-1){left_loc1=location_loc[0]+seed_len-loc-1;left_loc=location_loc[1]+seed_len-eit;}
		else if(u_k>0){left_loc1=location_loc[0]+seed_len-loc;left_loc=location_loc[1]+seed_len-eit+1;}
		else {left_loc1=location_loc[0];left_loc=location_loc[1]+1;}
		s_k=strlen(resultstore->right_store1);
		for(k=0,loc=0,eit=0;k<s_k;k++){
			if(resultstore->right_store1[k]!='-')loc++;
			if(resultstore->right_store2[k]!='-')eit++;				    
                }
		if(s_k>0){right_loc1=location_loc[0]+loc-1;right_loc=location_loc[1]+eit;}
		else     {right_loc1=location_loc[0]+seed_len-1;right_loc=location_loc[1]+seed_len;}
		if(s_k>=seed_len&&u_k>=seed_len){
			strcat(resultstore->out_store1,resultstore->right_store1+seed_len);
			strcat(resultstore->out_store2,resultstore->right_store2+seed_len);
                }
                else if(u_k<seed_len){
			strcpy(resultstore->out_store1,resultstore->right_store1);
			strcpy(resultstore->out_store2,resultstore->right_store2);
				}
		//strcat(resultstore->out_store1,resultstore->right_store1);
		//strcat(resultstore->out_store2,resultstore->right_store2);		
                FR=canidate_loc[i].chain;
				//output file 
		left_loc1=left_loc1-readstart1;
		right_loc1=right_loc1-readstart1;
		if(right_loc1-left_loc1>450){  
			
			numMatch=0;numMismatch=0;numIns=0;numDel=0;
			u_k=strlen(resultstore->out_store1);
			for(j=0;j<u_k;j++){
				if(resultstore->out_store1[j]==resultstore->out_store2[j]&&resultstore->out_store2[j]!='-'){
					numMatch++;
					matchPattern[j]='|';
				}
				else {
					numMismatch++;
					numIns++;
				    matchPattern[j]='*';
				}
			}
			matchPattern[u_k]='\0';
			/*m5 file
			strand1='+';
			if(FR=='F')strand2='+';
			else strand2='-';
		   if(strand2=='+')fprintf(outfile[threadint],"%d %d %d %d %c %d %d %d %d %c 100 %d %d %d %d 10 %s %s %s\n",indexread[readno].readno,indexread[readno].length,left_loc1-1,right_loc1,strand1,read_name,read_len,left_loc-1,right_loc,strand2,numMatch,numMismatch, numIns,numDel,resultstore->out_store1,matchPattern,resultstore->out_store2);
		   else fprintf(outfile[threadint],"%d %d %d %d %c %d %d %d %d %c 100 %d %d %d %d 10 %s %s %s\n",indexread[readno].readno,indexread[readno].length,left_loc1-1,right_loc1,strand1,read_name,read_len,read_len-right_loc,read_len-left_loc+1,strand2,numMatch,numMismatch, numIns,numDel,resultstore->out_store1,matchPattern,resultstore->out_store2);
			*/
			jscore=2*u_k-numMismatch;
			jscore=jscore*30*4/(u_k);
		   if(FR=='F')fprintf(outfile[threadint],"%d %d %.3f 100 0 %d %d %d 0 %d %d %d\n",indexread[readno].readno,read_name,jscore,left_loc1-1,right_loc1,indexread[readno].length,left_loc-1,right_loc,read_len);
		   else fprintf(outfile[threadint],"%d %d %.3f 100 0 %d %d %d 1 %d %d %d\n",indexread[readno].readno,read_name,jscore,left_loc1-1,right_loc1,indexread[readno].length,read_len-right_loc,read_len-left_loc+1,read_len);
			  //fprintf(fid,"%c\t%d\t%d\n",FR,indexread[readno].readno,read_name);
                }
 } 
      }
    
  }
   //fclose(fid);
   free(database);
   free(index_list);
   free(index_score);
   free(d_path);
   free(resultstore);
   //printf("xiao");
}






void *multithread(void *p) //
{
    int thread_arg,localnumber,localthreadno;
    sleep(1);
    pthread_mutex_lock(&mutilock);
    localthreadno=runthreadnum;
    runthreadnum++;
    //printf("threadno=%d\n",localthreadno);
    pthread_mutex_unlock(&mutilock); 
    sleep(1);
    //outputfile(localthreadno);
    pairwise_mapping(localthreadno);
}



int load_fastq(FILE *fq,int startno){
	int readlen,readno,sum=0,flag,i;
	char onedata[RM],*pre,tempstr[200];
	readcount=0;
	pre=savework;
	while(fscanf(fq,">%[^\n]s",tempstr)!=EOF&&fscanf(fq,"%s\n",pre)!=EOF&&readcount<SVM&&sum<MAXSTR){
		readno=startno+readcount;
		readinfo[readcount].seqloc=pre;
		readinfo[readcount].readno=readno;
        readlen=strlen(pre);
		for(i=0;i<readlen;i++)if(pre[i]>='a')pre[i]=toupper(pre[i]);
		readinfo[readcount].readlen=readlen;
		sum=sum+readlen+1;
	        pre=pre+readlen+1;
                readcount++;
	}
	flag=feof(fq);
        if(flag==0){
			readno=startno+readcount;
                readinfo[readcount].seqloc=pre;
		readinfo[readcount].readno=readno;
                readlen=strlen(pre);
		readinfo[readcount].readlen=readlen;
                readcount++;
                return(1);
        }
        else return(0);
}

int param_read(int argc1, char *argv1[], char *pathway, int *threadnum, int *starts, int *ende,int *isP, int *isT, int *isS,  int *isE)
{
    int i, j, k;
    char tempstr[300];
    *isP = *isT = *isS = *isE = 0;
    for (i = 1; i < argc1; i++)
    {
        k = strlen(argv1[i]);
        for (j = 2; j < k; j++)tempstr[j - 2] = argv1[i][j];
        tempstr[j - 2] = '\0';
        //printf("%s\n",tempstr);
        if (argv1[i][0] == '-')
        {
            switch (argv1[i][1])
            {
            case 'P':
            {
                strcpy(pathway, tempstr);
                *isP = 1;
                break;
            }
            case 'T':
            {
                *threadnum= str2num(tempstr);
                *isT = 1;
                break;
            }
            case 'S':
            {
                *starts= str2num(tempstr);
                *isS = 1;
                break;
            }
            case 'E':
            {
                *ende= str2num(tempstr);
                *isE= 1;
                break;
            }
            }
        }
        else return (-1);
    }
    if (strlen(pathway) < 1 || *starts <0 || *ende <0|| *threadnum<0)return (-1);
    return (1);
}

int filesize(char *stream){
	int curpos, length;
	FILE *fp;
	fp=fopen(stream,"r");
	curpos = ftell(fp);
	fseek(fp, 0L, SEEK_END);
	length = ftell(fp);
	fseek(fp, curpos, SEEK_SET);
	return length;
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
int main(int argc,char *argv[]){
	   char tempstr[300],workpath[300],tempstr1[50];
	    int i,k,flag,isP,isT,isS,isE,isF,Istart; 
        int splitsize,corenum,coreid,filecount,fileno,startid,endid;
        int fileflag,threadno,threadflag,*filestart,*fileend;
        FILE *fastq,*fp;

		flag = param_read(argc, argv, workpath,&threadnum, &startid,&endid, &isP, &isT, &isS, &isE);
		filecount=endid;

        //creat share memory
        filestart=(int *)malloc(filecount*sizeof(int));
		fileend=(int *)malloc(filecount*sizeof(int));
        seed_len=13;
	        //multiple thread running program 
		// -allreads -allbases -b 1 -e 20000
		sprintf(tempstr,"%s/ovlprep",workpath);
		fp=fopen(tempstr,"r");
        fileno=1;
        while(fscanf(fp," %s %s %s %d %s %d\n",tempstr,tempstr,tempstr,&k,tempstr,&curreadcount)!=EOF&&fileno<=endid){
					filestart[fileno-1]=k;
					fileend[fileno-1]=curreadcount;
					fileno++;
		}
		fclose(fp);

		curreadcount=fileend[startid-1]-filestart[startid-1]+1;
           //one file read load
       fileidconvert(tempstr1,startid); 
       sprintf(tempstr,"%s/%s.fasta",workpath,tempstr1);
	   splitsize=filesize(tempstr);
	   indexread=(readmemory *)malloc(sizeof(readmemory)*(curreadcount+1));
	   llocation=(int *)malloc(sizeof(int)*(curreadcount+1));
	   STRMEM=(char *)malloc(sizeof(char)*splitsize);   
	   savework=(char *)malloc((MAXSTR+RM)*sizeof(char));
        readinfo=(ReadFasta*)malloc((SVM+2)*sizeof(ReadFasta));
        thread=(pthread_t*)malloc(threadnum*sizeof(pthread_t));
        outfile=(FILE **)malloc(threadnum*sizeof(FILE *));
	    seqcount=load_read(curreadcount,STRMEM,tempstr,llocation,filestart[startid-1]);
           //one file creat share ref. index
        for(threadno=0;threadno<threadnum;threadno++){
           sprintf(tempstr,"%s/%d_%d.r",workpath,startid,threadno);
           outfile[threadno]=fopen(tempstr,"w");              
        }
         creat_ref_index(STRMEM,seqcount);
         
       for(i=startid;i<=filecount;i++){
       	   fileidconvert(tempstr1,i);
           sprintf(tempstr,"%s/%s.fasta",workpath,tempstr1);
           fastq=fopen(tempstr,"r");
           //multi process thread
          fileflag=1;
		   Istart=filestart[i-1];
		   readcount=0;
          while(fileflag){
			        Istart=Istart+readcount;
                    fileflag=load_fastq(fastq,Istart);
                   if(readcount%PLL==0)terminalnum=readcount/PLL;
                   else terminalnum=readcount/PLL+1;
				   if(readcount<=0)break;
                   runnumber=0;
                   runthreadnum=0;
                   pthread_mutex_init(&mutilock,NULL); 
                    //creat thread
                   if(readcount>0){
                       for(threadno=0;threadno<threadnum;threadno++){
                       threadflag= pthread_create(&thread[threadno], NULL, multithread, NULL);
                       if(threadflag){
                              printf("ERROR; return code is %d\n", threadflag);  
                              return EXIT_FAILURE;
                        }
                       }
                    //waiting thread
                      for(threadno=0;threadno<threadnum;threadno++)pthread_join(thread[threadno],NULL);
                   }
                   
                //   pairwise_mapping(1);
          }
         fclose(fastq);          
    }
    //clear creat index memory
    free(countin);free(databaseindex);
    free(allloc);free(indexread);free(llocation);
    free(STRMEM);
    for(threadno=0;threadno<threadnum;threadno++)fclose(outfile[threadno]);
    free(outfile);free(savework);free(readinfo);free(thread);
    return 0;
}

