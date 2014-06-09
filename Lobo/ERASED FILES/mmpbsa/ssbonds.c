/******************** Includes *********************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "/usr/local/src/num_rec_lib/nrutil.h"

/******************** Prototypes *********************************/

FILE *file_open(char *fname, char *acc);
void read_atom_pdb(char *buf,char *at_rec,int *p_n_at_pdb,
char *at_name,char *alt_loc,char *res_name,
char *chain, int *p_n_res, char *res_ins, float *p_x, float *p_y,
float *p_z, float *p_occ, float *p_temp_f, int *p_foo_num);

/******** General purpose utilities and subroutines **********/

void main(int argc, char *argv[]) {
      FILE *fp1, *fp2, *fp3, *fp4;
      char buf[121], buf1[6], fn[12], residue[50][121];
      char at_rec[5], at_name[5], alt_loc[2], c_null[2], res_name[4],
           chain[2], res_ins[2];
      int n_at_pdb, n_res, foo_num;
      float x, y, z, temp_f, occ;
      double *xsg, *ysg, *zsg;
      int res_counter = -999;
      int k, flag_nt = 0, flag_ct = 0, flag_crg = 0;
      int first = -999, is_first;
      int howmany = 0;
      int napr = 0;
      int *resnum, i, j;


if(argc != 3) 
{
printf("I need two filenames!\n"); 
exit(1);
}

/* Apro il file per vedere quante istidine ha */

fp1 = file_open(argv[1],"r");


        while(fgets(buf,120,fp1) != NULL)
            if(!strncmp("ATOM",buf,4)) 
            {
               read_atom_pdb(buf, at_rec, &n_at_pdb, at_name, alt_loc,
                       res_name, chain, &n_res, res_ins, &x, &y, &z, 
                       &occ, &temp_f, &foo_num);
                  napr++;
            }
xsg = vector(0,napr);
ysg = vector(0,napr);
zsg = vector(0,napr);
resnum = ivector(0,napr);
napr = 0;
fclose(fp1);
fp1 = file_open(argv[1],"r");
while(fgets(buf,120,fp1) != NULL)
	  {
            if(!strncmp("ATOM",buf,4)) 
            {
               read_atom_pdb(buf, at_rec, &n_at_pdb, at_name, alt_loc,
                       res_name, chain, &n_res, res_ins, &x, &y, &z, 
                       &occ, &temp_f, &foo_num);
            }
           napr++;
           xsg[napr] = x;
           ysg[napr] = y;
           zsg[napr] = z;
	   resnum[napr] = n_res;
	  }
                  fclose(fp1); 
fp2 = file_open(argv[2],"w");
fprintf(fp2,"* patch disulfide bonds\n*\n");

for(i=1; i<=napr; i++)
for(j=i+1; j<=napr; j++)
if ( ( (xsg[i] - xsg[j]) * (xsg[i] - xsg[j]) + 
       (ysg[i] - ysg[j]) * (ysg[i] - ysg[j]) + 
       (zsg[i] - zsg[j]) * (zsg[i] - zsg[j])) <= 5.0 ) 
	fprintf(fp2,"patch disu prot %5i prot %5i\n", resnum[i], resnum[j]);

}
FILE *file_open(char *fname,char *acc) {
        FILE *fp;
        fp =fopen(fname,acc);
        if (fp==NULL) {
                fprintf(stderr,"unable to open file %s\n",fname);
                exit(1);
        }
        return(fp);
}

void read_atom_pdb(char *buf, char *at_rec, int *p_n_at_pdb,
   char *at_name, char *alt_loc, char *res_name, char *chain,
   int *p_n_res, char *res_ins, float *p_x, float *p_y, float *p_z,
   float *p_occ, float *p_temp_f, int *p_foo_num) {

        char tok[10];

        strncpy(tok,buf,5);
        tok[5] = '\0';
        sscanf(tok,"%s", at_rec);

        strncpy(tok,buf + 6,6);
        tok[6] = '\0';
        sscanf(tok,"%i",p_n_at_pdb);

        strncpy(tok,buf + 12,5);
        tok[5] = '\0';
        sscanf(tok,"%s", at_name);

        strncpy(tok,buf + 16,1);
        tok[1] = '\0';
        if(sscanf(tok,"%s", alt_loc) == -1) strcpy(alt_loc," ");

        strncpy(tok,buf + 17,3);
        tok[3] = '\0';
        sscanf(tok,"%s", res_name);

        strncpy(tok,buf + 21,1);
        tok[1] = '\0';
        if(sscanf(tok,"%s", chain) == EOF) strcpy(chain," ");

        strncpy(tok,buf + 22,4);
        tok[4] = '\0';
        sscanf(tok,"%i", p_n_res);

        strncpy(tok,buf + 26,1);
        tok[1] = '\0';
        if (sscanf(tok,"%s", res_ins) == EOF) strcpy(res_ins," ");

        strncpy(tok,buf + 30,8);
        tok[8] = '\0';
        sscanf(tok,"%lf", p_x);

        strncpy(tok,buf + 38,8);
        tok[8] = '\0';
        sscanf(tok,"%lf", p_y);
        strncpy(tok,buf + 46,8);
        tok[8] = '\0';
        sscanf(tok,"%lf", p_z);

        strncpy(tok,buf + 54,6);
        tok[6] = '\0';
        sscanf(tok,"%lf", p_occ);

        strncpy(tok,buf + 60,6);
        tok[6] = '\0';
        sscanf(tok,"%lf", p_temp_f);

        strncpy(tok,buf + 67,3);
        tok[3] = '\0';
        sscanf(tok,"%i", p_foo_num);

}


