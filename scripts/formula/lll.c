#include <stdio.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
int main(int c,char**v){ fmpz_mat_t A; fmpz_mat_init(A,0,0); FILE*f=fopen(v[1],"r"); fmpz_mat_fread(f,A); fclose(f);
 fmpz_lll_t fl; fmpz_lll_context_init_default(fl); fmpz_lll(A,NULL,fl);
 f=fopen(v[2],"w"); fmpz_mat_fprint(f,A); fclose(f); fmpz_mat_clear(A); return 0;}
