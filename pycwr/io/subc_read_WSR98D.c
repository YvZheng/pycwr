/* subc_read.c */
# include <stdio.h>

void subc_read(float * dat_tmp, const int ndat, const int offset, const int scale){
    int iii;
    for(iii=0; iii<ndat; iii++){
        if(dat_tmp[iii] >= 5){
            dat_tmp[iii] = (dat_tmp[iii] - offset) / scale; }
        else{
            dat_tmp[iii] = -9999.0; }
    }
}


