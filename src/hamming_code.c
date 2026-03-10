
#include <stdio.h>                                                                             
#include <stdlib.h>                                                                            
#include <string.h>                                                                            
#include <ctype.h>                                                                             

                     

int **imatrix_alloc(int rows, int cols) {                                                      
    int **A = (int**)malloc(rows * sizeof *A);                                                 
    if (!A) return NULL;                                                                       
    for (int i = 0; i < rows; ++i) {                                                           
        A[i] = (int*)calloc(cols, sizeof **A);                                                 
        if (!A[i]) {                                                                           
            for (int t = 0; t < i; ++t) free(A[t]);                                            
            free(A);                                                                           
            return NULL;                                                                       
        }
    }
    return A;                                                                                  
}

void imatrix_free(int **A, int rows) {                                                         
    if (!A) return;                                                                            
    for (int i = 0; i < rows; ++i) free(A[i]);                                                 
    free(A);                                                                                   
}

                     

static void strip_comment(char *s) {                                                           
    char *p = strstr(s, "//");                                                                 
    if (p) *p = '\0';                                                                          
    for (p = s; *p; ++p) {                                                                     
        if (*p == '#' || *p == '%' || *p == ';') { *p = '\0'; break; }                         
    }
    int L = (int)strlen(s);                                                                    
    while (L > 0 && isspace((unsigned char)s[L-1])) s[--L] = '\0';                             
    for (int i = 0; s[i]; ++i) if (s[i] == ',' || s[i] == '=') s[i] = ' ';                     
}

static int read_bits_line(const char *line, int need, int *out) {                              
    int cnt = 0, v = 0, used = 0;                                                              
    const char *p = line;                                                                      
    while (*p && cnt < need) {                                                                 
        while (*p && !isdigit((unsigned char)*p)) ++p;                                         
        if (!*p) break;                                                                        
        if (sscanf(p, "%d%n", &v, &used) == 1) {                                               
            out[cnt++] = (v != 0) ? 1 : 0;                                                     
            p += used;                                                                         
        } else break;                                                                          
    }
    return cnt;                                                                                
}

                                                          
int read_sim(const char *path, int *n, int *r, int ***H_out, int *m, int ***ERR_out) {         
    FILE *fp = fopen(path, "r");                                                               
    if (!fp) { perror("open Sim.txt"); return 0; }                                             

    char line[4096];                                                                           
    int got_n = 0, got_r = 0, got_m = 0;                                                       

    while ((!got_n || !got_r) && fgets(line, sizeof line, fp)) {                               
        strip_comment(line);                                                                   
        int v;                                                                                 
        if (!got_n && sscanf(line, " %d", &v) == 1) { *n = v; got_n = 1; continue; }           
        if (!got_r && sscanf(line, " %d", &v) == 1) { *r = v; got_r = 1; continue; }           
    }
    if (!got_n || !got_r || *n <= 0 || *r <= 0 || *r >= *n) {                                  
        fprintf(stderr, "[read_sim] invalid n/r\n");                                           
        fclose(fp);                                                                            
        return 0;                                                                              
    }

    int **H = imatrix_alloc(*r, *n);                                                           
    if (!H) { fclose(fp); return 0; }                                                          

    int rcount = 0;                                                                            
    while (rcount < *r && fgets(line, sizeof line, fp)) {                                      
        strip_comment(line);                                                                   
        int *tmp = (int*)malloc(*n * sizeof *tmp);                                             
        if (!tmp) { imatrix_free(H, *r); fclose(fp); return 0; }                               
        int got = read_bits_line(line, *n, tmp);                                               
        if (got == *n) {                                                                       
            for (int j = 0; j < *n; ++j) H[rcount][j] = tmp[j];                                
            ++rcount;                                                                          
        }
        free(tmp);                                                                             
    }
    if (rcount != *r) {                                                                        
        fprintf(stderr, "[read_sim] H rows %d/%d\n", rcount, *r);                              
        imatrix_free(H, *r);                                                                   
        fclose(fp);                                                                            
        return 0;                                                                              
    }

    while (!got_m && fgets(line, sizeof line, fp)) {                                           
        strip_comment(line);                                                                   
        int v;                                                                                 
        if (sscanf(line, " %d", &v) == 1) { *m = v; got_m = 1; }                               
    }
    if (!got_m || *m < 0) {                                                                    
        fprintf(stderr, "[read_sim] invalid m\n");                                             
        imatrix_free(H, *r);                                                                   
        fclose(fp);                                                                            
        return 0;                                                                              
    }

    int **ERR = imatrix_alloc(*m, *n);                                                         
    if (!ERR) { imatrix_free(H, *r); fclose(fp); return 0; }                                   

    int ecnt = 0;                                                                              
    while (ecnt < *m && fgets(line, sizeof line, fp)) {                                        
        strip_comment(line);                                                                   
        int *tmp = (int*)malloc(*n * sizeof *tmp);                                             
        if (!tmp) { imatrix_free(ERR, *m); imatrix_free(H, *r); fclose(fp); return 0; }        
        int got = read_bits_line(line, *n, tmp);                                               
        if (got == *n) {                                                                       
            for (int j = 0; j < *n; ++j) ERR[ecnt][j] = tmp[j];                                
            ++ecnt;                                                                            
        }
        free(tmp);                                                                             
    }
    fclose(fp);                                                                                

    if (ecnt != *m) {                                                                          
        fprintf(stderr, "[read_sim] error patterns %d/%d\n", ecnt, *m);                        
        imatrix_free(ERR, *m);                                                                 
        imatrix_free(H, *r);                                                                   
        return 0;                                                                              
    }

    *H_out = H;                                                                                
    *ERR_out = ERR;                                                                            
    return 1;                                                                                  
}

                         

int **build_G_from_sysH(int **H, int r, int n) {                                               
    int k = n - r;                                                                             
    int **G = imatrix_alloc(k, n);                                                             
    if (!G) return NULL;                                                                       

    for (int i = 0; i < k; ++i) {                                                              
        for (int j = 0; j < n; ++j) G[i][j] = 0;                                               
        G[i][i] = 1;                                                                           
    }
    for (int rr = 0; rr < r; ++rr) {                                                           
        for (int cc = 0; cc < k; ++cc) {                                                       
            G[cc][k + rr] = H[rr][cc];                                                         
        }
    }
    return G;                                                                                  
}

                                               
int check_H_right_identity(int **H, int r, int n) {                                            
    int ok = 1;                                                                                
    int k = n - r;                                                                             
    for (int i = 0; i < r; ++i) {                                                              
        for (int j = 0; j < r; ++j) {                                                          
            int expect = (i == j) ? 1 : 0;                                                     
            if (H[i][k + j] != expect) { ok = 0; break; }                                      
        }
        if (!ok) break;                                                                        
    }
    return ok;                                                                                 
}

                           

void gen_u_lfsr(int total_bits, unsigned char *out) {                                          
    unsigned char s[6] = {1,0,0,0,0,0};                                                        
    for (int t = 0; t < total_bits; ++t) {                                                     
        out[t] = s[0];                                                                         
        unsigned char nxt = (unsigned char)(s[1] ^ s[0]);                                      
        s[0]=s[1]; s[1]=s[2]; s[2]=s[3]; s[3]=s[4]; s[4]=s[5]; s[5]=nxt;                       
    }
}

                        

void encode_u_to_x(const unsigned char *u, int **G, int k, int n, unsigned char *x) {          
    for (int j = 0; j < n; ++j) {                                                              
        int acc = 0;                                                                           
        for (int i = 0; i < k; ++i) acc ^= (u[i] & G[i][j]);                                   
        x[j] = (unsigned char)(acc & 1);                                                       
    }
}

void syndrome_vec(int **H, int r, int n, const unsigned char *y, int *s) {                     
    for (int i = 0; i < r; ++i) {                                                              
        int acc = 0;                                                                           
        for (int j = 0; j < n; ++j) acc ^= (H[i][j] & y[j]);                                   
        s[i] = acc & 1;                                                                        
    }
}

int s_is_zero(const int *s, int r) {                                                           
    for (int i = 0; i < r; ++i) if (s[i]) return 0;                                            
    return 1;                                                                                  
}

int find_col_equal_s(int **H, int r, int n, const int *s) {                                    
    for (int j = 0; j < n; ++j) {                                                              
        int ok = 1;                                                                            
        for (int i = 0; i < r; ++i) if (H[i][j] != s[i]) { ok = 0; break; }                   
        if (ok) return j;                                                                      
    }
    return -1;                                                                                 
}

int check_GH_T_zero(int **G, int k, int n, int **H, int r) {                                   
    for (int i = 0; i < k; ++i) {                                                              
        for (int j = 0; j < r; ++j) {                                                          
            int acc = 0;                                                                       
            for (int t = 0; t < n; ++t) acc ^= (G[i][t] & H[j][t]);                            
            if ((acc & 1) != 0) return 0;                                                      
        }
    }
    return 1;                                                                                  
}

            
 

int has_overall_parity(int **H, int r, int n) {                                                
    int rows = n;                                                                              
    int cols = r;                                                                              
    int **A = imatrix_alloc(rows, cols);                                                       
    if (!A) return 0;                                                                          
    int *b = (int*)malloc(rows * sizeof *b);                                                   
    if (!b) { imatrix_free(A, rows); return 0; }                                               

    for (int i = 0; i < rows; ++i) {                                                           
        for (int j = 0; j < cols; ++j) A[i][j] = H[j][i];                                      
        b[i] = 1;                                                                              
    }

    int row = 0;                                                                               
    for (int col = 0; col < cols && row < rows; ++col) {                                       
        int piv = -1;                                                                          
        for (int i = row; i < rows; ++i) if (A[i][col]) { piv = i; break; }                    
        if (piv == -1) continue;                                                               
        if (piv != row) { int *tA = A[piv]; A[piv] = A[row]; A[row] = tA; int tb = b[piv];     
                           b[piv] = b[row]; b[row] = tb; }                                     
        for (int i = 0; i < rows; ++i) if (i != row && A[i][col]) {                            
            for (int j = col; j < cols; ++j) A[i][j] ^= A[row][j];                             
            b[i] ^= b[row];                                                                    
        }
        ++row;                                                                                 
    }

    int ok = 1;                                                                                
    for (int i = 0; i < rows; ++i) {                                                           
        int all0 = 1;                                                                          
        for (int j = 0; j < cols; ++j) if (A[i][j]) { all0 = 0; break; }                       
        if (all0 && b[i]) { ok = 0; break; }                                                   
    }

    imatrix_free(A, rows);                                                                     
    free(b);                                                                                   
    return ok;                                                                                 
}

                                 

int main(void) {                                                                               
    const char *sim = "Sim.txt";                                                               

    int n = 0, r = 0, m = 0;                                                                   
    int **H = NULL, **ERR = NULL;                                                              

    if (!read_sim(sim, &n, &r, &H, &m, &ERR)) {                                                
        fprintf(stderr, "[FATAL] parse Sim.txt failed\n");                                     
        return EXIT_FAILURE;                                                                   
    }

    int k = n - r;                                                                             
    printf("[INFO] n=%d columns, r=%d rows of H, k=%d, m=%d\n", n, r, k, m);                  

    if (!check_H_right_identity(H, r, n)) {                                                    
        fprintf(stderr, "[WARN] H right half is NOT strictly I_r; proceeding anyway.\n");      
    }

    int **G = build_G_from_sysH(H, r, n);                                                      
    if (!G) {                                                                                  
        fprintf(stderr, "OOM G\n");                                                            
        imatrix_free(ERR, m); imatrix_free(H, r);                                              
        return EXIT_FAILURE;                                                                   
    }

    printf("[CHECK] G·H^T == 0 ? %s\n", check_GH_T_zero(G, k, n, H, r) ? "OK" : "FAIL");      

    int extended = has_overall_parity(H, r, n);                                                
    printf("[INFO] overall parity present (Extended)? %s\n", extended ? "YES" : "NO");         

    int M = m;                                                                                 
    size_t total = (size_t)M * (size_t)k;                                                      
    unsigned char *u_all = (unsigned char*)malloc(total);                                      
    if (!u_all) {                                                                              
        fprintf(stderr, "OOM u_all\n");                                                        
        imatrix_free(G, k); imatrix_free(ERR, m); imatrix_free(H, r);                          
        return EXIT_FAILURE;                                                                   
    }
    gen_u_lfsr((int)total, u_all);                                                             

    FILE *fo = fopen("u.txt", "w");                                                            
    if (!fo) {                                                                                 
        perror("u.txt");                                                                       
        free(u_all); imatrix_free(G, k); imatrix_free(ERR, m); imatrix_free(H, r);             
        return EXIT_FAILURE;                                                                   
    }

    unsigned char *x  = (unsigned char*)malloc(n);                                             
    unsigned char *y  = (unsigned char*)malloc(n);                                             
    unsigned char *xh = (unsigned char*)malloc(n);                                             
    unsigned char *z  = (unsigned char*)malloc(n);                                             
    int *s = (int*)malloc(r * sizeof *s);                                                      
    if (!x || !y || !xh || !z || !s) {                                                         
        fprintf(stderr, "OOM vecs\n");                                                         
        if (x) free(x); if (y) free(y); if (xh) free(xh); if (z) free(z); if (s) free(s);     
        fclose(fo); free(u_all); imatrix_free(G, k); imatrix_free(ERR, m); imatrix_free(H, r); 
        return EXIT_FAILURE;                                                                   
    }

    for (int a = 0; a < M; ++a) {                                                              
        encode_u_to_x(&u_all[a * k], G, k, n, x);                                              

        for (int j = 0; j < n; ++j) z[j] = (unsigned char)ERR[a][j];                           
        for (int j = 0; j < n; ++j) y[j] = x[j] ^ z[j];                                        

        syndrome_vec(H, r, n, y, s);                                                           
        int s_zero = s_is_zero(s, r);                                                          

        int parity = 0;                                                                         
        for (int j = 0; j < n; ++j) parity ^= y[j];                                            

        int err_col = -1;                                                                      
        int multi   = 0;                                                                       
        for (int j = 0; j < n; ++j) xh[j] = y[j];                                              

        if (extended) {                                                                        
            if (!s_zero && parity == 0) {                                                      
                multi = 1;                                                                     
            } else if (!s_zero && parity == 1) {                                               
                err_col = find_col_equal_s(H, r, n, s);                                        
                if (err_col >= 0) xh[err_col] ^= 1;                                            
                else multi = 1;                                                                
            } else if (s_zero && parity == 1) {                                                
                                              
            } else {                                                                           
                                                                               
            }
        } else {                                                                               
            if (!s_zero) {                                                                     
                err_col = find_col_equal_s(H, r, n, s);                                        
                if (err_col >= 0) xh[err_col] ^= 1;                                            
                else multi = 1;                                                                
            }                                                                                  
        }

        int orig_col1 = (err_col >= 0) ? (err_col + 1) : -1;                                   

        printf("[CASE %d] , syndrome=", a);                                    
        for (int i = 0; i < r; ++i) printf("%d", s[i]);                                        

        if (multi) {                                                                           
            printf(" | decision: DETECTED 2-errors\n");                           
        } else if (extended && s_zero && parity == 1) {                                        
            printf(" | decision: overall-parity-only error (u not affected)\n");               
        } else if (s_zero) {                                                                   
            printf(" | decision: no-error\n");                                                 
        } else {                                                                               
            printf(" | decision: single-error corrected @ column %d\n", orig_col1);            
        }

        printf("        û: ");                                                                 
        for (int i = 0; i < k; ++i) {                                                          
            int val = multi ? 2 : (int)xh[i];                                                  
            fprintf(fo, "%d%s", val, (i + 1 < k) ? " " : "");                                  
            printf("%d%s", val, (i + 1 < k) ? " " : "");                                       
        }
        fprintf(fo, " %%Est. u%d\n", a);                                                       
        printf("  %%Est. u%d\n", a);                                                           
    }

    fclose(fo);                                                                                
    printf("[OUT] wrote u.txt with %d rows × %d bits per row\n", m, k);                        
    printf("Press Enter to exit..."); fflush(stdout);                                          
    getchar();                                                                                 

    free(s); free(z); free(xh); free(y); free(x);                                              
    free(u_all);                                                                               
    imatrix_free(G, k);                                                                        
    imatrix_free(ERR, m);                                                                      
    imatrix_free(H, r);                                                                        
    return EXIT_SUCCESS;                                                                       
}
