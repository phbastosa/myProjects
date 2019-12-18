# include <stdio.h>
# include <stdlib.h>
# include <math.h>

void read_and_expand_model(int nx, int nz, int nxx, int nzz, int attenuation_layer, char model_vel_location[], char modelo_gerado[]) {
    int ii,jj,cont;
    float (*model)[nx] = malloc(sizeof(float[nz][nx]));   /* Original model */
    float (*aux_vp)[nxx] = malloc(sizeof(float[nzz][nxx])); /* Expanded model */

    FILE *rf,*wf;

    /* Reading the original model */
    rf = fopen(model_vel_location, "rb");
    for(ii=0; ii<nx; ii++) {
        for(jj=0; jj<nz; jj++) {
            fread(&model[jj][ii], sizeof(float), 1, rf);
        }
    }
    fclose(rf);

    /* Centering the original model in the expanded model*/
    for(ii=attenuation_layer; ii<nzz-attenuation_layer; ii++) {
        for(jj=attenuation_layer; jj<nxx-attenuation_layer; jj++) {
            aux_vp[ii][jj] = model[ii-attenuation_layer][jj-attenuation_layer]; // Centering
        }
    }

    /* Expanding the upper and bottom layer of the velocities array */
    for(ii=0; ii<attenuation_layer; ii++) {
        for(jj=attenuation_layer; jj<nxx-attenuation_layer; jj++) {
            aux_vp[ii][jj] = aux_vp[attenuation_layer][jj];              // Up            
            aux_vp[nzz-ii-1][jj] = aux_vp[nzz-attenuation_layer-1][jj];  // Botton
        }
    }

    /* Expanding the left and right layers of the velocities array */
    for(ii=0; ii<nzz; ii++) {
        for(jj=0; jj<attenuation_layer; jj++) {
            aux_vp[ii][nxx-jj-1] = aux_vp[ii][nxx-attenuation_layer-1]; // Right
            aux_vp[ii][jj] = aux_vp[ii][attenuation_layer];             // Left    
        }
    }

    /* Writing the expanded model */
    wf = fopen(modelo_gerado,"wb");
    for(ii=0; ii<nxx; ii++) {
        for(jj=0; jj<nzz; jj++) {
            fwrite(&aux_vp[jj][ii], sizeof(float), 1, wf);
        }
    }    
    fclose(wf);     
}

int main(int argc, char **argv) {

    int ii,jj;
    int nx = 2500;
    int nz = 270;
    int nxx = 2700;
    int nzz = 470;
    int attenuatiion_layer = 100;
    int *topo = (int *) malloc(nxx*sizeof(int));
    float (*vp)[nxx] = malloc(sizeof(float[nzz][nxx]));
    float (*vs)[nxx] = malloc(sizeof(float[nzz][nxx]));
    float (*rho)[nxx] = malloc(sizeof(float[nzz][nxx]));
    char * model = "modelo_teste.bin";
    char * expanded_model = "modelo_expandido_teste.bin";
    char * topography = "topografia_teste.txt";

    read_and_expand_model(nx,nz,nxx,nzz,attenuatiion_layer,model,expanded_model);

    FILE *rm = fopen(expanded_model,"rb");
    for(ii = 0; ii < nxx; ii++) {
        for(jj = 0; jj < nzz; jj++) {
            fread(&vp[jj][ii],sizeof(float),1,rm);
        } 
    }
    fclose(rm);

    FILE *rt = fopen(topography,"r");
    if(rt != NULL) {
        ii = 0;
        while(fscanf(rt,"%i",&topo[ii]) != EOF) {
            ++ii;                    
        }
    } 
    fclose(rt);

    // Parallel plane layer model    
    for(ii = 0; ii < nzz; ii++) 
    {
        for(jj = 0; jj < nxx; jj++)
        {
            if(ii < topo[jj]) 
            { 
            // if(ii < 100) {          // Vacuum layer
                rho[ii][jj] = 310*pow(vp[ii][jj],0.25);               
                vp[ii][jj]  = 0.0f;
                vs[ii][jj]  = 0.0f; 
            }  
            if(ii >= topo[jj])
            {        
                vs[ii][jj] = vp[ii][jj]/sqrt(3);
                rho[ii][jj] = 310*pow(vp[ii][jj],0.25);
            }
            // if((ii >= 100) && (ii < 150))  // Assembling the model itself 
            // {      
            //     vp[ii][jj] = 1800.0f;
            //     vs[ii][jj] = vp[ii][jj]/sqrt(3);
            //     rho[ii][jj] = 310*pow(vp[ii][jj],0.25);
            // }
            // if((ii >= 150) && (ii < 300))
            // {
            //     vp[ii][jj] = 2000.0f;
            //     vs[ii][jj] = vp[ii][jj]/sqrt(3);
            //     rho[ii][jj] = 310*pow(vp[ii][jj],0.25);
            // }
            // if(ii >= 300)
            // {
            //     vp[ii][jj] = 2900.0f;
            //     vs[ii][jj] = vp[ii][jj]/sqrt(3);
            //     rho[ii][jj] = 310*pow(vp[ii][jj],0.25);
            // }
        }
    }

    FILE *w_vp = fopen("vp.bin","wb");
    for(ii = 0; ii < nzz; ii++) {
        for(jj = 0; jj < nxx; jj++) {
            fwrite(&vp[ii][jj],sizeof(float),1,w_vp);
        } 
    }
    fclose(w_vp);

    FILE *w_vs = fopen("vs.bin","wb");
    for(ii = 0; ii < nzz; ii++) {
        for(jj = 0; jj < nxx; jj++) {
            fwrite(&vs[ii][jj],sizeof(float),1,w_vs);
        } 
    }
    fclose(w_vs);

    FILE *w_rho = fopen("rho.bin","wb");
    for(ii = 0; ii < nzz; ii++) {
        for(jj = 0; jj < nxx; jj++) {
            fwrite(&rho[ii][jj],sizeof(float),1,w_rho);
        } 
    }
    fclose(w_rho);

    return 0;
} 
