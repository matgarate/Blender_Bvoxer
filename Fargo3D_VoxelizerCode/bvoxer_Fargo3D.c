//Voxelizer code for Fargo3D simulations in spherical coordinates.
//By Pablo Benitez-Llambay  & Matías Gárate.


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define index i+j*nx+k*nx*ny
#define TRUE 1
#define FALSE 0

typedef double real;

int ns;
int nr;
int nc;

int logscale=FALSE;
int logvalues=FALSE;
int tri=FALSE;

//Voxel resolution
int nx = 400;
int ny = 400;
int nz = 100;

real *data;
float *vdata;



real *r;
real *theta;
real *phi;


real trilinear_interpolation(real x0, real x1,
			      real y0, real y1,
			      real z0, real z1,
			      real c000, real c100, real c110, real c010,
			      real c001, real c101, real c111, real c011,
			      real x, real y, real z) {

  //The cube is defined by the vertices (x0,y0,z0), (x1,y1,z1)
  real xd,yd,zd;
  real c00, c01, c10, c11;
  real c0, c1;
  real c;

  xd = (x-x0)/(x1-x0);
  yd = (y-y0)/(y1-y0);
  zd = (z-z0)/(z1-z0);

  // Linear interpolation along x
  c00 = c000*(1.0-xd) + c100*xd;
  c01 = c001*(1.0-xd) + c101*xd;
  c10 = c010*(1.0-xd) + c110*xd;
  c11 = c011*(1.0-xd) + c111*xd;
  // Now linear interpolation along y
  c0  = c00*(1.0-yd) + c10*yd;
  c1  = c01*(1.0-yd) + c11*yd;

  // Finally, linear interpolation along z
  return (c0*(1-zd) + c1*zd);
}



real GetCartesianValue(real r_cell, real phi_cell, real theta_cell){
//This function is returning the closest value in the data array to the real spherical coordinates.
  real value = 0;

  int i,j,k;

  i = (int)((phi_cell-phi[0])/(phi[ns]-phi[0])*ns);
  j = (int)((r_cell-r[0])/(r[nr]-r[0])*nr);

  if(logscale)
  	j = (int)(nr*log10(r_cell/r[0])/log10(r[nr]/r[0]));
  k = (int)((theta_cell-theta[0])/(theta[nc]-theta[0])*nc);

  if(j>=0 && k>=0 && i>=0 && j<nr && k<nc) {

    if(tri)
    {
        int ip=i+1;
        int jp=j+1;
        int kp=k+1;

        if(ip>=ns) ip=ip-ns;
        if(jp>=nr-1) jp=j;
        if(kp>=nc-1) kp=k;

        //For these values you are using the midpoints.
        //So go from i=0->ns-1  j=0->nr-1 k=0->nc-1, same for ip(cyclic),jp,kp
        real v000=data[i+j*ns+k*ns*nr];
        real v100=data[ip+j*ns+k*ns*nr];
        real v110=data[ip+jp*ns+k*ns*nr];
        real v010=data[i+jp*ns+k*ns*nr];
        real v001=data[i+j*ns+kp*ns*nr];
        real v101=data[ip+j*ns+kp*ns*nr];
        real v111=data[ip+jp*ns+kp*ns*nr];
        real v011=data[i+jp*ns+kp*ns*nr];

        //For here your are using the boundaries so i=0->ns-1  j=0->nr-1 k=0->nc-1, same for i+1,j+1,k+1
        value = trilinear_interpolation(phi[i],phi[i+1],
                                        r[j],r[j+1],
                                        theta[k],theta[k+1],
                                        v000,v100,v110,v010,
                                        v001,v101,v111,v011,
                                        phi_cell,r_cell,theta_cell);

    }
    else
        value = data[i+j*ns+k*ns*nr];


  } else
    value = 0.0;

  return value;
}

void usage()
{
    fprintf(stderr,"Welcome to the Fargo3D Voxelizer:\n");
    fprintf(stderr,"Developed by Pablo Benitez-Llambay  & Matias Garate:\n\n");


    fprintf(stderr,"USAGE:\n");
    fprintf(stderr,"[-in <file>] Input file from Fargo3d data field.\n");
    fprintf(stderr,"[-out <file>] Output file in blender voxel format.\n");
    fprintf(stderr,"[-ns <int>] Number of azimuthal sections in Fargo3D simulation.\n");
    fprintf(stderr,"[-nr <int>] Number of radial sections in Fargo3D simulation.\n");
    fprintf(stderr,"[-nc <int>] Number of colatitude sections in Fargo3D simulation.\n");

    fprintf(stderr,"[-vnx <int>] Blender voxel resolution in the x axis. Default: 400.\n");
    fprintf(stderr,"[-vny <int>] Blender voxel resolution in the y axis. Default: 400.\n");
    fprintf(stderr,"[-vnz <int>] Blender voxel resolution in the z axis. Default: 100.\n");

    fprintf(stderr,"[-logr] Fargo simulation radial coordinate is in logscale.\n");
    fprintf(stderr,"[-logv] Output fargo field values in logscale.\n");
    fprintf(stderr,"[-tri] Use Trilinear Interpolation.\n");


    exit(0);
}


int main(int argc, char **argv)
{


    char *filename;
    char *filename_out;

    int i,j,k;

    FILE* domain_x;
    FILE* domain_y;
    FILE* domain_z;
    FILE* fi;
    FILE* fo;



    if(argc==1) usage();
    i=1;
    while (i < argc) {
        if (!strcmp(argv[i],"-ns")) {
			++i;
			if (i >= argc) usage();
			ns = atoi(argv[i]);
			++i;
			}
	else if (!strcmp(argv[i],"-nr")) {
			++i;
			if (i >= argc) usage();
			nr = atoi(argv[i]);
			++i;
			}
        else if (!strcmp(argv[i],"-nc")) {
			++i;
			if (i >= argc) usage();
			nc = atoi(argv[i]);
			++i;
			}
        else if (!strcmp(argv[i],"-vnx")) {
			++i;
			if (i >= argc) usage();
			nx = atoi(argv[i]);
			++i;
			}
        else if (!strcmp(argv[i],"-vny")) {
			++i;
			if (i >= argc) usage();
			ny = atoi(argv[i]);
			++i;
			}
        else if (!strcmp(argv[i],"-vnz")) {
			++i;
			if (i >= argc) usage();
			nz = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-in")) {
			++i;
			if (i >= argc) usage();
			filename = argv[i];
			++i;
			}
		else if (!strcmp(argv[i],"-out")) {
			++i;
			if (i >= argc) usage();
			filename_out= argv[i];
			++i;
			}
		else if (!strcmp(argv[i],"-logr")) {
			logscale=TRUE;
			++i;
			}
		else if (!strcmp(argv[i],"-logv")) {
			logvalues=TRUE;
			++i;
			}
		else if (!strcmp(argv[i],"-tri")) {
			tri=TRUE;
			++i;
			}
        else usage();
    }


    //Initial min and max values. Ridiculously high and low to get the accurate min and max values.
    real min_value=1e30;
    real max_value=-1e30;

    real rr,ttheta,pphi;


    real *x;
    real *y;
    real *z;

    real *r_med;
    real *theta_med;
    real *phi_med;


    real dummy;

    real xmin,xmax,zmin,zmax,ymin,ymax;


    domain_x = fopen("domain_x.dat","r"); //phi
    domain_y = fopen("domain_y.dat","r"); //r
    domain_z = fopen("domain_z.dat","r"); //theta

    phi   = (real*)malloc((ns+1)*sizeof(real));
    r     = (real*)malloc((nr+1)*sizeof(real));
    theta = (real*)malloc((nc+1)*sizeof(real));

    phi_med   = (real*)malloc(ns*sizeof(real));
    r_med     = (real*)malloc(nr*sizeof(real));
    theta_med = (real*)malloc(nc*sizeof(real));

    x = (real*)malloc(nx*sizeof(real));
    y = (real*)malloc(ny*sizeof(real));
    z = (real*)malloc(nz*sizeof(real));


    //Avoiding ghost cells...
    for (i=0;i<3;i++) {
      fscanf(domain_y, "%lf", &dummy);
      fscanf(domain_z, "%lf", &dummy);
    }

    //Loading the domain data
    for (i=0;i<ns+1;i++) {
      fscanf(domain_x, "%lf", phi+i);
    }
    for (i=0;i<nr+1;i++) {
      fscanf(domain_y, "%lf", r+i);
    }
    for (i=0;i<nc+1;i++) {
      fscanf(domain_z, "%lf", theta+i);
    }

    for (i=0;i<ns;i++) {
      phi_med[i] = 0.5*(phi[i]+phi[i+1]);
    }
    for (i=0;i<nr;i++) {
      r_med[i] = 0.5*(r[i]+r[i+1]);
    }
    for (i=0;i<nc;i++) {
      theta_med[i] = 0.5*(theta[i]+theta[i+1]);
    }

    //Defining the voxel box limits
    xmin = -r_med[nr-1];
    xmax = r_med[nr-1];
    ymin = xmin;
    ymax = xmax;
    zmax = r_med[nr-1]*cos(theta_med[0]);
    zmin = r_med[nr-1]*cos(theta_med[nc-1]);

    for (i=0;i<nx;i++){
      x[i] = xmin+(xmax-xmin)/(nx-1)*i;
    }
    for (i=0;i<ny;i++){
      y[i] = ymin+(ymax-ymin)/(ny-1)*i;
    }
    for (i=0;i<nz;i++){
      z[i] = zmin+(zmax-zmin)/(nz-1)*i;
    }

    printf("xmin=%lf\n",xmin);
    printf("xmax=%lf\n",xmax);
    printf("ymin=%lf\n",ymin);
    printf("ymax=%lf\n",ymax);
    printf("zmin=%lf\n",zmin);
    printf("zmax=%lf\n",zmax);


    //Read Fargo data.
    data= (real*)malloc(nr*ns*nc*sizeof(real));
    fi = fopen(filename, "r");
    fread(data, sizeof(real), nr*ns*nc, fi);


    //Set logarithmic scale and find the max and min values.
    for (k=0; k<nc*ns*nr; k++) {
      if(logvalues) data[k] = log10(data[k]);

      if(data[k]>max_value) max_value = data[k];
      if(data[k]<min_value) min_value = data[k];

    }


    //min_value=-16.85;
    //max_value=-4.8;

    printf("MaxValue: %lf\n", max_value);
    printf("MinValue: %lf\n", min_value);

    //Normalize the data values to be between 0 - 1
    for (k=0; k<nc*ns*nr; k++)
      data[k] = (data[k] - min_value)/(max_value-min_value);


    //Set up the header
    //Store it in vdata array
    //Write the binary file that blender can read
    vdata= (float*)malloc(nx*ny*nz*sizeof(float));
    int header[]={nx,ny,nz,1};
    fo = fopen(filename_out,"w");
    fwrite(header,4*sizeof(int),1,fo);

    for (k=0;k<nz;k++) {
      for (j=0;j<ny;j++) {
        for (i=0;i<nx;i++) {
          rr = sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
          pphi = atan2(y[j],x[i]);
          ttheta = acos(z[k]/rr);
          vdata[index] = (float)GetCartesianValue(rr,pphi,ttheta);
        }
      }
    }

    fwrite(vdata, nx*ny*nz*sizeof(float), 1, fo);
    fclose(fo);
}
