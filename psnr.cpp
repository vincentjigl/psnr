#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
  #include "stdint_w32.h"
  #define alloca _alloca
#else
  #include <stdint.h>
#endif

typedef struct{
	double mean_ssim;
	double mean_psnr; 
	double	stdv_ssim;
	double	stdv_psnr;
}ssim_psnr_t;

/****************************************************************************
 * structural similarity metric [from x264]
 ****************************************************************************/

#define x264_alloca(x) (void *)(((intptr_t)alloca((x)+15)+15)&~15)
#define XCHG(type,a,b) { type t = a; a = b; b = t; }
#define X264_MIN(a,b) ( (a)<(b) ? (a) : (b) )

static void ssim_4x4x2_core( const uint8_t *pix1, int stride1,
                             const uint8_t *pix2, int stride2,
                             int sums[2][4]);


static float ssim_end1( int s1, int s2, int ss, int s12 );


static float ssim_end4( int sum0[5][4], int sum1[5][4], int width );

float x264_pixel_ssim_wxh(
                           uint8_t *pix1, int stride1,
                           uint8_t *pix2, int stride2,
                           int width, int height );

int calculate(int n, char *cl[], ssim_psnr_t *ssim_psnr );


typedef struct
{
  int             width;
  int             height;
  unsigned char*  data;
} ColorComponent;

typedef struct
{
  ColorComponent lum;
  ColorComponent cb;
  ColorComponent cr;
} YuvFrame;


static void ssim_4x4x2_core( const uint8_t *pix1, int stride1,
                             const uint8_t *pix2, int stride2,
                             int sums[2][4])
{
    int x, y, z;
    for(z=0; z<2; z++)
    {
        uint32_t s1=0, s2=0, ss=0, s12=0;
        for(y=0; y<4; y++)
            for(x=0; x<4; x++)
            {
                int a = pix1[x+y*stride1];
                int b = pix2[x+y*stride2];
                s1  += a;
                s2  += b;
                ss  += a*a;
                ss  += b*b;
                s12 += a*b;
            }
        sums[z][0] = s1;
        sums[z][1] = s2;
        sums[z][2] = ss;
        sums[z][3] = s12;
        pix1 += 4;
        pix2 += 4;
    }
}

static float ssim_end1( int s1, int s2, int ss, int s12 )
{
    static const int ssim_c1 = (int)(.01*.01*255*255*64 + .5);
    static const int ssim_c2 = (int)(.03*.03*255*255*64*63 + .5);
    int vars = ss*64 - s1*s1 - s2*s2;
    int covar = s12*64 - s1*s2;
    return (float)(2*s1*s2 + ssim_c1) * (float)(2*covar + ssim_c2)\
           / ((float)(s1*s1 + s2*s2 + ssim_c1) * (float)(vars + ssim_c2));
}

static float ssim_end4( int sum0[5][4], int sum1[5][4], int width )
{
    int i;
    float ssim = 0.0;
    for( i = 0; i < width; i++ )
        ssim += ssim_end1( sum0[i][0] + sum0[i+1][0] + sum1[i][0] + sum1[i+1][0],
                           sum0[i][1] + sum0[i+1][1] + sum1[i][1] + sum1[i+1][1],
                           sum0[i][2] + sum0[i+1][2] + sum1[i][2] + sum1[i+1][2],
                           sum0[i][3] + sum0[i+1][3] + sum1[i][3] + sum1[i+1][3] );
    return ssim;
}

float x264_pixel_ssim_wxh(
                           uint8_t *pix1, int stride1,
                           uint8_t *pix2, int stride2,
                           int width, int height )
{
    int x, y, z;
    float ssim = 0.0;
	int (*t)[4];
    int (*sum0)[4] = (int (*)[4])x264_alloca(4 * (width/4+3) * sizeof(int));
    int (*sum1)[4] = (int (*)[4])x264_alloca(4 * (width/4+3) * sizeof(int));
    width >>= 2;
    height >>= 2;
    z = 0;
    for( y = 1; y < height; y++ )
    {
        for( ; z <= y; z++ )
        {
            //XCHG( (int (*)[4]), sum0, sum1 );
            t =sum0; sum0= sum1; sum1=t;
            for( x = 0; x < width; x+=2 )
                ssim_4x4x2_core( &pix1[4*(x+z*stride1)], stride1, &pix2[4*(x+z*stride2)], stride2, &sum0[x] );
        }
        for( x = 0; x < width-1; x += 4 )
            ssim += ssim_end4( sum0+x, sum1+x, X264_MIN(4,width-x-1) );
    }
    return ssim / ((width-1) * (height-1));
}


int calculate(int n, char *cl[], ssim_psnr_t *ssim_psnr )
{
  FILE *f1, *f2;
  int ssim = 0, i, x, y, yuv, inc = 1, size = 0, N = 0, Y, F;
  double yrmse, diff, mean_ssim = 0, mean_psnr = 0, stdv_ssim = 0, stdv_psnr = 0,*ypsnr = 0, *yssim = 0;
  unsigned char *b1, *b2;
  clock_t t = clock();

  if (0) {
    puts("psnr x y <YUV format> <src.yuv> <dst.yuv> [multiplex] [ssim]");
    puts("  x\t\tframe width");
    puts("  y\t\tframe height");
    puts("  YUV format\t420, 422, etc.");
    puts("  src.yuv\tsource video");
    puts("  dst.yuv\tdistorted video");
    puts("  [multiplex]\toptional");

    return EXIT_FAILURE;
  }

  if ((f1 = fopen(cl[3], "rb")) == 0) goto A;
  if ((f2 = fopen(cl[4], "rb")) == 0) goto B;
  if (!(x = strtoul(cl[1], 0, 10)) ||
      !(y = strtoul(cl[2], 0, 10))) goto C; 
  //if ((yuv = strtoul(cl[3], 0, 10)) > 444) goto D;
  if (cl[6] && !strcmp(cl[6], "multiplex")) inc = 2;
  if (cl[6] && !strcmp(cl[6], "ssim")) ssim = 1;

  Y = x * y;
  switch (yuv) {
    case 400: F = Y; break;
    case 422: F = Y * 2; break;
    case 444: F = Y * 3; break;
    default :
    case 420: F = Y * 3 / 2; break;
  }
  F = Y * 3 / 2;

  if (!(b1 = (unsigned char*)malloc(F))) goto E;
  if (!(b2 = (unsigned char*)malloc(F))) goto E;

  for (;;) {
    if (1 != fread(b1, F, 1, f1) || 1 != fread(b2, F, 1, f2)) break;

    if (++N > size) {
      size += 0xffff;
	  if (!(yssim = (double *)realloc(yssim, size * sizeof *yssim))) 
	  	goto E;
      if (!(ypsnr = (double *)realloc(ypsnr, size * sizeof *ypsnr))) 
	  	goto E;
    }

    
      mean_ssim += yssim[N - 1] = x264_pixel_ssim_wxh(b1, x, b2, x, x, y);

      for (yrmse = 0, i = inc - 1; i < (inc == 1 ? Y : F); i += inc) {
        diff = b1[i] - b2[i];
        yrmse += diff * diff;
      }
      mean_psnr += ypsnr[N - 1] = yrmse ? 20 * (log10(255 / sqrt(yrmse / Y))) : 0;

    printf("\r ssim %.3f  \t psnr %.3f", yssim[N - 1], ypsnr[N - 1]);
  }

  printf("\n");
  if (N) {
    mean_ssim/= N;

    for (stdv_ssim = 0, i = 0; i < N; i++) {
      diff = yssim[i] - mean_ssim;
      stdv_ssim += diff * diff;
    }
    stdv_ssim = sqrt(stdv_ssim / (N - 1));

    free(yssim);
  }

  if (N) {
    mean_psnr /= N;

    for (stdv_psnr = 0, i = 0; i < N; i++) {
      diff = ypsnr[i] - mean_psnr;
      stdv_psnr += diff * diff;
    }
    stdv_psnr = sqrt(stdv_psnr / (N - 1));

    free(ypsnr);
  }

  fclose(f1);
  fclose(f2);

  ssim_psnr->mean_psnr = mean_psnr;
  ssim_psnr->mean_ssim = mean_ssim;
  ssim_psnr->stdv_ssim = stdv_ssim;
  ssim_psnr->stdv_psnr = stdv_psnr;
  
  fprintf(stderr, "\t%d frames (CPU: %lu s) psnr mean: %.2f psnr stdv: %.2f  ssim mean: %.2f ssim stdv: %.2f, \n",
     N, (unsigned long) ((clock() - t) / CLOCKS_PER_SEC), mean_psnr, stdv_psnr, mean_ssim, stdv_ssim);

  return 0;

A: fprintf(stderr, " Error opening source video file.\n"); goto X;
B: fprintf(stderr, " Error opening decoded video file.\n"); goto X;
C: fprintf(stderr, " Invalid width or height.\n"); goto X;
D: fprintf(stderr, " Invalid YUV format.\n"); goto X;
E: fprintf(stderr, " Not enough memory.\n");

X: return EXIT_FAILURE;

}


void createColorComponent( ColorComponent* cc )
{
  if( ! ( cc->data = (unsigned char *)malloc(cc->width * cc->height)))
  {
    fprintf(stderr, "\nERROR: memory allocation failed!\n\n");
    exit(-1);
  }
}

void deleteColorComponent( ColorComponent* cc )
{
  free(cc->data);
  cc->data = NULL;
}



void createFrame( YuvFrame* f, int width, int height )
{
  f->lum.width = width;    f->lum.height  = height;     createColorComponent( &f->lum );
  f->cb .width = width/2;  f->cb .height  = height/2;   createColorComponent( &f->cb  );
  f->cr .width = width/2;  f->cr .height  = height/2;   createColorComponent( &f->cr  );
}

void deleteFrame( YuvFrame* f )
{
  deleteColorComponent( &f->lum );
  deleteColorComponent( &f->cb  );
  deleteColorComponent( &f->cr  );
}

void readColorComponent( ColorComponent* cc, FILE* file )
{
  unsigned int size   = cc->width*cc->height;
  unsigned int rsize;

  rsize = (unsigned int)fread( cc->data, sizeof(unsigned char), size, file );

  if( size != rsize )
  {
    fprintf(stderr, "\nERROR: while reading from input file!\n\n");
    exit(-1);
  }
}

void writeColorComponent( ColorComponent* cc, FILE* file, int downScale )
{
  int outwidth  = cc->width   >> downScale;
  int outheight = cc->height  >> downScale;
  int wsize;

  for( int i = 0; i < outheight; i++ )
  {
    wsize = (int)fwrite( cc->data+i*cc->width, sizeof(unsigned char), outwidth, file );

    if( outwidth != wsize )
    {
      fprintf(stderr, "\nERROR: while writing to output file!\n\n");
      exit(-1);
    }
  }
}

double psnr( ColorComponent& rec, ColorComponent& org)
{
  unsigned char*  pOrg  = org.data;
  unsigned char*  pRec  = rec.data;
  double          ssd   = 0;
  int             diff;

  for  ( int r = 0; r < rec.height; r++ )
  {
    for( int c = 0; c < rec.width;  c++ )
    {
      diff  = pRec[c] - pOrg[c];
      ssd  += (double)( diff * diff );
    }
    pRec   += rec.width;
    pOrg   += org.width;
  }

  if( ssd == 0.0 )
  {
    return 99.99;
  }
  return ( 10.0 * log10( (double)rec.width * (double)rec.height * 65025.0 / ssd ) );
}

void getPSNR( double& psnrY, double& psnrU, double& psnrV, YuvFrame& rcFrameOrg, YuvFrame& rcFrameRec )
{
  psnrY = psnr( rcFrameRec.lum, rcFrameOrg.lum );
  psnrU = psnr( rcFrameRec.cb,  rcFrameOrg.cb  );
  psnrV = psnr( rcFrameRec.cr,  rcFrameOrg.cr  );
}

void readFrame( YuvFrame* f, FILE* file )
{
  readColorComponent( &f->lum, file );
  readColorComponent( &f->cb,  file );
  readColorComponent( &f->cr,  file );
}

void print_usage_and_exit( int test, const char* name, const char* message = 0 )
{
  if( test )
  {
    if( message )
    {
      fprintf ( stderr, "\nERROR: %s\n", message );
    }
    fprintf (   stderr, "\nUsage: %s <w> <h> <org> <rec> [<t> [<skip> [<strm> <fps> [strg]]]] [-r]\n\n", name );
    fprintf (   stderr, "\t    w : original width  (luma samples)\n" );
    fprintf (   stderr, "\t    h : original height (luma samples)\n" );
    fprintf (   stderr, "\t  org : original file\n" );
    fprintf (   stderr, "\t  rec : reconstructed file\n" );
    fprintf (   stderr, "\t    t : number of temporal downsampling stages (default: 0)\n" );
    fprintf (   stderr, "\t skip : number of frames to skip at start      (default: 0)\n" );
    fprintf (   stderr, "\t strm : coded stream\n" );
    fprintf (   stderr, "\t fps  : frames per second\n" );
    fprintf (   stderr, "\t strg : prefix string for summary output\n" );
	  fprintf (   stderr, "\t -r   : return Luma psnr (default: return -1 when failed and 0 otherwise)\n" );
    fprintf (   stderr, "\n" );
    exit    (   -1 );
  }
}



int main(int argc, char *argv[])
{
  int     acc = 10000;
#define   OUT "%d,%04d"

  //===== input parameters =====
  int           stream          = 0;
  unsigned int  width           = 0;
  unsigned int  height          = 0;
  unsigned int  temporal_stages = 0;
  unsigned int  skip_at_start   = 0;
  double        fps             = 0.0;
  FILE*         org_file        = 0;
  FILE*         rec_file        = 0;
  FILE*         str_file        = 0;
  char*         prefix_string   = 0;

  //===== variables =====
  unsigned int  index, skip, skip_between, sequence_length;
  int           py, pu, pv, br;
  double        bitrate = 0.0;
  double        psnrY, psnrU, psnrV;
  YuvFrame      cOrgFrame, cRecFrame;
  double        AveragePSNR_Y = 0.0;
  double        AveragePSNR_U = 0.0;
  double        AveragePSNR_V = 0.0;
  int		      	currarg = 5;
  int			      rpsnr   = 0;

  ssim_psnr_t ssim_psnr;
  memset(&ssim_psnr, 0, sizeof(ssim_psnr_t));


  //===== read input parameters =====
  print_usage_and_exit((argc < 5 || (argc > 11 )), argv[0]);
  width             = atoi  ( argv[1] );
  height            = atoi  ( argv[2] );
  org_file          = fopen ( argv[3], "rb" );
  rec_file          = fopen ( argv[4], "rb" );

  if(( argc >=  6 ) && strcmp( argv[5], "-r" ) )
  {
    temporal_stages = atoi  ( argv[5] );
    currarg++;
  }
  if(( argc >=  7 ) && strcmp( argv[6], "-r" ) )
  {
    skip_at_start   = atoi  ( argv[6] );
    currarg++;
  }
  if(( argc >= 9 ) && strcmp( argv[7], "-r" ) )
  {
    str_file        = fopen ( argv[7], "rb" );
	  print_usage_and_exit(!strcmp( argv[8], "-r" ), argv[0]);
	  fps             = atof  ( argv[8] );
    stream          = 1;
	  currarg+=2;
  }
  if(( argc >= 10 ) && strcmp( argv[9], "-r" ) )
  {
    prefix_string   = argv[9];
    currarg++;
  }

	if(currarg < argc )
	{
	  if(!strcmp( argv[currarg], "-r" ))
		  rpsnr=1;
	  else
      print_usage_and_exit (true,argv[0],"Wrong number of argument!" );
	}


  //===== check input parameters =====
  print_usage_and_exit  ( ! org_file,                                       argv[0], "Cannot open original file!" );
  print_usage_and_exit  ( ! rec_file,                                       argv[0], "Cannot open reconstructed file!" );
  print_usage_and_exit  ( ! str_file && stream,                             argv[0], "Cannot open stream!" );
  print_usage_and_exit  ( fps <= 0.0 && stream,                             argv[0], "Unvalid frames per second!" );

  //======= get number of frames and stream size =======
  fseek(    rec_file, 0, SEEK_END );
  fseek(    org_file, 0, SEEK_END );
  size_t rsize = ftell( rec_file );
  size_t osize = ftell( org_file );
  fseek(    rec_file, 0, SEEK_SET );
  fseek(    org_file, 0, SEEK_SET );

  if (rsize < osize)
    sequence_length = (unsigned int)((double)rsize/(double)((width*height*3)/2));
   else
    sequence_length = (unsigned int)((double)osize/(double)((width*height*3)/2));

  if( stream )
  {
    fseek(  str_file, 0, SEEK_END );
    bitrate       = (double)ftell(str_file) * 8.0 / 1000.0 / ( (double)(sequence_length << temporal_stages) / fps );
    fseek(  str_file, 0, SEEK_SET );
  }
  skip_between    = ( 1 << temporal_stages ) - 1;

  //===== initialization ======
  createFrame( &cOrgFrame, width, height );
  createFrame( &cRecFrame, width, height );

  //===== loop over frames =====
  for( skip = skip_at_start, index = 0; index < sequence_length; index++, skip = skip_between )
  {
    fseek( org_file, skip*width*height*3/2, SEEK_CUR);

    readFrame       ( &cOrgFrame, org_file );
    readFrame       ( &cRecFrame, rec_file );

    getPSNR         ( psnrY, psnrU, psnrV, cOrgFrame, cRecFrame);
    AveragePSNR_Y +=  psnrY;
    AveragePSNR_U +=  psnrU;
    AveragePSNR_V +=  psnrV;

    py = (int)floor( acc * psnrY + 0.5 );
    pu = (int)floor( acc * psnrU + 0.5 );
    pv = (int)floor( acc * psnrV + 0.5 );
    fprintf(stdout,"\r%d\t"OUT"\t"OUT"\t"OUT"",index,py/acc,py%acc,pu/acc,pu%acc,pv/acc,pv%acc);
  }
  fprintf(stdout,"\n");

	calculate(argc, argv, &ssim_psnr);
	
  py = (int)floor( acc * AveragePSNR_Y / (double)sequence_length + 0.5 );
  pu = (int)floor( acc * AveragePSNR_U / (double)sequence_length + 0.5 );
  pv = (int)floor( acc * AveragePSNR_V / (double)sequence_length + 0.5 );
  br = (int)floor( acc * bitrate                                 + 0.5 );
  if( stream )
  {
    if( prefix_string )
    {
      fprintf(stderr,"%s\t"OUT"\t"OUT"\t"OUT"\t"OUT"\n",prefix_string,br/acc,br%acc,py/acc,py%acc,pu/acc,pu%acc,pv/acc,pv%acc);
      fprintf(stdout,"%s\t"OUT"\t"OUT"\t"OUT"\t"OUT"\n",prefix_string,br/acc,br%acc,py/acc,py%acc,pu/acc,pu%acc,pv/acc,pv%acc);
    }
    else
    {
      fprintf(stderr,OUT"\t"OUT"\t"OUT"\t"OUT"\n",br/acc,br%acc,py/acc,py%acc,pu/acc,pu%acc,pv/acc,pv%acc);
      fprintf(stdout,OUT"\t"OUT"\t"OUT"\t"OUT"\n",br/acc,br%acc,py/acc,py%acc,pu/acc,pu%acc,pv/acc,pv%acc);
    }
  }
  else
  {
    //fprintf(stderr,"total\t"OUT"\t"OUT"\t"OUT"\n",py/acc,py%acc,pu/acc,pu%acc,pv/acc,pv%acc);
    fprintf(stdout,"total\t"OUT"\t"OUT"\t"OUT"  psnr stdv: %.3f  ssim mean: %.3f ssim stdv: %.3f\n",py/acc,py%acc,pu/acc,pu%acc,pv/acc,pv%acc, ssim_psnr.stdv_psnr, ssim_psnr.mean_ssim, ssim_psnr.stdv_ssim);
  }

  fprintf(stdout, "\n");


  //===== finish =====
  deleteFrame( &cOrgFrame );
  deleteFrame( &cRecFrame );
  fclose     ( org_file   );
  fclose     ( rec_file   );
  if( stream )
  {
    fclose   ( str_file   );
  }

  return (rpsnr*py);
}

