/**
 *
 *  This code is part of the publication:
 *
 *      "Contours, corners and T-junctions detection algorithm" by Antoni Buades,
 *      Rafael Grompone von Gioi, and Julia Navarro. Image Processing On Line, 2018.
 *
 *  Antoni Buades <toni.buades@uib.es> Universitat de les Illes Balears, Spain
 *  Rafael Grompone von Gioi <grompone@gmail.com> CMLA, ENS Cachan, France
 *  Julia Navarro <julia.navarro@uib.es> Universitat de les Illes Balears, Spain
 *
 *
 **/


#include "libDetectionIpol.h"




/*----------------------------------------------------------------------------*/
/** ln(10) */
#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif /* !M_LN10 */

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/*----------------------------------------------------------------------------*/
/** Fatal error, print a message to standard-error output and exit.
 */
static void error(char * msg)
{
    fprintf(stderr,"error: %s\n",msg);
    exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/** Doubles relative error factor
 */
#define RELATIVE_ERROR_FACTOR 100.0

/*----------------------------------------------------------------------------*/
/** Compare doubles by relative error.

 The resulting rounding error after floating point computations
 depend on the specific operations done. The same number computed by
 different algorithms could present different rounding errors. For a
 useful comparison, an estimation of the relative rounding error
 should be considered and compared to a factor times EPS. The factor
 should be related to the cumulated rounding error in the chain of
 computation. Here, as a simplification, a fixed factor is used.
 */
static int double_equal(double a, double b)
{
    double abs_diff,aa,bb,abs_max;

    /* trivial case */
    if( a == b ) return TRUE;

    abs_diff = fabs(a-b);  //求绝对值 。 abs fabs
    aa = fabs(a);
    bb = fabs(b);
    abs_max = aa > bb ? aa : bb;

    /* DBL_MIN is the smallest normalized number, thus, the smallest
     number whose relative error is bounded by DBL_EPSILON. For
     smaller numbers, the same quantization steps as for DBL_MIN
     are used. Then, for smaller numbers, a meaningful "relative"
     error should be computed by dividing the difference by DBL_MIN. */
    if( abs_max < DBL_MIN ) abs_max = DBL_MIN;

    /* equal if relative error <= factor x eps */
    return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);  //比较运算符
}

/*----------------------------------------------------------------------------*/
/** Computes the natural logarithm of the absolute value of
 the gamma function of x using the Lanczos approximation.
 See http://www.rskey.org/gamma.htm

 The formula used is
 @f[
 \Gamma(x) = \frac{ \sum_{n=0}^{N} q_n x^n }{ \Pi_{n=0}^{N} (x+n) }
 (x+5.5)^{x+0.5} e^{-(x+5.5)}
 @f]
 so
 @f[
 \log\Gamma(x) = \log\left( \sum_{n=0}^{N} q_n x^n \right)
 + (x+0.5) \log(x+5.5) - (x+5.5) - \sum_{n=0}^{N} \log(x+n)
 @f]
 and
 q0 = 75122.6331530,
 q1 = 80916.6278952,
 q2 = 36308.2951477,
 q3 = 8687.24529705,
 q4 = 1168.92649479,
 q5 = 83.8676043424,
 q6 = 2.50662827511.
 */
static double log_gamma_lanczos(double x)
{
    static double q[7] = { 75122.6331530, 80916.6278952, 36308.2951477,
        8687.24529705, 1168.92649479, 83.8676043424,
        2.50662827511 };
    double a = (x+0.5) * log(x+5.5) - (x+5.5);
    double b = 0.0;
    int n;

    for(n=0;n<7;n++)
    {
        a -= log( x + (double) n );
        b += q[n] * pow( x, (double) n );
    }
    return a + log(b);
}

/*----------------------------------------------------------------------------*/
/** Computes the natural logarithm of the absolute value of
 the gamma function of x using Windschitl method.
 See http://www.rskey.org/gamma.htm

 The formula used is
 @f[
 \Gamma(x) = \sqrt{\frac{2\pi}{x}} \left( \frac{x}{e}
 \sqrt{ x\sinh(1/x) + \frac{1}{810x^6} } \right)^x
 @f]
 so
 @f[
 \log\Gamma(x) = 0.5\log(2\pi) + (x-0.5)\log(x) - x
 + 0.5x\log\left( x\sinh(1/x) + \frac{1}{810x^6} \right).
 @f]
 This formula is a good approximation when x > 15.
 */
static double log_gamma_windschitl(double x)
{
    return 0.918938533204673 + (x-0.5)*log(x) - x
    + 0.5*x*log( x*sinh(1/x) + 1/(810.0*pow(x,6.0)) );
}

/*----------------------------------------------------------------------------*/
/** Computes the natural logarithm of the absolute value of
 the gamma function of x. When x>15 use log_gamma_windschitl(),
 otherwise use log_gamma_lanczos().
 */
#define log_gamma(x) ((x)>15.0?log_gamma_windschitl(x):log_gamma_lanczos(x))

/*----------------------------------------------------------------------------*/
/** Size of the table to store already computed inverse values.
 */
#define TABSIZE 100000

/*----------------------------------------------------------------------------*/
/** Computes -log10(NFA).

 NFA stands for Number of False Alarms:
 @f[
 \mathrm{NFA} = NT \cdot B(n,k,p)
 @f]

 - NT       - number of tests
 - B(n,k,p) - tail of binomial distribution with parameters n,k and p:
 @f[
 B(n,k,p) = \sum_{j=k}^n
 \left(\begin{array}{c}n\\j\end{array}\right)
 p^{j} (1-p)^{n-j}
 @f]

 The value -log10(NFA) is equivalent but more intuitive than NFA:
 - -1 corresponds to 10 mean false alarms
 -  0 corresponds to 1 mean false alarm
 -  1 corresponds to 0.1 mean false alarms
 -  2 corresponds to 0.01 mean false alarms
 -  ...

 Used this way, the bigger the value, better the detection,
 and a logarithmic scale is used.

 @param n,k,p binomial parameters.
 @param logNT logarithm of Number of Tests

 The computation is based in the gamma function by the following
 relation:
 @f[
 \left(\begin{array}{c}n\\k\end{array}\right)
 = \frac{ \Gamma(n+1) }{ \Gamma(k+1) \cdot \Gamma(n-k+1) }.
 @f]
 We use efficient algorithms to compute the logarithm of
 the gamma function.

 To make the computation faster, not all the sum is computed, part
 of the terms are neglected based on a bound to the error obtained
 (an error of 10% in the result is accepted).
 */
static double nfa(int n, int k, double p, double logNT)
{
    static double inv[TABSIZE];   /* table to keep computed inverse values */
    double tolerance = 0.1;       /* an error of 10% in the result is accepted */
    double log1term,term,bin_term,mult_term,bin_tail,err,p_term;
    int i;

    /* check parameters */
    if( n<0 || k<0 || k>n || p<=0.0 || p>=1.0 )
    error("nfa: wrong n, k or p values.");

    /* trivial cases */
    if( n==0 || k==0 ) return -logNT;
    if( n==k ) return -logNT - (double) n * log10(p);

    /* probability term */
    p_term = p / (1.0-p);

    /* compute the first term of the series */
    /*
     binomial_tail(n,k,p) = sum_{i=k}^n bincoef(n,i) * p^i * (1-p)^{n-i}
     where bincoef(n,i) are the binomial coefficients.
     But
     bincoef(n,k) = gamma(n+1) / ( gamma(k+1) * gamma(n-k+1) ).
     We use this to compute the first term. Actually the log of it.
     */
    log1term = log_gamma( (double) n + 1.0 ) - log_gamma( (double) k + 1.0 )
    - log_gamma( (double) (n-k) + 1.0 )
    + (double) k * log(p) + (double) (n-k) * log(1.0-p);
    term = exp(log1term);

    /* in some cases no more computations are needed */
    if( double_equal(term,0.0) )              /* the first term is almost zero */
    {
        if( (double) k > (double) n * p )     /* at begin or end of the tail?  */
        return -log1term / M_LN10 - logNT;  /* end: use just the first term  */
        else
        return -logNT;                      /* begin: the tail is roughly 1  */
    }

    /* compute more terms if needed */
    bin_tail = term;
    for(i=k+1;i<=n;i++)
    {
        /*
         As
         term_i = bincoef(n,i) * p^i * (1-p)^(n-i)
         and
         bincoef(n,i)/bincoef(n,i-1) = n-1+1 / i,
         then,
         term_i / term_i-1 = (n-i+1)/i * p/(1-p)
         and
         term_i = term_i-1 * (n-i+1)/i * p/(1-p).
         1/i is stored in a table as they are computed,
         because divisions are expensive.
         p/(1-p) is computed only once and stored in 'p_term'.
         */
        bin_term = (double) (n-i+1) * ( i<TABSIZE ?
                                       ( inv[i]!=0.0 ? inv[i] : ( inv[i] = 1.0 / (double) i ) ) :
                                       1.0 / (double) i );

        mult_term = bin_term * p_term;
        term *= mult_term;
        bin_tail += term;
        if(bin_term<1.0)
        {
            /* When bin_term<1 then mult_term_j<mult_term_i for j>i.
             Then, the error on the binomial tail when truncated at
             the i term can be bounded by a geometric series of form
             term_i * sum mult_term_i^j.                            */
            err = term * ( ( 1.0 - pow( mult_term, (double) (n-i+1) ) ) /
                          (1.0-mult_term) - 1.0 );

            /* One wants an error at most of tolerance*final_result, or:
             tolerance * abs(-log10(bin_tail)-logNT).
             Now, the error that can be accepted on bin_tail is
             given by tolerance*final_result divided by the derivative
             of -log10(x) when x=bin_tail. that is:
             tolerance * abs(-log10(bin_tail)-logNT) / (1/bin_tail)
             Finally, we truncate the tail if the error is less than:
             tolerance * abs(-log10(bin_tail)-logNT) * bin_tail        */
            if( err < tolerance * fabs(-log10(bin_tail)-logNT) * bin_tail ) break;
        }
    }
    return -log10(bin_tail) - logNT;
}

/*----------------------------------------------------------------------------*/
float probability(float den, float * val, float * acc, int N)
{
    float p = 1.0f;

    if (den <= val[0]) p = 1.0f;
    else if (den >= val[N-1]) p = 0.0f;
    else
    {
        int i = 1;
        for(; i < N && val[i] < den; i++);
        p = acc[i-1] + (acc[i] - acc[i-1]) * (den-val[i-1])/ (val[i]-val[i-1]);
        p = 1.0 - p;
    }

    return p;
}





//!
//!     This function computes the convolution of the image with several directional filters.
//!
void compute_directional_convolutions(float *input, float **dimages, int nDir, float fSigma, float anglePrecision, int width, int height)
{

    float sigma=fSigma;

    //! Build Gaussian filters for horizontal and vertical derivatives
    int ksize = 8 * (int) ceilf(sigma) + 1;   //取整     //或者说sigma 不能大于3  因为大于三 kernal 大于25 的部分就是乱序了
    int l2 = ksize/2;


    float *xkernelx = new float[ksize];  //横向响应的卷积。  ksize= 8*2+1 =17 sigma=2.0      但是初始化的时候要求大于24 所以改成了3
    float *xkernely = new float[ksize];  //纵向相应的卷积	 ksize= 8*3+1 =25

    float *xconvolved = new float[width*height];     //原图大小的空间 要存什么呢？？？？？
    float *yconvolved = new float[width*height];	//构建响应图



    float sigma2 = sigma*sigma;     //2 取平方

    for(int p = -l2; p <= l2; p++)      // p+12  所以好像 ksize 最少要24吧。。。。。   sigma》= 3.0    是kernel的初始化  
    {

        // \frac{\nabla}{\nabla x} e^{(x^2 + y^2)/(2 \sigma ^2)} = - (x / \sigma^2)  * e^{(x^2 + y^2)/(2 \sigma ^2)}
        xkernelx[p + l2] = - sigma * (1.0 / (2.0*PI*sigma2)) * (p/sigma2) * expf(- ((float)p*p)/(2.0f*sigma2));
        xkernely[p + l2] =  expf(- ((float)p*p)/(2.0f*sigma2));

        /* additional sigma, normalization to be comparable between scales,
         see p.38 of Mikolajczyk PhD thesis, 2002. */

    }

    libUSTG::fiSepConvol(input, xconvolved, width, height, xkernelx, ksize, xkernely, ksize, BOUNDARY_CONDITION_SYMMETRIC);
    libUSTG::fiSepConvol(input, yconvolved, width, height, xkernely, ksize, xkernelx, ksize, BOUNDARY_CONDITION_SYMMETRIC);


    //! Combine for computing directional convolutions
    for (int i=0; i < nDir; i++)
    {

        float a = (float) i * (float) anglePrecision  * PI / 180.0f;
        float sina = sin(a);
        float cosa = cos(a);

        for(int j = 0 ; j < width*height; j++)
        {
            dimages[i][j] = MAX(cosa * xconvolved[j] + sina * yconvolved[j] , 0.0f);
        }

    }



    //! delete memory
    delete[] xkernelx;
    delete[] xkernely;
    delete[] xconvolved;
    delete[] yconvolved;


}


//!
//! This function identifies meaningful points to be candidates for contours by means of lateral inhibition.
//!
//! X,Y : width and height
//! n: number of direction
//! fThreshold: threshold on magnitude of response to be taken into account
//! fAngleFiltPrec:  0 (all 90 degrees permit to inhibit), else  fAngleFiltPrec permits to inhibit
void lateral_inhibition_each_scale(float **g, float **gh, float sigma, int n, int X, int Y, float fAngleFiltPrec, float fThreshold, float anglePrec, int flagKeepValue)
{

    //! tabulate exp(-x), faster than using directly function expf
    int iLutLength = (int) rintf((float) LUTMAX * (float) LUTPRECISION); //30*29
    float *fpLut = new float[iLutLength];
    libUSTG::wxFillExpLut(fpLut,iLutLength);  //初始化wxFill

    int ksize = (int) rintf(4.0 * sigma);     // 定义ksize = 4*3 =12

    float fAux = 0.5f / (4.0*sigma*sigma);    // 很小的一个值 0.5/4*12*12

    //! Precompute spatial neighborhoods
    int *xMin = new int [X];      //x = 原图的宽   //序号与ksize 的min
    int *xMax = new int [X];

    for(int x=0; x<X; x++)
    {
        xMin[x] = MAX(0,x-ksize);
        xMax[x] = MIN(X-1,x+ksize);
    }


    int *yMin = new int [Y];
    int *yMax = new int [Y];

    for(int y=0; y<Y; y++)
    {
        yMin[y] = MAX(0,y-ksize);
        yMax[y] = MIN(Y-1,y+ksize);
    }


    //! Precompute exponential
    float *tableW = new float[2*ksize*ksize+1];       //2*12*12+1

    for (int iw=0; iw<2*ksize*ksize+1; iw++)
    {

        //!float w = exp( -0.5*(ox*ox + oy*oy) / (4.0*sigma*sigma) );
        float fValue=fAux * iw;    //0.5/4*12*12   *2*12*12

        tableW[iw] = libUSTG::wxSLUT(fValue,fpLut);

    }


    //! Precompute angular neighborhood
    int **tableJ = new int*[n];
    int *tableJcont = new int[n];
    for(int i=0; i<n; i++)
    {
        tableJ[i] = new int[n];
        int totalj = 0;

        for(int j=0; j<n; j++)
        {

            int ia=abs(i-j);
            if (abs(i-j+n) < ia) ia = abs(i-j+n);
            if (abs(i-j-n) < ia) ia = abs(i-j-n);

            float difA = (float) ia * anglePrec;

            if (difA <= fAngleFiltPrec)
            {
                tableJ[i][totalj] = j;
                totalj++;
            }
        }

        tableJcont[i] = totalj;
    }



//#pragma omp parallel for shared(fpLut, fAux, fThreshold, g, gh, xMin, xMax, yMin, yMax, tableJ, tableJcont, tableW)
    for(int i=0; i<n; i++)
    {

        float theta = i * 2.0 * M_PI / n;
        float dx = -sin(theta);
        float dy =  cos(theta);

        for(int x=0; x<X; x++)
        for(int y=0; y<Y; y++)
        if (g[i][x+y*X] > fThreshold)
        {

            if (!flagKeepValue)
                gh[i][x+y*X] = 1.0;
            else
                gh[i][x+y*X] = g[i][x+y*X];

            for(int xx=xMin[x]; xx<=xMax[x]; xx++)          //! and close pixels
            for(int yy=yMin[y]; yy<=yMax[y]; yy++)
            {
                float ox =  (xx-x) * dx + (yy-y) * dy;
                float oy = -(xx-x) * dy + (yy-y) * dx;

                //! Only directions at 45 degrees difference are able to inhibitate
                if( fabs(oy) >= fabs(ox) )
                {

                    float w = tableW[(int) rintf(ox*ox+oy*oy)];

                    int auxi = i;
                    int auxx=x+y*X;
                    int auxxx = xx+yy*X;

                    for(int jj=0; jj<tableJcont[i]; jj++)
                    {
                        int j = tableJ[i][jj];

                        if( g[auxi][auxx] < w*g[j][auxxx] )
                        gh[auxi][auxx] = 0.0;

                    }

                }
            }



        }
        else
        {
            gh[i][x+y*X] = 0.0;

        }


    }





    delete[] xMin;
    delete[] xMax;
    delete[] yMin;
    delete[] yMax;
    delete[] fpLut;


    delete[] tableW;
    delete[] tableJcont;
    for (int i=0; i < n; i++) delete[] tableJ[i];
    delete[]  tableJ;




}


//!
//! This function performs a filtering to reinforce points in the cube that are linked by the association field.
//!
//! Convolution of each slide with (sigma, 1.0) directional kernel
//!     float fAngleFiltPrec = 15.0f;
//!     float fAngleFiltStd  = 250.0f;  That's uniform average!
//!     float sigma  = Sigma;
//!     float sigmay = 1.0f;
void good_continuation_filter(float **input, float **output, float Sigma, float fAngleFiltPrec, int nDir, float anglePrec, int width, int height)
{

    float **tmpFiltered = new float*[nDir];
    for (int i=0; i < nDir; i++)  tmpFiltered[i] = new float[width*height];

    float sigma  = Sigma;
    float sigmay = 1.0f;


    #pragma omp parallel for shared(input, tmpFiltered, sigma, sigmay)
    for(int iDir=0; iDir < nDir; iDir++)
    {


        int kwidth = 2 * (int) rintf(3.0 * sigma) + 1;
        int klength = kwidth*kwidth;
        float *dirKernel = new float[klength];
        for (int ii=0; ii < klength; ii++) dirKernel[ii]=0.0f;

        float thetaDeg=0.0f;
        thetaDeg = (float) iDir * 360.0f / (float) nDir;

        libUSTG::fiFloatDirectionalGaussKernel(sigmay, sigma, thetaDeg, dirKernel,  kwidth, kwidth);
        libUSTG::fiConvol(input[iDir], tmpFiltered[iDir], width, height, dirKernel, kwidth, kwidth, BOUNDARY_CONDITION_SYMMETRIC);

        delete[] dirKernel;


    }

    //! Filter close orientations
    //! tabulate exp(-x), faster than using directly function expf
    int iLutLength = (int) rintf((float) LUTMAX * (float) LUTPRECISION);   //==30*29
    float *fpLut = new float[iLutLength];
    libUSTG::wxFillExpLut(fpLut,iLutLength);    //初始化wxF



    float fAngleFiltStd  = 250.0f;
    for (int i=0; i < nDir; i++)
    for(int ll=0; ll < width*height;ll++)
    output[i][ll] = 0.0f;



    #pragma omp parallel for shared(output, tmpFiltered, anglePrec, fAngleFiltPrec, fpLut)
    for(int iDir=0; iDir < nDir; iDir++)

    for(int jDir=0; jDir < nDir; jDir++)
    {
        float difA = fabsf((float) iDir-jDir) * anglePrec;
        difA = MIN(difA, fabsf((float) (iDir-jDir+nDir)) * anglePrec);
        difA = MIN(difA, fabsf((float) (iDir-jDir-nDir)) * anglePrec);

        if (difA <= fAngleFiltPrec)
        {
            float fWeight = libUSTG::wxSLUT(difA*difA/fAngleFiltStd, fpLut);

            for (int ll = 0; ll < width * height; ll++)
            {
                output[iDir][ll] += fWeight * tmpFiltered[jDir][ll];
            }

        }

    }



    for (int i=0; i < nDir; i++) delete[] tmpFiltered[i];
    delete[]  tmpFiltered;
    delete[]  fpLut;

}



//!
//! This function detects high curvature points to identify corners.
//!
//!  float minResponse = 0.01;
//!  float fCorThresholdv = fCorThreshold / sigma;
//!  libUSTG::fiImageDrawCircle(localMax, ii, jj, 2.0f*sigma, 128.0f, width, height);
void compute_corners(float **dimages, float *localMax, vector<str_point_descriptor> &cCorners, int nDir, float fSigma, float minRespCorner, float fCorThreshold, float radForZone, int width, int height)
{
	//用了dimages 来找roi的角点。 找到roi图后在做真正的角点处理。把结果 最后保存在cCorners里。 dimages没有变化。。

    float sigma = fSigma;

    //! corner response and convolved one
    float *corners  = new float[width*height];
    libUSTG::fpClear(corners,0.0f, width*height);


    //! sum number of responses for each pixel
    for (int kk=0; kk < nDir; kk++)
    for (int jj=0; jj < height; jj++)
    for (int ii=0; ii < width; ii++)
    if (dimages[kk][jj * width + ii] >  minRespCorner)      //minRespCorner=0.005
    {
        corners[jj * width + ii]++;
    }

    //! Convolve response
    libUSTG::fiGaussianConvol(corners, corners, width, height, sigma, BOUNDARY_CONDITION_SYMMETRIC);


    //! Look for local maxima
    float fCorThresholdv = fCorThreshold / sigma;  // 0.275*ndir /3

    libUSTG::fpClear(localMax, 255.0, width*height);  //其实之前已经做了这个工作
    for (int jj=1; jj < height-1; jj++)
    for (int ii=1; ii < width-1; ii++)    //像素点做遍历
    if (corners[ jj * width + ii] > fCorThresholdv) //corner已经又被激活了一次了 水平和垂直
    {

        int imax = 1;
        for (int ik=-1; ik<=1; ik++)
        for (int jk=-1; jk<=1 && imax==1; jk++)   //9邻域折腾
        if ( (abs(ik)>0 || abs(jk)>0) && corners[ (jj+jk) * width + ii + ik] >= corners[ jj * width + ii] )
        imax = 0;

        if (imax == 1)
        {

            libUSTG::fiImageDrawCircle(localMax, ii, jj, radForZone, 128.0f, width, height);

            str_point_descriptor pp;
            pp.px = ii;
            pp.py = jj;

            cCorners.push_back(pp);
        }

    }

    delete[] corners;




}




//!  通过反向验证确定边（曲线） 好像这之后就都不需要了。
//! This function performs the a-contrario validation to automatically select the relevant curves.
//!
void grouping_gestalt_binomial(float **dimages, float **curvesID, float *localMax, float fThreshold, std::vector < std::vector <str_point_descriptor> > &llista, int width, int height, int nDir, char *filename, float ntestsM, float lnfaThr, int connectivity)
{


    for (int ii=0; ii < nDir; ii++)
    libUSTG::fpClear(curvesID[ii], 0.0, width*height);

    //! Create mask image, 1 already visited
    float **auximages = new float*[nDir];
    for (int i=0; i < nDir; i++)
    {
        auximages[i]=new float[width*height];
        libUSTG::fpClear(auximages[i], 0.0f, width*height);
    }

    //! probability
    int nprob;
    std::ifstream f(filename);
    if (!f.good())
    {
        std::cout << "input file " << filename << "not found or impossible to open\n";
    }

    f >> nprob;

    float *fpx = new float[nprob];
    float *fpy = new float[nprob];
    float *afpy = new float[nprob];
    for(int ii = 0; ii < nprob; ii++)
    {
        f >> fpx[ii] >> fpy[ii];
    }

    float sum = 0.0f;
    for (int ii=0; ii < nprob; ii++)
    {
        sum += fpy[ii];
        afpy[ii] = sum;

    }

    //! Visit whole volume
    int ncurves = 0;
    for (int kk=0; kk < nDir; kk++)
    for (int jj=0; jj < height; jj++)
    for (int ii=0; ii < width; ii++)
    if (dimages[kk][jj * width + ii] > fThreshold && localMax[jj * width + ii] == 255.0f && auximages[kk][jj * width + ii] == 0.0f)
    {


        ncurves++;

        //! Begin new curve
        vector<str_point_descriptor> corba;

        //! Mark as visited
        str_point_descriptor pp;

        pp.po = kk;
        pp.px = ii;
        pp.py = jj;
        pp.den = dimages[kk][jj * width + ii];

        auximages[kk][jj * width + ii] = 1.0f;

        corba.push_back(pp);

        //! Fill corba while it can grow
        for (int ll=0; ll < (int) corba.size(); ll++)
        {
            //! Point in corbe
            int po = corba[ll].po;
            int px = corba[ll].px;
            int py = corba[ll].py;

            //! Check neighbors
            int ilength =  connectivity;
            for (int npo = -ilength; npo <= ilength; npo++)
            for (int npy = -ilength; npy <= ilength; npy++)
            for (int npx = -ilength; npx <= ilength; npx++)
            {
                int nnpx = px + npx;
                int nnpy = py + npy;
                int nnpo = (nDir + po + npo) % nDir;

                //! Add neighbor if non visited and magnited larger than t
                if (nnpx >= 0 && nnpx < width && nnpy >= 0 && nnpy < height && localMax[nnpy * width + nnpx] == 255.0f && dimages[nnpo][nnpy * width + nnpx] > fThreshold && auximages[nnpo][nnpy * width + nnpx] == 0.0f)
                {
                    auximages[nnpo][nnpy * width + nnpx] = 1.0f;

                    pp.po = nnpo;
                    pp.px = nnpx;
                    pp.py = nnpy;
                    pp.den = dimages[nnpo][nnpy * width + nnpx];

                    corba.push_back(pp);
                }
            }
        }


        double N = (double) width * (double) height;
        N *= N;
        N *= ntestsM;


        N = 5.0 * N; // count different values of p tested, about 10

        for( float p = 0.5; p >= 0.1; p -= 0.1 )
        {
            int n = corba.size();

            int k = 0;
            for(int i=0; i<(int)corba.size(); i++)
            if( probability(corba[i].den,fpx,afpy,nprob) <= p ) ++k;

            double lnfa = nfa( n, k, p, log10(N) );

            ///if( lnfa > 0.0 )
            if( lnfa > lnfaThr )
            {
                llista.push_back(corba);
                for (int i=0; i < (int) corba.size(); i++)
                {
                    curvesID[corba[i].po][corba[i].py * width + corba[i].px] = (float) llista.size();
                }
                break; // already validated, no need to test another value of p
            }
        }

    }

    for (int i=0; i < nDir; i++) delete[] auximages[i];
    delete[]  auximages;

    delete[]  fpx;
    delete[]  fpy;
    delete[]  afpy;
}







//!  只有至少两个不同方向的曲线经过此点，才被保留
//! This function retains detected corners are retained only when at least two detected curves of different orientation join them.
//!
void clear_corners(vector<str_point_descriptor> &cCorners, float *localMax, float **curvesID, float fSigma, int nDir, float anglePrecision, float rMult, float radForZone, float fCorValidationDegree, int width, int height)
{



    vector<str_point_descriptor> cCorners2;
    float sigma = fSigma;
    libUSTG::fpClear(localMax, 255.0f, width*height);


    //! for each corner
    for (int kk=0; kk < (int) cCorners.size(); kk++)
    {

        int px = cCorners[kk].px;
        int py = cCorners[kk].py;


        float radius = rMult *sigma;
        radius*=radius;

        int ksize = (int) ceilf(rMult*sigma);

        //! Add the index of different curves identified having points inside ball of squared radius (3*sigma)^2
        std::vector<int> indexCorbes;
        std::vector<int> indexDir;
        std::vector<float> indexDist;

        for (int oo=0; oo < nDir; oo++)
        for (int rr=-ksize; rr <= ksize; rr++)
        for (int ss=-ksize; ss <= ksize; ss++)
        if ( rr*rr+ss*ss <= radius && px+rr>=0 && px+rr < width && py+ss>=0 && py+ss<height)
        {
            if (curvesID[oo][(py+ss)*width+ px+rr]>0.0f)
            {

                int previous=0;
                for (int kkk=0; kkk < (int)indexCorbes.size() && previous==0; kkk++)
                if (indexCorbes[kkk] == (int) curvesID[oo][(py+ss)*width+ px+rr])
                {
                    previous=1;
                    if (rr*rr+ss*ss < indexDist[kkk])
                    {
                        indexDist[kkk]=rr*rr+ss*ss;
                        indexDir[kkk]=oo;
                    }

                }


                if (previous == 0)
                {
                    indexCorbes.push_back((int) curvesID[oo][(py+ss)*width+ px+rr]);
                    indexDir.push_back(oo);
                    indexDist.push_back(rr*rr+ss*ss);
                }

            }
        }


        float graus =  fCorValidationDegree / anglePrecision;
        int thrGraus = rintf( graus );


        if (indexCorbes.size() > 1)
        {

            int ncorbes=1;
            for (int ii=0; ii < (int) indexCorbes.size(); ii++)
            {
                int diferent = 0;
                for (int jj=ii+1; jj < (int) indexCorbes.size(); jj++)
                {
                    if ( abs(indexDir[ii] - indexDir[jj])>thrGraus &&  abs( (indexDir[ii] + nDir)%nDir - indexDir[jj])>thrGraus  )
                    if ( abs(indexDir[ii] - (indexDir[jj]+nDir/2)%nDir)>thrGraus &&  abs( (indexDir[ii] + nDir/2 + nDir)%nDir - indexDir[jj])>thrGraus )
                        diferent=1;

                }

                if (diferent==1) ncorbes++;

            }

            if (ncorbes > 1)
            {
                libUSTG::fiImageDrawCircle(localMax, px, py, radForZone, 128.0f, width, height);

                //! Mark as visited
                str_point_descriptor pp;
                pp.px = px;
                pp.py = py;

                cCorners2.push_back(pp);
            }

        }

    }

    cCorners = cCorners2;

}



//! 我并不要区分 T节点 和 角点
//! This function differentiates from corners and junctions based on the number of different edges arriving at a detected corner
//!
void differentiate_corners_junctions(vector<str_point_descriptor> &cCorners, float **curvesID, float fSigma, int nDir,  float rMult, int width, int height)
{


    float sigma = fSigma;

    //! for each corner
    for (int kk=0; kk < (int) cCorners.size(); kk++)
    {

        int px = cCorners[kk].px;
        int py = cCorners[kk].py;


        //! Identify different curves in neighborhood
        float radius = rMult *sigma;
        radius*=radius;

        int ksize = (int) ceilf(rMult*sigma);


        std::vector<int> indexCorbes;
        for (int oo=0; oo < nDir; oo++)
        for (int rr=-ksize; rr <= ksize; rr++)
        for (int ss=-ksize; ss <= ksize; ss++)
        if ( rr*rr+ss*ss <= radius && px+rr>=0 && px+rr < width && py+ss>=0 && py+ss<height)
        {
            if (curvesID[oo][(py+ss)*width+ px+rr]>0.0f)
            {
                int previous=0;
                for (int kkk=0; kkk < (int) indexCorbes.size() && previous==0; kkk++)
                if (indexCorbes[kkk] == (int) curvesID[oo][(py+ss)*width+ px+rr])
                {
                    previous=1;

                }


                if (previous == 0)
                {
                    indexCorbes.push_back((int) curvesID[oo][(py+ss)*width+ px+rr]);
                }

            }
        }

        cCorners[kk].ncorbes = (int) indexCorbes.size();
    }



}






