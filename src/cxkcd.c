#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


// upper bound of y: \frac{1}{\sqrt{2\pi} \sigma}
double maxy(int nsd){
    return 1 / (sqrt(2 * M_PI) * nsd);
}
    
// upper inverse of pdf of N(0, \sigma^2)
double invpdf(double y, double sd){
    return sqrt(-2 * pow(sd, 2) * log(y * sqrt(2 * M_PI) * sd));
}

// random number generator of xkcd distribution with C implementation
void crxkcd(int *nin, double *nsdin, int *swapin, double *x){
    int n = nin[0];
    double nsd = nsdin[0];
    int swap = swapin[0];

    int i;
    double rx;
    double ry;

    GetRNGstate();

    if (swap == 1) {
        for (i = 0; i < n; i++){
            rx = rnorm(0, nsd);
            ry = runif(0, dnorm(rx, 0, nsd, 0));
            x[i] = maxy(nsd) - ry;
        }
    }

    else {
        for (i = 0; i < n; i++){
            rx = rnorm(0, nsd);
            ry = runif(0, dnorm(rx, 0, nsd, 0));
            x[i] = ry;
        } 
    }
       
    PutRNGstate();
}

// density of xkcd distribution with C implementation
void cdxkcd(double *yin, double *nsdin, int *login, int *swapin, double *d){
    double y = yin[0];
    double nsd = nsdin[0];
    int logl = login[0];
    int swap = swapin[0];    

    double tempd;
    double w;

    if (swap == 0) {
        if (y < 0 || y > maxy(nsd))
            tempd = 0;
        else 
            tempd = invpdf(y, nsd);
        if (logl == 1)
            d[0] = log(2 * tempd);
        else
            d[0] = 2 * tempd;  
    }  
    else {
        w = maxy(nsd) - y;
        // the routine above with default setup
        if (w < 0 || w > maxy(nsd))
            tempd = 0;
        else 
            tempd = invpdf(w, nsd);
        tempd = tempd*2;

        if (logl == 1)
            d[0] = log(tempd);
        else
            d[0] = tempd;
    }
}


// cdf of xkcd distribution with C implementation
void cpxkcd(double *yin, double *nsdin, int *login, int *swapin, double *p) {
    double y = yin[0];
    double nsd = nsdin[0];
    int logp = login[0];
    int swap = swapin[0];    

    double Fy;
    double tempx;
    double dp;

    if (swap == 0){
        if (y <= 0)
            Fy = 0;
        else if (y >= maxy(nsd))
            Fy = 1;
        else {
            tempx = invpdf(y, nsd);
            Fy = 2 * pnorm(-tempx, 0, nsd, 1, 0) + 2 * y * tempx;
        }
        if (logp == 1)
            p[0] = log(Fy);
        else
            p[0] = Fy;
    }
    else {
        tempx = maxy(nsd) - y;

        // the routine above with default setup
        dp = invpdf(tempx, nsd);
        dp = 2 * pnorm(-dp, 0, nsd, 1, 0) + 2 * tempx * dp;

        dp = 1 - dp;
        if (logp == 1)
            p[0] = log(dp);
        else
            p[0] = dp;
    }
}

// quantile of xkcd distribution with C implementation
void cqxkcd(double *pin, double *nsdin, int *login, int *swapin, double *tolin, double *q) {
    double p = pin[0];
    double nsd = nsdin[0];
    int logp = login[0];
    int swap = swapin[0];       
    double tol = tolin[0];

    double left = 0, right = maxy(nsd), mid;
    double value;
    double tempx;

    if (swap == 1)
        p = 1 - p;

    // biniary search for root of a cdf
    while (right - left > tol) {
        mid = (left + right)/2;

        // naive call of cpxkcd
        tempx = invpdf(mid, nsd);
        value = 2 * pnorm(-tempx, 0, nsd, 1, 0) + 2 * mid * tempx;

        if (value > p)
            right = mid;
        else
            left = mid;
    }

    if (swap == 1)
        q[0] = maxy(nsd) - mid;
    else
        q[0] = mid;
}