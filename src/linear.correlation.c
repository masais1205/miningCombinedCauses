#include "common.h"

static SEXP covmat(SEXP data, SEXP length);

/* Linear Correlation, to be used in C code. */
double c_fast_cor(double *xx, double *yy, int *num) {

int i = 0;
double xm = 0, ym = 0;
long double xsd = 0, ysd = 0, sum = 0;
double tol = MACHINE_TOL;

  /* compute the mean values.  */
  for (i = 0 ; i < *num; i++) {

    xm += xx[i];
    ym += yy[i];

  }/*FOR*/

  xm /= (*num);
  ym /= (*num);

  /* compute the actual covariance. */
  for (i = 0; i < *num; i++) {

    sum += (xx[i] - xm) * (yy[i] - ym);
    xsd += (xx[i] - xm) * (xx[i] - xm);
    ysd += (yy[i] - ym) * (yy[i] - ym);

  }/*FOR*/

  /* safety check against "divide by zero" errors. */
  if ((xsd < tol) || (ysd < tol))
    sum = 0;
  else
    sum /= sqrt(xsd * ysd);

  /* double check that the coefficient is in the [-1, 1] range. */
  if (sum > 1) {

    warning("fixed correlation coefficient greater than 1, probably due to floating point errors.");

    sum = 1;

  }/*THEN*/
  else if (sum < -1) {

    warning("fixed correlation coefficient lesser than -1, probably due to floating point errors.");

    sum = -1;

  }/*ELSE*/

  return (double)sum;

}/*C_FAST_COR*/

/* Linear Correlation. */
SEXP fast_cor(SEXP x, SEXP y, SEXP length) {

int *n = INTEGER(length);
double *xx = REAL(x), *yy = REAL(y);
SEXP res;

  PROTECT(res = allocVector(REALSXP, 1));
  NUM(res) = c_fast_cor(xx, yy, n);
  UNPROTECT(1);

  return res;

}/*FAST_COR*/

double c_fast_pcor(double *covariance, int *ncols, double *u, double *d,
    double *vt, int *errcode) {

int i = 0, coord1 = 0, coord2 = 0;
double k11 = 0, k12 = 0, k22 = 0;
double res = 0, tol = MACHINE_TOL, sv_tol = 0;

  c_svd(covariance, u, d, vt, ncols, ncols, ncols, FALSE, errcode);

  if (*errcode != 0)
    return 0;

  /* set the threshold for the singular values as in corpcor. */
  sv_tol = (*ncols) * d[0] * tol * tol;

  /* compute the three elements of the pseudoinverse needed
   * for the partial correlation coefficient. */
  for (i = 0; i < *ncols; i++) {

    if (d[i] > sv_tol) {

      coord1 = CMC(0, i, *ncols);
      coord2 = CMC(i, 1, *ncols);

      k11 += u[coord1] * vt[CMC(i, 0, *ncols)] / d[i];
      k12 += u[coord1] * vt[coord2] / d[i];
      k22 += u[CMC(1, i, *ncols)] * vt[coord2] / d[i];

    }/*THEN*/

  }/*FOR*/

  /* safety check against "divide by zero" errors. */
  if ((k11 < tol) || (k22 < tol))
    res = 0;
  else
    res = -k12 / sqrt(k11 * k22);

  return res;

}/*C_FAST_PCOR*/

/* Partial Linear Correlation. */
SEXP fast_pcor(SEXP data, SEXP length, SEXP shrinkage, SEXP strict) {

int i = 0, ncols = length(data), errcode = 0;
int *shrink = LOGICAL(shrinkage);
double *u = NULL, *d = NULL, *vt = NULL, *res = NULL;
double k11 = 0, k12 = 0, k22 = 0;
double tol = MACHINE_TOL, sv_tol = 0;
SEXP result, cov;

  /* compute the covariance matrix. */
  if (*shrink > 0)
    PROTECT(cov = cov_lambda(data, length));
  else
    PROTECT(cov = covmat(data, length));

  /* allocate the matrices needed for the SVD decomposition. */
  u = alloc1dreal(ncols * ncols);
  d = alloc1dreal(ncols);
  vt = alloc1dreal(ncols * ncols);

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);

  /* compute the singular value decomposition of the covariance matrix. */
  c_svd(REAL(cov), u, d, vt, &ncols, &ncols, &ncols, FALSE, &errcode);

  if (errcode != 0) {

    /* unprotect everything and return. */
    UNPROTECT(2);

    if (isTRUE(strict)) {

      error("failed to compute the pseudoinverse of the covariance matrix.");

    }/*THEN*/
    else {

      /* if computing SVD decomposition fails, assume a null correlation. */
      *res = 0;
      /* warn the user that something went wrong.*/
      warning("failed to compute the pseudoinverse of the covariance matrix, assuming independence.");

      return result;

    }/*ELSE*/

  }/*THEN*/

  /* set the threshold for the singular values as in corpcor. */
  sv_tol = ncols * d[0] * tol * tol;

  /* compute the three elements of the pseudoinverse needed
   * for the partial correlation coefficient. */
  for (i = 0; i < ncols; i++) {

    if (d[i] > sv_tol) {

      k11 += u[CMC(0, i, ncols)] * vt[CMC(i, 0, ncols)] / d[i];
      k12 += u[CMC(0, i, ncols)] * vt[CMC(i, 1, ncols)] / d[i];
      k22 += u[CMC(1, i, ncols)] * vt[CMC(i, 1, ncols)] / d[i];

    }/*THEN*/

  }/*FOR*/

  /* safety check against "divide by zero" errors and negative variances. */
  if ((k11 < tol) || (k22 < tol))
    *res = 0;
  else
    *res = -k12 / sqrt(k11 * k22);

  /* double check that partial correlation really is in [-1, 1], ill
   * conditioned matrices and numeric errors in SVD can result in
   * invalid partial correlation coefficients. */
  if (*res > 1) {

    warning("fixed partial correlation coefficient greater than 1, probably due to floating point errors.");

    *res = 1;

  }/*THEN*/
  else if (*res < -1) {

    warning("fixed partial correlation coefficient lesser than -1, probably due to floating point errors.");

    *res = -1;

  }/*ELSE*/

  UNPROTECT(2);
  return result;

}/*FAST_PCOR*/

void c_covmat(double **data, double *mean, int *ncols, int *nrows,
    double *mat) {

int i = 0, j = 0, k = 0;
long double temp = 0;

  /* compute the actual covariance. */
  for (i = 0; i < *ncols; i++) {

    for (j = i; j < *ncols; j++) {

      for (k = 0, temp = 0; k < *nrows; k++)
        temp += (data[i][k] - mean[i]) * (data[j][k] - mean[j]);

      /* fill in the symmetric element of the matrix. */
      mat[CMC(j, i, *ncols)] = mat[CMC(i, j, *ncols)] =
        (double)(temp / (*nrows - 1));

    }/*FOR*/

  }/*FOR*/

}/*C_COVMAT*/

static SEXP covmat(SEXP data, SEXP length) {

int i = 0, j = 0, *n = INTEGER(length), ncols = length(data);
double *var = NULL, *mean = NULL, **column = NULL;
SEXP res;

  /* allocate the covariance matrix. */
  PROTECT(res = allocMatrix(REALSXP, ncols, ncols));
  var = REAL(res);
  memset(var, '\0', ncols * ncols * sizeof(double));

  /* allocate and initialize an array of pointers for the variables. */
  column = (double **) alloc1dpointer(ncols);
  for (i = 0; i < ncols; i++)
    column[i] = REAL(VECTOR_ELT(data, i));

  /* allocate an array to store the mean values. */
  mean = alloc1dreal(ncols);

  /* compute the mean values  */
  for (i = 0; i < ncols; i++) {

    for (j = 0 ; j < *n; j++)
      mean[i] += column[i][j];

    mean[i] /= (*n);

  }/*FOR*/

  /* call the C backend that does the actual work. */
  c_covmat(column, mean, &ncols, n, var);

  UNPROTECT(1);
  return res;

}/*COVMAT*/

/* update only a single row/column in a covariance matrix. */
void c_update_covmat(double **data, double *mean, int update, int *ncols,
    int *nrows, double *mat) {

int j = 0, k = 0;
long double temp = 0;

  /* compute the actual covariance. */
  for (j = 0; j < *ncols; j++) {

    for (k = 0, temp = 0; k < *nrows; k++)
      temp += (data[update][k] - mean[update]) * (data[j][k] - mean[j]);

    /* fill the symmetric elements of the matrix. */
    mat[CMC(j, update, *ncols)] = mat[CMC(update, j, *ncols)] =
      (double)(temp / (*nrows - 1));

  }/*FOR*/

}/*C_UPDATE_COVMAT*/

