#include "Matrix.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <cstdio>
//#include <omp.h>

  // LAPACK/BLAS functions
  extern"C"{

    /*! 
     * SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     *               LDVR, WORK, LWORK, INFO )
     *   for info, see: http://www.netlib.org/lapack/double/dgeev.f              
     */
    void dgeev_(char* jobvl, char *jobvr, int *n, double *a, int *lda,
        double *wr, double *wi, double *vl, int *ldvl, double *vr,
        int *ldvr, double *work, int *lwork, int *info);

  /*! 
   * SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
   *   for info, see: http://www.netlib.org/blas/dgemm.f
   */
  void dgemm_(char *transa, char* transb, int *m, int *n, int *k, double *alpha, 
      double *a, int *lda, double *b, int *ldb, double *beta, double *c,
      int *ldc);

  /*! 
   *  SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
   *  for info, see http://www.netlib.org/lapack/double/dgetri.f
   */
  void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

  /*!
   * SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
   * for info, see http://www.netlib.org/lapack/double/dgetrf.f
   */
  void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
}
/*!
 *  Default constructor, constructs a matrix with size (0,0)
 */
Matrix::Matrix(): nr_rows(0), nr_cols(0), m_data(0){
}
/*!
 *  Constructs a square matrix, filled with 0
 *  @param size The side of the square matrix
 */
Matrix::Matrix(std::size_t size): nr_rows(size), nr_cols(size), m_data(size*size, 0){
}

/*!
 *  Constructs a matrix with size (rows, cols), filled with 0
 *  @param rows Amount of rows 
 *  @param cols Amount of columns
 */
Matrix::Matrix(std::size_t rows, std::size_t cols): nr_rows(rows), nr_cols(cols), m_data(cols*rows, 0){
}

/*!
 *  Constructs a matrix with size (rows, cols) and fills it with elements from 
 *  an array. The array needs to be of size rows*cols
 *  @param array[] The elements in row-major order
 *  @param rows The amount of rows
 *  @param cols The amount of cols
 */
Matrix::Matrix(const double array[], std::size_t rows, std::size_t cols): m_data(rows*cols, 0), nr_rows(rows), nr_cols(cols){
  for (int i=0; i<cols; i++) {
    for (int j=0; j<rows; j++) {
      m_data[i*rows + j] = array[j*cols + i];
    }
  }
}
/*! Constructs a diagonal matrix from a vector
 *  @param d A reference to a vector containing the diagonal elements
 */
Matrix::Matrix(const std::vector<double> &d): nr_rows(d.size()), nr_cols(d.size()), m_data(d.size()*d.size()){
  for (int i=0; i<d.size(); i++)
    m_data[i*d.size()+i] = d[i];
}

/*! Creates a copy of a matrix
 * @param m The matrix to be copied
 */
Matrix::Matrix(const Matrix &m) : nr_rows(m.get_rows()), nr_cols(m.get_cols()), m_data(m.m_data){
}
/*! Assignment operator
 *  @param m The matrix to be 
 *  */
Matrix &Matrix::operator=(const Matrix &m){
  if (&m != this){
    nr_rows = m.get_rows();
    nr_cols = m.get_cols();
    m_data.clear();
    m_data.reserve(m.get_rows()*m.get_cols());
    for (int i=0; i<m.get_rows()*m.get_cols(); i++)
      m_data.push_back(m(i));
  }
  return *this;
}
/*! Function that multiplies two matrices and returns a copy of the product
 *  @param lhv The left matrix
 *  @param rhv The right matrix
 *  @return a const copy with the product of lhv*rhv
 */
Matrix Matrix::mult(const Matrix &lhv,const Matrix &rhv) {
  return mult(lhv, rhv, false, false);
}

/*! 
 *  Static function that multiplies to matrices and possiby transposes them,
 *  returns a copy of the product
 *  @param lhv The left matrix 
 *  @param rhv The right matrix
 *  @param tr_left true if the left matrix should be transposed, false otherwise
 *  @param tr_right true if the right matrix should be transposed, false otherwise
 *  @return A const copy with the product lhv*rhv
 */
Matrix Matrix::mult(const Matrix &lhv, const Matrix &rhv, bool tr_left, bool tr_right) {
    if (lhv.get_cols() != rhv.get_rows())
      throw std::out_of_range("Matrix dimensions doesn't agree");

    char transl, transr;
    int left_op_rows = lhv.get_rows();
    int right_op_cols = rhv.get_cols();
    int left_op_cols = lhv.get_cols();
    double dummy_one = 1;
    double dummy_zero = 0;
    // Matrix to store the result
    Matrix result(left_op_rows, right_op_cols);

    // T - the matrix should be transposed
    // N - the matrix should not be transposed
    tr_left ? transl = 'T' : transl = 'N';
    tr_right ? transr = 'T' : transr = 'N';

    // dgemm_ is a BLAS subroutine that multiplies two matrices
    // const_cast<double*> because dgemm_ will not change those arguments 
    dgemm_(&transl, &transr, &left_op_rows, &right_op_cols, &left_op_cols, 
        &dummy_one, const_cast<double *>(&lhv.m_data[0]), &left_op_cols, 
        const_cast<double *>(&rhv.m_data[0]), &right_op_cols,
        &dummy_zero, &result.m_data[0], &left_op_rows); 
    return result;
  }

//Should this function be turned into a non-member too?
/*! 
 * 
 * @return A copy of the product
 */
Matrix Matrix::diag_mult(const std::vector<double> &diag) const{
  Matrix temp(*this);

  for (int i=0; i<get_cols(); i++)
    for (int j=0; j<get_rows(); j++)
      temp(i,j) = (*this)(i*get_rows()+j)*diag[i];  

  return temp;
}
/*!
 * A function that applies ln() on every element in the matrix.
 */
void Matrix::mlog(){
  int size = get_rows()*get_cols();
  for (int i=0; i<size; i++)
    (*this)(i) = log((*this)(i));
}

/*! 
 *  A const function that calculates the matrix exponential and returns a copy.
 *  The exponential is calculated by finding the eigenvalues and eigenvectors
 *  of the matrix, exponentiating the eigenvalues. The eigenvalues is stored in
 *  a matrix V, eigenvectors is stored in a matrix A, inv(A) is calculated.
 *  The product A*V*inv(A) is returned.
 *  @return A copy of the exponentiated matrix
 */
Matrix Matrix::expm() const{
  return expm(DblVec(1,1))[0];
}

/*! 
 *  A const function that for every value in a vector calculates the matrix 
 *  exponential of the matrix multiplied with that value
 *  The exponential is calculated by finding the eigenvalues and eigenvectors
 *  of the matrix, exponentiating the eigenvalues. The eigenvalues is stored in
 *  a matrix V, eigenvectors is stored in a matrix A, inv(A) is calculated.
 *  The product A*V*inv(A) is returned.
 *  @param s A vector with values to be multiplied with the matrix before 
 *            the exponent is calculated.
 *  @return A vector with the exponential of the matrix multiplied with every
 *            value in s
 */
MatVec Matrix::expm(const DblVec &s) const {

  // Can only calculate eigenvalues and vectors of square matrices
  if (get_rows() != get_cols())
    throw std::out_of_range("Matrix needs to be square");

  int size = get_rows();
  DblVec eg_val_real(size, 0); // Real part of eigenvalues
  DblVec eg_val_im(size, 0);   // Imaginary part of eigenvalues
                               // should be zero
  double dummy[1];
  int dummy_size = 1;
  double dummy_one = 1;
  int info[1];
  char n = 'N';   // Do not want to use this argument
  char v = 'V';   // Want to use this argument
  double workspace_size[1];
  int w_query = -1;

  // Need to make a copy of the data in Q to send into dgeev_ because
  // the data sent in is overwritten
  int data_size = get_rows()*get_cols();
  DblVec data(m_data);

  // Matrix for the eigenvectors  
  Matrix t_mat = Matrix(size, size);

  //workspace-query
  // SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
  //               LDVR, WORK, LWORK, INFO )
  dgeev_(&n, &v, &size, &data[0], &size, &eg_val_real[0], &eg_val_im[0], dummy, 
      &dummy_size, &t_mat.m_data[0], &size, workspace_size, &w_query, info);

  DblVec workspace_vec(static_cast<int>(workspace_size[0]), 0);
  int w_size = static_cast<int>(workspace_size[0]);

  // Real calculation of eigenvalues and eigenvectors for Q
  dgeev_(&n, &v, &size, &data[0], &size, &eg_val_real[0], &eg_val_im[0], dummy, 
      &dummy_size, &t_mat.m_data[0], &size, &workspace_vec[0], &w_size, info);

  // Calculating inverse of matrix with eigenvectors
  Matrix t_mat_inv(t_mat);
  int ipiv[size];

  // LU factorization, t_mat_inv.m_data is overwritten with the LU factorization
  dgetrf_(&size, &size, &t_mat_inv.m_data[0], &size, ipiv, info);

  //workspace-query, nothing happens with t_mat_inv.m_data
  dgetri_(&size, &t_mat_inv.m_data[0], &size, ipiv, workspace_size, &w_query, info);

  double workspace_vec2[static_cast<int>(workspace_size[0])];
  w_size = static_cast<int>(workspace_size[0]);

  // Inverse calculation from LU values, the inverse is stored in t_mat_inv.m_data
  dgetri_(&size, &t_mat_inv.m_data[0], &size, ipiv, workspace_vec2, &w_size, info);

  MatVec result;
  result.reserve(s.size());

  // e^(this) = T*D*T^-1
  // T = matrix with eigenvectors (t_mat), D = matrix with exponentiated eigenvalues
  // Calculate for every value in incoming vector s
  DblVec eg_val_exp; 
  eg_val_exp.reserve(size);
  for (DblVec::const_iterator it=s.begin(); it != s.end(); it++){
    for (int i=0; i<size; i++)
      eg_val_exp.push_back(exp(eg_val_real[i]*(*it)));
    Matrix left = Matrix::mult(t_mat, Matrix(eg_val_exp));
    Matrix res = Matrix::mult( left, t_mat_inv);
    result.push_back(res);
    eg_val_exp.clear();
  }
  return result;
}

/*! 
 *  A function for printing the matrix
 */
void Matrix::printm() const{
  for (int i=0; i<get_rows(); i++) {
    for (int j=0; j<get_cols(); j++)
      printf(" %.20e ", (*this)(i,j));
    std::cout << std::endl;
  }
}
void Matrix::printmi() const{
  for (int i=0; i<get_rows(); i++) {
    for (int j=0; j<get_cols(); j++)
      printf(" %.3d ", (int) (*this)(i,j));
    std::cout << std::endl;
  }
}
/*! 
 *  A const function that sums all the elements in the matrix.
 *  @return The sum of all the elements in the matrix
 */
double Matrix::sum() const{
  return std::accumulate(m_data.begin(), m_data.end(), 0.0);
}

/*! 
 *  A const function that sums all the elements on the diagonal of the matrix
 *  @return The sum of all the diagonal elements of the matrix
 */
double Matrix::sum_diag() const{
  double sum=0.0;
  if (get_rows() != get_cols())
    throw std::out_of_range("Matrix needs to be square");

  for (int i=0;i<get_rows();i++)
    sum += (*this)(i,i);

  return sum;
}
//--------------------------------------------------------------
// Non-member functions 

/*! 
 *  Non-member function that multiplies two matrices elementwise. The product is 
 *  C(i,j) = A(i,j)*B(i,j). The matrices must have the same dimensions.
 *  @param lhv The matrix to the right
 *  @param rhv The matrix to the right
 *  @param A copy of the product
 */
  Matrix elem_mult(const Matrix &lhv, const Matrix &rhv){
    if (lhv.get_rows() != rhv.get_rows() || lhv.get_cols() != rhv.get_cols())
      throw std::out_of_range("Matrix dimensions must agree");

    Matrix temp(lhv.get_rows(), rhv.get_cols());

    for (int i=0; i<lhv.get_rows()*lhv.get_cols(); i++)
      temp(i) = lhv(i)*rhv(i);

    return temp;

  }

/*! 
 *  Non-member function that divides a matrix with another elementwise. A copy 
 *  of the quotient C(i,j) = A(i,j)/B(i,j) is returned.
 *  @param lhv The dividend
 *  @param rhv The divisor
 *  @return The quotient
 */ 
  Matrix elem_div(const Matrix &lhv, const Matrix &rhv){
    if (lhv.get_rows() != rhv.get_rows() || lhv.get_cols() != rhv.get_cols())
      throw std::out_of_range("Matrix dimensions must agree");

    Matrix temp(lhv.get_rows(), lhv.get_cols());

    for (int i=0; i<lhv.get_rows()*lhv.get_cols(); i++)
      temp(i) = lhv(i)/rhv(i);

    return temp;
  }
