#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#include <stdexcept>
#include <vector>

  /*
   * Class for handling matrices.
   * Internal data is stored in column major format for easier use with LAPACK.
   * Element at matrix[1][2] is accessed with matrix(1,2);
   */
  class Matrix;
  typedef std::vector<double> DblVec;
  typedef std::vector<Matrix> MatVec;

  class Matrix {
    public:
      //! Default constructor
      Matrix();
      //! Constructor for square matrix
      explicit Matrix(std::size_t size);
      //! Constructor for rectangular matrix
      Matrix(std::size_t rows, std::size_t cols);
      //! Creates full matrix from array
      Matrix(const double[], std::size_t, std::size_t);
      //! Creates diagonal matrix from std::vector
      Matrix(const DblVec &);
      //! Copy constructor
      Matrix(const Matrix &);
      //! Assignment operator
      Matrix &operator=(const Matrix&);
      //! Returns the row dimension
      std::size_t get_rows() const { return nr_rows;}
      //! Returns the column dimension
      std::size_t get_cols() const { return nr_cols;}
      //! Index operator (row, col)
      double &operator()(const std::size_t, const std::size_t);
      //! Const index operator (row, col)
      const double operator()(const std::size_t,const std::size_t) const;
      //! Index operator for raw data
      double &operator()(const std::size_t);
      //! Const index operator for raw data
      const double operator()(const std::size_t) const;
      //! Matrix multiplication, with possibility to transpose the matrices
      static Matrix mult(const Matrix &lhv,const Matrix &rhv, bool tr_left, bool tr_right);
      //! Matrix multiplication
      static Matrix mult(const Matrix &lhv,const Matrix &rhv);
      //! Multiplication with a diagonal matrix 
      Matrix diag_mult(const DblVec &) const;
      //! Elementwise logarithm of the matrix
      void mlog();
      //! Calculates the matrix exponential
      Matrix expm() const;
      //! Multiplies the matrix with a scalar, and then calculates the matrix exponential
      MatVec expm(const DblVec &) const;
      //! Prints the matrix
      void printm() const;
      void printmi() const;
      //! Sums all the elements
      double sum() const;
      //! Sums the diagonal elements
      double sum_diag() const;

    private:
      //! Stores the matrix in column-major order
      DblVec m_data;
      //! The amount of rows
      std::size_t nr_rows;
      //! The columns
      std::size_t nr_cols;
      
  };

  //Non-member functions
  //! Elementwise matrix multiplication
  Matrix elem_mult(const Matrix &lhv, const Matrix &rhv);
  //! Elementwise matrix division
  Matrix elem_div(const Matrix &lhv, const Matrix &rhv);

  /*!
   *  Inlined index operator
   *  @param row The row 
   *  @param col The column
   *  @return A reference to the element at position (row, col)
   */
  inline double &Matrix::operator() (const std::size_t row, const std::size_t col){
    if (row >= get_rows() || col >= get_cols())
      throw std::out_of_range("Indexing outside of matrix");
    return m_data[col*nr_rows + row];   
  }

  /*!
   *  Inlined const index operator
   *  @param row The row 
   *  @param col The column
   *  @return A const copy of the element at position (row, col)
   */
  inline const double Matrix::operator() (const std::size_t row, const std::size_t col) const{
    if (row >= get_rows() || col >= get_cols())
      throw std::out_of_range("Indexing outside of matrix");
    return m_data[col*nr_rows + row];   
  }

  /*!
   *  Inlined index operator
   *  @param ind Index in the representation of the matrix
   *  @return A reference to the element at position (ind) in the datastructure for the matrix
   */
  inline double &Matrix::operator() (const std::size_t ind){
    if (ind >= get_cols()*get_rows())
      throw std::out_of_range("Indexing outside of matrix");
    return m_data[ind];   
  }

  /*!
   *  Inlined const index operator
   *  @param ind Index in the representation of the matrix
   *  @return A const copy of the element at position (ind) in the datastructure for the matrix
   */
  inline const double Matrix::operator() (const std::size_t ind) const{
    if (ind >= get_cols()*get_rows())
      throw std::out_of_range("Indexing outside of matrix");
    return m_data[ind];   
  }

#endif
