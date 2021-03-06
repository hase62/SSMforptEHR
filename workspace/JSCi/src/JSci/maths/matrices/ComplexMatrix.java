/* AUTO-GENERATED */
package JSci.maths.matrices;

import JSci.maths.Complex;
import JSci.maths.ComplexMapping;
import JSci.maths.DimensionException;
import JSci.maths.vectors.AbstractComplexVector;
import JSci.maths.vectors.ComplexVector;
import JSci.maths.groups.AbelianGroup;
import JSci.maths.algebras.*;
import JSci.maths.fields.*;

/**
* The ComplexMatrix class provides an object for encapsulating matrices containing complex numbers.
* @version 2.2
* @author Mark Hale
*/
public class ComplexMatrix extends AbstractComplexMatrix {
        /**
        * Arrays containing the elements of the matrix.
        */
        protected double matrixRe[][],matrixIm[][];
        /**
        * Constructs an empty matrix.
        * @param rows the number of rows
        * @param cols the number of columns
        */
        public ComplexMatrix(final int rows,final int cols) {
                super(rows,cols);
                matrixRe=new double[rows][cols];
                matrixIm=new double[rows][cols];
        }
        /**
        * Constructs a matrix by wrapping two arrays.
        * @param arrayRe an array of real values
        * @param arrayIm an array of imaginary values
        */
        public ComplexMatrix(final double arrayRe[][],final double arrayIm[][]) {
                this(arrayRe.length,arrayRe[0].length);
                matrixRe=arrayRe;
                matrixIm=arrayIm;
        }
        /**
        * Constructs a matrix from an array.
        * @param array an assigned value
        */
        public ComplexMatrix(final Complex array[][]) {
                this(array.length,array[0].length);
                for(int j,i=0;i<numRows;i++) {
                        for(j=0;j<numCols;j++) {
                                matrixRe[i][j]=array[i][j].real();
                                matrixIm[i][j]=array[i][j].imag();
                        }
                }
        }
        /**
        * Constructs a matrix from an array of vectors (columns).
        * @param array an assigned value
        */
        public ComplexMatrix(ComplexVector array[]) {
                this(array[0].dimension(),array.length);
                for(int j,i=0;i<numRows;i++) {
                        for(j=0;j<numCols;j++) {
                                matrixRe[i][j]=array[j].getComponent(i).real();
                                matrixIm[i][j]=array[j].getComponent(i).imag();
                        }
                }
        }
        /**
        * Compares two complex matrices for equality.
        * @param m a complex matrix
        */
        public boolean equals(AbstractComplexMatrix m, double tol) {
                if(m != null && numRows == m.rows() && numCols == m.columns()) {
			double sumSqr = 0.0;
                        for(int i=0;i<numRows;i++) {
                                for(int j=0;j<numCols;j++) {
					double deltaRe = matrixRe[i][j]-m.getRealElement(i,j);
					double deltaIm = matrixIm[i][j]-m.getImagElement(i,j);
					sumSqr += deltaRe*deltaRe + deltaIm*deltaIm;
                                }
                        }
                        return (sumSqr <= tol*tol);
                } else {
                        return false;
                }
        }
        /**
        * Returns a string representing this matrix.
        */
        public String toString() {
                final StringBuffer buf=new StringBuffer(5*numRows*numCols);
                for(int j,i=0;i<numRows;i++) {
                        for(j=0;j<numCols;j++) {
                                buf.append(Complex.toString(matrixRe[i][j],matrixIm[i][j]));
                                buf.append(' ');
                        }
                        buf.append('\n');
                }
                return buf.toString();
        }
        /**
        * Returns a hashcode for this matrix.
        */
        public int hashCode() {
                return (int)Math.exp(infNorm());
        }
        /**
        * Returns the real part of this complex matrix.
        * @return a double matrix
        */
        public AbstractDoubleMatrix real() {
                return new DoubleMatrix(matrixRe);
        }
        /**
        * Returns the imaginary part of this complex matrix.
        * @return a double matrix
        */
        public AbstractDoubleMatrix imag() {
                return new DoubleMatrix(matrixIm);
        }
        /**
        * Returns an element of the matrix.
        * @param i row index of the element
        * @param j column index of the element
        * @exception MatrixDimensionException If attempting to access an invalid element.
        */
        public Complex getElement(final int i, final int j) {
                if(i>=0 && i<numRows && j>=0 && j<numCols)
                        return new Complex(matrixRe[i][j],matrixIm[i][j]);
                else
                        throw new MatrixDimensionException(getInvalidElementMsg(i,j));
        }
        public double getRealElement(final int i, final int j) {
                if(i>=0 && i<numRows && j>=0 && j<numCols)
                        return matrixRe[i][j];
                else
                        throw new MatrixDimensionException(getInvalidElementMsg(i,j));
        }
        public double getImagElement(final int i, final int j) {
                if(i>=0 && i<numRows && j>=0 && j<numCols)
                        return matrixIm[i][j];
                else
                        throw new MatrixDimensionException(getInvalidElementMsg(i,j));
        }
        /**
        * Sets the value of an element of the matrix.
        * Should only be used to initialise this matrix.
        * @param i row index of the element
        * @param j column index of the element
        * @param z a complex number
        * @exception MatrixDimensionException If attempting to access an invalid element.
        */
        public void setElement(final int i, final int j, final Complex z) {
                if(i>=0 && i<numRows && j>=0 && j<numCols) {
                        matrixRe[i][j]=z.real();
                        matrixIm[i][j]=z.imag();
                } else
                        throw new MatrixDimensionException(getInvalidElementMsg(i,j));
        }
        /**
        * Sets the value of an element of the matrix.
        * Should only be used to initialise this matrix.
        * @param i row index of the element
        * @param j column index of the element
        * @param x the real part of a complex number
        * @param y the imaginary part of a complex number
        * @exception MatrixDimensionException If attempting to access an invalid element.
        */
        public void setElement(final int i, final int j, final double x, final double y) {
                if(i>=0 && i<numRows && j>=0 && j<numCols) {
                        matrixRe[i][j]=x;
                        matrixIm[i][j]=y;
                } else
                        throw new MatrixDimensionException(getInvalidElementMsg(i,j));
        }
        /**
        * Returns the l<sup><img border=0 alt="infinity" src="doc-files/infinity.gif"></sup>-norm.
        * @author Taber Smith
        */
        public double infNorm() {
                double result=0.0,tmpResult;
                for(int i=0;i<numRows;i++) {
                        tmpResult=0.0;
                        for(int j=0;j<numCols;j++)
                                tmpResult += Math.sqrt((matrixRe[i][j]*matrixRe[i][j] + matrixIm[i][j]*matrixIm[i][j]));
                        if(tmpResult>result)
                                result=tmpResult;
                }
                return result;
        }
        /**
        * Returns the Frobenius or Hilbert-Schmidt (l<sup>2</sup>) norm.
        * @jsci.planetmath FrobeniusMatrixNorm
        * @author Taber Smith
        */
        public double frobeniusNorm() {
                double result=0.0;
                for(int j,i=0;i<numRows;i++)
                        for(j=0;j<numCols;j++)
                                result+=matrixRe[i][j]*matrixRe[i][j]+matrixIm[i][j]*matrixIm[i][j];
                return Math.sqrt(result);
        }

//============
// OPERATIONS
//============

        /**
        * Returns the negative of this matrix.
        */
        public AbelianGroup.Member negate() {
                final double arrayRe[][]=new double[numRows][numCols];
                final double arrayIm[][]=new double[numRows][numCols];
                for(int j,i=0;i<numRows;i++) {
                        arrayRe[i][0]=-matrixRe[i][0];
                        arrayIm[i][0]=-matrixIm[i][0];
                        for(j=1;j<numCols;j++) {
                                arrayRe[i][j]=-matrixRe[i][j];
                                arrayIm[i][j]=-matrixIm[i][j];
                        }
                }
                return new ComplexMatrix(arrayRe,arrayIm);
        }

// ADDITION

        /**
        * Returns the addition of this matrix and another.
        * @param m a complex matrix
        * @exception MatrixDimensionException If the matrices are different sizes.
        */
        public AbstractComplexMatrix add(final AbstractComplexMatrix m) {
                if(m instanceof ComplexMatrix) {
                        return add((ComplexMatrix)m);
                } else {
                        if(numRows==m.rows() && numCols==m.columns()) {
                                final double arrayRe[][]=new double[numRows][numCols];
                                final double arrayIm[][]=new double[numRows][numCols];
                                for(int j,i=0;i<numRows;i++) {
                                        arrayRe[i][0]=matrixRe[i][0]+m.getElement(i,0).real();
                                        arrayIm[i][0]=matrixIm[i][0]+m.getElement(i,0).imag();
                                        for(j=1;j<numCols;j++) {
                                                arrayRe[i][j]=matrixRe[i][j]+m.getElement(i,j).real();
                                                arrayIm[i][j]=matrixIm[i][j]+m.getElement(i,j).imag();
                                        }
                                }
                                return new ComplexMatrix(arrayRe,arrayIm);
                        } else
                                throw new MatrixDimensionException("Matrices are different sizes.");
                }
        }
        public ComplexMatrix add(final ComplexMatrix m) {
                if(numRows==m.numRows && numCols==m.numCols) {
                        final double arrayRe[][]=new double[numRows][numCols];
                        final double arrayIm[][]=new double[numRows][numCols];
                        for(int j,i=0;i<numRows;i++) {
                                arrayRe[i][0]=matrixRe[i][0]+m.matrixRe[i][0];
                                arrayIm[i][0]=matrixIm[i][0]+m.matrixIm[i][0];
                                for(j=1;j<numCols;j++) {
                                        arrayRe[i][j]=matrixRe[i][j]+m.matrixRe[i][j];
                                        arrayIm[i][j]=matrixIm[i][j]+m.matrixIm[i][j];
                                }
                        }
                        return new ComplexMatrix(arrayRe,arrayIm);
                } else
                        throw new MatrixDimensionException("Matrices are different sizes.");
        }

// SUBTRACTION

        /**
        * Returns the subtraction of this matrix by another.
        * @param m a complex matrix
        * @exception MatrixDimensionException If the matrices are different sizes.
        */
        public AbstractComplexMatrix subtract(final AbstractComplexMatrix m) {
                if(m instanceof ComplexMatrix) {
                        return subtract((ComplexMatrix)m);
                } else {
                        if(numRows==m.rows() && numCols==m.columns()) {
                                final double arrayRe[][]=new double[numRows][numCols];
                                final double arrayIm[][]=new double[numRows][numCols];
                                for(int j,i=0;i<numRows;i++) {
                                        arrayRe[i][0]=matrixRe[i][0]-m.getElement(i,0).real();
                                        arrayIm[i][0]=matrixIm[i][0]-m.getElement(i,0).imag();
                                        for(j=1;j<numCols;j++) {
                                                arrayRe[i][j]=matrixRe[i][j]-m.getElement(i,j).real();
                                                arrayIm[i][j]=matrixIm[i][j]-m.getElement(i,j).imag();
                                        }
                                }
                                return new ComplexMatrix(arrayRe,arrayIm);
                        } else
                                throw new MatrixDimensionException("Matrices are different sizes.");
                }
        }
        public ComplexMatrix subtract(final ComplexMatrix m) {
                if(numRows==m.numRows && numCols==m.numCols) {
                        final double arrayRe[][]=new double[numRows][numCols];
                        final double arrayIm[][]=new double[numRows][numCols];
                        for(int j,i=0;i<numRows;i++) {
                                arrayRe[i][0]=matrixRe[i][0]-m.matrixRe[i][0];
                                arrayIm[i][0]=matrixIm[i][0]-m.matrixIm[i][0];
                                for(j=1;j<numCols;j++) {
                                        arrayRe[i][j]=matrixRe[i][j]-m.matrixRe[i][j];
                                        arrayIm[i][j]=matrixIm[i][j]-m.matrixIm[i][j];
                                }
                        }
                        return new ComplexMatrix(arrayRe,arrayIm);
                } else
                        throw new MatrixDimensionException("Matrices are different sizes.");
        }

// SCALAR MULTIPLICATION

        /**
        * Returns the multiplication of this matrix by a scalar.
        * @param z a complex number
        * @return a complex matrix
        */
        public AbstractComplexMatrix scalarMultiply(final Complex z) {
                final double real=z.real();
                final double imag=z.imag();
                final double arrayRe[][]=new double[numRows][numCols];
                final double arrayIm[][]=new double[numRows][numCols];
                for(int j,i=0;i<numRows;i++) {
                        arrayRe[i][0]=real*matrixRe[i][0]-imag*matrixIm[i][0];
                        arrayIm[i][0]=imag*matrixRe[i][0]+real*matrixIm[i][0];
                        for(j=1;j<numCols;j++) {
                                arrayRe[i][j]=real*matrixRe[i][j]-imag*matrixIm[i][j];
                                arrayIm[i][j]=imag*matrixRe[i][j]+real*matrixIm[i][j];
                        }
                }
                return new ComplexMatrix(arrayRe,arrayIm);
        }
        /**
        * Returns the multiplication of this matrix by a scalar.
        * @param x a double
        * @return a complex matrix
        */
        public AbstractComplexMatrix scalarMultiply(final double x) {
                final double arrayRe[][]=new double[numRows][numCols];
                final double arrayIm[][]=new double[numRows][numCols];
                for(int j,i=0;i<numRows;i++) {
                        arrayRe[i][0]=x*matrixRe[i][0];
                        arrayIm[i][0]=x*matrixIm[i][0];
                        for(j=1;j<numCols;j++) {
                                arrayRe[i][j]=x*matrixRe[i][j];
                                arrayIm[i][j]=x*matrixIm[i][j];
                        }
                }
                return new ComplexMatrix(arrayRe,arrayIm);
        }

// SCALAR DIVISON

        /**
        * Returns the division of this matrix by a scalar.
        * @param z a complex number
        * @return a complex matrix
        */
        public AbstractComplexMatrix scalarDivide(final Complex z) {
                final Complex array[][]=new Complex[numRows][numCols];
                for(int j,i=0;i<numRows;i++) {
                        array[i][0]=new Complex(matrixRe[i][0],matrixIm[i][0]).divide(z);
                        for(j=1;j<numCols;j++)
                                array[i][j]=new Complex(matrixRe[i][j],matrixIm[i][j]).divide(z);
                }
                return new ComplexMatrix(array);
        }
        /**
        * Returns the division of this matrix by a scalar.
        * @param x a double
        * @return a complex matrix
        */
        public AbstractComplexMatrix scalarDivide(final double x) {
                final double arrayRe[][]=new double[numRows][numCols];
                final double arrayIm[][]=new double[numRows][numCols];
                for(int j,i=0;i<numRows;i++) {
                        arrayRe[i][0]=matrixRe[i][0]/x;
                        arrayIm[i][0]=matrixIm[i][0]/x;
                        for(j=1;j<numCols;j++) {
                                arrayRe[i][j]=matrixRe[i][j]/x;
                                arrayIm[i][j]=matrixIm[i][j]/x;
                        }
                }
                return new ComplexMatrix(arrayRe,arrayIm);
        }

// MATRIX MULTIPLICATION

        /**
        * Returns the multiplication of a vector by this matrix.
        * @param v a complex vector
        * @exception DimensionException If the matrix and vector are incompatible.
        */
        public AbstractComplexVector multiply(final AbstractComplexVector v) {
                if(numCols==v.dimension()) {
                        final double arrayRe[]=new double[numRows];
                        final double arrayIm[]=new double[numRows];
                        Complex comp;
                        for(int j,i=0;i<numRows;i++) {
                                comp=v.getComponent(0);
                                arrayRe[i]=(matrixRe[i][0]*comp.real() - matrixIm[i][0]*comp.imag());
                                arrayIm[i]=(matrixIm[i][0]*comp.real() + matrixRe[i][0]*comp.imag());
                                for(j=1;j<numCols;j++) {
                                        comp=v.getComponent(j);
                                        arrayRe[i]+=(matrixRe[i][j]*comp.real() - matrixIm[i][j]*comp.imag());
                                        arrayIm[i]+=(matrixIm[i][j]*comp.real() + matrixRe[i][j]*comp.imag());
                                }
                        }
                        return new ComplexVector(arrayRe,arrayIm);
                } else
                        throw new DimensionException("Matrix and vector are incompatible.");
        }
        /**
        * Returns the multiplication of this matrix and another.
        * @param m a complex matrix
        * @return an AbstractComplexMatrix or an AbstractComplexSquareMatrix as appropriate
        * @exception MatrixDimensionException If the matrices are incompatible.
        */
        public AbstractComplexMatrix multiply(final AbstractComplexMatrix m) {
                if(m instanceof ComplexMatrix) {
                        return multiply((ComplexMatrix)m);
                } else {
                        if(numCols==m.rows()) {
                                final double arrayRe[][]=new double[numRows][m.columns()];
                                final double arrayIm[][]=new double[numRows][m.columns()];
                                int n,k;
                                Complex elem;
                                for(int j=0;j<numRows;j++) {
                                        for(k=0;k<m.columns();k++) {
                                                elem=m.getElement(0,k);
                                                arrayRe[j][k]=(matrixRe[j][0]*elem.real() - matrixIm[j][0]*elem.imag());
                                                arrayIm[j][k]=(matrixIm[j][0]*elem.real() + matrixRe[j][0]*elem.imag());
                                                for(n=1;n<numCols;n++) {
                                                        elem=m.getElement(n,k);
                                                        arrayRe[j][k]+=(matrixRe[j][n]*elem.real() - matrixIm[j][n]*elem.imag());
                                                        arrayIm[j][k]+=(matrixIm[j][n]*elem.real() + matrixRe[j][n]*elem.imag());
                                                }
                                        }
                                }
                                if(numRows==m.columns())
                                        return new ComplexSquareMatrix(arrayRe,arrayIm);
                                else
                                        return new ComplexMatrix(arrayRe,arrayIm);
                        } else
                                throw new MatrixDimensionException("Incompatible matrices.");
                }
        }
        public AbstractComplexMatrix multiply(final ComplexMatrix m) {
                if(numCols==m.numRows) {
                        final double arrayRe[][]=new double[numRows][m.numCols];
                        final double arrayIm[][]=new double[numRows][m.numCols];
                        int n,k;
                        for(int j=0;j<numRows;j++) {
                                for(k=0;k<m.numCols;k++) {
                                        arrayRe[j][k]=(matrixRe[j][0]*m.matrixRe[0][k] - matrixIm[j][0]*m.matrixIm[0][k]);
                                        arrayIm[j][k]=(matrixIm[j][0]*m.matrixRe[0][k] + matrixRe[j][0]*m.matrixIm[0][k]);
                                        for(n=1;n<numCols;n++) {
                                                arrayRe[j][k]+=(matrixRe[j][n]*m.matrixRe[n][k] - matrixIm[j][n]*m.matrixIm[n][k]);
                                                arrayIm[j][k]+=(matrixIm[j][n]*m.matrixRe[n][k] + matrixRe[j][n]*m.matrixIm[n][k]);
                                        }
                                }
                        }
                        if(numRows==m.numCols)
                                return new ComplexSquareMatrix(arrayRe,arrayIm);
                        else
                                return new ComplexMatrix(arrayRe,arrayIm);
                } else
                        throw new MatrixDimensionException("Incompatible matrices.");
        }

// DIRECT SUM

        /**
        * Returns the direct sum of this matrix and another.
        */
        public AbstractComplexMatrix directSum(final AbstractComplexMatrix m) {
                final double arrayRe[][]=new double[numRows+m.numRows][numCols+m.numCols];
                final double arrayIm[][]=new double[numRows+m.numRows][numCols+m.numCols];
                for(int j,i=0;i<numRows;i++) {
                        for(j=0;j<numCols;j++) {
                                arrayRe[i][j]=matrixRe[i][j];
                                arrayIm[i][j]=matrixIm[i][j];
                        }
                }
                for(int j,i=0;i<m.numRows;i++) {
                        for(j=0;j<m.numCols;j++) {
                                Complex elem=m.getElement(i,j);
                                arrayRe[i+numRows][j+numCols]=elem.real();
                                arrayIm[i+numRows][j+numCols]=elem.imag();
                        }
                }
                return new ComplexMatrix(arrayRe,arrayIm);
        }

// TENSOR PRODUCT

        /**
        * Returns the tensor product of this matrix and another.
        */
        public AbstractComplexMatrix tensor(final AbstractComplexMatrix m) {
                final double arrayRe[][]=new double[numRows*m.numRows][numCols*m.numCols];
                final double arrayIm[][]=new double[numRows*m.numRows][numCols*m.numCols];
                for(int i=0;i<numRows;i++) {
                        for(int j=0;j<numCols;j++) {
                                for(int k=0;k<m.numRows;j++) {
                                        for(int l=0;l<m.numCols;l++) {
                                                Complex elem=m.getElement(k,l);
                                                arrayRe[i*m.numRows+k][j*m.numCols+l]=(matrixRe[i][j]*elem.real() - matrixIm[i][j]*elem.imag());
                                                arrayIm[i*m.numRows+k][j*m.numCols+l]=(matrixIm[i][j]*elem.real() + matrixRe[i][j]*elem.imag());
                                        }
                                }
                        }
                }
                return new ComplexMatrix(arrayRe,arrayIm);
        }

// HERMITIAN ADJOINT

        /**
        * Returns the hermitian adjoint of this matrix.
        * @return a complex matrix
        */
        public AbstractComplexMatrix hermitianAdjoint() {
                final double arrayRe[][]=new double[numCols][numRows];
                final double arrayIm[][]=new double[numCols][numRows];
                for(int j,i=0;i<numRows;i++) {
                        arrayRe[0][i]=matrixRe[i][0];
                        arrayIm[0][i]=-matrixIm[i][0];
                        for(j=1;j<numCols;j++) {
                                arrayRe[j][i]=matrixRe[i][j];
                                arrayIm[j][i]=-matrixIm[i][j];
                        }
                }
                return new ComplexMatrix(arrayRe,arrayIm);
        }

// CONJUGATE

        /**
        * Returns the complex conjugate of this matrix.
        * @return a complex matrix
        */
        public AbstractComplexMatrix conjugate() {
                final double arrayIm[][]=new double[numRows][numCols];
                for(int j,i=0;i<numRows;i++) {
                        arrayIm[i][0]=-matrixIm[i][0];
                        for(j=1;j<numCols;j++)
                                arrayIm[i][j]=-matrixIm[i][j];
                }
                return new ComplexMatrix(matrixRe,arrayIm);
        }

// TRANSPOSE

        /**
        * Returns the transpose of this matrix.
        * @return a complex matrix
        */
        public Matrix transpose() {
                final double arrayRe[][]=new double[numCols][numRows];
                final double arrayIm[][]=new double[numCols][numRows];
                for(int j,i=0;i<numRows;i++) {
                        arrayRe[0][i]=matrixRe[i][0];
                        arrayIm[0][i]=matrixIm[i][0];
                        for(j=1;j<numCols;j++) {
                                arrayRe[j][i]=matrixRe[i][j];
                                arrayIm[j][i]=matrixIm[i][j];
                        }
                }
                return new ComplexMatrix(arrayRe,arrayIm);
        }

// MAP ELEMENTS

        /**
        * Applies a function on all the matrix elements.
        * @param f a user-defined function
        * @return a complex matrix
        */
        public AbstractComplexMatrix mapElements(final ComplexMapping f) {
                final Complex array[][]=new Complex[numRows][numCols];
                for(int j,i=0;i<numRows;i++) {
                        array[i][0]=f.map(matrixRe[i][0],matrixIm[i][0]);
                        for(j=1;j<numCols;j++)
                                array[i][j]=f.map(matrixRe[i][j],matrixIm[i][j]);
                }
                return new ComplexMatrix(array);
        }
}
