package Hasegawa.matrix;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import RandomGenerator.Sfmt;
import Jama.SingularValueDecomposition;

public class Matrix {
	
	private  final int MAX_VECTOR_LENGTH = 1000;
	private  final int MAX_MATRIX_SIZE = 500;

	private  double[] vector = new double[MAX_VECTOR_LENGTH];
	private  double[][] matrix = new double[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];

	public Matrix(){
		
	}
	
	public Matrix(final int i){
		this.vector = new double[i];
		this.matrix = new double[i][i];
	}
	
	public ArrayList<Double> transStringListToDoubeList(ArrayList<String> str){
		ArrayList<Double> data = new ArrayList<Double>();
		for (int i = 0; i < str.size(); i++) {
			data.add(Double.parseDouble(str.get(i)));
		}
		return data;
	}
	
	public  void setvalue(double[][] A, double value) {
		for (int i = 0; i < A.length; i++)
			for (int j = 0; j < A[0].length; j++)
				A[i][j] = value;
	}

	public  void setDiagvalue(double[][] A, double value) {
		for (int i = 0; i < A.length; i++)
			A[i][i] = value;
	}

	public  void setvalue(double[] x, double value) {
		for (int i = 0; i < x.length; i++)
			x[i] = value;
	}

	public  void setvalue(int[] x, int value) {
		for (int i = 0; i < x.length; i++)
			x[i] = value;
	}
	
	public  void setGaussian(double[][] A, double mean, double sd) {
		Random r = new Random();
		for (int i = 0; i < A.length; i++)
			for (int j = 0; j < A[0].length; j++)
				A[i][j] = sd * (r.nextGaussian() - mean);
	}

	public  void setDiagGaussian(double[][] A, double mean, double sd) {
		Random r = new Random();
		for (int i = 0; i < A.length; i++)
			A[i][i] = sd * (r.nextGaussian() - mean);
	}

	public  void setGaussian(double[] x, double mean, double sd) {
		Random r = new Random();
		for (int i = 0; i < x.length; i++)
			x[i] = sd * (r.nextGaussian() - mean);
	}

	public  void setGaussian(double[][] A, double mean, double sd, Sfmt sf) {
		for (int i = 0; i < A.length; i++)
			for (int j = 0; j < A[0].length; j++)
				A[i][j] = sd * (sf.NextNormal() - mean);
	}

	public  void setDiagGaussian(double[][] A, double mean, double sd,
			Sfmt sf) {
		for (int i = 0; i < A.length; i++)
			A[i][i] = sd * (sf.NextNormal() - mean);
	}

	public  void setGaussian(double[] x, double mean, double sd, Sfmt sf) {
		for (int i = 0; i < x.length; i++)
			x[i] = sd * (sf.NextNormal() - mean);
	}

	public  void setRandomvalue(double[][] A, double upper, double lower) {
		for (int i = 0; i < A.length; i++)
			for (int j = 0; j < A[0].length; j++)
				A[i][j] = Math.random() * (upper - lower) + lower;
	}

	public  void setRandomvalue(double[] x, double upper, double lower) {
		for (int i = 0; i < x.length; i++)
			x[i] = Math.random() * (upper - lower) + lower;
	}

	public  void setRandomvalue(double[][] A, double upper, double lower,
			Sfmt sf) {
		for (int i = 0; i < A.length; i++)
			for (int j = 0; j < A[0].length; j++)
				A[i][j] = sf.NextUnif() * (upper - lower) + lower;
	}

	public  void setRow(double[][] A, double[] x, int r) {
		int n = A[r].length;
		for (int i = 0; i < n; i++) {
			A[r][i] = x[i];
		}
	}

	public  void zero(double[][] A) {
		int m = A.length;
		int n = A[0].length;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = 0;
			}
		}
	}

	public  void zero(int[][] A) {
		int m = A.length;
		int n = A[0].length;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = 0;
			}
		}
	}

	public  void zero(double[] x) {
		int length = x.length;
		for (int i = 0; i < length; i++) {
			x[i] = 0;
		}
	}

	public  void zero(int[] x) {
		int length = x.length;
		for (int i = 0; i < length; i++) {
			x[i] = 0;
		}
	}

	public  void setIdentityMatrix(double[][] I) {
		int n = I.length;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				I[i][j] = 0;
			}
			I[i][i] = 1;
		}
	}
	
	/**
	 * y = x;
	 */
	public void copy(double[] y, double[] x) {
		int length = y.length;
		if(x.length != y.length) System.err.println("Vector Length not Equal!!!");
		for (int i = 0; i < length; i++) {
			y[i] = x[i];
		}
	}

	/**
	 * y = A[, col];
	 */
	public void copy_A_col(double[] y, double[][] A, int col) {
		int length = y.length;
		if(A[0].length != y.length) System.err.println("Vector Length not Equal!!!");
		for (int i = 0; i < length; i++) {
			y[i] = A[i][col];
		}
	}

	/**
	 * y = A[row, ];
	 */
	public void copy_A_row(double[] y, double[][] A, int row) {
		int length = y.length;
		if(A.length != y.length) System.err.println("Vector Length not Equal!!!");
		for (int i = 0; i < length; i++) {
			y[i] = A[row][i];
		}
	}

	/**
	 * y = x;
	 */
	public void copy(int[] y, int[] x) {
		int length = y.length;
		if(x.length != y.length) System.err.println("Vector Length not Equal!!!");
		for (int i = 0; i < length; i++) {
			y[i] = x[i];
		}		
	}

	/**
	 * y = x;
	 */
	public void copy(String[] name, String[] name_) {
		int length = name.length;
		if(name.length != name_.length) System.err.println("Vector Length not Equal!!!");
		for (int i = 0; i < length; i++) {
			name[i] = name_[i];
		}
	}

	/**
	 * trace(A);
	 */
	public  double trace(double[][] A) {
		int m = A.length;
		double trace = 0;
		for (int i = 0; i < m; i++) {
			trace += A[i][i];
		}
		return trace;
	}

	/**
	 * A = B;
	 */
	public  void copy(double[][] A, double[][] B) {
		int m = A.length;
		if(m == 0) return;
		int n = A[0].length;
		if(A.length != B.length || A[0].length != B[0].length){
			System.err.println("Matrix Dimension not Equal!!!");
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = B[i][j];
			}
		}
	}

	/**
	 * A = B^T;
	 */
	public  void copyt(double[][] A, double[][] B) {
		int m = A.length;
		if(m == 0) return;
		int n = A[0].length;
		if(A.length != B[0].length || A[0].length != B.length){
			System.err.println("Matrix Dimension not Equal!!!");
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = B[j][i];
			}
		}
	}

	/**
	 * A = B;
	 */
	public void copy(double[][] A, double[] b) {
		int m = A.length;
		if(A.length != b.length){
			System.err.println("Matrix Dimension not Equal!!!");
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				A[i][j] = 0;
			}
			A[i][i] = b[i];
		}
	}

	/**
	 * A = B;
	 */
	public void copy_part(double[][] A, double[][] B, int m, int n) {
		if(A.length != B.length || A[0].length != B[0].length){
			System.err.println("Matrix Dimension not Equal!!!");
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = B[i][j];
			}
		}
	}

	/**
	 * y = A * x
	 */
	public void multAx(double[][] A, double[] x, double[] y) {
		for (int i = 0; i < A.length; i++) {
			y[i] = 0;
			for (int j = 0; j < A[0].length; j++) {
				y[i] += A[i][j] * x[j];
			}
		}
	}
	
	/**
	 * (return) = A * x
	 * 
	 */
	public double[] multAxReturn(double[][] A, double[] x) {
		double[] y = new double[A.length];
		int m = y.length;
		int n = x.length;

		for (int i = 0; i < m; i++) {
			double[] a = A[i];
			double value = 0;
			for (int j = 0; j < n; j++) {
				value += a[j] * x[j];
			}
			y[i] = value;
		}
		return y;
	}
	
	/**
	 * y = A(diag) * x
	 */
	public void multAx(double[] A, double[] x, double[] y) {
		int m = A.length;
		for (int i = 0; i < m; i++) {
			y[i]=A[i]*x[i];
		}
	}
	
	/**
	 * y = A^T * x
	 * 
	 */
	public  void multAtx(double[][] A, double[] x, double[] y) {
		int m = y.length;
		int n = x.length;

		for (int i = 0; i < m; i++) {
			double value = 0;
			for (int j = 0; j < n; j++) {
				value += A[j][i] * x[j];
			}
			y[i] = value;
		}
	}

	/**
	 * y = y + A * x
	 * 
	 */
	public  void multAddAx(double[][] A, double[] x, double[] y) {
		int m = y.length;
		int n = x.length;

		for (int i = 0; i < m; i++) {
			double[] a = A[i];
			double value = 0;
			for (int j = 0; j < n; j++) {
				value += a[j] * x[j];
			}
			y[i] += value;
		}
	}

	/**
	 * y = y + A^T * x
	 * 
	 */
	public  void multAddAtx(double[][] A, double[] x, double[] y) {
		int m = y.length;
		int n = x.length;

		for (int i = 0; i < m; i++) {
			double value = 0;
			for (int j = 0; j < n; j++) {
				value += A[j][i] * x[j];
			}
			y[i] += value;
		}
	}

	/**
	 * y = y + alpha * A * x
	 * 
	 */
	public  void multAddAx(double alpha, double[][] A, double[] x,
			double[] y) {
		int m = y.length;
		int n = x.length;

		for (int i = 0; i < m; i++) {
			double[] a = A[i];
			double value = 0;
			for (int j = 0; j < n; j++) {
				value += a[j] * x[j];
			}
			y[i] += alpha * value;
		}
	}

	/**
	 * y = y - A * x
	 * 
	 */
	public  void multSubAx(double[][] A, double[] x, double[] y) {
		int m = y.length;
		int n = x.length;

		for (int i = 0; i < m; i++) {
			double[] a = A[i];
			double value = 0;
			for (int j = 0; j < n; j++) {
				value += a[j] * x[j];
			}
			y[i] -= value;
		}
	}

	/**
	 * A = D - A
	 * 
	 */
	public  void SubDA(double[] D, double[][] A) {
		int m = A.length;
		int n = A[0].length;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] *= -1;
			}
		}
		for (int i = 0; i < m; i++) {
			A[i][i] += D[i];
		}
	}
	
	/**
	 * A = A - D
	 * 
	 */
	public void SubAD(double[][] A, double[] D) {
		int m = A.length;
		for (int i = 0; i < m; i++) {
			A[i][i]-=D[i];
		}
	}

	/**
	 * C = A * B, where A is diagonal
	 * 
	 */
	public  void multAB(double[] A, double[][] B, double[][] C) {
		int m = A.length;
		int n = B[0].length;

		for (int i = 0; i < m; i++) {
			double a = A[i];
			for (int j = 0; j < n; j++) {
				C[i][j] = a * B[i][j];
			}
		}
	}

	/**
	 * B = A * B, where A is diagonal
	 * 
	 */
	public  void multAB(double[] A, double[][] B) {
		int m = A.length;
		int n = B[0].length;

		for (int i = 0; i < m; i++) {
			double a = A[i];
			for (int j = 0; j < n; j++) {
				B[i][j] *= a;
			}
		}
	}

	/**
	 * A = A * B, where B is diagonal
	 * 
	 */
	public  void multAB(double[][] A, double[] B) {
		int m = B.length;
		int n = A.length;

		for (int i = 0; i < m; i++) {
			double b = B[i];
			for (int j = 0; j < n; j++) {
				A[j][i] *= b;
			}
		}
	}
	
	/**
	 * C = A * B, where B is diagonal
	 * 
	 */
	public  void multAB(double[][] A, double[] B, double[][] C) {
		int m = B.length;
		int n = A.length;
				
		for (int i = 0; i < m; i++) {
			double b = B[i];
			for (int j = 0; j < n; j++) {
				C[j][i]= A[j][i] * b;
			}
		}
	}
	
	/**
	 * C = diag(D) * x, where B is diagonal
	 * 
	 */
	public  void multDx(double[] D, double[] x, double[] C) {
		int n = D.length;
				
		for (int i = 0; i < n; i++) {
			C[i] = D[i] * x[i];
		}
	}

	/**
	 * temp = A * B
	 * C = temp
	 * 
	 */
	/*public  void multAB(double[][] A, double[][] B, double[][] C) {
		int m = A.length;
		int n = B.length;
		int o = B[0].length;
		double[][] D = new double[C.length][C[0].length];

		for (int i = 0; i < m; i++) {
			double[] a = A[i];
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += a[k] * B[k][j];
				}
				D[i][j] = value;
			}
		}
		this.copy(C, D);
	}*/
	
	/**
	 * C = A * B
	 * 
	 */
	public void multAB(double[][] A, double[][] B, double[][] C) {
		int cash = 32;
		int m = A.length;
		int n = B.length;
		int o = B[0].length;
		this.setvalue(C, 0);
		for (int ii = 0; ii < m; ii += cash) {
			for (int jj = 0; jj < o; jj += cash) {
				for (int kk = 0; kk < n; kk += cash) {
					int th_1 = cash;
					if(ii + th_1 > m) th_1 = m%cash;
					for (int i = 0; i < th_1; i++) {
						double[] ai = A[ii + i];
						double[] ci = C[ii + i];
						int th_2 = cash;
						if(jj + th_2 > o) th_2 = o%cash;
						for (int j = 0; j < th_2; j++) {
							double cij = ci[jj + j];
							int th_3 = cash;
							if(kk + th_3 > n) th_3 = n%cash;
							for (int k = 0; k < th_3; k++) {
								cij += ai[kk + k] * B[kk + k][jj + j];
							}
							ci[jj + j] = cij;
						}
					}
				}
			}
		}
	}

	/**
	 * C = C - A * B
	 * 
	 */
	public  void multSubAB(double[][] A, double[][] B, double[][] C) {
		int m = A.length;
		int n = B.length;
		int o = B[0].length;

		for (int i = 0; i < m; i++) {
			double[] a = A[i];
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += a[k] * B[k][j];
				}
				C[i][j] -= value;
			}
		}
	}

	/**
	 * C = C - D * A
	 * 
	 */
	public  void multSubDA(double[][] A, double[] D, double[][] C) {
		int m = A.length;
		int n = A[0].length;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] -= A[i][j] * D[i];
			}
		}
	}

	/**
	 * C = C - x * y^t
	 * 
	 */
	public  void multSubxyt(double[] x, double[] y, double[][] C) {
		int m = x.length;
		int n = y.length;

		for (int i = 0; i < m; i++) {
			double x_ = x[i];
			for (int j = 0; j < n; j++) {
				double value=x_*y[j];
				C[i][j] -= value;
			}
		}
	}

	/**
	 * C = A^T * B
	 * 
	 */
	public  void multAtB(double[][] A, double[][] B, double[][] C) {
		int m = A[0].length;
		int n = B.length;
		int o = B[0].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[k][i] * B[k][j];
				}
				C[i][j] = value;
			}
		}
	}

	/**
	 * C = A^T * B + C
	 * 
	 */
	public  void multAddAtB(double[][] A, double[][] B, double[][] C) {
		int m = A[0].length;
		int n = B.length;
		int o = B[0].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[k][i] * B[k][j];
				}
				C[i][j] += value;
			}
		}
	}

	/**
	 * C = A * B^T
	 * 
	 */
	public  void multABt(double[][] A, double[][] B, double[][] C) {
		int m = A.length;
		int n = B[0].length;
		int o = B.length;

		for (int i = 0; i < m; i++) {
			double[] a = A[i];
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += a[k] * B[j][k];
				}
				C[i][j] = value;
			}
		}
	}

	/**
	 * C = C + A * B^T
	 * 
	 */
	public  void multAddABt(double[][] A, double[][] B, double[][] C) {
		int m = A.length;
		int n = B[0].length;
		int o = B.length;

		for (int i = 0; i < m; i++) {
			double[] a = A[i];
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += a[k] * B[j][k];
				}
				C[i][j] += value;
			}
		}
	}

	/**
	 * C = C - A * B^T
	 * 
	 */
	public  void multSubABt(double[][] A, double[][] B, double[][] C) {
		int m = A.length;
		int n = B[0].length;
		int o = B.length;

		for (int i = 0; i < m; i++) {
			double[] a = A[i];
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += a[k] * B[j][k];
				}
				C[i][j] -= value;
			}
		}
	}

	/**
	 * y = y + A * B * x
	 * 
	 */
	public  void multAddABx(double[][] A, double[][] B, double[] x,	double[] y) {
		int m = A.length;
		int n = B.length;
		int o = B[1].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[i][k] * B[k][j];
				}
				y[i] += value * x[j];
			}
		}
	}

	/**
	 * y = A * B * x
	 * 
	 */
	public  void multABx(double[][] A, double[][] B, double[] x, double[] y) {
		int m = A.length;
		int n = B.length;
		int o = B[1].length;

		for (int i = 0; i < m; i++) {
			y[i] = 0;
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[i][k] * B[k][j];
				}
				y[i] += value * x[j];
			}
		}
	}

	/**
	 * y = y - A * B * x
	 * 
	 */
	public  void multSubABx(double[][] A, double[][] B, double[] x,	double[] y) {
		int m = A.length;
		int n = B.length;
		int o = B[0].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[i][k] * B[k][j];
				}
				y[i] -= value * x[j];
			}
		}
	}

	/**
	 * y = A * B^T * x + y
	 * 
	 */
	public  void multAddABtx(double[][] A, double[][] B, double[] x, double[] y) {
		int m = A.length;
		int n = B[1].length;
		int o = B.length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[i][k] * B[j][k];
				}
				y[i] += value * x[j];
			}
		}
	}

	/**
	 * D = A * B * C
	 * 
	 */
	public  void multABC(double[][] A, double[][] B, double[][] C, double[][] D) {
		int m = A.length;
		int n = B.length;
		int o = B[0].length;
		int p = C[0].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < p; j++) {
				D[i][j] = 0;
			}
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[i][k] * B[k][j];
				}
				for (int k = 0; k < p; k++) {
					D[i][k] += value * C[j][k];
				}
			}
		}
	}
	
	/**
	 * D = A * B * C
	 * 
	 */
	public  void multABC(double[][] A, double[] B, double[] C, double[][] D) {
		int m = B.length;
		int n = A.length;
				
		for (int i = 0; i < m; i++) {
			double b = B[i];
			double c = C[i];
			for (int j = 0; j < n; j++) {
				D[j][i]= A[j][i] * b * c;
			}
		}
	}

	/**
	 * D = D - A^T * B * C
	 * 
	 */
	public  void multSubAtBC(double[][] A, double[][] B, double[][] C, double[][] D) {
		int m = A[0].length;
		int n = B.length;
		int o = B[0].length;
		int p = C[0].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[k][i] * B[k][j];
				}
				for (int k = 0; k < p; k++) {
					D[i][k] -= value * C[j][k];
				}
			}
		}
	}

	/**
	 * D = A * B^T * C
	 * 
	 */
	public  void multABtC(double[][] A, double[][] B, double[][] C, double[][] D) {
		int m = A.length;
		int n = B[0].length;
		int o = B.length;
		int p = C[0].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < p; j++) {
				D[i][j] = 0;
			}
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[i][k] * B[j][k];
				}
				for (int k = 0; k < p; k++) {
					D[i][k] += value * C[j][k];
				}
			}
		}
	}

	/**
	 * D = A * B * C^T
	 * 
	 */
	public  void multABCt(double[][] A, double[][] B, double[][] C, double[][] D) {
		int m = A.length;
		int n = B.length;
		int o = B[0].length;
		int p = C.length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < p; j++) {
				D[i][j] = 0;
			}
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[i][k] * B[k][j];
				}
				for (int k = 0; k < p; k++) {
					D[i][k] += value * C[k][j];
				}
			}
		}
	}

	/**
	 * D = D - A * B * C^T
	 * 
	 */
	public  void multSubABCt(double[][] A, double[][] B, double[][] C, double[][] D) {
		int m = A.length;
		int n = B.length;
		int o = B[0].length;
		int p = C.length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[i][k] * B[k][j];
				}
				for (int k = 0; k < p; k++) {
					D[i][k] -= value * C[k][j];
				}
			}
		}
	}

	/**
	 * D = D + A * B * C^T
	 * 
	 */
	public  void multAddABCt(double[][] A, double[][] B, double[][] C, double[][] D) {
		int m = A.length;
		int n = B.length;
		int o = B[0].length;
		int p = C.length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < o; j++) {
				double value = 0;
				for (int k = 0; k < n; k++) {
					value += A[i][k] * B[k][j];
				}
				for (int k = 0; k < p; k++) {
					D[i][k] += value * C[k][j];
				}
			}
		}
	}

	/**
	 * A = A + B
	 * 
	 */
	public  void add(double[][] A, double[][] B) {
		int m = A.length;
		int n = A[0].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] += B[i][j];
			}
		}
	}

	/**
	 * A = A - B
	 * 
	 */
	public  void sub(double[][] A, double[][] B) {
		int m = A.length;
		int n = A[0].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] -= B[i][j];
			}
		}
	}

	/**
	 * x = x - y
	 * 
	 */
	public  void sub(double[] x, double[] y) {
		int m = x.length;

		for (int i = 0; i < m; i++) {
			x[i] -= y[i];
		}
	}

	/**
	 * z = x - y
	 * 
	 */
	public  void sub(double[] x, double[] y, double[] z) {
		int m = x.length;

		for (int i = 0; i < m; i++) {
			z[i] = x[i] - y[i];
		}
	}


	/**
	 * z = x - y
	 * 
	 */
	public  void sub(double[] x, double[] y, double[] w, double[] z) {
		int m = x.length;

		for (int i = 0; i < m; i++) {
			z[i] = x[i] - y[i] - w[i];
		}
	}


	/**
	 * A = A + d
	 * 
	 */
	public  void add(double[][] A, double[] d) {
		int m = A.length;

		for (int i = 0; i < m; i++) {
			A[i][i] += d[i];
		}
	}

	/**
	 * A = A + d
	 * 
	 */
	public  double[][] add_return(double[][] A, double[] d) {
		double[][] B = new double[A.length][A[0].length];
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				B[i][j] = A[i][j] + d[i];
			}
		}
		return(B);
	}
	
	/**
	 * x = x + y
	 * 
	 */
	public  void add(double[] x, double[] y) {
		int length = x.length;
		for (int i = 0; i < length; i++) {
			x[i] += y[i];
		}
	}

	/**
	 * x = x + alpha * y
	 * 
	 */
	public  void add(double alpha, double[] x, double[] y) {
		int length = x.length;
		for (int i = 0; i < length; i++) {
			x[i] += alpha * y[i];
		}
	}

	public  double dotProduct(double[] x, double[] y) {
		int length = x.length;
		double value = 0;
		for (int i = 0; i < length; i++) {
			value += x[i] * y[i];
		}
		return value;
	}

	public  double dotProduct(double[] x, int x_pos, double[] y, int y_pos, int n, int max) {
		double value = 0;
		for (int i = 0; i < n && i < max; i++) {
			value += x[x_pos + i] * y[y_pos + i];
		}
		return value;
	}

	
	/**
	 * D = (A^{-1} + B^{t} C B)^{-1}, where C is diagonal matrix. Note that A is
	 * regular matrix.
	 * 
	 */
	public  void inversionTheorem(double[][] A, double[][] B, double[] C, double[][] D) {
		int m = A.length;
		int n = C.length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				D[i][j] = A[i][j];
			}
		}
		for (int i = 0; i < n; i++) {
			double[] b = B[i];
			// double coef = 1.0/(1.0/C[i]+1+multxtAx(D,b));
			double coef = 1.0 / (1.0 / C[i] + multxtAx(D, b));
			for (int j = 0; j < m; j++) {
				double[] d = D[j];
				double value = 0;
				for (int k = 0; k < m; k++) {
					value += d[k] * b[k];
				}
				vector[j] = value;
			}
			for (int j = 0; j < m; j++) {
				for (int k = 0; k < m; k++) {
					D[j][k] -= coef * vector[j] * vector[k];
				}
			}
		}
	}
	
	/**
	 * D = (A^{-1} + C )^{-1}, where C is diagonal matrix. 
	 * Note that A is regular matrix.
	 * 
	 */
	public  void inversionTheorem(double[][] A, double[] C, double[][] D) {
		int m = A.length;
		int n = C.length;
		
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				D[i][j] = A[i][j];
			}
		}
		for (int i = 0; i < n; i++) {
			double coef = 1.0 / (1.0 / C[i] + D[i][i]);
			for (int j = 0; j < m; j++) {
				vector[j] = D[j][i];
			}
			for (int j = 0; j < m; j++) {
				for (int k = 0; k < m; k++) {
					D[j][k] -= coef * vector[j] * vector[k];
				}
			}
		}
	}

	/**
	 * x^T * A * x
	 */
	public  double multxtAx(double[][] A, double[] x) {
		int length = x.length;

		double value = 0;
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < length; j++) {
				value += x[i] * x[j] * A[i][j];
			}
		}
		return value;
	}
	
	/**
	 * A = x * y^T
	 * 
	 */
	public  void multxyt(double[] x, double[] y, double[][] A) {
		int xLength = x.length;
		int yLength = y.length;
		for (int i = 0; i < xLength; i++) {
			for (int j = 0; j < yLength; j++) {
				A[i][j] = x[i] * y[j];
			}
		}
	}

	/**
	 * A += x * y^T
	 * 
	 */
	public  void multAddxyt(double[] x, double[] y, double[][] A) {
		int xLength = x.length;
		int yLength = y.length;
		for (int i = 0; i < xLength; i++) {
			for (int j = 0; j < yLength; j++) {
				A[i][j] += x[i] * y[j];
			}
		}
	}

	/**
	 * A += D * x * y^T
	 * 
	 */
	public  void multAddDxyt(double[] x, double[] y, double[] D, double[][] A) {
		int xLength = x.length;
		int yLength = y.length;
		for (int i = 0; i < xLength; i++) {
			for (int j = 0; j < yLength; j++) {
				A[i][j] += D[i] * x[i] * y[j];
			}
		}
	}

	/**
	 * A += D1 * D2 * x * y^T
	 * 
	 */
	public  void multAddD2xyt(double[] x, double[] y, double[] D1, double[] D2, double[][] A) {
		int xLength = x.length;
		int yLength = y.length;
		for (int i = 0; i < xLength; i++) {
			for (int j = 0; j < yLength; j++) {
				A[i][j] += D1[i] * D2[i] * x[i] * y[j];
			}
		}
	}

	/**
	 * B += D * A
	 */
	public  void multAddDA(double[][] A, double[] D, double[][] B) {
		if(A.length != B.length) System.err.println("Error at Clculating multAddDA");
		if(A.length != D.length) System.err.println("Error at Clculating multAddDA");
		if(A[0].length != B[0].length) System.err.println("Error at Clculating multAddDA");
		for (int i = 0; i < B.length; i++) {
			for (int j = 0; j < B[0].length; j++) {
				B[i][j] += D[i] * A[i][j];				
			}
		}
	}
	
	/**
	 * A -= x * y^T
	 * 
	 */
	public  void multSubDxyt(double[] x, double[] y, double[] D, double[][] A) {
		int xLength = x.length;
		int yLength = y.length;
		for (int i = 0; i < xLength; i++) {
			for (int j = 0; j < yLength; j++) {
				A[i][j] -= D[i] * x[i] * y[j];
			}
		}
	}

	
	public void absolute(double[] v) {
		for (int i = 0; i < v.length; i++) {
			v[i] = Math.abs(v[i]);
		}
	}

	public  double Max(double[] v) {
		double max = -1 * Double.MAX_VALUE;
		for (double value : v) {
			if (value > max) {
				max = value;
			}
		}
		return max;
	}

	public  double absoluteMax(double[] v) {
		double max = 0;
		for (double value : v) {
			if (Math.abs(value) > max) {
				max = Math.abs(value);
			}
		}
		return max;
	}

	public  double absoluteMax(double[] v, double[] w) {
		double max = 0;
		int length = v.length;
		for (int i = 0; i < length; i++) {
			double value = Math.abs(v[i] / w[i]);
			if (value > max) {
				max = value;
			}
		}
		return max;
	}

	public  double l1norm(double[] v) {
		double norm = 0;
		for (double value : v) {
			norm += Math.abs(value);
		}
		return norm;
	}

	public  void rescale(double[] vector, double alpha) {
		int length = vector.length;
		for (int i = 0; i < length; i++) {
			vector[i] *= alpha;
		}
	}

	/**
	 * dest = alpha*vector
	 */
	public  void rescale(double[] vector, double alpha, double[] dest) {
		int length = vector.length;
		for (int i = 0; i < length; i++) {
			dest[i] = alpha * vector[i];
		}
	}

	public  void rescale(double[][] A, double alpha) {
		int m = A.length;
		if(m==0) return;
		int n = A[0].length;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] *= alpha;
			}
		}
	}

	public  void rescale(double[] x, double[] alpha) {
		int m = x.length;

		for (int i = 0; i < m; i++) {
			x[i] *= alpha[i];
		}
	}

	public  void symmetricSimultaneousEquation(double[][] A, double[] y,
			double[][] LD, int length) {
		if (length == 0) {
			return;
		}
		modifiedCholeskyDecomposition(A, LD, length);
		double value;
		for (int r = 0; r < length; r++) {
			value = y[r];
			for (int i = 0; i < r; i++) {
				value -= LD[r][i] * y[i];
			}
			y[r] = value;
		}
		for (int r = 0; r < length; r++) {
			if (LD[r][r] == 0) {
				y[r] = 0;
			} else {
				y[r] = y[r] / LD[r][r];
			}
		}
		for (int r = length - 1; r >= 0; r--) {
			value = y[r];
			for (int i = length - 1; i > r; i--) {
				value -= LD[i][r] * y[i];
			}
			y[r] = value;
		}
	}
	
	public  void symmetricInverseWithWS(double[][] A, double[][] ws) {
		int length = A.length;
		modifiedCholeskyDecomposition(A, ws, length);
		double value;

		setIdentityMatrix(A);
		for (double[] y : A) {
			for (int r = 0; r < length; r++) {
				value = y[r];
				for (int i = 0; i < r; i++) {
					value -= ws[r][i] * y[i];
				}
				y[r] = value;
			}
			for (int r = 0; r < length; r++) {
				if (ws[r][r] == 0) {
					y[r] = 0;
				} else {
					y[r] = y[r] / ws[r][r];
				}
			}
			for (int r = length - 1; r >= 0; r--) {
				value = y[r];
				for (int i = length - 1; i > r; i--) {
					value -= ws[i][r] * y[i];
				}
				y[r] = value;
			}
		}
		changesymmetric(A);
	}


	public  void symmetricInverse(double[][] A, double[][] Ai) {
		int length = A.length;
		double[][] LD = matrix;
		modifiedCholeskyDecomposition(A, LD, length);
		double value;

		setIdentityMatrix(Ai);
		for (double[] y : Ai) {
			for (int r = 0; r < length; r++) {
				value = y[r];
				for (int i = 0; i < r; i++) {
					value -= LD[r][i] * y[i];
				}
				y[r] = value;
			}
			for (int r = 0; r < length; r++) {
				if (LD[r][r] == 0) {
					y[r] = 0;
				} else {
					y[r] = y[r] / LD[r][r];
				}
			}
			for (int r = length - 1; r >= 0; r--) {
				value = y[r];
				for (int i = length - 1; i > r; i--) {
					value -= LD[i][r] * y[i];
				}
				y[r] = value;
			}
		}
		changesymmetric(Ai);
	}

	public  void symmetricInverse(double[][] A) {
		int length = A.length;
		double[][] LD = matrix;
		modifiedCholeskyDecomposition(A, LD, length);
		double value;

		setIdentityMatrix(A);
		for (double[] y : A) {
			for (int r = 0; r < length; r++) {
				value = y[r];
				for (int i = 0; i < r; i++) {
					value -= LD[r][i] * y[i];
				}
				y[r] = value;
			}
			for (int r = 0; r < length; r++) {
				if (LD[r][r] == 0) {
					y[r] = 0;
				} else {
					y[r] = y[r] / LD[r][r];
				}
			}
			for (int r = length - 1; r >= 0; r--) {
				value = y[r];
				for (int i = length - 1; i > r; i--) {
					value -= LD[i][r] * y[i];
				}
				y[r] = value;
			}
		}
		changesymmetric(A);
	}

	public  void symmetricInverse(double[][] A, int length) {
		double[][] LD = matrix;
		modifiedCholeskyDecomposition(A, LD, length);
		double value;

		setIdentityMatrix(A);
		for (double[] y : A) {
			for (int r = 0; r < length; r++) {
				value = y[r];
				for (int i = 0; i < r; i++) {
					value -= LD[r][i] * y[i];
				}
				y[r] = value;
			}
			for (int r = 0; r < length; r++) {
				if (LD[r][r] == 0) {
					y[r] = 0;
				} else {
					y[r] = y[r] / LD[r][r];
				}
			}
			for (int r = length - 1; r >= 0; r--) {
				value = y[r];
				for (int i = length - 1; i > r; i--) {
					value -= LD[i][r] * y[i];
				}
				y[r] = value;
			}
		}
		changesymmetric(A);
	}

	public  double symmetricDeterminant(double[][] A, int length) {
		double[][] LD = matrix;
		modifiedCholeskyDecomposition(A, LD, length);
		double determinant = 1;
		for (int i = 0; i < length; i++) {
			determinant *= LD[i][i];
		}
		return determinant;
	}

	public  double symmetricLogDeterminant(double[][] A, int length) {
		double[][] LD = matrix;
		modifiedCholeskyDecomposition(A, LD, length);
		double determinant = 0;
		for (int i = 0; i < length; i++) {
			determinant += Math.log(LD[i][i]);
		}
		return determinant;
	}

	public  double symmetricLogDeterminant(double[][] A) {
		double[][] LD = matrix;
		int length = A.length;
		modifiedCholeskyDecomposition(A, LD, length);
		double determinant = 0;
		for (int i = 0; i < length; i++) {
			determinant += Math.log(LD[i][i]);
		}
		return determinant;
	}

	public  double symmetricDeterminant(double[][] A) {
		int length = A.length;
		double[][] LD = matrix;
		modifiedCholeskyDecomposition(A, LD, length);
		double determinant = 1;
		for (int i = 0; i < length; i++) {
			determinant *= LD[i][i];
		}
		return determinant;
	}

	public  void modifiedCholeskyDecomposition(double[][] A, double[][] LD) {
		modifiedCholeskyDecomposition(A, LD, A.length);
	}

	public  void modifiedCholeskyDecomposition(double[][] A,
			double[][] LD, int length) {
		if (length == 0) {
			return;
		}
		double value;
		for (int i = 0; i < length; i++) {
			for (int j = 0; j <= i; j++) {
				value = A[i][j];
				for (int k = 0; k < j; k++) {
					value -= LD[i][k] * LD[j][k] * LD[k][k];
				}
				if (j != i) {
					value /= LD[j][j];
				}
				LD[i][j] = value;
			}
		}
	}

	public  boolean isIdentity(double[][] A) {
		int m = A.length;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				if (i == j) {
					if (Math.abs(A[i][i] - 1.0) >= epsilon) {
						return false;
					}
				} else {
					if (Math.abs(A[i][j]) >= epsilon) {
						return false;
					}
				}
			}
		}
		return true;
	}

	private  double epsilon = 0.001;

	public  boolean isEqual(double[][] A, double[][] B) {
		if (A.length != B.length) {
			return false;
		}
		if (A[0].length != B[0].length) {
			return false;
		}
		int m = A.length;
		int n = A[0].length;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Math.abs(A[i][j] - B[i][j]) >= epsilon) {
					return false;
				}
			}
		}
		return true;
	}

	public  void print(double[][] A) {
		for (double[] a : A) {
			for (int i = 0; i < a.length; i++) {
				// System.err.print(((int)(a[i]+0.000001)));
				System.err.print(a[i]);
				if (i < a.length - 1) {
					System.err.print(" ");
				}
			}
			System.err.println();
		}
	}

	public  void print(double[] x) {
		for (int i = 0; i < x.length; i++) {
			// System.err.print(((int)(a[i]+0.000001)));
			System.err.print(x[i]);
			if (i < x.length - 1) {
				System.err.print(" ");
			}
		}
		System.err.println();
	}

	public  void print(int[] x) {
		for (int i = 0; i < x.length; i++) {
			System.err.print(x[i]);
			if (i < x.length - 1) {
				System.err.print(" ");
			}
		}
		System.err.println();
	}

	public  void print(int[][] A) {
		for (int[] a : A) {
			for (int i = 0; i < a.length; i++) {
				System.err.print(a[i]);
				if (i < a.length - 1) {
					System.err.print(" ");
				}
			}
			System.err.println();
		}
	}

	/**
	 * Mean: 0, Variance: 1;
	 */
	public void normalizeRow(double[][] A) {
		int row = A.length;
		int col = A[0].length;
		
		for (int i = 0; i < row; i++) {
			double mean = 0;
			double mean2 = 0;
			for (int j = 0; j < col; j++) {
				mean += A[i][j];
				mean2 += A[i][j] * A[i][j];
			}
			double std = Math.sqrt((mean2/col) - ((mean/col)*(mean/col)));
			for (int j = 0; j < col; j++) {
				A[i][j] = (A[i][j] - (mean/col)) / std;
			}
		}
	}

	/**
	 * Mean: 0, Variance: 1;
	 */
	public void normalizeCol(double[][] A) {
		int row = A.length;
		int col = A[0].length;
		
		for (int j = 0; j < col; j++) {
			double mean = 0;
			double mean2 = 0;
			for (int i = 0; i < row; i++) {
				mean += A[i][j];
				mean2 += A[i][j] * A[i][j];
			}
			double std = Math.sqrt((mean2/row) - ((mean/row)*(mean/row)));
			for (int i = 0; i < row; i++) {
				A[i][j] = (A[i][j] - (mean / row)) / std;
			}
		}
	}

	/**
	 * Mean: 0, Variance: 1;
	 */
	public  void normalize(double[] x) {
		int n = x.length;

		double mean = 0;
		double mean2 = 0;
		for (int d = 0; d < n; d++) {
			mean += x[d];
			mean2 += x[d] * x[d]; 
		}
		
		double std = Math.sqrt((mean2 / n) - ((mean / n) * (mean / n)));
		for (int d = 0; d < n; d++) {
			x[d] = (x[d] - mean) / std;
		}
	}

	public  void normalizeTimeSerise(double[][] A) {
		int p = A[0].length;
		int T = A.length;

		for (int i = 0; i < p; i++) {
			double mean = 0;
			double var = 0;
			int S = 0;
			for (int d = 0; d < T; d++) {
				if (A[d] == null) {
					continue;
				}
				mean += A[d][i];
				S++;
			}
			mean /= S;
			for (int d = 0; d < T; d++) {
				if (A[d] == null) {
					continue;
				}
				A[d][i] -= mean;
				var += A[d][i] * A[d][i];
			}
			var /= S;

			double std = Math.sqrt(var);
			for (int d = 0; d < T; d++) {
				if (A[d] == null) {
					continue;
				}
				A[d][i] /= std;
			}
		}
	}

	/*
	 * only symmetric
	 */
	public  void transpose(double[][] A) {
		int m = A.length;
		for (int i = 0; i < m; i++) {
			for (int j = i + 1; j < m; j++) {
				double swap = A[i][j];
				A[i][j] = A[j][i];
				A[j][i] = swap;
			}
		}
	}

	/**
	 * B=At
	 */
	public  void transpose(double[][] A, double[][] B) {
		int m = A.length;
		int n = A[0].length;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				B[j][i] = A[i][j];
			}
		}
	}

	public  boolean checker(double[][] A) {

		int m = A.length;

		double epsilon = 0.0001;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				if (i == j) {
					if (A[i][i] <= 0) {
						return true;
					}
				} else {
					if (Math.abs(A[i][j] - A[j][i]) > epsilon) {
						System.err.println(A[i][j] + ", " + A[j][i]);
						return true;
					}
				}

			}
		}
		return false;
	}

	public  void changesymmetric(double[][] A) {
		for (int i = 0; i < A.length - 1; i++) {
			for (int j = i + 1; j < A[0].length; j++) {
				double value = (A[i][j] + A[j][i]) / 2;
				A[i][j] = A[j][i] = value;
			}
		}
	}

	public void DiagInvese(double[] x) {
		for (int i = 0; i < x.length; i++)
			x[i] = 1 / x[i];
	}

	/*
	 * only implemented about A=[N][N]
	 */
	public  void LUDecomposition(double[][] A, double[][] L, double[][] U) {
		int i, j, k;
		int N = A.length;
		double T;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
				L[i][j] = U[i][j] = 0;
		L[0][0] = 1;
		for (j = 0; j < N; j++)
			U[0][j] = A[0][j];
		for (i = 1; i < N; i++) {
			U[i][0] = L[0][i] = 0;
			L[i][0] = A[i][0] / U[0][0];
		}
		for (i = 1; i < N; i++) {
			L[i][i] = 1;
			T = A[i][i];
			for (k = 0; k <= i; k++)
				T -= L[i][k] * U[k][i];
			U[i][i] = T;
			for (j = i + 1; j < N; j++) {
				U[j][i] = L[i][j] = 0;
				T = A[j][i];
				for (k = 0; k <= i; k++)
					T -= L[j][k] * U[k][i];
				L[j][i] = T / U[i][i];
				T = A[i][j];
				for (k = 0; k <= i; k++)
					T -= L[i][k] * U[k][j];
				U[i][j] = T;
			}
		}
	}

	public  double diagDeterminant(double D[]) {
		double value = 1;
		for (int i = 0; i < D.length; i++)
			value *= D[i];
		return value;
	}

	public  double diagLogDeterminant(double D[]) {
		double value = 0;
		for (int i = 0; i < D.length; i++)
			value += Math.log(Math.abs(D[i]));
		return value;
	}

	public  double diagLogDeterminant(double D1[], double D2[]) {
		double value = 0;
		for (int i = 0; i < D1.length; i++)
			value += Math.log(Math.abs(D1[i]));
		for (int i = 0; i < D2.length; i++)
			value += Math.log(Math.abs(D2[i]));
		return value;
	}

	public  double diagInvLogDeterminant(double D1[]) {
		double value = 0;
		for (int i = 0; i < D1.length; i++)
			value += Math.log(Math.abs(1.0 / D1[i]));
		return value;
	}

	public  double diagLogDeterminant(double D[][]) {
		double value = 0;
		for (int i = 0; i < D.length; i++)
			value += Math.log(Math.abs(D[i][i]));
		return value;
	}

	public  void QRDecomposition(double[][] A, double[][] Q, double[][] R, int N) {
		int i, j, k;
		double Rjk, T;
		double[] X = new double[N];
		for (k = 0; k < N; k++) {
			for (i = 0; i < N; i++)
				X[i] = A[i][k];
			for (j = 0; j <= k - 1; j++) {
				Rjk = 0;
				for (i = 0; i < N; i++)
					Rjk += A[i][k] * Q[i][j];
				R[j][k] = Rjk;
				R[k][j] = 0;
				for (i = 0; i < N; i++)
					X[i] -= Rjk * Q[i][j];
			}
			T = 0;
			for (i = 0; i < N; i++)
				T += X[i] * X[i];
			R[k][k] = Math.sqrt(T);
			for (j = 0; j < N; j++)
				Q[j][k] = X[j] / R[k][k];
		}
	}
	
	/*
	 * A = US(=diag)V';
	 */
	public  void SVDecomposition(double[][] A, double[][] U, double[][] V, double[] S){
		Jama.Matrix mA = new Jama.Matrix(A);
		Jama.Matrix mU = new Jama.Matrix(U);
		Jama.Matrix mV = new Jama.Matrix(V);
		double[] mS = this.copy_generate(S);
		Jama.SingularValueDecomposition SVD = new SingularValueDecomposition(mA);
		mU = SVD.getU();
		mV = SVD.getV();
		mS = SVD.getSingularValues();
		for (int i = 0; i < mU.getRowDimension(); i++) {
			for (int j = 0; j < mU.getColumnDimension(); j++) {
			U[i][j]	= mU.get(i, j);
			}			
		}
		for (int i = 0; i < mV.getRowDimension(); i++) {
			for (int j = 0; j < mV.getColumnDimension(); j++) {
			V[i][j]	= mV.get(i, j);
			}			
		}
		for (int i = 0; i < mS.length; i++) {
			S[i] = mS[i];
		}		
	}
	
	/*
	 * A = LL';
	 */
	public  void CholeskyDecomposition(double[][] A, double[][] L){
		Jama.Matrix mA = new Jama.Matrix(A);
		Jama.Matrix mL = new Jama.Matrix(L);
		Jama.CholeskyDecomposition CD = new Jama.CholeskyDecomposition(mA);
		mL = CD.getL();
		for (int i = 0; i < mL.getRowDimension(); i++) {
			for (int j = 0; j < mL.getColumnDimension(); j++) {
			L[i][j]	= mL.get(i, j);
			}			
		}		
	}
	
	/*
	 * A = LDL';
	 */
	public  void ImprovedCholeskyDecomposition(double[][] A, double[][] L, double[] D){
		Jama.Matrix mA = new Jama.Matrix(A);
		Jama.CholeskyDecomposition CD = new Jama.CholeskyDecomposition(mA);
		for (int i = 0; i < CD.getL().getRowDimension(); i++) {
			double eigen = CD.getL().get(i, i);
			D[i] = eigen * eigen;
			for (int j = 0; j < CD.getL().getColumnDimension(); j++) {
				L[j][i] = CD.getL().get(j, i) / eigen;
			}
		}
	}
	
	/*
	 * sqrtA * sqrtA = A;
	 */
	public int getSquareRoot(double[][] A, double[][] sqrtA){
		if(this.checkNaN(A) || this.checkAbsValues(A, 1.0e100)){
			System.err.println("A includes NaN. System terminated.");
			return 1;
			//System.exit(0);
		}
		Jama.Matrix mA = new Jama.Matrix(A);
		Jama.SingularValueDecomposition SVD = new SingularValueDecomposition(mA);
		double[][] U = new double[A.length][A[0].length];
		double[][] V = new double[A.length][A[0].length];
		double[][] S = new double[A.length][A[0].length];
		for (int i = 0; i < SVD.getU().getRowDimension(); i++) {
			for (int j = 0; j < SVD.getU().getColumnDimension(); j++) {
				U[i][j]	= SVD.getU().get(i, j);
			}			
		}
		for (int i = 0; i < SVD.getV().getRowDimension(); i++) {
			for (int j = 0; j < SVD.getV().getColumnDimension(); j++) {
				V[i][j]	= SVD.getV().get(i, j);
			}			
		}
		for (int i = 0; i < SVD.getSingularValues().length; i++) {
			S[i][i] = Math.sqrt(SVD.getSingularValues()[i]);
		}
		this.multABCt(U, S, V, sqrtA);
		return 0;
	}
	
	/*
	 * x^T = x^T*D
	 */
	public  void multxtD(double[] x, double[] D){
		for (int i = 0; i < D.length; i++) 
			x[i]=x[i]*D[i];
	}
	
	/*
	 * y^T = x^T*A
	 */
	public void multxtA(double[] x, double[][] A, double[] y){
		for (int i = 0; i < A[0].length; i++) {
			y[i]=0;
			for (int j = 0; j < A.length; j++) {
				y[i]+=x[j]*A[j][i];
			}
		}
	}
	
    public  void write(String fileName, double[][] A) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);

            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A[0].length; j++) {
                	writer.write(String.format("%.15f", A[i][j]));
                	if(j+1 < A[0].length)writer.write("\t");
                }
                writer.write(String.format("\n"));
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    } 
    
    public  void write(String fileName, double[] v) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);

            for (int i = 0; i < v.length; i++) {
                writer.write(String.format("%.15f", v[i]));
                writer.write(String.format("\n"));
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    }
    
    public  void write(String fileName, double v) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);
            writer.write(String.format("%.15f", v));
            writer.write(String.format("\n"));
            
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    }
    
    public double sumofVector(double[] v){
    	double s = 0;
    	for (double i : v) 
    		s+=i;
    	return s;
    }
    
    public void sumofRow(double[][] A, double[] s){
    	for (int i = 0; i < s.length; i++) {
        	s[i] = 0;
    		for (double x : A[i]) {
        		s[i]+=x;
        	}
		}
    }
    
    public double[] sumofColumn(double[][] A){
    	double[] s = new double[A[0].length];
    	for (int i = 0; i < A[0].length; i++) {
			double sum = 0;
			for (int j = 0; j < A.length; j++) {
				sum += A[j][i]; 
			}
			s[i] = sum;
		}
    	return s;
    }
    
    public void prune(double[][] A, double th){
    	for (int i = 0; i < A.length; i++) {
			this.prune(A[i], th);
		}
    }
    
    public void prune(double[] v, double th){
    	for (int i = 0; i < v.length; i++) {
			if(Math.abs(v[i])<th) v[i]=0;
		}
    }
    
    /*
	 * Get and Return sub-matrix of matrix A at Index's row and column. 
	 * (ex) Index = {1,2,4}
	 * B = A[Index][Index]
	 */
	public void getSubMatrix(double[][] A, ArrayList<Integer> Index, double[][] B){
		for (int i = 0; i < B.length; i++) {
			for (int j = 0; j < B[0].length; j++) {
				B[i][j]=A[Index.get(i)][Index.get(j)];
			}
		}
	}	
	
	/*
	 * Get and Return sub-matrix of matrix A at Index's row and column. 
	 * (ex) Index = {1,2,4}
	 * B = A[Index][Index]
	 */
	public void getSubMatrix(double[][] A, int[] Index, double[][] sub){
		for (int i = 0; i < sub.length; i++) {
			for (int j = 0; j < sub[0].length; j++) {
				sub[i][j]=A[Index[i]][Index[j]];
			}
		}
	}	
	
	
	/*
	 * Get and Return sub-matrix of matrix A at Index's column. 
	 * (ex) Index = {1,2,4}
	 */
	public double[] getColVector(double[][] A, int Index){
		double[] re = new double[A.length];
		for (int i = 0; i < re.length; i++) {
			re[i] = A[i][Index];
		}
		return re;
	}
	
	/*
	 * Get and Return sub-matrix of matrix A at Index's column. 
	 * (ex) Index = {1,2,4}
	 */
	public double[][] getPartOfColVectors(double[][] A, ArrayList<Integer> Index){
		double[][] re = new double[A.length][Index.size()];
		for (int i = 0; i < re.length; i++) {
			for (int j = 0; j < re[0].length; j++) {
				re[i][j]=A[i][Index.get(j)];
			}
		}
		return re;
	}
	
	/*
	 * Get and Return sub-matrix of matrix A at Index's row. 
	 * (ex) Index = {1,2,4}
	 */
	public double[][] getPartOfRowVectors(double[][] A, ArrayList<Integer> Index){
		double[][] re = new double[Index.size()][A[0].length];
		for (int i = 0; i < re.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				re[i][j]=A[Index.get(i)][j];
			}
		}
		return re;
	}	
	
	/**
	 * x^T * A * x
	 */
	public  double multxtDx(double[] D, double[] x) {
		int length = x.length;

		double value = 0;
		for (int i = 0; i < length; i++) {
				value += x[i] * x[i] * D[i];
		}
		return value;
	}
	
	/**
	 * x^T * (A + epsilon)* x
	 */
	public  double multxtDxEp(double[] D, double[] x, double epsilon) {
		int length = x.length;

		double value = 0;
		for (int i = 0; i < length; i++) {
				value += x[i] * x[i] * (D[i] + epsilon);
		}
		return value;
	}
	
	public void multxx(double[] x, double[] xx){
		for (int i = 0; i < x.length; i++) {
			xx[i] = x[i] * x[i];
		}
	}
	
	public double[][] copy_generate(double[][] M){
		if(M.length ==0) {
			//System.err.println("Earning(Calculator.copy): Matrix is null!!");
			return new double[0][0];
		}
		double[][] C = new double[M.length][M[0].length];
		for (int i = 0; i < C.length; i++) {
			for (int j = 0; j < C[0].length; j++) {
				C[i][j] = M[i][j];
			}
		}
		return C;
	}
	
	public double[] copy_generate(double[] x){
		double[] v = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			v[i] = x[i];
		}
		return v;
	}
	
	public int[] copy_generate(int[] x){
		int[] v = new int[x.length];
		for (int i = 0; i < x.length; i++) {
			v[i] = x[i];
		}
		return v;
	}
	
	public double[][] ArrayListToMatrix(ArrayList<double[]> array){
		double[][] M = new double[array.size()][array.get(0).length];
		for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				M[i][j] = array.get(i)[j];
			}
		}
		return M;
	}
	
	public double[] ArrayListToArrayDouble(ArrayList<Double> array){
		double[] v = new double[array.size()];
		for (int i = 0; i < v.length; i++) {
			v[i] = array.get(i);
		}
		return v;
	}
	
	public int[] ArrayListToArrayInt(ArrayList<Integer> array){
		int[] v = new int[array.size()];
		for (int i = 0; i < v.length; i++) {
			v[i] = array.get(i);
		}
		return v;
	}
		
	public boolean checkAbsValues(double[][] A, double th){
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				if(Math.abs(A[i][j]) > th) 
					return true;
			}
		}
		return false;
	}
	
	public boolean checkAbsValues(double[][] A, double up_th, double low_th){
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				if(Math.abs(A[i][j]) > up_th || (Math.abs(A[i][j]) < low_th)) return true;
			}
		}
		return false;
	}
	
	public boolean checkAbsValues(double[] x, double th){
		for (int i = 0; i < x.length; i++) {
			if(Math.abs(x[i]) > th) return true;
		}
		return false;
	}
	
	public boolean checkAbsValues(double[] x, double up_th, double low_th){
		for (int i = 0; i < x.length; i++) {
			if(Math.abs(x[i]) > up_th || (Math.abs(x[i]) < low_th)) return true;
		}
		return false;
	}

	public boolean checkValues(double[] x, double up_th, double low_th){
		for (int i = 0; i < x.length; i++) {
			if(x[i] > up_th || x[i] < low_th) return true;
		}
		return false;
	}
	
	/*
	 * Check whether it includes NaN or not. 
	 */
	public boolean checkNaN(double[] v){
		for (int i = 0; i < v.length; i++) {
			if(Double.isNaN(v[i])) return true;
		}
		return false;
	}
	
	/*
	 * Check whether it includes NaN or not. 
	 */
	public boolean checkNaN(double[][] A){
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				if(Double.isNaN(A[i][j])) 
					return true;
			}
		}
		return false;
	}
	
	public int getNumberOfEn(int time, int row, int particle, int dimension, double maxparticle){
		return (int)((time * maxparticle * dimension) + (particle * dimension) + row);
	}

	public void sqrt(double[] D, double[] Original) {
		for (int i = 0; i < D.length; i++) {
			D[i] = Math.sqrt(Original[i]);
		}
	}
		
	public void display(double[][] A){
		for (int i = 0; i < A.length; i++) {
			System.out.println();
			for (int j = 0; j < A[0].length; j++) {
				if(j < A[0].length - 1) System.out.print(A[i][j] + "\t");
				else System.out.print(A[i][j]);
			}
		}
		System.out.println();
	}
	
	public void display(double[] v){
		for (int i = 0; i < v.length; i++) {
			System.out.print(v[i] + ", ");
		}
		System.out.println();
	}

	public String[] copy_generate(String[] str) {
		String[] temp =  new String[str.length];
		this.copy(temp, str);
		return temp;
	}
	
	/**
	 * 
	 * @param R Original Triangular Matrix
	 * @param Rinv Inverse Matrix (Must be Initialized)
	 */
	public void triangularInverse(double[][] R, double[][] Rinv){
	    double bij;
	    for (int k = 0; k < R[0].length; k++) {
			for (int i = 0; i < R.length; i++) {
				for (int j = 0; j <= i; j++) {
					if (i != j + k) continue;
					/* Diagonal Element */
				    if (k == 0){
				        if (i == j) Rinv[i][j] = 1.0 / R[i][j];
				    } else {
				        bij = 0;
				        for (int l = j; l < i; l++){
				            bij += R[i][l] * Rinv[l][j];
				        }
				        bij *= -1.0 / R[i][i];
				        Rinv[i][j] = bij;
				    }
			    }
			}
		}
	}
	
	public void replace(double[][] A, int i1, int i2){
		if(i1 == i2) return;
		double temp;
		for (int j = 0; j < A[0].length; j++) {
			temp = A[i1][j];
			A[i1][j] = A[i2][j];
			A[i2][j] = temp;
		}
		for (int j = 0; j < A.length; j++) {
			temp = A[j][i1];
			A[j][i1] = A[j][i2];
			A[j][i2] = temp;
		}
	}
	
	public void replace2(double[][] A, int i1, int i2){
		if(i1 == i2) return;
		double temp;
		for (int j = 0; j < A[0].length; j++) {
			temp = A[i1][j];
			A[i1][j] = A[i2][j];
			A[i2][j] = temp;
		}
	}
	
	public void kernelABC(final double[] observedData, final double[] simulatedData, int NumberOfParticles, double variance, int max_reduced_dimension, 
			double epsilon, double sigma, int multiple_particle, int system_dimension,double eps_N, int repSize, 
			double[] copy, int[] p, double[] r_vec, double[] rii, double[] current_kernel, double[] current_weight){
		/* Initialize */
		for (int i = 0; i < p.length; i++) {
			p[i] = i;
		}
		
		/* 
		 * Note that in incompleteCholeskyDecompositionUsingGramMatrix(), 
		 * input vector will be overwritten. Copy Observation. 
		 */
		this.copy(copy, simulatedData);
		
		/* Incomplete Cholesky Decomposition */
		this.setvalue(r_vec, 0);
		this.setvalue(rii, 0);
		int fin_m = this.incompleteCholeskyDecompositionUsingGramMatrix
					(copy, r_vec, p, rii, epsilon, NumberOfParticles, max_reduced_dimension, sigma, 
							system_dimension * multiple_particle, system_dimension);
		final double[][] fin_r = sMatrix.sMP.returnICDecomposedMatrix(r_vec, NumberOfParticles, fin_m, max_reduced_dimension);
		System.out.println("Row-Rank: " + fin_m);
		
		/* Sorting */
		for (int i = 0; i < p.length; i++) {
			if(p[i] != i){
				for (int j = i + 1; j < p.length; j++) {
					if(i == p[j]){
						if(i != j){
							double temp;
							for (int j2 = 0; j2 < fin_r[0].length; j2++) {
								temp = fin_r[i][j2];
								fin_r[i][j2] = fin_r[j][j2];
								fin_r[j][j2] = temp;
							}
						}
						p[j] = p[i];
						p[i] = i;
						break;
					}
				}
			}
		}
		
		/* Wood-Bury and Kernel Bayes */
		double[][] RRt = new double[fin_m][fin_m];
		this.multtAA(fin_r, RRt);
		for (int i = 0; i < RRt.length; i++) {
			RRt[i][i] += eps_N;
		}
		this.symmetricInverse(RRt);
		double[][] RGram = new double[NumberOfParticles][fin_m];
		this.multAB(fin_r, RRt, RGram);
		for (int j = 0; j < NumberOfParticles; j++) {
			double total_sum = 0;
			double sum = 0;
			for (int rep = 0; rep < repSize; rep++) {
				for (int i = 0; i < system_dimension * multiple_particle; i++) {
					sum += (simulatedData[j * system_dimension * multiple_particle  + i] - observedData[rep * system_dimension + i % system_dimension])
							* (simulatedData[j * system_dimension * multiple_particle  + i] - observedData[rep * system_dimension + i % system_dimension]);
					if(i % system_dimension == system_dimension - 1) {
						total_sum += Math.exp(-1 * sum / (2 * sigma));
						sum= 0;
					}
				}
			}
			current_kernel[j] = total_sum / repSize * multiple_particle;
			current_weight[j] = current_kernel[j] * (1 - dotProduct(RGram[j], fin_r[j]));
		}
		
	    for (int i = 0; i < system_dimension; i++) {
	    	for (int j = i + 1; j < system_dimension; j++) {
	    		double temp = this.dotProduct(RGram[i], fin_r[j]);
	    		current_weight[NumberOfParticles + i] -=  temp * current_kernel[j];
	    		current_weight[NumberOfParticles + j] -=  temp * current_kernel[i];
	    	}
	    }
	
	    /* Calculate Current Weight */
		int minus = 0;
		for (int i = 0; i < NumberOfParticles; i++) {
			current_weight[i] = current_weight[i] / eps_N;
			if(current_weight[i] <= 0) {
				current_weight[i] = 0;
				minus++;
			}
		}
		System.out.println(minus);
		
		/* Update Weight */
		this.rescale(current_weight, 1.0 / this.sumofVector(current_weight));
	}
	
	
	/**
	 * 
	 * @param A Vector for Gram Matrix
	 * @param R Vector for Row-rank Matrix
	 * @param P Transformation Vector
	 * @param Rii Diagonal Elements of Row-rank Matrix
	 * @param eps Threshold
	 * @param d a + diagonal(d) = Original Matrix
	 * @param n Number of Row in a
	 * @param m Number of Column of r
	 */
	public int incompleteCholeskyDecompositionUsingGramMatrix(double[] a, double[] R,
			int[] P, double[] Rii, double eps, final int n, final int m, double sigma, int dim, int elementNum) {
		/* Initialize */
		int index;
		double sum, temp;
		for (int j = 0; j < m; j++) {
			R[j * m + j] = Rii[j] = this.forKernelBayes(a, j * dim, j * dim, dim, sigma, elementNum);
		}
		for (int j = m; j < n; j++) {
			Rii[j] = this.forKernelBayes(a, j * dim, j * dim, dim, sigma, elementNum);
		}
		for (int j = n; j < Rii.length; j++) {
			Rii[j] = 0;
		}
		
		/* Main Routine */
		for (int i = 0; i < m; i++) {
			/* Check Tolerance */
			sum = Rii[i];
			for (int j = i + 1; j < n; j++) {
				sum += Rii[j];
			}
			//System.out.println(sum);
			if(sum < eps) {
				for (int j = i; j < m; j++) {
					R[j * m + j] = 0;
				}
				return i;
			}
			//System.err.println(sum);
			
			/* Next Evaluation Row */
			index = i;
			for (int j = i; j < n; j++) {
				if (Math.abs(Rii[j]) > Math.abs(Rii[index])) index = j;
			}
			
			/* Replace a */
			for (int j = 0; j < dim; j++) {
				temp = a[i * dim + j];
				a[i * dim + j] = a[index * dim + j];
				a[index * dim + j] = temp;
			}
			
			/* Replace R */
			temp = Rii[index];
			Rii[index] = Rii[i];
			if(m > index) R[index * m + index] = Rii[i];
			Rii[i] = R[i * m + i] = temp;
			for (int j = 0; j < i; j++) {
				temp = R[index * m + j];
				R[index * m + j] = R[i * m + j];
				R[i * m + j] = temp;
			}
			
			/* Replace P */
			temp = P[index];
			P[index] = P[i];
			P[i] = (int)temp;
			
			/* Calculate ith Column */
			Rii[i] = R[i * m + i] = Math.sqrt(Rii[i]);
			for (int j = i + 1; j < n; j++) {
				R[j * m + i] = this.forKernelBayes(a, j * dim, i * dim, dim, sigma, elementNum);
				for (int k = 0; k < i; k++) {
					R[j * m + i] -= R[j * m + k] * R[i * m + k];
				}
				R[j * m + i] /= R[i * m + i];
				Rii[j] = this.forKernelBayes(a, j * dim, j * dim, dim, sigma, elementNum);
				for (int k = 0; k < i + 1; k++) {
					Rii[j] -= R[j * m + k] * R[j * m + k];
				}
				if(j < m) R[j * m + j] = Rii[j];
			}
			
			/* Update Diagonal Element */
			for (int j = i + 1; j < m; j++) {
				R[j * m + j] = this.forKernelBayes(a, j * dim, j * dim, dim, sigma, elementNum);
				for (int k = 0; k < i + 1; k++) {
					R[j * m + j] -= R[j * m + k] * R[j * m + k];
				}
				Rii[j] = R[j * m + j];
			}
			for (int j = m; j < n; j++) {
				Rii[j] = this.forKernelBayes(a, j * dim, j * dim, dim, sigma, elementNum);
				for (int k = 0; k < i + 1; k++) {
					Rii[j] -= R[j * m + k] * R[j * m + k];
				}
			}
		}
		return m;
	}
	
	/**
	 * 
	 * @param A Original Matrix (n * n)
	 * @param R A = RR' (n * m)
	 * @param P Transformation Vector (n)
	 * @param Rii Diagonal Elements of R (n)
	 * @param eps Tolerance Parameter
	 */
	public int incompleteCholeskyDecomposition(final double[][] A, double[][] R, 
			int[] P, double[] Rii, final double eps) {
		/* Initialize */
		int n = A.length;
		int m = R[0].length;
		int index;
		double sum, temp;
		for (int j = 0; j < m; j++) {
			R[j][j] = Rii[j] = A[j][j];
		}
		for (int j = m; j < n; j++) {
			Rii[j] = A[j][j];
		}
		
		/* Main Routine */
		for (int i = 0; i < m; i++) {
			/* Check Tolerance */
			sum = Rii[i];
			for (int j = i + 1; j < n; j++) {
				sum += Rii[j];
			}
			//System.out.println(sum);
			if(sum < eps) {
				for (int j = i; j < m; j++) {
					R[j][j] = 0;
				}
				return i;
			}
			
			/* Next Evaluation Row */
			index = i;
			for (int j = i; j < n; j++) {
				if (Math.abs(Rii[j]) > Math.abs(Rii[index])) index = j;
			}
			
			/* Replace A */
			this.replace(A, index, i);
			
			/* Replace R */
			temp = Rii[index];
			Rii[index] = Rii[i];
			if(m > index) R[index][index] = Rii[i];
			Rii[i] = R[i][i] = temp;
			for (int j = 0; j < i; j++) {
				temp = R[index][j];
				R[index][j] = R[i][j];
				R[i][j] = temp;
			}
			
			/* Replace P */
			temp = P[index];
			P[index] = P[i];
			P[i] = (int)temp;
			
			/* Calculate ith Column */
			R[i][i] = Rii[i] = Math.sqrt(Rii[i]);
			for (int j = i + 1; j < n; j++) {
				R[j][i] = A[j][i];
				for (int k = 0; k < i; k++) {
					R[j][i] -= R[j][k] * R[i][k];
				}
				R[j][i] /= R[i][i];
			}
			
			/* Update Diagonal Element */
			for (int j = i + 1; j < m; j++) {
				R[j][j] = A[j][j];
				for (int k = 0; k < i + 1; k++) {
					R[j][j] -= R[j][k] * R[j][k];
				}
				Rii[j] = R[j][j];
			}
			for (int j = m; j < n; j++) {
				Rii[j] = A[j][j];
				for (int k = 0; k < i + 1; k++) {
					Rii[j] -= R[j][k] * R[j][k];
				}
			}
		}
		return m;
	}
	

	
	public double[][] returnICDecomposedMatrix(double[] R, int row, int col, int ColOfOriginalR){
		double[][] d = new double[row][col];
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				d[i][j] = R[i * ColOfOriginalR + j];
			}
		}
		return d;
	}
	
	public double forKernelBayes(double[] a, int x, int y, int dim, double sigma, int elementNum){
		return this.gaussKernel(a, x, y, dim, sigma, elementNum);
	}
	
	public double l2Kernel(double[] a, int x, int y, int dim, double sigma, int elementNum){
		double sum = 0;
		for (int i = 0; i < dim; i++) {
			sum += (a[x + i] - a[y + i]) * (a[x + i] - a[y + i]);
		}
		return Math.exp(-1 * Math.sqrt(sum) / (2 * sigma));
	}

	public double expKernel(double[] a, int x, int y, int dim, double sigma, int elementNum){
		double sum = 0;
		for (int i = 0; i < dim; i++) {
			sum += Math.abs(a[x + i] - a[y + i]);
		}
		return Math.exp(-1 * sum / (2 * sigma));
	}
	
	public double gaussKernel(double[] a, int x, int y, int dim, double sigma, int elementNum){
		double total_sum = 0;
		for (int m1 = 0; m1 < dim / elementNum; m1++) {
			for (int m2 = 0; m2 < dim / elementNum; m2++) {
				double sum = 0;
				for (int i = 0; i < elementNum; i++) {
					sum += (a[x + elementNum * m1 + i] - a[y + elementNum * m2  + i])
						* (a[x + elementNum * m1 + i] - a[y + elementNum * m2  + i]);
				}
				total_sum += Math.exp(-1 * sum / (2 * sigma));
			}
		}
		return total_sum / (dim / elementNum * dim / elementNum);
	}

	public void multtAA(double[][] A, double[][] C) {
		int n = A[0].length;
		int m = A.length;

		for (int i = 0; i 
				< n; i++) {
			for (int j = i; j < n; j++) {
				double value = 0;
				for (int k = 0; k < m; k++) {
					value += A[k][i] * A[k][j];
				}
				C[i][j] = C[j][i] = value;
			}
		}
	}
	
	public void multtAAscale(double[][] A, double eps, double[][] C) {
		int n = A[0].length;
		int m = A.length;

		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				double value = 0;
				for (int k = 0; k < m; k++) {
					value += A[k][i] * A[k][j];
				}
				C[i][j] = C[j][i] = value * eps;
			}
		}
	}

	public double getDiagSum(double[][] A) {
		double x = 0;
		for (int i = 0; i < A.length; i++) {
			x += A[i][i];
		}
		return x;
	}

	public double sumofVector(int[] v) {
		int s = 0;
	    	for (double i : v) s += i;
	    	return s;
	}

    public double sumofVector(Double[] v){
	    	double s = 0;
	    	for (double i : v) 
	    		s+=i;
	    	return s;
    }

	public int sumofVector(Integer[] v) {
	    	int s = 0;
	    	for (int i : v) 
	    		s += i;
	    	return s;
	}

	public void copy(double[][][] A, double[][][] B) {
		for (int t = 0; t < A.length; t++) {
			for (int i = 0; i < A[0].length; i++) {
				for (int j = 0; j < A[0][0].length; j++) {
					A[t][i][j]= B[t][i][j]; 
				}
			}
		}
	}
}
