package Hasegawa.stat;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import Hasegawa.matrix.Matrix;
import RandomGenerator.Sfmt;

public class simpleMath {

	public simpleMath() {

	}
	
	public <T> List<T> getUniqueElement(List<T> arg0) {
		List<T> ret = new ArrayList<T>();
		for (int i = 0; i < arg0.size(); i++) {
			T x = arg0.get(i);
			if (!ret.contains(x)) {
				ret.add(x);
			}
		}
		return ret;
	}
	
	public int[] match(String[] x, String[] y){
		int[] data = new int[x.length];
		for (int i = 0; i < data.length; i++) {
			data[i] = -1;
			for (int j = 0; j < data.length; j++) {
				if(x[i]==y[j]) {data[i] = j;break;}
			}
		}
		return data;
	}
	
	public int[] remove(int remove, int[] from){
		ArrayList<Integer> data = new ArrayList<Integer>();
		for (int i = 0; i < from.length; i++) {
			if(from[i]!=remove)data.add(from[i]);
		}
		int[] data2 = new int[data.size()];
		for (int i = 0; i < data2.length; i++) {
			data2[i] = data.get(i);
		}
		return data2;
	}

	public int getSumOfLength(double[] a, double[] b) {
		int length = 0;
		if (a != null)
			length += a.length;
		if (b != null)
			length += b.length;
		return length;
	}

	@SuppressWarnings("rawtypes")
	public int getSumOfLength(ArrayList a, ArrayList b) {
		int length = 0;
		if (a != null)
			length += a.size();
		if (b != null)
			length += b.size();
		return length;
	}
	
	public int getSumOfRowLength(double[][] A, double[][] B) {
		int length = 0;
		if (A != null)
			length += A.length;
		if (B != null)
			length += B.length;
		return length;
	}
	
	public int getSumOfColLength(double[][] A, double[][] B) {
		int length = 0;
		if (A != null)
			length += A[0].length;
		if (B != null)
			length += B[0].length;
		return length;
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public ArrayList getMergedList(ArrayList a, ArrayList b) {
		ArrayList re = new ArrayList();
		if (a != null){
			for (Object o1 : a) {
				re.add(o1);
			}
		}
		if (b != null){
			for (Object o2 : b) {
				re.add(o2);
			}
		}
		return re;
	}
	
	public double getMean(ArrayList<Double> v) {
		double mean = 0;
		for (int i = 0; i < v.size(); i++) {
			mean += v.get(i);
		}
		return mean/v.size();
	}
	
	public double getMean(double[] v) {
		return getMean(v, 0, v.length - 1);
	}
		
	public double getMean(double[] v, int min_index, int max_index) {
		double sum = 0;
		for (int i = min_index; i <= max_index; i++) {
			sum += v[i];
		}
		sum /= (max_index - min_index + 1);
		return sum;
	}
	
	public double[] getMeanOfRow(double[][] A) {
		double[] ave = new double[A.length];
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				ave[i] += A[i][j] / A.length;
			}
		}
		return ave;
	}
	
	public double getMeanOfRow(double[][] A, int row) {
		double m = 0;
		for (int j = 0; j < A[0].length; j++) {
			m += A[row][j];
		}
		return m/A[0].length;
	}
	
	public double[] getMeanOfCol(double[][] A) {
		double[] ave = new double[A[0].length];
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				ave[j] += A[i][j] / A[0].length;
			}
		}
		return ave;
	}
	
	public double getMeanOfCol(double[][] A, int col) {
		double m = 0;
		for (int i = 0; i < A.length; i++) {
			m += A[i][col];
		}
		return m/A.length;
	}
	
	public double getMedianFromSamples(double[] array, int sampling, Sfmt Sf){
		double[] sample = new double[sampling];
		for (int i = 0; i < sample.length; i++) {
			sample[i] = array[(int)(Sf.NextUnif() * 0.9999 * array.length)];
		}
		this.sort(sample, 0, sample.length - 1);
		return sample[sample.length / 2];
	}
	
	public double getMedianFromSamples(double[] array, int sampling, Sfmt Sf, int point, int dimension){
		double[] sample = new double[sampling];
		for (int i = 0; i < sample.length; i++) {
			sample[i] = array[((int)(Sf.NextUnif() * 0.9999 * (array.length/dimension)) * dimension) + point];
		}
		this.sort(sample, 0, sample.length - 1);
		return sample[sample.length / 2];
	}
	
	public double getWidth(double[] v, int point, int dimension, double maxWidth) {
		double min = Double.MAX_VALUE;
		double max = (-1) * Double.MAX_VALUE;
		for (int i = 0; i < v.length/dimension; i++) {
			if(min > v[(dimension * i) + point]) min = v[(dimension * i) + point];
			if(max < v[(dimension * i) + point]) max = v[(dimension * i) + point];
		}
		if(max -min > maxWidth) return maxWidth;
		else return (max - min);
	}
	
	public double getUnbiased(int length){
		return (double)length/((double)length - 1.0);
	}

	/**
	 * variance = (x - mean)*(x - mean)/n, n = the amount of elements in v
	 */
	public double getVariance(double[] v, double mean) {
		double var = 0;
		for (int i = 0; i < v.length; i++) {
			var += (v[i] - mean) * (v[i] - mean);
		}
		return var / v.length;
	}
	
	public double getVariance(ArrayList<Double> v) {
		double m1 = 0;
		double m2 = 0;
		for (int i = 0; i < v.size(); i++) {
			m1 += v.get(i);
			m2 += v.get(i) * v.get(i);
		}
		return (m2 / v.size()) - (m1 * m1 /(v.size() * v.size()));
	}
	
	/**
	 * variance = (x - mean)*(x - mean)/n n = the amount of elements in v
	 */
	public double getVariance(ArrayList<Double> v, double mean) {
		double var = 0;
		for (int i = 0; i < v.size(); i++) {
			var += (v.get(i) - mean) * (v.get(i) - mean);
		}
		return var / v.size();
	}
	
	public double[] getVarianceOfCol(double[][] A){
		double[] variance = new double[A[0].length];
		for (int j = 0; j < variance.length; j++) {
			double mean = 0;
			double mean2 = 0;
			for (int i = 0; i < A.length; i++) {
				mean += A[i][j];
				mean2 += A[i][j] * A[i][j];
			}
			variance[j] = (mean2/A.length - ((mean/A.length)*(mean/A.length)));
		}		
		return variance;
	}
	
	public double getVarianceOfCol(double[][] A, int col){
		double mean = 0;
		double mean2 = 0;
		for (int i = 0; i < A.length; i++) {
			mean += A[i][col];
			mean2 += A[i][col] * A[i][col];
		}
		return (mean2/A.length - ((mean/A.length)*(mean/A.length)));
	}
	
	public double[] getVarianceOfRow(double[][] A){
		double[] variance = new double[A.length];
		for (int i = 0; i < variance.length; i++) {
			double mean = 0;
			double mean2 = 0;
			for (int j = 0; j < A[0].length; j++) {
				mean += A[i][j];
				mean2 += A[i][j] * A[i][j];
			}
			variance[i] = (mean2/A.length - ((mean/A.length)*(mean/A.length)));
		}		
		return variance;
	}
	
	public double getVarianceOfRow(double[][] A, int row){
		double mean = 0;
		double mean2 = 0;
		for (int j = 0; j < A[0].length; j++) {
			mean += A[row][j];
			mean2 += A[row][j] * A[row][j];
		}
		return (mean2/A.length - ((mean/A.length)*(mean/A.length)));
	}

	public double[] normlization(double[] y) {
		double[] y2 = new double[y.length];
		double mean = 0;
		double mean2 = 0;
		for (int i = 0; i < y.length; i++) {
			mean += y[i];
			mean2 += y[i] * y[i];
		}
		mean /= y.length;
		mean = mean * mean;
		mean2 /= y.length;

		
		double var = mean2 - mean;
		if (var < 1.0e-3) {
			mean = 0;
			var = 1;
		}
		for (int i = 0; i < y.length; i++) {
			y2[i] = (y[i] - mean) / Math.sqrt(var);
		}
		return y2;
	}
	
	/**
	 * mean(A[][j])=0, variance(A[][j])=1,
	 */
	public double[][] normlizationColumn(double[][] A) {
		double[][] A2 = new double[A.length][A[0].length];
		for (int j = 0; j < A[0].length; j++) {
			double mean = 0;
			double mean2 = 0;
			for (int i = 0; i < A.length; i++) {
				mean += A[i][j];
				mean2 += A[i][j] * A[i][j];
			}
			mean /= A.length;
			mean = mean * mean;
			mean2 /= A.length;
			double var = mean2 - mean;
			if (var == 0) {
				mean = 0;
				var = 1;
			}
			for (int i = 0; i < A.length; i++) {
				A2[i][j] = (A[i][j] - mean) / Math.sqrt(var);
			}
		}
		return A2;
	}
	
	/**
	 * mean(A[i][])=0, variance(A[i][])=1,
	 */
	public double[][] normlizationRow(double[][] A) {
		double[][] A2 = new double[A.length][A[0].length];
		for (int i = 0; i < A.length; i++) {
			A2[i] = this.normlization(A[i]);
		}
		return A2;
	}

	public int getMaxAbsValueIndex(double[] v) {
		int max = 0;
		for (int i = 0; i < v.length; i++) {
			if (Math.abs(v[i]) > Math.abs(v[max]))
				max = i;
		}
		return max;
	}

	public int getMinAbsValueIndex(double[] v) {
		int min = 0;
		for (int i = 0; i < v.length; i++) {
			if (Math.abs(v[i]) < Math.abs(v[min]))
				min = i;
		}
		return min;
	}

	public int getMaxValueIndex(double[] v) {
		int max = 0;
		for (int i = 0; i < v.length; i++) {
			if (v[i] > v[max])
				max = i;
		}
		return max;
	}
	
	public int getMaxValueIndex(double[] v, double th) {
		int max = 0;
		for (int i = 0; i < v.length; i++) {
			if (v[i] > v[max] && v[i] < th)
				max = i;
		}
		return max;
	}
	
	public int getMaxValueIndex(ArrayList<Double> v) {
		int max = 0;
		for (int i = 0; i < v.size(); i++) {
			if (v.get(i) > v.get(max))
				max = i;
		}
		return max;
	}
	
	public int getMaxValueIndex(ArrayList<double[]> v, int col) {
		int max = 0;
		for (int i = 0; i < v.size(); i++) {
			if (v.get(i)[col] > v.get(max)[col])
				max = i;
		}
		return max;
	}
	
	public int getMaxValueIndex(double[] v, int end) {
		double max = (-1) * Double.MAX_VALUE;
		int maxIndex = -1;
		for (int i = 0; i < end; i++) {
			if (v[i] > max){
				maxIndex = i;
				max = v[i];
			}
		}
		return maxIndex;
	}
	
	public double getMaxValue(double[] v) {
		double max = (-1) * Double.MAX_VALUE;
		for (int i = 0; i < v.length; i++) {
			if (v[i] > max){
				max = v[i];
			}
		}
		return max;
	}
	
	public double getMaxValue(double[] v, int end) {
		double max = (-1) * Double.MAX_VALUE;
		for (int i = 0; i < end; i++) {
			if (v[i] > max){
				max = v[i];
			}
		}
		return max;
	}
	
	public int getMinValueIndex(ArrayList<Double> v) {
		int min = 0;
		for (int i = 0; i < v.size(); i++) {
			if (v.get(i) < v.get(min))
				min = i;
		}
		return min;
	}
	
	public int getMinValueIndex(double[] v) {
		int min = 0;
		for (int i = 0; i < v.length; i++) {
			if (v[i] < v[min])
				min = i;
		}
		return min;
	}

	/*
	 * Get nth max absolute value's element
	 * setting 1 returns the highest value index
	 */
	public int getNthMaxAbsValueIndexA(double[] v, int n) {
		int max = 0;
		if(n == 0) {
			n = 1;
			System.err.println("n must be more than 0");
		}
		double[] vv = new double[v.length];
		for (int i = 0; i < v.length; i++) {
			vv[i] = Math.abs(v[i]);
		}
		for (int nn = 0; nn < n; nn++) {
			max = 0;
			for (int i = 0; i < v.length; i++) {
				if (vv[i] > vv[max])
					max = i;
			}
			vv[max] = 0;
		}
		return max;
	}
	
	/*
	 * Get nth min absolute value's element
	 * setting 1 returns the highest value index
	 */
	public int getNthMinAbsValueIndex(double[] v, int n) {
		int min = 0;
		if(n == 0) {
			n = 1;
			System.err.println("n must be more than 0");
		}
		double[] vv = new double[v.length];
		for (int i = 0; i < v.length; i++) {
			vv[i] = Math.abs(v[i]);
		}
		for (int nn = 0; nn < n; nn++) {
			min = 0;
			for (int i = 0; i < v.length; i++) {
				if (vv[i] < vv[min])
					min = i;
			}
			vv[min] = Double.MAX_VALUE;
		}
		return min;
	}

	/* 
	 * Get nth max value's element
	 * setting 1 returns the highest value index
	 */
	public int getNthMinAbsValueIndex(double[] v, int n, double without) {
		int min = 0;
		if(n == 0) {
			n = 1;
			System.err.println("n must be more than 0");
		}
		double[] vv = new double[v.length];
		for (int i = 0; i < v.length; i++) {
			vv[i] = Math.abs(v[i]);
		}
		for (int nn = 0; nn < n; nn++) {
			min = 0;
			for (int i = 0; i < v.length; i++) {
				if ((vv[i] < vv[min] && vv[i] != without) | vv[min] == without)
					min = i;
			}
			vv[min] = Double.MAX_VALUE;
		}
		return min;
	}

	/*
	 * Get nth max value's element
	 * setting 1 returns the highest value index
	 */
	public int getNthMaxValueIndex(double[] v, int n) {
		int max = 0;
		if(n == 0) {
			n = 1;
			System.err.println("n must be more than 0");
		}
		double[] vv = new double[v.length];
		for (int i = 0; i < v.length; i++) {
			vv[i] = v[i];
		}
		for (int nn = 0; nn < n; nn++) {
			max = 0;
			for (int i = 0; i < v.length; i++) {
				if (vv[i] > vv[max])
					max = i;
			}
			vv[max] = (-1) * Double.MAX_VALUE;
		}
		return max;
	}
	
	/*
	 * Get nth max value's element
	 * setting 1 returns the highest value index
	 */
	public int getNthMinValueIndex(double[] v, int n) {
		int min = 0;
		if(n == 0) {
			n = 1;
			System.err.println("n must be more than 0");
		}
		double[] vv = new double[v.length];
		for (int i = 0; i < v.length; i++) {
			vv[i] = v[i];
		}
		for (int nn = 0; nn < n; nn++) {
			min = 0;
			for (int i = 0; i < v.length; i++) {
				if (vv[i] < vv[min])
					min = i;
			}
			vv[min] = (1) * Double.MAX_VALUE;
		}
		return min;
	}
	
	public double sumAbsValueVector(double[] v) {
		double sum = 0;
		for (int i = 0; i < v.length; i++) {
			sum += Math.abs(v[i]);
		}
		return sum;
	}

	public int sumAbsValueVector(int[] v) {
		int sum = 0;
		for (int i = 0; i < v.length; i++) {
			sum += Math.abs(v[i]);
		}
		return sum;
	}

	public double[] sumAbsValueMatrix(double[][] A) {
		double[] sum = new double[A.length];
		for (int i = 0; i < A.length; i++) {
			sum[i] = this.sumAbsValueVector(A[i]);
		}
		return sum;
	}

	public ArrayList<Integer> getActiveRowElements(int row, double[][] A) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i < A[0].length; i++) {
			// if(i==row) continue;
			if (A[row][i] != 0)
				list.add(i);
		}
		return list;
	}

	public ArrayList<Integer> getNonActiveRowElements(int row, double[][] A) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i < A[0].length; i++) {
			// if(i==row) continue;
			if (A[row][i] == 0)
				list.add(i);
		}
		return list;
	}

	public int getIndexOf(double val, double[] array) {
		for (int i = 0; i < array.length; i++) {
			if (val == array[i])
				return i;
		}
		return -1;
	}

	public int getIndexOfAbs(double val, double[] array) {
		for (int i = 0; i < array.length; i++) {
			if (Math.abs(val) == Math.abs(array[i]))
				return i;
		}
		return -1;
	}

	public int[] getRandomOrder(int min, int max, Sfmt s) {
		if (min == max)
			return new int[] { min };
		int[] array = new int[max - min + 1];
		double[] array1 = new double[max - min + 1];
		double[] array2 = new double[max - min + 1];
		for (int i = 0; i < array1.length; i++) {
			array1[i] = s.NextUnif();
			array2[i] = array1[i];
		}
		this.sort(array1, 0, max - min);
		for (int i = 0; i < array1.length; i++) {
			for (int j = 0; j < array2.length; j++) {
				if (array1[i] == array2[j]) {
					array[i] = j;
					break;
				}
			}
		}
		return array;
	}

	/*
	 * Get the number of active elements
	 */
	public int getActiveCount(double[][] M) {
		int c = 0;
		if(M == null) return c;
		for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				if (M[i][j] != 0)
					c++;
			}
		}
		return c;
	}
	
	public int getActiveCount(double[] v) {
		int c = 0;
		for (int i = 0; i < v.length; i++) {
			if (v[i] != 0)
				c++;
		}
		return c;
	}
	
	public int getActiveRowCount(double[][] M1, double[][] M2, int row) {
		int c = 0;
		if(M1 != null && M1.length >0){
			for (int j = 0; j < M1[0].length; j++) {
				if (M1[row][j] != 0)
					c++;
			}
		}
		if(M2 != null && M2.length >0){
			for (int j = 0; j < M2[0].length; j++) {
				if (M2[row][j] != 0)
					c++;
			}
		}
		return c;
	}

	/**
	 * Must be Positive
	 */
	public int getEncodedIntegerList(ArrayList<Integer> list) {
		if(list.size() == 0) {
			System.err.println("List Size is zero!!");
			return -1;
		}
		int re = list.get(0);
		if (list.size() == 1)
			return re;
		int count0 = 0;
		if (re == 0)
			count0++;
		for (int i = 1; i < list.size(); i++) {
			if (re == 0 && list.get(i) == 0) {
				count0++;
			} else if (re == 0 && list.get(i) > 0) {
				re += Math.pow(10, count0) * list.get(i);
			} else {
				re += Math.pow(10, Math.ceil(Math.log10(re + 0.1)))
						* list.get(i);
			}
		}
		return re;
	}
	
	public double squrePredictionError(double[][] X, double[] y, double[] beta){
		double error = 0;
		double[] y_pre = new double[y.length];
		for (int i = 0; i < y.length; i++) {
			for (int j = 0; j < beta.length; j++) {
				y_pre[i] += X[i][j] * beta[j];
			}
			error += (y[i] - y_pre[i]) * (y[i] - y_pre[i]);
		}
		return error;
	}
	
	public double squrePredictionError(double[] y, double[] y_pre){
		double error = 0;
		for (int i = 0; i < y.length; i++) {
			error += (y[i] - y_pre[i]) * (y[i] - y_pre[i]);
		}
		return error;
	}

	
	public double squreSum(double[] y){
		double error = 0;
		for (int i = 0; i < y.length; i++) {
			error += y[i] * y[i];
		}
		return error;
	}

	/*
	 * MLE
	 */
	public double calculateLogLikelihoodOfMLE(double[] residual, Matrix MP) {
		double v2 = MP.dotProduct(residual, residual) / residual.length;
		return (-0.5 * MP.dotProduct(residual, residual) / v2) - (0.5 * residual.length * Math.log(2 * Math.PI * v2));
	}

	/*
	 * Ridge Regression
	 */
	public double Ridge(double[][] X, double[] y, double[] beta,
			double[] y_MLE, Matrix Calculator, double Ridge) {
		double[][] M = new double[X[0].length][X[0].length];
		double[][] M2 = new double[X[0].length][X.length];
		// X'X
		Calculator.multAtB(X, X, M);
		// Ridge
		for (int i = 0; i < M.length; i++) {
			M[i][i] += Ridge;
		}
		// (X'X)-1
		Calculator.symmetricInverse(M);
		// (X'X)-1X'
		Calculator.multABt(M, X, M2);
		// (X'X)-1X'Y=beta
		Calculator.multAx(M2, y, beta);
		// X(X'X)-1X'Y=y_MLE
		Calculator.multAx(X, beta, y_MLE);
		return this.squrePredictionError(y, y_MLE);
	}

	/*
	 * LAR Least Angle Regression of n-folds cross validation
	 */
	public void LARS_NfoldsCV(double[][] origX, double[] origY, double[] beta, Matrix Calculator, int nfolds, boolean printprogress, int eval_amount) {

		double[] y = this.normlization(origY);
		double[][] X = this.normlizationColumn(origX);
		
		int[] cv_order = this.getRandomOrder(0, origY.length-1, new Sfmt(0));
		double[][] Lambda = new double[eval_amount][3];
		
		double[] residual_correlation = new double[X[0].length];
		Calculator.multxtA(y, X, residual_correlation);
		double max = Math.abs(residual_correlation[this.getMaxAbsValueIndex(residual_correlation)]);
		for (int i = 0; i < Lambda.length; i++) {
			Lambda[i][0] = max * (Lambda.length - i) / Lambda.length;
		}
		
		double active = 0;		
		int length = 0;
		int start = 0;
		int end = 0;
		
		for (int k = 0; k < nfolds + 1; k++) {

			/* Make Cross Validation Samples */
			if(nfolds > k){
				length = (origY.length - end)/(nfolds - k);
				start = end;
				end += length;
			} else {
				active = (int) Math.round((-1.0) * active/nfolds);
				length = 0;
				start = -1;
				end = -1;
			}
			
			/* Initialize */
			Calculator.setvalue(beta, 0);			
			for (int i = 0; i < Lambda.length; i++) {
				Lambda[i][1] = Lambda[i][2] = 0;
			}
			
			/* CV Sample */
			double[] y_cv = new double[length];
			double[][] X_cv =  new double[length][X[0].length];
					
			/* Learning Sample */
			y = new double[origY.length - length];
			X =  new double[origX.length - length][X[0].length];
			
			for (int c = 0; c < origY.length; c++) {
				if(c < start){
					y[c] = origY[cv_order[c]];
					Calculator.copy(X[c], origX[cv_order[c]]);
				} else if(c >=start && c < end){
					y_cv[c - start] = origY[cv_order[c]];
					Calculator.copy(X_cv[c - start], origX[cv_order[c]]);
				} else if(c >= end){
					y[c - length] = origY[cv_order[c]];
					Calculator.copy(X[c - length], origX[cv_order[c]]);
				}
			}
			
			/* Normalization */
			y = this.normlization(y);
			X = this.normlizationColumn(X);
			if(nfolds > k){
				y_cv = this.normlization(y_cv);
				X_cv = this.normlizationColumn(X_cv);
			}
			
			/* Main LARS */
			ArrayList<Integer> Active = new ArrayList<Integer>();
			double[] residual = new double[y.length];
			double[] next_additional_beta_direction = new double[X[0].length];
			double[] next_additional_y = new double[y.length];
			double[] next_additional_cor = new double[X[0].length];
			double[] rate = new double[X[0].length];
			
			Calculator.multxtA(y, X, residual_correlation);
			Active.add(this.getMaxAbsValueIndex(residual_correlation));
			if(printprogress) System.out.println("\nInitial Residual: " + Calculator.dotProduct(y, y) + "\n");

			int counter = 0;
			while (counter < Lambda.length) {
				
				/* Calculate Next OLS Estimation Value for y */
				double[][] aaM = new double[Active.size()][Active.size()];
				double[][] aM = Calculator.getPartOfColVectors(X, Active);
				Calculator.multAtB(aM, aM, aaM);
				Calculator.symmetricInverse(aaM);

				Calculator.setvalue(next_additional_beta_direction, 0);
				double[] next_directionActive = new double[Active.size()];
				double[] memo = new double[Active.size()];
				for (int i = 0; i < Active.size(); i++) {
					memo[i] = residual_correlation[Active.get(i)];
				}
				Calculator.multAx(aaM, memo, next_directionActive);
				for (int i = 0; i < Active.size(); i++) {
					next_additional_beta_direction[Active.get(i)] = next_directionActive[i];
				}

				/* Get Starting Point (maxcol = Lambda) */
				Calculator.copy(residual, y);
				Calculator.multSubAx(X, beta, residual);
				Calculator.multxtA(residual, X, residual_correlation);
				double maxcol = Math.abs(residual_correlation[this.getMaxAbsValueIndex(residual_correlation)]);
				if(printprogress) System.out.println("Maximum Correlation: " + maxcol);
				
				/* Get Direction */
				Calculator.multAx(X, next_additional_beta_direction, next_additional_y);
				Calculator.multxtA(next_additional_y, X, next_additional_cor);

				for (int i = 0; i < rate.length; i++) {
					if (Active.contains(i)) {
						rate[i] = 1.0;
						continue;
					}
					/* maxcol = Math.abs(next_additional_cor[Active.get(0)]) */
					boolean plus = true;
					if (residual_correlation[i] > 0 && next_additional_cor[i] > 0){
						if(Math.abs(residual_correlation[i]) > Math.abs(next_additional_cor[i])) plus = false;
						else plus = true;		
					} else if (residual_correlation[i] < 0 && next_additional_cor[i] < 0){
						if(Math.abs(residual_correlation[i]) > Math.abs(next_additional_cor[i])) plus = true;
						else plus = false;
					} else if (residual_correlation[i] > 0 && next_additional_cor[i] < 0) {
						plus = false;
					} else if (residual_correlation[i] < 0 && next_additional_cor[i] > 0) {
						plus = true;
					}					
					if(plus) rate[i] = (maxcol + residual_correlation[i]) / (maxcol + next_additional_cor[i]);
					else rate[i] = (maxcol - residual_correlation[i]) / (maxcol - next_additional_cor[i]);
				}
				
				/* Delete Element */
				for (int i = 0; i < rate.length; i++) {
					if (!Active.contains(i) || this.getSign(beta[i]) == this.getSign(next_additional_beta_direction[i]) || this.getSign(beta[i]) == 0){
						continue;
					}
					rate[i] = (-1) * beta[i] / next_additional_beta_direction[i];
				}

				/* Update Beta */
				int next = this.getMinValueIndex(rate);
				boolean updateLambda=false;
				while(maxcol < Lambda[counter][0]){
					Lambda[counter][1] = this.squrePredictionError(X_cv, y_cv, beta);
					Lambda[counter][2] = 0;
					counter++;
				}
				double s = (maxcol - Lambda[counter][0]) / maxcol;
				if(s < rate[next]){
					rate[next] = s;
					updateLambda=true;
				} else if(Active.contains(next)){
					Active.remove(Active.indexOf(next));
				} else {
					Active.add(next);
				}
				Collections.sort(Active);
				Calculator.rescale(next_additional_beta_direction, rate[next]);
				Calculator.add(beta, next_additional_beta_direction);

				/* Calculate Residual */
				Calculator.copy(residual, y);
				Calculator.multSubAx(X, beta, residual);
				if(printprogress) System.out.println("Residual: " + Calculator.dotProduct(residual, residual));
				if(updateLambda) {
					Lambda[counter][1] = this.squrePredictionError(X_cv, y_cv, beta);
					Lambda[counter][2] = this.getActiveCount(beta);
					if(active == this.getActiveCount(beta)) return;
					if(printprogress) System.out.println("CV error: " + Lambda[counter][1] + "\n");
					counter++;
				}
				
				/* Update Correlation */
				Calculator.multxtA(residual, X, residual_correlation);
			}
			
			int index = -1;
			double min = Double.MAX_VALUE;
			for (int i = 0; i < Lambda.length; i++) {
				if(min > Lambda[i][1]) {
					index = i;
					min = Lambda[i][1];
				}
			}
			active -= Lambda[index][2];
		}
	}
	
	/*
	 * LAR Least Angle Regression of n-folds cross validation
	 */
	public void LARS_BIC(double[][] origX, double[] origY, double[] beta, Matrix Calculator, 
			boolean printprogress, int eval_amount) {
		
		/* Main LARS */
		int index=-1;
		for (int loop = 0; loop < 2; loop++) {
			double[] y = this.normlization(origY);
			double[][] X = this.normlizationColumn(origX);
			
			double[][] Lambda = new double[eval_amount][3];
			
			double[] residual_correlation = new double[X[0].length];
			Calculator.multxtA(y, X, residual_correlation);
			double max = Math.abs(residual_correlation[this.getMaxAbsValueIndex(residual_correlation)]);
			for (int i = 0; i < Lambda.length; i++) {
				Lambda[i][0] = max * (Lambda.length - i) / Lambda.length * 1.1;
			}
			
			/* Initialize */
			Calculator.setvalue(beta, 0);			
			for (int i = 0; i < Lambda.length; i++) {
				Lambda[i][1] = Lambda[i][2] = 0;
			}
			
			int counter = 0;
			ArrayList<Integer> Active = new ArrayList<Integer>();
			double[] residual = new double[y.length];
			double[] next_additional_beta_direction = new double[X[0].length];
			double[] next_additional_y = new double[y.length];
			double[] next_additional_cor = new double[X[0].length];
			double[] rate = new double[X[0].length];
			
			Calculator.multxtA(y, X, residual_correlation);
			Active.add(this.getMaxAbsValueIndex(residual_correlation));
			if(printprogress && loop==0) System.out.println("\nInitial Residual: " + Calculator.dotProduct(y, y) + "\n");
	
			while (counter < Lambda.length) {
					
				/* Calculate Next OLS Estimation Value for y */
				double[][] aaM = new double[Active.size()][Active.size()];
				double[][] aM = Calculator.getPartOfColVectors(X, Active);
				Calculator.multAtB(aM, aM, aaM);
				Calculator.symmetricInverse(aaM);
	
				Calculator.setvalue(next_additional_beta_direction, 0);
				double[] next_directionActive = new double[Active.size()];
				double[] memo = new double[Active.size()];
				for (int i = 0; i < Active.size(); i++) {
					memo[i] = residual_correlation[Active.get(i)];
				}
				Calculator.multAx(aaM, memo, next_directionActive);
				for (int i = 0; i < Active.size(); i++) {
					next_additional_beta_direction[Active.get(i)] = next_directionActive[i];
				}
	
				/* Get Starting Point (maxcol = Lambda) */
				Calculator.copy(residual, y);
				Calculator.multSubAx(X, beta, residual);
				Calculator.multxtA(residual, X, residual_correlation);
				double maxcol = Math.abs(residual_correlation[this.getMaxAbsValueIndex(residual_correlation)]);
				if(printprogress && loop==0) System.out.println("Maximum Correlation: " + maxcol);
					
				/* Get Direction */
				Calculator.multAx(X, next_additional_beta_direction, next_additional_y);
				Calculator.multxtA(next_additional_y, X, next_additional_cor);
	
				for (int i = 0; i < rate.length; i++) {
					if (Active.contains(i)) {
						rate[i] = 1.0;
						continue;
					}
					/* maxcol = Math.abs(next_additional_cor[Active.get(0)]) */
					boolean plus = true;
					if (residual_correlation[i] > 0 && next_additional_cor[i] > 0){
						if(Math.abs(residual_correlation[i]) > Math.abs(next_additional_cor[i])) plus = false;
						else plus = true;		
					} else if (residual_correlation[i] < 0 && next_additional_cor[i] < 0){
						if(Math.abs(residual_correlation[i]) > Math.abs(next_additional_cor[i])) plus = true;
						else plus = false;
					} else if (residual_correlation[i] > 0 && next_additional_cor[i] < 0) {
						plus = false;
					} else if (residual_correlation[i] < 0 && next_additional_cor[i] > 0) {
						plus = true;
					}					
					if(plus) rate[i] = (maxcol + residual_correlation[i]) / (maxcol + next_additional_cor[i]);
					else rate[i] = (maxcol - residual_correlation[i]) / (maxcol - next_additional_cor[i]);
				}
				
				/* Delete Element */
				for (int i = 0; i < rate.length; i++) {
					if (!Active.contains(i) || this.getSign(beta[i]) == this.getSign(next_additional_beta_direction[i]) || this.getSign(beta[i]) == 0){
						continue;
					}
					rate[i] = (-1) * beta[i] / next_additional_beta_direction[i];
				}
	
				/* Update Beta */
				int next = this.getMinValueIndex(rate);
				boolean updateLambda=false;
				Calculator.copy(residual, y);
				Calculator.multSubAx(X, beta, residual);
				while(maxcol < Lambda[counter][0]){
					Lambda[counter][1] = this.calculateBIC(this.calculateLogLikelihoodOfMLE(y, Calculator), y.length, beta);
					Lambda[counter][2] = 0;
					counter++;
				}
				double s = (maxcol - Lambda[counter][0]) / maxcol;
				if(s < rate[next]){
					rate[next] = s;
					updateLambda=true;
				} else if(Active.contains(next)){
					Active.remove(Active.indexOf(next));
				} else {
					Active.add(next);
				}
				Collections.sort(Active);
				Calculator.rescale(next_additional_beta_direction, rate[next]);
				Calculator.add(beta, next_additional_beta_direction);
	
				/* Calculate Residual */
				Calculator.copy(residual, y);
				Calculator.multSubAx(X, beta, residual);
				if(printprogress && loop==0) System.out.println("Residual: " + Calculator.dotProduct(residual, residual));
				if(updateLambda) {
					Lambda[counter][1] = this.calculateBIC(this.calculateLogLikelihoodOfMLE(residual, Calculator), y.length, beta);
					Lambda[counter][2] = this.getActiveCount(beta);
					if(counter==index) return;
					if(printprogress) System.out.println("BIC: " + Lambda[counter][1] + "\n");
					counter++;
				}
					
				/* Update Correlation */
				Calculator.multxtA(residual, X, residual_correlation);
			}
			double min = Double.MAX_VALUE;
			for (int i = 0; i < Lambda.length; i++) {
				if(min > Lambda[i][1]) {
					index = i;
					min = Lambda[i][1];
				}
			}
		}
	}

	public double calculateBIC(double logLike, int dataLength, double[] beta){
		return (-2.0 * logLike) + this.getActiveCount(beta) * Math.log(dataLength);
	}

	/*
	 * Output n-folds Cross Validation error
	 */
	public double OLS_withCV(double[][] origX, double[] origY, int nfolds, Sfmt Sf, Matrix Calculator) {

		int[] cv_order = this.getRandomOrder(0, origY.length-1, Sf);

		int start = 0;
		int end = 0;
		int length = 0;
		double error = 0;
		for (int k = 0; k < nfolds; k++) {
			
			/* Make Cross Validation Samples */
			length = (origY.length - end)/(nfolds - k);
			start = end;
			end += length;
			
			/* CV Sample */
			double[] y_cv = new double[length];
			double[][] X_cv =  new double[length][origX[0].length];
					
			/* Learning Sample */
			double[] y = new double[origY.length - length];
			double[][] X =  new double[origX.length - length][origX[0].length];
			
			for (int c = 0; c < origY.length; c++) {
				if(c < start){
					y[c] = origY[cv_order[c]];
					Calculator.copy(X[c], origX[cv_order[c]]);
				} else if(c >=start && c < end){
					y_cv[c - start] = origY[cv_order[c]];
					Calculator.copy(X_cv[c - start], origX[cv_order[c]]);
				} else if(c >= end){
					y[c - length] = origY[cv_order[c]];
					Calculator.copy(X[c - length], origX[cv_order[c]]);
				}
			}
			
			/* Normalization */
			/*y = this.normlization(y);
			X = this.normlizationColumn(X);
			y_cv = this.normlization(y_cv);
			X_cv = this.normlizationColumn(X_cv);
			*/
			
			double[] beta = new double[X[0].length];
			error += this.Ridge(X, y, beta, new double[y.length], Calculator, 1.0e-5);
		}
		return error / nfolds;
	}

	public void sort(double[] array, int min, int max) {
		if (array.length == 1 || array.length == 0 || min > max || min < 0
				|| max > array.length - 1 || array[0] == Double.NaN) {
			// induce error
			array[Integer.MAX_VALUE] = 0;
		}

		double[] tmp = null;
		// Bubble sort if the number of elements is enough small
		if ((max - min) < 5) {
			bubbleSort(array, min, max);
			return;
		}

		// Get heap space
		if (tmp == null)
			tmp = new double[array.length];

		// Divide space into large one and small one
		double m = getMean(array, min, max);
		int largeCount = 0;
		int smallCount = 0;
		for (int i = min; i <= max; i++) {
			if (array[i] > m) {
				tmp[min + largeCount] = array[i];
				largeCount++;
			} else {
				tmp[max - smallCount] = array[i];
				smallCount++;
			} 
		}
		for (int i = min; i <= max; i++) {
			array[i] = tmp[i];
		}
		
		if(largeCount==0){
			largeCount = (max - min + 1)/2 ;
			smallCount = (max - min + 1) - largeCount;
		}
		if(smallCount==0){
			largeCount = (max - min + 1)/2 ;
			smallCount = (max - min + 1) - largeCount;
		}

		// Sorting for each space
		sort(array, min, min + largeCount - 1);
		sort(array, max - smallCount + 1, max);
	}
	
	public void sort(double[] array, double order[], int min, int max) {
		if(array.length == 1 || array.length == 0){
			if(array.length == 0) System.err.println("Caution, Array.length=0!!");
			return;
		}
		if (min > max || min < 0 || max > array.length - 1 || array[0] == Double.NaN) {
			System.err.println("Sorting ERROR!!");
			array[Integer.MAX_VALUE] = 0;
		}
		
		for (int i = 0; i < array.length; i++) {
			if(array[i] > 1.0e20) {
				array[i] = 1.0e20;
				System.out.println("simpleMath - Too large value is replaced");
			}
		}

		double[] tmp = null;
		double[] tmp_ = null;
		// Bubble sort if the number of elements is enough small
		if ((max - min) < 5) {
			bubbleSort(array, order, min, max);
			return;
		}
		
		// Get heap space
		if (tmp == null){
			tmp = new double[array.length];
			tmp_ = new double[order.length];
		}

		// Divide space into large one and small one
		double m = getMean(array, min, max);
		int largeCount = 0;
		int smallCount = 0;
		for (int i = min; i <= max; i++) {
			if (array[i] > m) {
				tmp[min + largeCount] = array[i];
				tmp_[min + largeCount] = order[i];
				largeCount++;
			} else {
				tmp[max - smallCount] = array[i];
				tmp_[max - smallCount] = order[i];
				smallCount++;
			}
		}
		if(largeCount==0 || smallCount==0) return;
		for (int i = min; i <= max; i++) {
			array[i] = tmp[i];
			order[i] = tmp_[i];
		}

		// Sorting for each space
		sort(array, order, min, min + largeCount - 1);
		sort(array, order, max - smallCount + 1, max);
	}


	public void bubbleSort(double[] array, int min, int max) {

		for (int i = min; i < max; i++) {
			for (int j = min; j < max - (i - min); j++) {
				if (array[j] < array[j + 1]) {
					double tmp = array[j];
					array[j] = array[j + 1];
					array[j + 1] = tmp;
				}
			}
		}
	}
	
	public void bubbleSort(double[] array, double[] order, int min, int max) {

		for (int i = min; i < max; i++) {
			for (int j = min; j < max - (i - min); j++) {
				if (array[j] < array[j + 1]) {
					double tmp = array[j];
					double tmp_ = order[j];
					array[j] = array[j + 1];
					order[j] = order[j+1];
					array[j + 1] = tmp;
					order[j + 1] = tmp_;
				}
			}
		}
	}

	public ArrayList<Integer> getDownStreamNodes(ArrayList<Integer> array,
			double[][] A, int parent) {
		for (int i = 0; i < A.length; i++) {
			if (A[i][parent] != 0 && array.indexOf(i) == -1) {
				array.add(i);
				this.getDownStreamNodes(array, A, i);
			}
		}
		return array;
	}
	
	public int[] getIthOrder(int length){
		int[] re = new int[length];
		for (int i = 0; i < re.length; i++) {
			re[i] = i;
		}
		return re;
	}
	
	public double[] getIthOrder(double length){
		double[] re = new double[(int)length];
		for (int i = 0; i < re.length; i++) {
			re[i] = i;
		}
		return re;
	}
	
	public int getSign(double v) {
		return (int) (v / (Math.abs(v)));
	}

	public ArrayList<Integer> getRandomExtract(int max_size, int extract_size) {
		Sfmt sf = new Sfmt((int) (max_size * 100000 * Math.random()));
		return getRandomExtract(max_size, extract_size, sf);
	}

	public ArrayList<Integer> getRandomExtract(int max_size, int extract_size, Sfmt sf) {
		ArrayList<Integer> result = new ArrayList<Integer>();
		ArrayList<Integer> test = new ArrayList<Integer>();
		for (int i = 0; i < max_size; i++)
			test.add(i);
		if (max_size <= extract_size)
			return test;
		for (int i = 0; i < extract_size; i++) {
			int r = getRandomInt(test.size(), sf);
			result.add(test.get(r));
			test.remove(r);
		}
		Collections.sort(result);
		return result;
	}

	public ArrayList<Integer> getSwap(int size) {
		Sfmt sf = new Sfmt((int) (size * 100000 * Math.random()));
		return getSwap(size, sf);
	}

	public ArrayList<Integer> getSwap(int size, Sfmt sf) {
		ArrayList<Integer> result = new ArrayList<Integer>();
		ArrayList<Integer> test = new ArrayList<Integer>();
		for (int i = 0; i < size; i++)
			test.add(i);
		while (test.size() > 0) {
			int r = getRandomInt(test.size(), sf);
			result.add(test.get(r));
			test.remove(r);
		}
		return result;
	}

	public int getRandomInt(int size) {
		Sfmt sf = new Sfmt((int) (size * 100000 * Math.random()));
		return getRandomInt(size, sf);
	}

	public int getRandomInt(int size, Sfmt sf) {
		int x = size;
		while (x == size)
			x = (int) (sf.NextUnif() * size);
		return x;
	}

	/**
	 * x = [node1, 2, , 3,..., N][observed 1, 2, ..., e.g., 1000]
	 */
	public double[][] getCovarianceRow(double[][] x) {
		double[][] cov = new double[x.length][x.length];
		for (int j = 0; j < x[0].length; j++) {
			for (int i = 0; i < x.length; i++) {
				for (int i2 = i; i2 < x.length; i2++) {
					double val = x[i][j] * x[i2][j];
					cov[i][i2] += val;
					if(i!=i2) cov[i2][i] += val;
				}
			}
		}
		return cov;
	}
	
	/**
	 * x = [observed 1, 2, ..., e.g., 1000][node1, 2, , 3,..., N]
	 */
	public double[][] getCovarianceCol(double[][] x) {
		double[][] cov = new double[x[0].length][x[0].length];
		for (int i = 0; i < x.length; i++) {
			for (int j = 0; j < x[0].length; j++) {
				for (int j2 = j; j2 < x[0].length; j2++) {
					double val = x[i][j] * x[i][j2];
					cov[j][j2] += val;
					if(j!=j2) cov[j2][j] += val;
				}
			}
		}
		return cov;
	}
	
	public double getLagPartialCorrelation(double[][] lagCor, double[][] Cor, int target_i, int target_j, ArrayList<Integer> Given){
		if(Given.size()==0 || Given == null){
			return lagCor[target_i][target_j];
		} else if(Given.size()==1){
			return (lagCor[target_i][target_j] - (lagCor[target_i][Given.get(0)] * Cor[target_j][Given.get(0)])) 
					/ (Math.sqrt(1 - (lagCor[target_i][Given.get(0)] * lagCor[target_i][Given.get(0)])) 
							* Math.sqrt(1 - (Cor[target_j][Given.get(0)] * Cor[target_j][Given.get(0)])));
		} else {
			ArrayList<Integer> newGiven = new ArrayList<Integer>(Given);
			int rem = Given.get(0);
			newGiven.remove(0);
			double r1 = this.getLagPartialCorrelation(lagCor, Cor, target_i, rem, newGiven);
			double r2 = this.getPartialCorrelation(Cor, rem, target_j, newGiven);
			return (this.getLagPartialCorrelation(lagCor, Cor, target_i, target_j, newGiven) - ( r1 * r2)) 
					/ (Math.sqrt(1 - (r1 * r1)) * Math.sqrt(1 - (r2 * r2)));
		}
	}
	
	public double getPartialCorrelation(double[][] Cor, int target_i, int target_j, ArrayList<Integer> Given){
		
		if(Given.size()==0 || Given == null){
			return Cor[target_i][target_j];
		} else if(Given.size()==1){
			return (Cor[target_i][target_j] - (Cor[target_i][Given.get(0)] * Cor[target_j][Given.get(0)])) 
					/ (Math.sqrt(1 - (Cor[target_i][Given.get(0)] * Cor[target_i][Given.get(0)])) 
							* Math.sqrt(1 - (Cor[target_j][Given.get(0)] * Cor[target_j][Given.get(0)])));
		} else {
			ArrayList<Integer> newGiven = new ArrayList<Integer>(Given);
			int rem = Given.get(0);
			newGiven.remove(0);
			double r1 = this.getPartialCorrelation(Cor, target_i, rem, newGiven);
			double r2 = this.getPartialCorrelation(Cor, rem, target_j, newGiven);
			return (this.getPartialCorrelation(Cor, target_i, target_j, newGiven) - ( r1 * r2)) 
					/ (Math.sqrt(1 - (r1 * r1)) * Math.sqrt(1 - (r2 * r2)));
		}
	}
	
	/**
	 * Not Yet
	 */
	public int GrubbsTest(double[] list, double alpha){
		double mu = this.getMean(list);
		double var = this.getVariance(list, mu) * list.length / (list.length - 1);
		double[] tau = new double[list.length];
		for (int i = 0; i < list.length; i++) {
			tau[i] = (list[i] - mu) / Math.sqrt(var);
		}
		if(this.getMaxValue(tau) < alpha) return this.getMaxValueIndex(tau);
		else return -1;		
	}

	public double getVariance(ArrayList<Double> currentLikelihoodList, ArrayList<Double> previousLikelihoodList) {
		double mean = 0;
		double mean2 = 0;
		int length = currentLikelihoodList.size() + previousLikelihoodList.size();
		for (int i = 0; i < currentLikelihoodList.size(); i++) {
			mean += currentLikelihoodList.get(i);
			mean2 += currentLikelihoodList.get(i) * currentLikelihoodList.get(i);
		}
		for (int i = 0; i < previousLikelihoodList.size(); i++) {
			mean += previousLikelihoodList.get(i);
			mean2 += previousLikelihoodList.get(i) * previousLikelihoodList.get(i);
		}
		return (mean2 / length) - ((mean / length) * (mean / length));
	}

	public void getIthOrder(ArrayList<Integer> order, int n) {
		// TODO Auto-generated method stub
		for (int i = 0; i < n; i++) {
			order.add(i);
		}
	}

	public boolean match(double[][] timeSeriesData, double[][] currentNet) {
		boolean match = true;
		for (int i = 0; i < currentNet.length; i++) {
			for (int j = 0; j < currentNet.length; j++) {
				if(timeSeriesData[i][j] != currentNet[i][j]) match = false;
			}
		}
		return match;
	}
	
	public double getNormalDense(double x, double mean, double var){
		return (-1)/(Math.sqrt(2*Math.PI*var)*Math.exp(-0.5*(x - mean)*(x - mean)/var));
	}

	public int countContain(ArrayList<Integer> up, int target) {
		int count = 0;
		for (int i = 0; i < up.size(); i++) {
			if(up.get(i)==target) count++;
		}
		return count;
	}

	public double[][] getVarianceMatrix(double[] v, int dim) {
		int N = (int)Math.round((double)v.length/dim);
		double[] mean = new double[dim];
		double[] mean2 = new double[dim];
		double[][] cor = new double[dim][dim];
		for (int n = 0; n < N; n++) {
			for (int i = 0; i < dim; i++) {
				double val = v[dim * n + i];
				mean[i] += val;
				mean2[i] += val * val;
			}
		}
		for (int i = 0; i < dim; i++) {
			mean[i] /= N;
			mean2[i] /= N;
		}
		for (int i = 0; i < dim; i++) {
			cor[i][i] = mean2[i] - mean[i] * mean[i];
			for (int n = 0; n < N; n++) {
				for (int j = i + 1; j < dim; j++) {
					cor[i][j] += (v[dim * n + i] - mean[i])
							* (v[dim * n + j] - mean[j]);
				}
			}
			for (int j = i + 1; j < dim; j++) {
				cor[i][j] /= N;
				cor[j][i] = cor[i][j];
			}
		}
		return cor;
	}
	
	public double[][] getCor(double[] v, int dim) {
		double[][] cor = this.getVarianceMatrix(v, dim);
		double[] dev = new double[dim];
		for (int i = 0; i < dim; i++) {
			dev[i] = Math.sqrt(cor[i][i]);
		}
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				cor[i][j] /= dev[i] * dev[j];
			}
		}
		return cor;
	}
	
	public ArrayList<Integer> convertIntToBinary(double v) {
		ArrayList<Integer> binary= new ArrayList<Integer>();
		while(v > 0){
			if(v % 2 == 1) {
				v = (v - 1) / 2;
				binary.add(1);
			} else{
				v = v / 2;
				binary.add(0);
			}
		}
		return binary;
	}

	public double[][] getCor(double[] v, int dim, ArrayList<Integer> updatingParameter, Matrix MP) {
		int valid_thetaNum = (int)MP.sumofVector(updatingParameter.toArray(new Integer[updatingParameter.size()]));
		
		double[] valid_v = new double[v.length / dim * valid_thetaNum];
		int count = 0;
		for (int i = 0; i < v.length; i++) {
			if(updatingParameter.get(i % dim) == 1) {
				valid_v[count] = v[i];
				count++;
			}
		}
		
		double[][] cor = this.getVarianceMatrix(valid_v, valid_thetaNum);
		double[] dev = new double[valid_thetaNum];
		for (int i = 0; i < valid_thetaNum; i++) {
			dev[i] = Math.sqrt(cor[i][i]);
		}
		for (int i = 0; i < valid_thetaNum; i++) {
			for (int j = 0; j < valid_thetaNum; j++) {
				cor[i][j] /= dev[i] * dev[j];
			}
		}
		return cor;
	}
}







