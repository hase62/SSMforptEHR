package Hasegawa.TimeSeries;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.matrix.Matrix;

public class Storage {
	
	public double Criterion = Double.MAX_VALUE;
	public double logLikelihood = (-1)*Double.MAX_VALUE;
	public double predictionError[];
	public double predictionAbility[];
	
	/* Simulation Expression Profiles*/
	protected double[][] x_p;
	protected double[][] x_f;
	protected double[][] x_s;
	protected double[][] x_s0;
	protected double[][] y_p;
	
	/* Expectations */
	protected double[][] Txx_obs;
	protected double[][] Txx;
	protected double[][] Tyx;
	protected double[][] Txx_m;
	protected double[][] Txx_mm;
	protected double[][] Txz_m;
	protected double[][] Txz_mm;
	protected double[] sum_x;
	protected double[] sum_xm;
	protected double[] sum_zm;
	
	protected double[][] Skew_mean_t;
	protected double[][] Skew_sq_mean_t;
	protected double[][] Delta_skew_mean_t;

	protected double[][] Kurtosis_mean_t;
	
	protected double[][] Kurtosis_mean_Rinv_t;
	protected double[][][] Htrans_Kurtosis_mean_Rinv_H_t;
	protected double[][][] Kurtosis_mean_Rinv_H_t;
	
	protected double[] Kurtosis_skew_sq_mean_sum;
	protected double[] Dif_y_x_skew_mean_sum;
	protected double[] Dif_y_x_sq_sum;
	protected double[][] Skew_mean_x_sum;
	protected double[] Kurtosis_H_Sigma_Htrans_sum;
	
	/*
	 * Initialize
	 */
	public void Initialize(int sysDim, TimeSeriesDataArray tsda, Setting set){
		this.predictionError = new double[sysDim];
		this.predictionAbility = new double[sysDim];
		
		this.x_p = new double[tsda.allTime][sysDim];
		this.x_f = new double[tsda.allTime][sysDim];
		this.x_s = new double[tsda.allTime][sysDim];
		this.x_s0 = new double[tsda.repSize][sysDim];
		this.y_p = new double[tsda.allTime][tsda.elementNum];
		
		this.Txx_obs = new double[sysDim][sysDim];
		this.Txx = new double[sysDim][sysDim];
		this.Tyx = new double[tsda.elementNum][sysDim];
		this.Txx_m = new double[sysDim][sysDim];
		this.Txx_mm = new double[sysDim][sysDim];
		
		//this.WeightA = new double[sysDim][sysDim];
		
		if(set.Drug){
			this.Txz_m = new double[sysDim][tsda.drugMulRepSize[0].length];
			this.Txz_mm = new double[sysDim][tsda.drugMulRepSize[0].length];
			//this.WeightG = new double[sysDim][fTimeSeries.drugMulRepSize[0].length];
			if(set.Input){
				this.sum_zm = new double[tsda.drugMulRepSize[0].length];
			}
		}
		if(set.Input){
			this.sum_x = new double[sysDim];
			this.sum_xm = new double[sysDim];
		}
		
		this.Skew_mean_t = new double[tsda.allTime][tsda.elementNum];
		this.Skew_sq_mean_t = new double[tsda.allTime][tsda.elementNum];
		this.Delta_skew_mean_t = new double[tsda.allTime][tsda.elementNum];
		
		this.Kurtosis_mean_t = new double[tsda.allTime][tsda.elementNum];

		this.Kurtosis_mean_Rinv_t = new double[tsda.allTime][tsda.elementNum];
		this.Htrans_Kurtosis_mean_Rinv_H_t = new double[tsda.allTime][sysDim][sysDim];
		this.Kurtosis_mean_Rinv_H_t = new double[tsda.allTime][tsda.elementNum][sysDim];
		
		this.Kurtosis_skew_sq_mean_sum = new double[tsda.elementNum];
		this.Dif_y_x_skew_mean_sum = new double[tsda.elementNum];
		this.Dif_y_x_sq_sum = new double[tsda.elementNum];
		this.Skew_mean_x_sum = new double[tsda.elementNum][sysDim];
		this.Kurtosis_H_Sigma_Htrans_sum = new double[tsda.elementNum];
	}
	
	/*
	 * Set all Profiles
	 */
	public void setProfiles(double[][] x_s0_, double[][] x_p_, double[][] x_f_, double[][] x_s_, double[][] y_p_){
		this.setx_s0(x_s0_);
		this.setx_p(x_p_);
		this.setx_f(x_f_);
		this.setx_s(x_s_);
		if(y_p_ != null) this.sety_p(y_p_);
	}
	
	/*
	 * Get all Profiles
	 */
	public void getProfiles(double[][] x_s0_, double[][] x_p_, double[][] x_f_, double[][] x_s_, double[][] y_p_, Matrix MP){
		MP.copy(x_s0_, this.getx_s0());
		MP.copy(x_p_, this.getx_p());
		MP.copy(x_f_, this.getx_f());
		MP.copy(x_s_, this.getx_s());
		if(y_p_ != null) MP.copy(y_p_, this.gety_p());
	}
	
	/*
	 * x_p
	 */
	private void setx_p(double[][] M) {
		this.x_p = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.x_p[i][j] = M[i][j];
			}
    	}	
	}
	
    public double[][] getx_p(){
    	assert(this.x_p!= null);
    	return this.x_p;
    }
	
    /*
     * x_f
     */
    private void setx_f(double[][] M) {
		this.x_f = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.x_f[i][j] = M[i][j];
			}
    	}
	}
	
    public double[][] getx_f(){
    	assert(this.x_f!= null);
    	return this.x_f;
    }
	
    /*
     * x_s
     */
    private void setx_s(double[][] M) {
		this.x_s = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.x_s[i][j] = M[i][j];
			}
    	}	
	}
	
    public double[][] getx_s(){
    	assert(this.x_s!= null);
    	return this.x_s;
    }
	
    /*
     * x_s0
     */
    private void setx_s0(double[][] M) {
		this.x_s0 = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.x_s0[i][j] = M[i][j];
			}
    	}	
	}
	
    public double[][] getx_s0(){
    	assert(this.x_s0!= null);
    	return this.x_s0;
    }
    
	/*
	 * y_p
	 */
	private void sety_p(double[][] M) {
		this.y_p = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.y_p[i][j] = M[i][j];
			}
    	}	
	}
	
    public double[][] gety_p(){
    	assert(this.y_p!= null);
    	return this.y_p;
    }
    
    /* Set Expectation Value*/
	public void setExpectations(double[][] Txx_obs, double[][] Txx, double[][] Tyx, double[][] Txx_m, 
				double[][] Txx_mm, double[][] Txz_mm, double[][] Txz_m, double[] sum_x, 
				double[] sum_xm, double[] sum_zm, double[][] Skew_mean_t, double[][] Skew_sq_mean_t, 
				double[][] Delta_skew_mean_t, double[][] Kurtosis_mean_t, 
				double[][] Kurtosis_mean_Rinv_t, double[][][] Htrans_Kurtosis_mean_Rinv_H_t, 
				double[][][] Kurtosis_mean_Rinv_H_t, 
				double[] Kurtosis_skew_sq_mean_sum, 
				double[] Dif_y_x_skew_mean_sum, 
				double[] Dif_y_x_sq_sum, 
				double[][] Skew_mean_x_sum, 
				double[] Kurtosis_H_Sigma_Htrans_sum, 
				Matrix Calculator){
		if(Txx_obs!=null) Calculator.copy(this.Txx_obs, Txx_obs);
		if(Txx!=null) Calculator.copy(this.Txx, Txx);
		if(Tyx!=null) Calculator.copy(this.Tyx, Tyx);
		if(Txx_m!=null) Calculator.copy(this.Txx_m, Txx_m);
		if(Txx_mm!=null) Calculator.copy(this.Txx_mm, Txx_mm);
		if(Txz_mm!=null) Calculator.copy(this.Txz_mm, Txz_mm);
		if(Txz_m!=null) Calculator.copy(this.Txz_m, Txz_m);
		if(sum_x!=null) Calculator.copy(this.sum_x, sum_x);
		if(sum_xm!=null) Calculator.copy(this.sum_xm, sum_xm);
		if(sum_zm!=null) Calculator.copy(this.sum_zm, sum_zm);
		if(Skew_mean_t!=null) Calculator.copy(this.Skew_mean_t, Skew_mean_t);
		if(Skew_sq_mean_t!=null) Calculator.copy(this.Skew_sq_mean_t, Skew_sq_mean_t);
		if(Delta_skew_mean_t!=null) Calculator.copy(this.Delta_skew_mean_t, Delta_skew_mean_t);
		if(Kurtosis_mean_t!=null) Calculator.copy(this.Kurtosis_mean_t, Kurtosis_mean_t);
		if(Kurtosis_mean_Rinv_t!=null) Calculator.copy(this.Kurtosis_mean_Rinv_t, Kurtosis_mean_Rinv_t);
		if(Htrans_Kurtosis_mean_Rinv_H_t!=null) Calculator.copy(this.Htrans_Kurtosis_mean_Rinv_H_t, Htrans_Kurtosis_mean_Rinv_H_t);
		if(Kurtosis_mean_Rinv_H_t!=null) Calculator.copy(this.Kurtosis_mean_Rinv_H_t, Kurtosis_mean_Rinv_H_t);
		if(Kurtosis_skew_sq_mean_sum!=null) Calculator.copy(this.Kurtosis_skew_sq_mean_sum, Kurtosis_skew_sq_mean_sum);
		if(Dif_y_x_skew_mean_sum!=null) Calculator.copy(this.Dif_y_x_skew_mean_sum, Dif_y_x_skew_mean_sum);
		if(Dif_y_x_sq_sum!=null) Calculator.copy(this.Dif_y_x_sq_sum, Dif_y_x_sq_sum);
		if(Skew_mean_x_sum!=null) Calculator.copy(this.Skew_mean_x_sum, Skew_mean_x_sum);
		if(Kurtosis_H_Sigma_Htrans_sum!=null) Calculator.copy(this.Kurtosis_H_Sigma_Htrans_sum, Kurtosis_H_Sigma_Htrans_sum);
	}
	
	/* Get Expectation Values*/
	public void getExpectations(double[][] Txx_obs, double[][] Txx, double[][] Tyx, double[][] Txx_m, 
							   double[][] Txx_mm, double[][] Txz_mm, double[][] Txz_m, double[] sum_x, 
							   double[] sum_xm, double[] sum_zm, double[][] Skew_mean_t, double[][] Skew_sq_mean_t, 
							   double[][] Delta_skew_mean_t, double[][] Kurtosis_mean_t, 
							   double[][] Kurtosis_mean_Rinv_t, double[][][] Htrans_Kurtosis_mean_Rinv_H_t, 
							   double[][][] Kurtosis_mean_Rinv_H_t,
							   double[] Kurtosis_skew_sq_mean_sum, 
							   double[] Dif_y_x_skew_mean_sum, 
							   double[] Dif_y_x_sq_sum, 
							   double[][] Dif_skeSkew_mean_x_sumw_mean_x_sum, 
							   double[] Kurtosis_H_Sigma_Htrans_sum, 
							   Matrix Calculator){
		if(this.Txx_obs!=null) Calculator.copy(Txx_obs, this.Txx_obs);
		if(this.Txx!=null) Calculator.copy(Txx, this.Txx);
		if(this.Tyx!=null) Calculator.copy(Tyx, this.Tyx);
		if(this.Txx_m!=null) Calculator.copy(Txx_m, this.Txx_m);
		if(this.Txx_mm!=null) Calculator.copy(Txx_mm, this.Txx_mm);
		if(this.Txz_mm!=null) Calculator.copy(Txz_mm, this.Txz_mm);
		if(this.Txz_m!=null) Calculator.copy(Txz_m, this.Txz_m);
		if(this.sum_x!=null) Calculator.copy(sum_x, this.sum_x);
		if(this.sum_xm!=null) Calculator.copy(sum_xm, this.sum_xm);
		if(this.sum_zm!=null) Calculator.copy(sum_zm, this.sum_zm);
		if(this.Skew_mean_t!=null && Skew_mean_t!=null) Calculator.copy(Skew_mean_t, this.Skew_mean_t);
		if(this.Skew_sq_mean_t!=null && Skew_sq_mean_t!=null) Calculator.copy(Skew_sq_mean_t, this.Skew_sq_mean_t);
		if(this.Delta_skew_mean_t!=null && Delta_skew_mean_t!=null) Calculator.copy(Delta_skew_mean_t, this.Delta_skew_mean_t);
		if(this.Kurtosis_mean_t!=null && Kurtosis_mean_t!=null) Calculator.copy(Kurtosis_mean_t, this.Kurtosis_mean_t);
		if(this.Kurtosis_mean_Rinv_t!=null && Kurtosis_mean_Rinv_t!=null) Calculator.copy(Kurtosis_mean_Rinv_t, this.Kurtosis_mean_Rinv_t);
		if(this.Htrans_Kurtosis_mean_Rinv_H_t!=null && Htrans_Kurtosis_mean_Rinv_H_t!=null) Calculator.copy(Htrans_Kurtosis_mean_Rinv_H_t, this.Htrans_Kurtosis_mean_Rinv_H_t);
		if(this.Kurtosis_mean_Rinv_H_t!=null && Kurtosis_mean_Rinv_H_t!=null) Calculator.copy(Kurtosis_mean_Rinv_H_t, this.Kurtosis_mean_Rinv_H_t);
		if(this.Kurtosis_skew_sq_mean_sum!=null && Kurtosis_skew_sq_mean_sum!=null) Calculator.copy(Kurtosis_skew_sq_mean_sum, this.Kurtosis_skew_sq_mean_sum);
		if(this.Dif_y_x_skew_mean_sum!=null && Dif_y_x_skew_mean_sum!=null) Calculator.copy(Dif_y_x_skew_mean_sum, this.Dif_y_x_skew_mean_sum);
		if(this.Dif_y_x_sq_sum!=null && Dif_y_x_sq_sum!=null) Calculator.copy(Dif_y_x_sq_sum, this.Dif_y_x_sq_sum);
		if(this.Skew_mean_x_sum!=null && Skew_mean_x_sum!=null) Calculator.copy(Skew_mean_x_sum, this.Skew_mean_x_sum);
		if(this.Kurtosis_H_Sigma_Htrans_sum!=null && Kurtosis_H_Sigma_Htrans_sum!=null) Calculator.copy(Kurtosis_H_Sigma_Htrans_sum, this.Kurtosis_H_Sigma_Htrans_sum);
	}
}
