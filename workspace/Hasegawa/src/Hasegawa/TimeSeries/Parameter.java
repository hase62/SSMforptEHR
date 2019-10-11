package Hasegawa.TimeSeries;

import Hasegawa.matrix.Matrix;

public class Parameter {
	
	/*
	 * Parameters
	 */
	protected double[][] x0_;
	protected double[][] v0;
	protected double[][] A;
	protected double[][] F;
	protected double[][] H;
	protected double[][] G;
	protected double[] D;
	protected double[] Q;
	protected double[] R;
	protected double[] U;
	protected double[] I;
	protected double[][] Effect1;
	protected double[][] Effect2;
	protected Matrix MP = new Matrix();
	
	protected double[] Delta;
	protected double[] Nu;
	
	protected double[][] Skew_mean_t;
	protected double[][] Skew_sq_mean_t;
	
	protected double[][] Kurtosis_mean_t;
	
	/*
	 * Get all Parameters
	 */
	public double getParameters(double[][] x0_to, double[][] v0_to, double[][] A_to, double[][] F_to, 
			double[][] Ftrans_to, double[][] H_to, double[][] G_to, double[][] RinvH_to, 
			double[][] HtransRinvH_to, double[] D_to, double[] Q_to, double[] Qinv_to, double[] R_to, 
			double[] Rinv_to, double[][] RinvMatrix_to, double[] U_to, double[] I_to, 
			double[] Delta_to, double[] Nu_to, double[][] Skew_mean_t_to, double[][] Skew_sq_mean_t_to, 
			double[][] Kurtosis_mean_t_to, double[][] Delta_skew_mean_t_to, double[][] Kurtosis_mean_Rinv_t_to, 
			double[][][] Htrans_Kurtosis_mean_Rinv_H_t_to, double[][][] Kurtosis_mean_Rinv_H_t_to, 
			double[][] Effect1_to, double[][] Effect2_to){
			
			if(this.x0_!=null) MP.copy(x0_to, this.getx0());
			if(this.v0!=null) MP.copy(v0_to, this.getv0());
		    	if(this.A!=null) MP.copy(A_to, this.getA());
		    	if(this.F!=null) {
				MP.copy(F_to, this.getF());
				MP.copy(Ftrans_to, F_to);
				MP.transpose(Ftrans_to);
		    	}
			if(this.G!=null) MP.copy(G_to, this.getG());
			if(this.D!=null) MP.copy(D_to, this.getD());
			if(this.Q!=null) {
				MP.copy(Q_to, this.getQ());
			    	MP.copy(Qinv_to, Q_to);
			    	MP.DiagInvese(Qinv_to);
			    }
			if(this.R!=null) {
				MP.copy(R_to, this.getR());
				MP.copy(Rinv_to, R_to);
				MP.DiagInvese(Rinv_to);
				MP.setvalue(RinvMatrix_to, 0);
				MP.add(RinvMatrix_to, Rinv_to);
			}
			if(this.H!=null) {
				MP.copy(H_to, this.getH());
				MP.multAB(Rinv_to, H_to, RinvH_to);
				MP.multAtB(H_to, RinvH_to, HtransRinvH_to);
			}
			if(this.U!=null) MP.copy(U_to, this.getU());
			if(this.I!=null) MP.copy(I_to, this.getI());
			if(this.Effect1!=null) MP.copy(Effect1_to, this.getEffect1());
			if(this.Effect2!=null) MP.copy(Effect2_to, this.getEffect2());
			
			if(Delta_to != null) MP.copy(Delta_to, this.getDelta());
			if(Nu_to != null) MP.copy(Nu_to, this.getNu());
			if(Skew_mean_t_to != null) MP.copy(Skew_mean_t_to, this.getSkew_mean_t());
			if(Skew_sq_mean_t_to != null) MP.copy(Skew_sq_mean_t_to, this.getSkew_sq_mean_t());
			if(Kurtosis_mean_t_to != null) MP.copy(Kurtosis_mean_t_to, this.getsetKurtosis_mean_t());
			
			if(this.H!=null && this.Delta!=null && Skew_mean_t_to!=null) {
				for (int t = 0; t < Skew_mean_t_to.length; t++) {
					MP.multDx(Skew_mean_t_to[t], Delta_to, Delta_skew_mean_t_to[t]);					
				}
				for (int t = 0; t < Kurtosis_mean_t_to.length; t++) {
					MP.multDx(Kurtosis_mean_t_to[t], Rinv_to, Kurtosis_mean_Rinv_t_to[t]);
					MP.multAB(Kurtosis_mean_Rinv_t_to[t], H_to, Kurtosis_mean_Rinv_H_t_to[t]);
					MP.multAtB(H_to, Kurtosis_mean_Rinv_H_t_to[t], Htrans_Kurtosis_mean_Rinv_H_t_to[t]);
				}
			}

			return MP.diagLogDeterminant(R_to);
		}
	
	/*
	 * Set all Parameters
	 */
	public void setParameters(double[][] x0_to, double[][] v0_to, double[][] A_to, double[][] F_to, 
			double[][] H_to, double[][] G_to, double[] D_to, double[] Q_to, double[] R_to, double[] U_to, 
			double[] I_to, double[] Delta_to, double[] Nu_to, double[][] Skew_mean_t_to, 
			double[][] Skew_sq_mean_t_to, double[][] Kurtosis_mean_t_to, double[][] Effect1_to, double[][] Effect2_to){
			if(x0_to!=null) this.setx0(x0_to);
			if(v0_to!=null) this.setv0(v0_to);
		    	if(A_to!=null) this.setA(A_to);
		    	if(F_to!=null) this.setF(F_to);
		    	if(H_to!=null) this.setH(H_to);
		    	if(G_to!=null) this.setG(G_to);
		    	if(D_to!=null) this.setD(D_to);
		    	if(Q_to!=null) this.setQ(Q_to);
		    	if(R_to!=null) this.setR(R_to);
		    	if(U_to!=null) this.setU(U_to);
		    	if(I_to!=null) this.setI(I_to);
		    	if(Delta_to!=null) this.setDelta(Delta_to);
		    	if(Nu_to!=null) this.setNu(Nu_to);
		    	if(Skew_mean_t_to!=null) this.setSkew_mean_t(Skew_mean_t_to);
		    	if(Skew_sq_mean_t_to!=null) this.setSkew_sq_mean_t(Skew_sq_mean_t_to);
		    	if(Kurtosis_mean_t_to!=null) this.setKurtosis_mean_t(Kurtosis_mean_t_to);
		    	if(Effect1_to!=null) this.setEffect1(Effect1_to);
		    	if(Effect2_to!=null) this.setEffect2(Effect2_to);
		}
		
	/*
	 * x0
	 */
    public double[][] getx0(){
    	if(this.x0_==null) return null;
    	assert(this.x0_ != null);
    	return this.x0_;
    }
    
    public void setx0(double[][] x){
    	this.x0_ = new double[x.length][x[0].length];
    	for (int i = 0; i < x.length; i++) {
    		for (int j = 0; j < x[0].length; j++) {
    			this.x0_[i][j] = x[i][j];
    		}
		}    	
    }
    
    public void setx0(double x, int rep, int i){
    	this.x0_[rep][i] = x;
	}

	/*
	 * v0
	 */
	public double[][] getv0(){
		if(this.v0==null) return null;
    	assert(this.v0 != null);
    	return this.v0;
    }
    
    public void setv0(double[][] A){
    	this.v0 = new double[A.length][A[0].length];
    	for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				this.v0[i][j] = A[i][j];
			}
    	}
    }
    
    /*
     * A
     */
	public double[][] getA() {
		if(this.A==null) return null;
		assert(this.A !=null);
		return this.A;
	}
        
	public void setA(double[][] M) {
		this.A = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.A[i][j] = M[i][j];
			}
    	}
	}	
    
	/*
	 * F
	 */
    public double[][] getF(){
    	if(this.F==null) return null;
    	assert(this.F != null);
    	return this.F;
    }
    
    public void setF(double[][] M){
    	this.F = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.F[i][j] = M[i][j];
			}
    	}
    }
    
	/*
     * H
     */
    public double[][] getH(){
    	if(this.H==null) return null;
    	assert(this.H!= null);
    	return this.H;
    }
    
    public void setH(double[][] A){
    	this.H = new double[A.length][A[0].length];
    	for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				this.H[i][j] = A[i][j];
			}
    	}
    }
    
    /*
     * G
     */
    public double[][] getG(){
    	if(this.G==null) return null;
    	assert(this.G != null);
    	return this.G;
    }
    
    public void setG(double[][] M){
    	this.G = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.G[i][j] = M[i][j];
			}
    	}
    }
    
	public void setG(int i, int j, double g) {
		this.G[i][j] = g;
	}
    
    /*
     * D
     */
	public double[] getD() {
		if(this.D==null) return null;
		assert(this.D !=null);
		return this.D;
	}
	
	public void setD(double[] v) {
		this.D = new double[v.length];
    		for (int i = 0; i < v.length; i++) {
			this.D[i] = v[i];
		}    	
	}
	
	public void setD(double v, int row) {
		this.D[row] = v;  	
	}
	
	 /*
     * Q
     */
	public double[] getQ(){
    	if(this.Q==null) return null;
    	assert(this.Q!= null);
    	return this.Q;
    }
    
    public void setQ(double[] x){
	    	this.Q = new double[x.length];
	    	for (int i = 0; i < x.length; i++) {
					this.Q[i] = x[i];
	    	}
    }
	
    public void setQ(double v, int row) {
		this.Q[row] = v;  	
	}
    
	 /*
     * R
     */
    public double[] getR(){
    	if(this.R==null) return null;
    	assert(this.R != null);
    	return this.R;
    }
    
    public void setR(double[] x){
    	this.R = new double[x.length];
    	for (int i = 0; i < x.length; i++) {
			this.R[i] = x[i];
    	}
    }
	
    public void setR(double v, int row) {
		this.R[row] = v;  	
	}
    
    /*
     * U
     */
    public double[] getU(){
    	if(this.U==null) return null;
    	assert(this.U !=null);
    	return this.U;    	
    }
    
    public void setU(double[] v){
    	this.U = new double[v.length];
    	for (int i = 0; i < v.length; i++) {
			this.U[i] = v[i];
    	}
    }
	
    public void setU(double v, int row) {
		this.U[row] = v;  	
	}
    
    /*
     * I
     */
    public double[] getI(){
    	if(this.I==null) return null;
    	assert(this.I !=null);
    	return this.I;    	
    }
    
    public void setI(double[] v){
    	this.I = new double[v.length];
    	for (int i = 0; i < v.length; i++) {
			this.I[i] = v[i];
    	}
    }
    
    /*
     * Effect1
     */
    public double[][] getEffect1() {
		if(this.Effect1 ==null) return null;
		assert(this.Effect1 !=null);
		return this.Effect1;
	}
	
	public void setEffect1(double[][] M) {
		if(M.length==0) {
			this.Effect1 = new double[0][];
			return;
		}
		this.Effect1 = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.Effect1[i][j] = M[i][j];
			}
    	}	
    }
	
	public void setEffect1(double eff, int count, int from) {
		this.Effect1[count][1] = from;
		this.Effect1[count][2] = eff;
	}
	
	/*
	 * Effect2
	 */
	public double[][] getEffect2() {
		if(this.Effect2 ==null) return null;
		assert(this.Effect2 !=null);
		return this.Effect2;
	}
	
	public void setEffect2(double[][] M) {
		if(M.length == 0) {
			this.Effect2 = new double[0][];
			return;
		}
		this.Effect2 = new double[M.length][M[0].length];
    	for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				this.Effect2[i][j] = M[i][j];
			}
    	}	
    }
	
	public void setEffect2(double eff, int count, int from1, int from2) {
		this.Effect2[count][1] = from1;
		this.Effect2[count][2] = from2;
	    this.Effect2[count][3] = eff;
    }
	
	public void setDelta(double[] x){
		this.Delta = new double[x.length];
	    	for (int i = 0; i < x.length; i++) {
			this.Delta[i] = x[i];
	    	}
	}
	
	public double[] getDelta() {
		if(this.Delta==null) return null;
		assert(this.Delta !=null);
		return this.Delta;
	}
	
	public void setNu(double[] x){
		this.Nu = new double[x.length];
	    	for (int i = 0; i < x.length; i++) {
	    		this.Nu[i] = x[i];
	    	}
	}

	public double[] getNu() {
		if(this.Nu==null) return null;
		assert(this.Nu !=null);
		return this.Nu;
	}
	
	public void setSkew_mean_t(double[][] x){
		this.Skew_mean_t = new double[x.length][x[0].length];
	    	for (int i = 0; i < x.length; i++) {
			this.Skew_mean_t[i] = x[i];
	    	}
	}

	public double[][] getSkew_mean_t() {
		if(this.Skew_mean_t==null) return null;
		assert(this.Skew_mean_t !=null);
		return this.Skew_mean_t;
	}
	
	public void setSkew_sq_mean_t(double[][] x){
		this.Skew_sq_mean_t = new double[x.length][x[0].length];
	    	for (int i = 0; i < x.length; i++) {
			this.Skew_sq_mean_t[i] = x[i];
	    	}
	}

	public double[][] getSkew_sq_mean_t() {
		if(this.Skew_sq_mean_t==null) return null;
		assert(this.Skew_sq_mean_t !=null);
		return this.Skew_sq_mean_t;
	}
	
	public void setKurtosis_mean_t(double[][] x){
		this.Kurtosis_mean_t = new double[x.length][x[0].length];
	    for (int i = 0; i < x.length; i++) {
	    	for (int j = 0; j < x[0].length; j++) {
		    	this.Kurtosis_mean_t[i][j] = x[i][j];
			}
	    }
	}

	public double[][] getsetKurtosis_mean_t() {
		if(this.Kurtosis_mean_t==null) return null;
		assert(this.Kurtosis_mean_t !=null);
		return this.Kurtosis_mean_t;
	}
	
	public double[][] getKurtosis_inv_R_t() {
		double[][] Kurtosis_inv_R = new double[this.Kurtosis_mean_t.length][this.R.length];
		for (int t = 0; t < this.Kurtosis_mean_t.length; t++) {
			for (int i = 0; i < this.R.length; i++) {
				Kurtosis_inv_R[t][i] = this.R[i] / this.Kurtosis_mean_t[t][i];
			}
		}
		return Kurtosis_inv_R;
	}
}
