package Hasegawa.TimeSeries.Linear.SSM;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.TimeSeries.Parameter;
import Hasegawa.matrix.Matrix;
import Hasegawa.stat.simpleMath;
import RandomGenerator.Sfmt;

public class ssmParameter extends Parameter{
	
	protected double[][] ssmD;
	protected double[][] Phi;
	protected double[][] rootRinvHG;
    
    /*
     * ssmD
     */
    public double[][] getssmD(){
    	if(this.ssmD==null) return null;
    	assert(this.ssmD!= null);
    	return this.ssmD.clone();
    }
    
    public void setssmD(double[][] A){
    	this.ssmD = new double[A.length][A[0].length];
    	for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				this.ssmD[i][j] = A[i][j];
			}
    	}
    }
    
    /*
     * Phi
     */
    public double[][] getPhi(){
    	if(this.Phi==null) return null;
    	assert(this.Phi !=null);
    	return this.Phi.clone();    	
    }
    
    public void setPh(double[][] A){
    	this.Phi = new double[A.length][A[0].length];
    	for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				this.Phi[i][j] = A[i][j];
			}
    	}
    }
    
    /*
     * rootRinvHG
     */
    public double[][] getrootRinvHG(){
    	if(this.rootRinvHG==null) return null;
    	assert(this.rootRinvHG !=null);
    	return this.rootRinvHG.clone();    	
    }
    
    public void setrootRinvHG(double[][] A){
    	this.rootRinvHG = new double[A.length][A[0].length];
    	for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[0].length; j++) {
				this.rootRinvHG[i][j] = A[i][j];
			}
    	}
    }
	
	public ssmParameter(final int syd, Sfmt sf, ssmSetting setting, TimeSeriesDataArray tsd, Matrix Calculator){

		/* Mu and v0 */
		this.x0_ = new double[tsd.repSize][syd];
		this.v0 = new double[syd][syd];
		Calculator.setDiagGaussian(this.v0, setting.Mean_of_Mu0, setting.Var_of_Mu0, sf);
		for (int i = 0; i < syd; i++) this.v0[i][i]=Math.abs(this.v0[i][i]);
				
		/* F */
		this.F = new double[syd][syd];
		Calculator.setRandomvalue(this.F, setting.UPPer_bound_of_F, setting.Lower_bound_of_F, sf);
		
		/* H */
		this.H = new double[tsd.elementNum][syd];
		Calculator.setGaussian(this.H, setting.Mean_of_H , setting.SD_of_H, sf);
		
		/* G */
		if(setting.Drug) {
			this.G = new double[syd][tsd.drugMulRepSize[0].length];
			Calculator.setRandomvalue(this.G, setting.UPPer_bound_of_F, setting.Lower_bound_of_F, sf);
			/* rootRinvHG */
			this.rootRinvHG = new double[tsd.elementNum][tsd.drugMulRepSize[0].length];
		}
		
		/* Q */
		this.Q = new double[syd];
		Calculator.setvalue(this.Q, 1);
		
        /* R  */
		this.R = new double[tsd.elementNum];
		
		for (int i = 0; i < tsd.elementNum; i++) {
			double mean = 0;
			double val =0;
			double var = 0;
			for (int j = 0; j < tsd.observationalTimeNum; j++) {
				val = tsd.ObservationData[j][i];
            	mean += val;
                var += val * val;
            }
			mean /= (double)tsd.observationalTimeNum;
			var /= (double)tsd.observationalTimeNum;
			var -= mean * mean;
				
			this.R[i]=var;	

		}
		if(setting.R_rI <= 0.0) {
			Calculator.setvalue(this.R, new simpleMath().getMean(this.R));	
		}
		
		/* I */
		this.I = new double[syd];
		Calculator.setvalue(I, 1.0);

		/* Phi */
		this.Phi = new double[tsd.elementNum][tsd.elementNum];
		
		/* ssmD */
		this.ssmD = new double[syd][tsd.elementNum];
		
	}
    
    public void storeParameters(ssmParameter Pa){
    	this.setParameters(Pa.getx0(), Pa.getv0(), Pa.getA(), Pa.getF(), Pa.getH(), Pa.getG(), Pa.getD(), 
    			Pa.getQ(), Pa.getR(), Pa.getU(), Pa.getI(), null, null, null, null, null, null, null);
    	if(Pa.getssmD()!=null) this.setssmD(Pa.getssmD());
    	if(Pa.getPhi()!=null) this.setPh(Pa.getPhi());
    	if(Pa.getrootRinvHG()!=null) this.setrootRinvHG(Pa.getrootRinvHG());    	
    }
}