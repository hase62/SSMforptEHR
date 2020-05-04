package Hasegawa.TimeSeries.Linear.VARSSM;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.TimeSeries.Parameter;
import Hasegawa.matrix.Matrix;
import RandomGenerator.Sfmt;

public class vsParameter extends Parameter{

	public double[] L1;
	public double[] L1h;
	
	public vsParameter(final int sysDim, Sfmt sf, vsSetting set, TimeSeriesDataArray tsda, Matrix Calculator, double initial){
		
		/* Mu0 */
		this.x0_ = new double[tsda.repSize][sysDim];
		if(set.Mean_of_Mu0==-1.0){
			for (int rep = 0; rep < tsda.repSize; rep++) {
				for (int i = 0; i < sysDim; i++) {
					/* Need to Change */
					this.x0_[rep][i] = tsda.ObservationData[0][i];
				}
			}
		}
		this.v0 = new double[sysDim][sysDim];
		Calculator.setDiagvalue(this.v0, Math.abs(set.Var_of_Mu0));
				
		/* F, A and D */
		this.F = new double[sysDim][sysDim];
		this.A = new double[sysDim][sysDim];
		if(set.Degradation < 0.0){ 
			Calculator.setDiagvalue(this.F, 1 + set.Degradation);
			Calculator.setDiagvalue(this.A, initial);
		} else if(set.Degradation == 0.0){ 
			Calculator.setDiagvalue(this.F, initial);
			Calculator.copy(this.A, this.F);
		} else {
			this.D = new double[this.F.length];
			Calculator.setvalue(this.D, set.Degradation);
			Calculator.setDiagvalue(this.F, 1 - set.Degradation + initial);
			Calculator.setDiagvalue(this.A, initial);
		}
		
		/* G */
		if(set.Drug) {
			this.G = new double[sysDim][tsda.drugMulRepSize[0].length];
			Calculator.setvalue(this.G, 0.00);
		}
		
		/* U */
		if(set.Input){
			this.U = new double[sysDim];
			Calculator.setvalue(this.U, 0.0);
			if(set.Degradation!=0.0) {
				for (int i = 0; i < this.U.length; i++) {
					this.U[i] = tsda.ObservationData[0][i]*set.Degradation;
				}
			}
		}
		
		/* Q */
		this.Q = new double[sysDim];
		if(set.upDateQ > 0) Calculator.setvalue(this.Q, set.upDateQ);
		else Calculator.setvalue(this.Q, set.Var_of_Mu0);
		
		/* H */
		if(sysDim != tsda.elementNum){
			this.H = new double[tsda.elementNum][sysDim];
			Calculator.setGaussian(this.H, 0, 0.1, sf);
		}
		
        /* R  */
		this.R = new double[tsda.elementNum];
		double sum = 0;
		for (int i = 0; i < tsda.elementNum; i++) {
			double mean = 0;
			double val =0;
			double vs = 0;
			for (int j = 0; j < tsda.observationalTimeNum; j++) {
				val = tsda.ObservationData[j][i];
				mean += val;
				vs += val * val;
            }
			mean /= (double)tsda.observationalTimeNum;
			vs /= (double)tsda.observationalTimeNum;
			vs -= mean * mean;
			this.R[i]=vs / 5.0;
			if(this.R[i] > Math.abs(set.R_rI)) this.R[i] = Math.abs(set.R_rI);
			sum += this.R[i];
		}
       	if(set.R_rI <= 0.0) Calculator.setvalue(this.R, sum/this.R.length);
       	
       	/* I */
       	this.I = new double[sysDim];
		Calculator.setvalue(this.I, 1.0);
		
		/* LASSO */
		this.L1 = new double[sysDim];
		Calculator.setvalue(this.L1, 0.01);
		if(this.H != null) {
			this.L1h = new double[tsda.elementNum];
			Calculator.setvalue(this.L1h, 0.01);
		}
		
		/* st */
		this.Delta = new double[tsda.elementNum];
		//Calculator.setRandomvalue(this.Delta, -0.1, 0.1);
		//Calculator.setvalue(this.Delta, 0.1);

		this.Nu = new double[tsda.elementNum];
		if(set.nu_kurtosis != null) Calculator.copy(this.Nu, set.nu_kurtosis);
		
		this.Skew_mean_t = new double[tsda.allTime][tsda.elementNum];
		this.Skew_sq_mean_t = new double[tsda.allTime][tsda.elementNum];
		Calculator.setvalue(this.Skew_mean_t, 0.0);
		Calculator.setvalue(this.Skew_sq_mean_t, 1.0);
		
		this.Kurtosis_mean_t = new double[tsda.allTime][tsda.elementNum];
		Calculator.setvalue(this.Kurtosis_mean_t, 1.0);
	}

	public void setParameters(vsParameter stPa, Matrix Calculator){
		this.setParameters(stPa.x0_, stPa.v0, stPa.A, stPa.F, stPa.H, stPa.G, stPa.D, stPa.Q, stPa.R, stPa.U, stPa.I, 
						   stPa.Delta, stPa.Nu, stPa.Skew_mean_t, stPa.Skew_sq_mean_t, stPa.Kurtosis_mean_t, null, null);
	}
}
