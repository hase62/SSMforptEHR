package Hasegawa.TimeSeries.Linear.stFilter;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.TimeSeries.Linear.VARSSM.vsSetting;
import Hasegawa.TimeSeries.Linear.VARSSM.vsStorage;
import Hasegawa.TimeSeries.Linear.VARSSMEx.vseInference;
import Hasegawa.matrix.Matrix;
import RandomGenerator.Sfmt;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import java.util.ArrayList;

public class stInference extends vseInference{

	protected org.apache.commons.math3.distribution.NormalDistribution normal= new NormalDistribution();
	protected org.apache.commons.math3.distribution.ChiSquaredDistribution chi = new ChiSquaredDistribution(3);
	
	public stInference(final int SYD, vsSetting SET, TimeSeriesDataArray TSDA, Matrix MP, Sfmt Sf, vsStorage vStor) {
		super(SYD, SET, TSDA, MP, Sf, vStor);
	}
	
	protected void getParameters() {
		this.LogRDet = this.vPa.getParameters(this.x0_, this.v0, this.A, this.F, this.Ftrans, this.H, this.G, 
											 this.RinvH, this.HtransRinvH, this.D, this.Q, this.Qinv, this.R, 
											 this.Rinv, this.RinvMatrix, this.U, this.I, 
											 this.Delta, 
											 this.Nu, 
											 this.Skew_mean_t, 
											 this.Skew_sq_mean_t, 
											 this.Kurtosis_mean_t, 
											 this.Delta_skew_mean_t, 
											 this.Kurtosis_mean_Rinv_t, 
											 this.Htrans_Kurtosis_mean_Rinv_H_t, 
											 this.Kurtosis_mean_Rinv_H_t, null, null);
		for (int rep = 0; rep < this.tsda.repSize; rep++) {
			this.Calculator.copy(this.x0_s[rep], this.x0_[rep]); 
			this.Calculator.copy(this.v0_s[rep], this.v0); 
		}
		this.Calculator.setvalue(this.predictionError, 0);
		this.Calculator.setvalue(this.predictionAbility, 0);
	    if(this.H != null) {
			/* Not Required? */
			this.Calculator.multAB(this.Rinv, this.H, this.RinvH);
			this.Calculator.multAtB(this.H, this.RinvH, this.HtransRinvH);
			for (int time = 0; time < this.tsda.allTime; time++) {
				this.Calculator.multDx(this.Delta, this.Skew_mean_t[time], this.Delta_skew_mean_t[time]);				
				this.Calculator.multDx(this.Kurtosis_mean_t[time], this.Rinv, this.Kurtosis_mean_Rinv_t[time]);
				this.Calculator.multAB(this.Kurtosis_mean_Rinv_t[time], this.H, this.Kurtosis_mean_Rinv_H_t[time]);
				this.Calculator.multAtB(this.H, this.Kurtosis_mean_Rinv_H_t[time], this.Htrans_Kurtosis_mean_Rinv_H_t[time]);
			}
			/* Not Required? */
	    }
	}
	
	public void setParameters(){
		if(vSet.R_rI <= 0) this.Calculator.setvalue(this.R, this.Calculator.sumofVector(this.R)/this.R.length);
		if(vSet.upDateQ < 0) this.Calculator.setvalue(this.Q, this.Calculator.sumofVector(this.Q)/this.Q.length);
		this.vPa.setParameters(this.x0_, this.v0, this.A, this.F, this.H, this.G, this.D, this.Q, 
				   this.R, this.U, this.I, this.Delta, this.Nu, this.Skew_mean_t, 
				   this.Skew_sq_mean_t, this.Kurtosis_mean_t, null, null);
	}
	
	/*
	 * Store current parameters to vSto
	 */
	protected void storeCurrentSettings(vsStorage vSto){
		/* Parameters */
		vSto.vPa.setParameters(this.x0_, this.v0, this.A, this.F, this.H, this.G, this.D, this.Q, 
							   this.R, this.U, this.I, this.Delta, this.Nu, this.Skew_mean_t, 
							   this.Skew_sq_mean_t, this.Kurtosis_mean_t, null, null);
		
		/* Profiles */
		vSto.setExpectations(this.Txx_obs, this.Txx, this.Tyx, this.Txx_m, this.Txx_mm, this.Txz_mm, 
							 this.Txz_m, this.sum_x, this.sum_xm, this.sum_zm, this.Skew_mean_t, 
							 this.Skew_sq_mean_t, this.Delta_skew_mean_t, this.Kurtosis_mean_t, 
							 this.Kurtosis_mean_Rinv_t, this.Htrans_Kurtosis_mean_Rinv_H_t, 
							 this.Kurtosis_mean_Rinv_H_t, 
							 this.Kurtosis_skew_sq_mean_sum, 
							 this.Dif_y_x_delta_skew_mean_sum, 
							 this.Dif_y_x_skew_mean_sum, 
							 this.Dif_y_x_sq_sum, this.Skew_mean_x_sum, this.Kurtosis_H_Sigma_Htrans_sum, 
							 this.Calculator);
		vSto.setProfiles(this.x0_s, this.x_p, this.x_f, this.x_s, this.y_p);
		
		/* likelihood */
		vSto.logLikelihood = this.currentLogLikelihood;
		vSto.Criterion = this.Criterion;
		this.Calculator.copy(vSto.predictionError, this.predictionError);
		this.Calculator.copy(vSto.predictionAbility, this.predictionAbility);
				
		/* Penalty */
		this.Calculator.copy(vSto.vPa.L1, this.vPa.L1);
		if(this.H != null) this.Calculator.copy(vSto.vPa.L1h, this.vPa.L1h);
	}

	/*
	 * Recall vSto setting to current setting
	 */
	public void recallPreviousSettings(vsStorage vSto){
		/* Parameters */
		this.LogRDet = vSto.vPa.getParameters(this.x0_, this.v0, this.A, this.F, this.Ftrans, this.H, this.G, 
											  this.RinvH, this.HtransRinvH, this.D, this.Q, this.Qinv, this.R, 
											  this.Rinv, this.RinvMatrix, this.U, this.I, 
											  this.Delta, 
											  this.Nu, 
											  this.Skew_mean_t, 
											  this.Skew_sq_mean_t, 
											  this.Kurtosis_mean_t, 
											  this.Delta_skew_mean_t, 
											  this.Kurtosis_mean_Rinv_t, 
											  this.Htrans_Kurtosis_mean_Rinv_H_t, 
											  this.Kurtosis_mean_Rinv_H_t, null, null);

		/* Profiles */
		/* Not Required? */
		vSto.getExpectations(this.Txx_obs, this.Txx, this.Tyx, this.Txx_m, this.Txx_mm, this.Txz_mm, this.Txz_m, 
							this.sum_x, this.sum_xm, this.sum_zm, this.Skew_mean_t, 
							this.Skew_sq_mean_t, this.Delta_skew_mean_t, this.Kurtosis_mean_t, 
							this.Kurtosis_mean_Rinv_t, this.Htrans_Kurtosis_mean_Rinv_H_t, 
							this.Kurtosis_mean_Rinv_H_t, 
							this.Kurtosis_skew_sq_mean_sum, 
							this.Dif_y_x_delta_skew_mean_sum, 
							this.Dif_y_x_skew_mean_sum, 
							this.Dif_y_x_sq_sum, this.Skew_mean_x_sum, 
							this.Kurtosis_H_Sigma_Htrans_sum, 
							this.Calculator);
		vSto.getProfiles(this.x0_s, this.x_p, this.x_f, this.x_s, this.y_p, this.Calculator);
		
		/* likelihood */
		this.currentLogLikelihood = this.previousLogLikelihood = vSto.logLikelihood;
		this.Criterion = vSto.Criterion;
		this.Calculator.copy(this.predictionError, vSto.predictionError);
		this.Calculator.copy(this.predictionAbility, vSto.predictionAbility);
		
		/* Penalty */
		this.Calculator.copy(this.vPa.L1, vSto.vPa.L1);
		if(this.H != null) this.Calculator.copy(this.vPa.L1h, vSto.vPa.L1h);
	}
	
	protected int Filter(int time, boolean observed, boolean validity, int rep_num, int count){
		if(this.H != null) this.Calculator.multAx(this.H, this.x_p[time], this.y_p[time]);
		else this.Calculator.copy(this.y_p[time], this.x_p[time]);
		if (observed && validity) {
			if (this.H == null){
				/* V_t */
				this.Calculator.inversionTheorem(this.v_p[time], this.Kurtosis_mean_Rinv_t[time], this.v_f[time]);
				this.Calculator.changesymmetric(this.v_f[time]);
				/* x_t */
				this.Calculator.copy(this.x_f[time], this.x_p[time]);
				this.Calculator.multAB(this.v_f[time], this.Kurtosis_mean_Rinv_t[time], this.v_fRinv);
				this.Calculator.sub(this.y[this.tsda.repIDSum.get(rep_num) + count], this.x_p[time], this.Delta_skew_mean_t[time], this.y_ws);
				this.Calculator.multAddAx(this.v_fRinv, this.y_ws, this.x_f[time]);
			} else {
				/* V_t */
				double[][] IKH = new double[this.sysDim][this.sysDim];
				this.Calculator.symmetricInverse(this.v_p[time], IKH);
				this.Calculator.changesymmetric(IKH);
				this.Calculator.add(IKH, this.Htrans_Kurtosis_mean_Rinv_H_t[time]);
				this.Calculator.symmetricInverse(IKH, this.v_f[time]);
				this.Calculator.changesymmetric(this.v_f[time]);

				/* x_t */
				this.Calculator.copy(this.x_f[time], this.x_p[time]);
				this.Calculator.sub(this.y[this.tsda.repIDSum.get(rep_num) + count], this.Delta_skew_mean_t[time], this.y_ws);				
				this.Calculator.multAddABtx(this.v_f[time], this.Kurtosis_mean_Rinv_H_t[time], this.y_ws, this.x_f[time]);
				this.Calculator.multSubABx(this.v_f[time], this.Htrans_Kurtosis_mean_Rinv_H_t[time], this.x_p[time], this.x_f[time]);
			}
		} else {
			this.Calculator.copy(this.x_f[time], this.x_p[time]);
			this.Calculator.copy(this.v_f[time], this.v_p[time]);
		}
		if(observed) count++;
		this.Calculator.copy(this.x_s[time], this.x_f[time]);
		this.Calculator.copy(this.v_s[time], this.v_f[time]);
		return count;
	}
	
	protected void CalculateComplement() {
		/* Clear Sums */
		this.Calculator.setvalue(this.Kurtosis_skew_sq_mean_sum, 0);
		this.Calculator.setvalue(this.Kurtosis_H_Sigma_Htrans_sum, 0);

		this.Calculator.setvalue(this.Dif_y_x_sq_sum, 0);
		this.Calculator.setvalue(this.Dif_y_x_abs_max_sign, 0);
		this.Calculator.setvalue(this.Dif_y_x_delta_skew_mean_sum, 0);
		this.Calculator.setvalue(this.Dif_y_x_skew_mean_sum, 0);
		this.Calculator.setvalue(this.Skew_mean_x_sum, 0);
		this.Calculator.setvalue(this.Tyx, 0);
		
		double[] Skew_parental_mean = new double[this.tsda.elementNum];
		double[] Skew_parental_var = new double[this.tsda.elementNum];
		double[] K_u = new double[this.tsda.elementNum];
		for (int i = 0; i < this.tsda.elementNum; i++) {
    		K_u[i] = this.Delta[i] / (this.Delta[i] * this.Delta[i] + this.R[i]);
    	}
		
		for (int rep = 0; rep < this.tsda.repSize; rep++) {
			int count = 0;
			int time = 0;
			int ytime = 0;
			int MaxTimeOfRep = 0;
			for (int j = this.tsda.TimeArray.get(rep).size() - 1; j >= 0 ; j--) {
				if(this.tsda.Validity.get(rep).get(j)) {
					MaxTimeOfRep = this.tsda.TimeArray.get(rep).get(j);
					break;
				}
			}
			for (int j = 0; j < MaxTimeOfRep; j++) {
				time = this.tsda.maxTime * rep + j;
				if (this.tsda.TimeArray.get(rep).get(count) == j + 1 && this.tsda.Validity.get(rep).get(count)) {
					ytime = this.tsda.repIDSum.get(rep) + count;
					
					/* Parental Skew Mean */
					this.Calculator.copy(this.y_ws, this.y[ytime]);
					this.Calculator.multSubAx(this.H, this.x_s[time], this.y_ws);
					this.Calculator.multAx(K_u, this.y_ws, Skew_parental_mean);
					for (int i = 0; i < this.Dif_y_x_abs_max_sign.length; i++) {
						//if(Math.abs(this.Dif_y_x_abs_max_sign[i]) < Math.abs(this.y_ws[i]))
							this.Dif_y_x_abs_max_sign[i] += this.y_ws[i] * this.y_ws[i] * this.y_ws[i]; 
					}
					
					/* Parental Skew Variance */
					for (int i = 0; i < Skew_parental_var.length; i++) {
						Skew_parental_var[i] = (1 - K_u[i] * this.Delta[i]) / this.Kurtosis_mean_t[time][i];
					}
					
					/* Skew */
					for (int i = 0; i < Skew_mean_t[time].length; i++) {
						double t = -1 * Skew_parental_mean[i] / Math.sqrt(Skew_parental_var[i]);
						if(t >  7) {
							t =  7;
						}
						if(t < -7) {
							t = -7;
						}
						double c_t = 1.0 / (Math.sqrt(2.0 * Math.PI) * (1 - normal.cumulativeProbability(t)));
						this.Skew_mean_t[time][i] = Skew_parental_mean[i] + Math.sqrt(Skew_parental_var[i]) * c_t * Math.pow(Math.E, -1.0 * t * t / 2.0);
						this.Delta_skew_mean_t[time][i] = this.Delta[i] * this.Skew_mean_t[time][i];
						if(t > 0) { 
							this.Skew_sq_mean_t[time][i] = Skew_parental_var[i] * c_t * Math.sqrt(Math.PI / 2.0) * (1 - chi.cumulativeProbability(t * t))
									+ 2 * Skew_parental_mean[i] * this.Skew_mean_t[time][i] + Skew_parental_mean[i] * Skew_parental_mean[i];
						} else {
							this.Skew_sq_mean_t[time][i] = Skew_parental_var[i] * c_t * Math.sqrt(Math.PI / 2.0) * (1 + chi.cumulativeProbability(t * t))
									+ 2 * Skew_parental_mean[i] * this.Skew_mean_t[time][i] + Skew_parental_mean[i] * Skew_parental_mean[i];
						}
					}
					
					/* Kurtosis */
					// y_ws = y_{t} - H x_{t}
					for (int i = 0; i < this.tsda.elementNum; i++) {
						double hsh = this.Calculator.multxtAx(this.v_s[time], this.H[i]);
						double phi = this.y_ws[i] * this.y_ws[i] + hsh; 
						phi -=  2 * this.Delta[i] * this.y_ws[i] * this.Skew_mean_t[time][i];
						phi *= this.Rinv[i];
						phi += (this.Delta[i] * this.Rinv[i] * this.Delta[i] + 1) * this.Skew_sq_mean_t[time][i];
						if(phi > 5 * this.Nu[i]) phi = 5 * this.Nu[i];
						
						this.Kurtosis_mean_t[time][i] = (this.Nu[i] + 2) / (this.Nu[i] + phi);
						this.Kurtosis_skew_sq_mean_sum[i] += this.Kurtosis_mean_t[time][i] * this.Skew_sq_mean_t[time][i];
						
						this.Kurtosis_H_Sigma_Htrans_sum[i] += this.Kurtosis_mean_t[time][i] * hsh;
						this.Dif_y_x_sq_sum[i] += this.y_ws[i] * this.y_ws[i] * this.Kurtosis_mean_t[time][i];
						this.Dif_y_x_delta_skew_mean_sum[i] += this.y_ws[i] * this.Delta[i] * this.Skew_mean_t[time][i] * this.Kurtosis_mean_t[time][i];
						this.Dif_y_x_skew_mean_sum[i] += this.y_ws[i] * this.Skew_mean_t[time][i] * this.Kurtosis_mean_t[time][i];
					}
					
					this.Calculator.multDx(this.Kurtosis_mean_t[time], this.Rinv, this.Kurtosis_mean_Rinv_t[time]);
					this.Calculator.multAB(this.Kurtosis_mean_Rinv_t[time], this.H, this.Kurtosis_mean_Rinv_H_t[time]);
					this.Calculator.multAtB(this.H, this.Kurtosis_mean_Rinv_H_t[time], this.Htrans_Kurtosis_mean_Rinv_H_t[time]);
					
					/* For Update */
					this.Calculator.multAddD2xyt(this.Skew_mean_t[time], this.x_s[time], this.Delta, this.Kurtosis_mean_t[time], this.Skew_mean_x_sum);
					this.Calculator.multAddDxyt(this.y[ytime], this.x_s[time], this.Kurtosis_mean_t[time], this.Tyx);

					count++;
					if (count == this.tsda.TimeArray.get(rep).size()) count = 0;
				}
			}
		}

		/* Save Skew */
		if(this.miss_update_eternal < 150) {
			this.vPa.setSkew_mean_t(this.Skew_mean_t);
			this.vPa.setSkew_sq_mean_t(this.Skew_sq_mean_t);
		}
		
		if(this.miss_update_eternal < 300) {
			this.vPa.setKurtosis_mean_t(this.Kurtosis_mean_t);	
		}
	}
	
	protected double logLikelihood(double[][] HtRiH, double[][] ssmH, double[][] RiH) {
		double logLikelihood = 0;
		double[][] LU= new double[this.sysDim][this.sysDim];
		double[] yx= new double[this.tsda.elementNum];
		double[][] L_d= new double[this.sysDim][this.sysDim];
		double[][] U_d= new double[this.sysDim][this.sysDim];
		for (int i = 0; i < this.tsda.repSize; i++) {
			int count = 0;
			for (int k = 0; k < this.tsda.maxTime; k++) {
				if (this.tsda.TimeArray.get(i).get(count) == k + 1) {
					int time = this.tsda.maxTime * i + k;
					int num = this.tsda.repIDSum.get(i) + count;
					double[] ws = new double[this.tsda.elementNum];
					double[] ws2 = new double[this.sysDim];
					
					// yx = y - H * xp
					this.Calculator.copy(yx, this.y[num]);
					if(ssmH == null) Calculator.sub(yx, this.x_p[time]);
					else this.Calculator.multSubAx(ssmH, this.x_p[time], yx);
					this.Calculator.sub(yx, this.Delta_skew_mean_t[time]);
					this.Calculator.copy(ws, yx);
					
					if(ssmH == null){
						if(!this.tsda.Validity.get(i).get(count)){
							for (int j = 0; j < yx.length; j++) {
								this.predictionAbility[j] += yx[j] * yx[j];
							}
							count++;
							if (count == this.tsda.TimeArray.get(i).size())	count = 0;
							continue;
						}
						for (int j = 0; j < yx.length; j++) {
							this.predictionError[j] += yx[j] * yx[j];
						}
					}

					/**
					 * Wood-burry Equation -> (A + UCV)^{-1} = A^{-1} - A^{-1}*U*(C^{-1} + V*A^{-1}*U)^{-1}*V*A^{-1}
					 * (i) : (H*Vp*Ht + R)^{-1} = Ri - Ri*H*(Vp^{-1} + HtRiH)^{-1}*Ht*Ri
					 * Here, (Vp^{-1} + HtRiH)^{-1} = Vf, 
					 * (i) = Ri - RiH*Vf*HtRi
					 * (ii) (y - Hxp)*(i)*(y - Hxp)t = (iii)ws*Ri*ws - (iv)(ws*RiH)*Vf*(HtRi*ws)
					 */
					
					// (iii) 
					this.Calculator.multxtD(ws, this.Kurtosis_mean_Rinv_t[time]);
					logLikelihood -= 0.5 *this.Calculator.dotProduct(ws, yx);
					
					// (iv)
					if(RiH == null) this.Calculator.multAx(this.Kurtosis_mean_Rinv_t[time], yx, ws2);
					else this.Calculator.multAtx(this.Kurtosis_mean_Rinv_H_t[time], yx, ws2);
					logLikelihood += 0.5 * this.Calculator.multxtAx(this.v_f[time], ws2);
					
					/**
					 * |I +UVt| = |I + UtV|
					 *  |(H*Vp*Ht + K^{-1}R)| = |(H*Vp*Ht*KRi + I) * K^{-1}R| = |(H*Vp*Ht*KRi + I)| * |K^{-1}R|
					 *  	= |(Vp*Ht*Ri*H + I)| * |K^{-1}R|
					 */
					// LU = I + Vp*HtRiH
					this.Calculator.copy(LU, this.I);
					/* must set HtRiH */
					this.Calculator.multAddAtB(this.v_p[time], this.Htrans_Kurtosis_mean_Rinv_H_t[time], LU);
					this.Calculator.LUDecomposition(LU, L_d, U_d);
					logLikelihood -= 0.5 * this.Calculator.diagLogDeterminant(U_d);
					logLikelihood -= 0.5 * this.LogRDet;
					logLikelihood -= 0.5 * this.Calculator.diagInvLogDeterminant(this.Kurtosis_mean_t[time]);
					
					logLikelihood -= 0.5 * this.tsda.elementNum * Math.log(2 * Math.PI);
					
					count++;
					if (count == this.tsda.TimeArray.get(i).size())	count = 0;
				}
			}
		}
		return logLikelihood;
	}
	
	protected void Update(boolean variable_cutting, boolean allUpdate, boolean update_H, 
						  boolean initialLoop, boolean without_L1_fix_variables) {
		CalculateComplement();
		
		/*
		 * update A, F and G
		 */
		this.Calculator.setvalue(this.A_up_ws, 0);
		this.Calculator.setvalue(this.F_up_ws, 0);
		if(this.vSet.Drug) this.Calculator.setvalue(this.G_up_ws, 0);
		double[] aws = new double[this.sysDim];
		double[] gws = null;
		int min_search_edges_A = Math.max(Math.min((int)(this.th_EdgeNum / 2) - 1, this.sysDim), 3);
		int min_search_edges_G = 0;
		if(this.vSet.Drug) {
			min_search_edges_G = Math.max(Math.min((int)(this.th_EdgeNum / 2) - 1, this.z0_[0].length), 3);
			if(this.G[0].length > 30) min_search_edges_G--;
		}
		
		if(this.vSet.Drug) gws = new double[this.tsda.drugMulRepSize[0].length];
		for (int i = 0; i < this.F.length; i++) {
			/* A (= F + diag(1 - d)) */
			this.Calculator.copy(aws, this.Txx_m[i]);
			if(this.vSet.Drug)	this.Calculator.multSubAx(this.Txz_mm, this.G[i], aws);
			if(this.vSet.Input) {
				this.Calculator.copy(this.x_rescale_ws, this.sum_xm);
				this.Calculator.rescale(this.x_rescale_ws, this.U[i]);
				this.Calculator.sub(aws, this.x_rescale_ws);
			}
			if(this.vSet.Degradation != 0.0) {
				/* Set D[i] 0 */
				this.Calculator.copy(this.x_rescale_ws, this.Txx_mm[i]);
				if(this.vSet.Degradation < 0) {
					this.Calculator.rescale(this.x_rescale_ws, (1 + this.vSet.Degradation));
				} else if(this.vSet.Degradation > 0){
					this.Calculator.rescale(this.x_rescale_ws, (1 - this.vSet.Degradation) - this.D[i]);
				}
				this.Calculator.sub(aws, this.x_rescale_ws);
			}
			
			/* G */
			if(this.vSet.Drug){
				this.Calculator.copy(gws, this.Txz_m[i]);
				this.Calculator.transpose(this.Txz_mm, this.Txz_mm_trans);
				this.Calculator.multSubAx(this.Txz_mm_trans, this.F[i], gws);
				if(this.vSet.Input) {
					this.Calculator.copy(this.z_rescale_ws, this.sum_zm);
					this.Calculator.rescale(this.z_rescale_ws, this.U[i]);
					this.Calculator.sub(gws, this.z_rescale_ws);
				}
			}
			
			/* Determine Next Active Set */
			this.activeSetA.clear();
			this.argMinA = 0;
			if(this.vSet.Drug) {
				this.activeSetG.clear();
				this.argMinG = 0;
			}
			
			if(allUpdate || (this.calculationOrder[(int)this.updatingRow]==i && !variable_cutting && initialLoop)){
				ArrayList<Integer> tempActiveSetA = new ArrayList<Integer>();
				for (int j = 0; j < this.F[0].length; j++) {
					if(this.vSet.Degradation > 0 && i == j) continue;	
					tempActiveSetA.add(j);
				}
				this.calcArgMin(tempActiveSetA, aws, i, 0, this.Q[i], this.Qinv[i]);
				this.activeSetAList[i] = new ArrayList<Integer>(activeSetA);
				
				if(this.vSet.Drug){
					ArrayList<Integer> tempActiveSetG = new ArrayList<Integer>();
					for (int j = 0; j < this.G[0].length; j++) {
							tempActiveSetG.add(j);
					}
					this.calcArgMin(tempActiveSetG, gws, i, 1, this.Q[i], this.Qinv[i]);
					this.activeSetGList[i] = new ArrayList<Integer>(activeSetG);
				}
			} else if((this.calculationOrder[(int)this.updatingRow]==i && !variable_cutting) || this.activeSetAList[i].size() ==0) {
				/* Test Existing Active Set A */
				this.Finished.clear();
				ArrayList<Integer> tempActiveSetA = new ArrayList<Integer>();
				ArrayList<Integer> nextNonActiveSetA = new ArrayList<Integer>(this.activeSetAList[i]);
				//this.calcArgMin(tempActiveSetA, nextNonActiveSetA, aws, i, 0, (int)this.th_EdgeNum, 0);
				
				/* Test Active Set A +-1 */
				this.Finished.clear();
				for (int remove = -1; remove < this.sm.getActiveRowElements(i, this.A).size(); remove++) {
					tempActiveSetA = new ArrayList<Integer>(this.sm.getActiveRowElements(i, this.A));
					if(remove > -1) tempActiveSetA.remove(remove);
					this.calcArgMin(tempActiveSetA, aws, i, 0, this.Q[i], this.Qinv[i]);
					if(without_L1_fix_variables) break;
					nextNonActiveSetA = new ArrayList<Integer>();
					for (int j = 0; j < aws.length; j++) {
						if(this.vSet.Degradation > 0 && i == j) continue;
						if(!tempActiveSetA.contains(j)) nextNonActiveSetA.add(j);
					}
					this.calcArgMin(tempActiveSetA, nextNonActiveSetA, aws, i, 0, tempActiveSetA.size() + 1, 0, this.Q[i], this.Qinv[i]);
				}
				if(this.activeSetA.size()==0) {
					this.calcArgMin(new ArrayList<Integer>(), new ArrayList<Integer>(this.activeSetAList[i]), 
							aws, i, 0, min_search_edges_A, 0, this.Q[i], this.Qinv[i]);
				}
				this.activeSetAList[i] = new ArrayList<Integer>(this.activeSetA);
				
				/* 
				 * Active Set of G
				 */
				if(this.vSet.Drug){
					/* Test Existing Active Set A */
					this.Finished.clear();
					ArrayList<Integer> tempActiveSetG = new ArrayList<Integer>();
					ArrayList<Integer> nextNonActiveSetG = new ArrayList<Integer>(this.activeSetGList[i]);
					//this.calcArgMin(tempActiveSetG, nextNonActiveSetG, gws, i, 0, (int)this.th_EdgeNum, 1);
					
					/* Test Active Set G +-1 */
					this.Finished.clear();
					for (int remove = -1; remove < this.sm.getActiveRowElements(i, this.G).size(); remove++) {
						tempActiveSetG = new ArrayList<Integer>(this.sm.getActiveRowElements(i, this.G));
						if(remove > -1) tempActiveSetG.remove(remove);
						this.calcArgMin(tempActiveSetG, gws, i, 1, this.Q[i], this.Qinv[i]);
						if(without_L1_fix_variables) break;
						nextNonActiveSetG = new ArrayList<Integer>();
						for (int j = 0; j < gws.length; j++) {
							if(!tempActiveSetG.contains(j)) nextNonActiveSetG.add(j);
						}
						this.calcArgMin(tempActiveSetG, nextNonActiveSetG, gws, i, 0, tempActiveSetG.size() + 1, 1, this.Q[i], this.Qinv[i]);
					}
					if(this.activeSetG.size()==0) {
						this.calcArgMin(new ArrayList<Integer>(), new ArrayList<Integer>(this.activeSetGList[i]), 
								gws, i, 0, min_search_edges_G, 1, this.Q[i], this.Qinv[i]);
					}
					this.activeSetGList[i] = new ArrayList<Integer>(this.activeSetG);
				}
			} else {
				/* A */
				/* Test Existing Active Set A */
				this.Finished.clear();
				ArrayList<Integer> tempActiveSetA = new ArrayList<Integer>();
				ArrayList<Integer> nextNonActiveSetA = new ArrayList<Integer>(this.activeSetAList[i]);
				//this.calcArgMin(tempActiveSetA, nextNonActiveSetA, aws, i, 0, nextNonActiveSetA.size(), 0);
				
				/* Test Active Set A 0-1 */
				this.Finished.clear();
				for (int remove = -1; remove < this.sm.getActiveRowElements(i, this.A).size(); remove++) {
					tempActiveSetA = new ArrayList<Integer>(this.sm.getActiveRowElements(i, this.A));
					if(remove > -1) tempActiveSetA.remove(remove);
					this.calcArgMin(tempActiveSetA, aws, i, 0, this.Q[i], this.Qinv[i]);
					if(without_L1_fix_variables) break;
					for (int j = 0; j < this.activeSetAList[i].size(); j++) {
						if(remove > -1) break;
						nextNonActiveSetA.clear();
						//Search {{{Temporal} and {Temporal -1}} + 1 (from activeSetAList but not used)}
						if(!tempActiveSetA.contains(this.activeSetAList[i].get(j))) {
							nextNonActiveSetA.add(this.activeSetAList[i].get(j));
							this.calcArgMin(tempActiveSetA, nextNonActiveSetA, aws, i, 0, tempActiveSetA.size() + 1, 0, this.Q[i], this.Qinv[i]);							
						}
					}
				}
				
				/* G */
				if(this.vSet.Drug){
					/* Test All Active-Sets*/
					this.Finished.clear();
					ArrayList<Integer> tempActiveSetG = new ArrayList<Integer>();
					ArrayList<Integer> nextNonActiveSetG = new ArrayList<Integer>(this.activeSetGList[i]);
					//this.calcArgMin(tempActiveSetG, nextNonActiveSetG, gws, i, 0, nextNonActiveSetG.size(), 1);
					
					this.Finished.clear();
					for (int remove = -1; remove < this.sm.getActiveRowElements(i, this.G).size(); remove++) {
						tempActiveSetG = new ArrayList<Integer>(this.sm.getActiveRowElements(i, this.G));
						if(remove > -1) tempActiveSetG.remove(remove);
						this.calcArgMin(tempActiveSetG, gws, i, 1, this.Q[i], this.Qinv[i]);
						if(without_L1_fix_variables) break;
						for (int j = 0; j < this.activeSetGList[i].size(); j++) {
							if(remove > -1) break;
							nextNonActiveSetG.clear();
							//Search {{{Temporal} and {Temporal -1}} + 1 (from activeSetAList but not used)}
							if(!tempActiveSetG.contains(this.activeSetGList[i].get(j))) {
								nextNonActiveSetG.add(this.activeSetGList[i].get(j));
								this.calcArgMin(tempActiveSetG, nextNonActiveSetG, gws, i, 0, tempActiveSetG.size() + 1, 1, this.Q[i], this.Qinv[i]);
							}
						}
					}
				}
			}
			
			/* Set Next Parameters */
			if(this.activeSetA.size()!=0) {
				double sum_of_row = 0;
				if(this.vSet.Degradation < 0) sum_of_row += 1 + this.vSet.Degradation;
				for (int j = 0; j < activeSetA.size(); j++) {
					A_up_ws[i][activeSetA.get(j)] = this.nextA.get(j);
					sum_of_row += this.nextA.get(j);
				}
				if(sum_of_row > 1.5){
					for (int j2 = 0; j2 < this.A[i].length; j2++) {
						A_up_ws[i][j2] = this.A[i][j2];
					}
				}
				if(this.Calculator.checkAbsValues(A_up_ws[i], 1.1) || 
					(this.vSet.Degradation < 0 && 1 + this.vSet.Degradation + A_up_ws[i][i] > 1.1)) {
					for (int j2 = 0; j2 < this.A[i].length; j2++) {
						A_up_ws[i][j2] = this.A[i][j2];
					}
				}
				
			}
			
			if(this.vSet.Drug && this.activeSetG.size()!=0){
				for (int j = 0; j < activeSetG.size(); j++) {
					G_up_ws[i][activeSetG.get(j)] = this.nextG.get(j);
				}
			}
		}
		
		/*
		 * update D
		 */
		if(this.vSet.Degradation > 0.0){
			for (int i = 0; i < this.D.length; i++) {
				D_up_ws[i] = this.Txx_m[i][i] - this.Calculator.dotProduct(A_up_ws[i], this.Txx_mm[i]);
				if(this.vSet.Drug) D_up_ws[i] -= this.Calculator.dotProduct(G_up_ws[i], this.Txz_mm[i]);
				if(this.vSet.Input) D_up_ws[i] -= this.U[i] * this.sum_xm[i];
				D_up_ws[i] = (1 - this.vSet.Degradation) - (D_up_ws[i] / this.Txx_mm[i][i]);

				double sum_of_row = (1 - this.vSet.Degradation) - D_up_ws[i];
				for (int j = 0; j < A_up_ws[i].length; j++) {
					sum_of_row += A_up_ws[i][j];
				}
				if(sum_of_row > 1.5){
					for (int j2 = 0; j2 < this.A[i].length; j2++) {
						A_up_ws[i][j2] = this.A[i][j2];						
					}
					D_up_ws[i] = this.D[i];
				}
				if(D_up_ws[i] < 0) D_up_ws[i] = this.D[i];
			}
			this.vPa.setD(D_up_ws);
		}
		if(this.miss_update_eternal < 200) {
			this.vPa.setA(A_up_ws);
			if(this.vSet.Drug) this.vPa.setG(G_up_ws);
		}
		/* 
		 * update F
		 */
		this.Calculator.copy(F_up_ws, A_up_ws);
		if(this.vSet.Degradation < 0.0){
			for (int i = 0; i < this.F.length; i++) {
				F_up_ws[i][i] += (1 + this.vSet.Degradation);
			}
		} else if(this.vSet.Degradation > 0.0){
			for (int i = 0; i < this.F.length; i++) {
				F_up_ws[i][i] += (1 - this.vSet.Degradation) - D_up_ws[i];
			}
		}
		if(this.miss_update_eternal < 200) {
			this.vPa.setF(F_up_ws);
		}
		
		/*
		 * update H
		 */
		this.H_up_ws = null;
		if((this.H != null && !update_H)) {
			this.H_up_ws = new double[this.tsda.elementNum][this.sysDim];
			double[] kurtosis_element = new double[this.sysDim];
			double[] hws = new double[this.sysDim];
			double[][] Txx_obsInv = new double[this.sysDim][this.sysDim];
			for (int i = 0; i < this.H.length; i++) {
				this.Calculator.setvalue(this.Txx_obs, 0);
				for (int rep = 0; rep < this.tsda.repSize; rep++) {
					int count = 0;
					int time = 0;
					int MaxTimeOfRep = 0;
					for (int j = this.tsda.TimeArray.get(rep).size() - 1; j >= 0 ; j--) {
						if(this.tsda.Validity.get(rep).get(j)) {
							MaxTimeOfRep = this.tsda.TimeArray.get(rep).get(j);
							break;
						}
					}
					for (int j = 0; j < MaxTimeOfRep; j++) {
						time = this.tsda.maxTime * rep + j;
						if (this.tsda.TimeArray.get(rep).get(count) == j + 1 && this.tsda.Validity.get(rep).get(count)) {
							this.Calculator.setvalue(kurtosis_element, this.Kurtosis_mean_t[time][i]);
							this.Calculator.multAddDxyt(this.x_s[time], this.x_s[time], kurtosis_element, this.Txx_obs);
							this.Calculator.multAddDA(this.v_s[time], kurtosis_element, this.Txx_obs);
							count++;
							if (count == this.tsda.TimeArray.get(rep).size()) count = 0;
						}
					}
				}
				this.Calculator.symmetricInverse(this.Txx_obs, Txx_obsInv);
				this.Calculator.changesymmetric(Txx_obsInv);
				this.Calculator.copy(hws, this.Tyx[i]);
				this.Calculator.sub(hws, this.Skew_mean_x_sum[i]);
				this.Calculator.multAx(Txx_obsInv, hws, this.H_up_ws[i]);
			}
			
			if(this.miss_update_eternal < 200) {
				this.vPa.setH(H_up_ws);
			}
		} else if(this.H != null && update_H) {
			this.H_up_ws = new double[this.tsda.elementNum][this.sysDim];
			double[] kurtosis_element = new double[this.sysDim];
			double[] hws = new double[this.sysDim];
			for (int i = 0; i < this.H.length; i++) {
				this.Calculator.setvalue(this.Txx_obs, 0);
				for (int rep = 0; rep < this.tsda.repSize; rep++) {
					int count = 0;
					int time = 0;
					int MaxTimeOfRep = 0;
					for (int j = this.tsda.TimeArray.get(rep).size() - 1; j >= 0 ; j--) {
						if(this.tsda.Validity.get(rep).get(j)) {
							MaxTimeOfRep = this.tsda.TimeArray.get(rep).get(j);
							break;
						}
					}
					for (int j = 0; j < MaxTimeOfRep; j++) {
						time = this.tsda.maxTime * rep + j;
						if (this.tsda.TimeArray.get(rep).get(count) == j + 1 && this.tsda.Validity.get(rep).get(count)) {
							this.Calculator.setvalue(kurtosis_element, this.Kurtosis_mean_t[time][i]);
							this.Calculator.multAddDxyt(this.x_s[time], this.x_s[time], kurtosis_element, this.Txx_obs);
							this.Calculator.multAddDA(this.v_s[time], kurtosis_element, this.Txx_obs);
							count++;
							if (count == this.tsda.TimeArray.get(rep).size()) count = 0;
						}
					}
				}
				
				this.Calculator.copy(hws, this.Tyx[i]);
				this.Calculator.sub(hws, this.Skew_mean_x_sum[i]);
				
				/* Determine Next Active Set */
				this.activeSetH.clear();
				this.argMinH = 0;
				
				/* Test Existing Active Set H */
				this.Finished.clear();
				ArrayList<Integer> tempActiveSetH = new ArrayList<Integer>();
				ArrayList<Integer> nextNonActiveSetH = new ArrayList<Integer>(this.activeSetHList[i]);
				//this.calcArgMin(tempActiveSetH, nextNonActiveSetH, hws, i, 0, nextNonActiveSetH.size(), 2);
				
				/* Test Active Set H +-1 */
				this.Finished.clear();
				for (int remove = -1; remove < this.sm.getActiveRowElements(i, this.H).size(); remove++) {
					tempActiveSetH = new ArrayList<Integer>(this.sm.getActiveRowElements(i, this.H));
					if(remove > -1) tempActiveSetH.remove(remove);
					this.calcArgMin(tempActiveSetH, hws, i, 2, this.R[i], this.Rinv[i]);
					if(without_L1_fix_variables) break;
					nextNonActiveSetH = new ArrayList<Integer>();
					for (int j = 0; j < hws.length; j++) {
						if(!tempActiveSetH.contains(j) && !variable_cutting) nextNonActiveSetH.add(j);
					}
					this.calcArgMin(tempActiveSetH, nextNonActiveSetH, hws, i, 0, tempActiveSetH.size() + 1, 2, this.R[i], this.Rinv[i]);
				}
				this.activeSetHList[i] = new ArrayList<Integer>(this.activeSetH);

				/* Set Next Parameters */
				if(this.activeSetH.size()!=0) {
					for (int j = 0; j < activeSetH.size(); j++) {
						H_up_ws[i][activeSetH.get(j)] = this.nextH.get(j);
					}
				}
			}
			if(this.miss_update_eternal < 250) {
				this.vPa.setH(H_up_ws);
			}
		}
		
		/* 
		 * update R 
		 */
		for (int i = 0; i < this.tsda.elementNum; i++) {
			if(this.H == null) this.R_up_ws[i] = this.Txx_obs[i][i] - 2 * this.Tyx[i][i] + this.Tyy[i];
			else R_up_ws[i] = this.Dif_y_x_sq_sum[i] - 2 * this.Dif_y_x_delta_skew_mean_sum[i] + this.Kurtosis_H_Sigma_Htrans_sum[i];
			this.R_up_ws[i] += this.Delta[i] * this.Delta[i] * this.Kurtosis_skew_sq_mean_sum[i];
			this.R_up_ws[i] /= this.usedObservationalTimeNum;
			if(R_up_ws[i] < 0) 
				this.R_up_ws[i] = this.R[i];
			if(R_up_ws[i] > Math.abs(this.vSet.R_rI)) 
				this.R_up_ws[i] = Math.abs(this.vSet.R_rI);
		}
		if (this.vSet.R_rI<=0.0) {
			this.Calculator.setvalue(this.R_up_ws, this.Calculator.sumofVector(this.R_up_ws) / (double) this.R_up_ws.length);
		}
		if(this.miss_update_eternal < 100) {
			this.vPa.setR(this.R_up_ws);
		}
		
		/*
		 * update U
		 */
		if(this.vSet.Input){
			this.Calculator.copy(U_up_ws, this.sum_x);
			Calculator.sub(U_up_ws, Calculator.multAxReturn(F_up_ws, sum_xm));
			if(this.vSet.Drug) this.Calculator.sub(U_up_ws, this.Calculator.multAxReturn(G_up_ws, this.sum_zm));
			this.Calculator.rescale(U_up_ws, 1.0 / this.usedAllHiddenTime);
			this.vPa.setU(U_up_ws);
		}
		
		/* 
		 * update Mu
		 */
		this.Calculator.setvalue(this.x_up_ws, 0);
		double[] x_sum = new double[this.sysDim];
		for (int rep = 0; rep < this.tsda.repSize; rep++) {
			this.Calculator.add(this.x_up_ws[rep], this.x0_s[rep]);
			this.Calculator.add(x_sum, this.x0_s[rep]);			
		}
		this.Calculator.rescale(x_sum, (1.0/this.tsda.repSize));
		if (this.vSet.Mu0_Update > 0) {
			if (this.vSet.Mu0_Update == 1) {
				for (int rep = 0; rep < this.tsda.repSize; rep++) {
					for (int i = 0; i < this.sysDim; i++) {
						this.x_up_ws[rep][i] = x_sum[i];
					}
				}
			}
			if(this.miss_update_eternal < 350) {
				this.vPa.setx0(this.x_up_ws);
			}
		}
		
		/*
		 * update Q
		 */
		if(this.vSet.upDateQ <= 0){
			for (int i = 0; i < this.Q.length; i++) {
				this.Q_up_ws[i] = this.Txx[i][i] - 2 * this.Calculator.dotProduct(this.Txx_m[i], F_up_ws[i]) 
						+ this.Calculator.multxtAx(this.Txx_mm, F_up_ws[i]);
				if(this.vSet.Drug) {
					double[][] Txz_mmws = new double[this.Txz_mm.length][this.Txz_mm[0].length];
					this.Calculator.multAB(F_up_ws, this.Txz_mm, Txz_mmws);
					this.Q_up_ws[i] += 2 * this.Calculator.dotProduct(Txz_mmws[i], G_up_ws[i]) 
							- 2 * this.Calculator.dotProduct(this.Txz_m[i], G_up_ws[i]) 
							+ this.Calculator.multxtAx(this.Tzz_mm, G_up_ws[i]);
				} 
				if(this.vSet.Input) {
					this.Q_up_ws[i] += 2 * this.Calculator.dotProduct(F_up_ws[i], this.sum_xm) * U_up_ws[i] 
							- 2 * this.sum_x[i] * U_up_ws[i] + U_up_ws[i] * U_up_ws[i] * this.usedAllHiddenTime;
				}
				if(this.vSet.Input && this.vSet.Drug) {
					this.Q_up_ws[i] += 2 * this.Calculator.dotProduct(G_up_ws[i], this.sum_zm) * U_up_ws[i];
				}
			}
			this.Calculator.rescale(this.Q_up_ws, 1.0 / this.usedAllHiddenTime);
			for (int i = 0; i < this.Q_up_ws.length; i++) {
				if(this.Q_up_ws[i] > Math.abs(this.vSet.R_rI)) this.Q_up_ws[i] = Math.abs(this.vSet.R_rI);
				if(this.Q_up_ws[i] > Math.abs(this.vSet.upDateQ)) this.Q_up_ws[i] = Math.abs(this.vSet.upDateQ);
				if(this.Q_up_ws[i] < 0) this.Q_up_ws[i] = this.Q[i];
			}
			if (this.vSet.upDateQ < 0.0) {
				this.Calculator.setvalue(this.Q_up_ws, this.Calculator.sumofVector(this.Q_up_ws) / (double) this.Q_up_ws.length);
			}
			if(this.miss_update_eternal < 50) {
				this.vPa.setQ(this.Q_up_ws);
			}
		}
		
		this.Calculator.setvalue(this.Delta_up_ws, 0);
		for (int i = 0; i < this.Delta_up_ws.length; i++) {
			this.Delta_up_ws[i] = this.Dif_y_x_skew_mean_sum[i] / this.Kurtosis_skew_sq_mean_sum[i];
		}
		if(this.miss_update_eternal < 50) {
			this.vPa.setDelta(this.Delta_up_ws);
		}
	}
}