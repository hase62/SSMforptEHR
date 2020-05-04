package Hasegawa.TimeSeries;

import java.util.ArrayList;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.matrix.Matrix;
import Hasegawa.stat.simpleMath;
import RandomGenerator.Sfmt;

public class Inference {
	
	protected TimeSeriesDataArray tsda;
	protected Matrix Calculator;
	protected Sfmt sfmt;
	protected simpleMath sm = new simpleMath();
	protected int sysDim;
	
	/* Observation */
	protected double[][] y;
	protected double[] Tyy;
	
	/* Parameters */
	protected double[][] x0_;
	protected double[][] v0;
	protected double[][] A;
	protected double[][] F;
	protected double[][] H;
	protected double[][] G;
	protected double[] D;
	protected double[] U;
	protected double[] Q;
	protected double[] R;
	protected double[] I;
	
	/* Parameters for Monte-Carlo */
	protected int NumberOfParticle;
	protected int[] point1;
	protected int[] point2;
	protected double[][] Effect1;
	protected double[][] Effect2;
	
	/* Auxiliary Parameters */
	protected double[][] Ftrans;
	protected double[] Rinv;
	protected double[][] RinvMatrix;
	protected double[] Qinv;
	protected double[] sqR;
	protected double[] sqQ;
	protected double LogRDet;
	protected double LogQDet;
	protected double[][] Jtrans;
	protected double[][] HtransRinvH;
	protected double[][] RinvH;

	/* Likelihood */
	protected double previousLogLikelihood = (-1) * Double.POSITIVE_INFINITY;
	protected double currentLogLikelihood = (-1) * Double.POSITIVE_INFINITY;
	protected double Criterion = Double.POSITIVE_INFINITY;
	public double[] predictionError;
	public double[] predictionAbility;
	volatile protected double[] En_logLikelihood;
	volatile protected double[] En_temporal_logLikelihood;
	protected int usedObservationalTimeNum;
	protected int usedAllHiddenTime;
	
	/* Auxiliary Indicators */
	protected Boolean convergence = false;

	/* Workspace */
	protected double[][] sm_ws;
	protected double[][] sm_ws2;
	protected double[] sm_ws_v;
	protected double[][] v_pInv_ws;
	
	/* Simulation Values */
	protected double[][] x0_s;
	protected double[][][] v0_s;
	protected double[][] x_p;
	protected double[][] x_f;
	protected double[][] x_s;
	protected double[][] y_p;
	protected double[][][] v_p;
	protected double[][][] v_f;
	protected double[][][] v_s;
	protected double[][][] covariance;
	protected double[][][] J;
	protected double[][][] J0;
	
	/* Ensembles of Simulation Values */
	volatile protected double[] En_xs;
	volatile protected double[] En_xp;
	volatile protected double[] En_xp_no_noise;
	volatile protected double[] En_xf;
	volatile protected double[] En_x0_;
	
	protected double[] wi_x0;
	volatile protected double[] storedResidual;
	
	/* Smoothed Expectations */
	protected double[][] Txx;
	protected double[][] Txx_obs;
	protected double[][] Tyx;
	//protected double[][] Txy;
	protected double[][] Txx_m;
	protected double[][] Txx_mm;
	protected double[][] Txz_mm;
	protected double[][] Txz_mm_trans;
	protected double[][] Txz_m;
	protected double[] sum_x;
	protected double[] sum_xm;
	
	/* Inverse Matrix*/
	volatile protected double[][][] Vxx_inv;
	volatile protected double[][][] Vxxxx_inv;
	
	/* Extended Expectation Only Active Set */
	protected double[][] Txxxx_mmmm;
	protected double[][] Txxx_mm;
	protected double[][] Txxx_mmm;
	protected double[] vec;
	protected double[][] Txxz_mmm;
	
	/* Drug */
	protected double[][] z0_;
	protected double[][] z;
	protected double[][] Tzz_mm;
	protected double[][] Tzz_mmInv;
	protected double[] sum_zm;
	
	/* Correlations */
	protected double[][] Cor;
	protected double[][] lagCor;
	
	/* Test Prediction Ability */
	protected boolean testPredictionAbility = false;
	
	/* initialize Commons */
	protected void initilizeCommons(boolean DEGRADATION, boolean DRUG, boolean INPUT, boolean H){
		
		this.previousLogLikelihood = (-1) * Double.POSITIVE_INFINITY;
		this.currentLogLikelihood = (-1) * Double.POSITIVE_INFINITY;
		this.Criterion = Double.POSITIVE_INFINITY;
		this.convergence = false;
		this.predictionError = new double[this.sysDim];
		this.predictionAbility = new double[this.sysDim];
		
		/*
		 * Observation
		 */
		this.y = new double[this.tsda.observationalTimeNum][this.tsda.elementNum];
		this.Calculator.copy(this.y, this.tsda.ObservationData);		
		
		this.Tyy = new double[this.tsda.elementNum];
		this.usedObservationalTimeNum = 0;
		for (int r = 0; r < this.tsda.Validity.size(); r++) {
			for (int i = 0; i < this.tsda.Validity.get(r).size(); i++) {
				if(this.tsda.Validity.get(r).get(i)) {
					this.usedObservationalTimeNum++;
					int ytime = this.tsda.repIDSum.get(r) + i;
					for (int j = 0; j < this.tsda.elementNum; j++) {
						this.Tyy[j]+=this.y[ytime][j]*this.y[ytime][j];
					}
				}
			}
		}
		
		/*
		 * Parameters
		 */
		this.x0_ = new double[this.tsda.repSize][this.sysDim];
		this.v0 = new double[this.sysDim][this.sysDim];
	    this.A=new double[this.sysDim][this.sysDim];
	    this.F=new double[this.sysDim][this.sysDim];
	    if(H) {
	    	this.H = new double[this.tsda.elementNum][this.sysDim];
			this.RinvH = new double[this.tsda.elementNum][this.sysDim];
			this.HtransRinvH = new double[this.sysDim][this.sysDim];
	    }
		if(DRUG) this.G = new double[this.sysDim][tsda.drug0_[0].length];
	    if(DEGRADATION) this.D = new double[this.sysDim];
	    if(INPUT) this.U=new double[this.sysDim];
		this.Q=new double[this.sysDim];
	    this.R=new double[this.tsda.elementNum];
	    this.I = new double[this.sysDim];
		
	    /*
		 * Auxiliary Parameters
		 */
		this.Ftrans = new double[this.sysDim][this.sysDim];
	    this.Qinv = new double[this.sysDim];
	    this.Rinv = new double[this.tsda.elementNum];
	    this.RinvMatrix = new double[this.tsda.elementNum][this.tsda.elementNum];
		this.LogRDet = 0;
		this.LogQDet = 0;
		
		/* Workspace */
		this.sm_ws = new double[this.sysDim][this.sysDim];
		this.sm_ws2 = new double[this.sysDim][this.sysDim];
		this.sm_ws_v = new double[this.sysDim];
		this.v_pInv_ws = new double[this.sysDim][this.sysDim];
		
		/*
		 * Simulation Values
		 */
		this.x0_s = new double[this.tsda.repSize][this.sysDim];
		this.v0_s = new double[this.tsda.repSize][this.sysDim][this.sysDim];
		this.x_p = new double[this.tsda.allTime][this.sysDim];
		this.x_f = new double[this.tsda.allTime][this.sysDim];
		this.x_s = new double[this.tsda.allTime][this.sysDim];
		this.y_p = new double[this.tsda.allTime][this.tsda.elementNum];
		this.v_p = new double[this.tsda.allTime][this.sysDim][this.sysDim];
		this.v_f = new double[this.tsda.allTime][this.sysDim][this.sysDim];
		this.v_s = new double[this.tsda.allTime][this.sysDim][this.sysDim];
		this.covariance = new double[this.tsda.allTime][this.sysDim][this.sysDim];
		this.J = new double[this.tsda.allTime][this.sysDim][this.sysDim];
		this.J0 = new double[this.tsda.repSize][this.sysDim][this.sysDim];
		this.Jtrans = new double[this.sysDim][this.sysDim];
	}
	
	/*
	 * initialize Drug
	 */
	protected void initilizeDrug(boolean INPUT){
		/*
		 * Drug
		 */
		this.z = this.Calculator.copy_generate(this.tsda.drugMulRepSize);
		this.z0_ = this.Calculator.copy_generate(this.tsda.drug0_);
		this.Tzz_mm = new double[this.tsda.drugMulRepSize[0].length][this.tsda.drugMulRepSize[0].length];
		this.Tzz_mmInv = new double[tsda.drugMulRepSize[0].length][this.tsda.drugMulRepSize[0].length];
		for (int rep = 0; rep < tsda.repSize; rep++) {
			this.Calculator.multAddxyt(this.z0_[rep], this.z0_[rep], this.Tzz_mm);
			for (int t = 0; t < this.tsda.maxTime - 1; t++) {
				this.Calculator.multAddxyt(this.z[t + rep * this.tsda.maxTime], 
						this.z[t + rep * this.tsda.maxTime], this.Tzz_mm);
			}
		}
		this.Calculator.symmetricInverse(this.Tzz_mm, this.Tzz_mmInv);
		this.Calculator.changesymmetric(this.Tzz_mmInv);
		
		if(INPUT){
			this.sum_zm = new double[this.tsda.drugMulRepSize[0].length];
			for (int rep = 0; rep < this.tsda.repSize; rep++) {
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
					if (j == 0) {
						this.Calculator.add(this.sum_zm, this.z0_[rep]);
					} else {
						this.Calculator.add(this.sum_zm, this.z[time-1]);
					}
				}
			}
		}
	}
	
	/*
	 * initialize for updating Expectation of Smoothers
	 */
	protected void initializeExp(boolean DRUG, boolean INPUT) {
		this.Txx_obs = new double[this.sysDim][this.sysDim];
		this.Txx = new double[this.sysDim][this.sysDim];
		this.Tyx = new double[this.tsda.elementNum][this.sysDim];
		this.Txx_m = new double[this.sysDim][this.sysDim];
		this.Txx_mm = new double[this.sysDim][this.sysDim];
		if(DRUG){
			this.Txz_m = new double[this.sysDim][this.tsda.drugMulRepSize[0].length];
			this.Txz_mm = new double[this.sysDim][this.tsda.drugMulRepSize[0].length];
			this.Txz_mm_trans = new double[this.tsda.drugMulRepSize[0].length][this.sysDim];
		}
		if(INPUT){
			this.sum_x = new double[this.sysDim];
			this.sum_xm = new double[sysDim];
		}
	}
	
	protected void initializeEn(double[][] eff1, double eff2[][], double var_mu) {	
		this.En_logLikelihood = new double[this.NumberOfParticle];
		this.En_temporal_logLikelihood = new double[this.NumberOfParticle * this.tsda.observationalTimeNum];
    		this.y = new double[this.tsda.observationalTimeNum][this.tsda.elementNum];
		this.Calculator.copy(this.y, this.tsda.ObservationData);
		this.Tyy = new double[this.tsda.elementNum];
		this.usedObservationalTimeNum = 0;
		for (int r = 0; r < this.tsda.Validity.size(); r++) {
			for (int i = 0; i < this.tsda.Validity.get(r).size(); i++) {
				if(this.tsda.Validity.get(r).get(i)) {
					this.usedObservationalTimeNum++;
					int ytime = this.tsda.repIDSum.get(r) + i;
					for (int j = 0; j < this.tsda.elementNum; j++) {
						this.Tyy[j] += this.y[ytime][j] * this.y[ytime][j];
					}
				}
			}
		}
		
		this.x0_ = new double[this.tsda.repSize][this.sysDim];
		this.D = new double[this.sysDim];
	    if(this.sysDim != this.tsda.elementNum) {
	    	this.H = new double[this.tsda.elementNum][this.sysDim];
			this.RinvH = new double[this.tsda.elementNum][this.sysDim];
			this.HtransRinvH = new double[this.sysDim][this.sysDim];
	    }
		this.Q = new double[this.sysDim];
		this.Qinv = new double[this.sysDim];
		this.sqQ = new double[this.sysDim];
		this.R = new double[this.tsda.elementNum];
		this.Rinv = new double[this.tsda.elementNum];
		this.RinvMatrix = new double[this.tsda.elementNum][this.tsda.elementNum];
		this.sqR = new double[this.tsda.elementNum];
		this.U = new double[this.sysDim];
		this.I = new double[this.sysDim];
		if(this.tsda.givenG != null) this.G = new double[this.tsda.givenG.length][this.tsda.givenG[0].length];
		
		this.convergence=true;
		this.Effect1 = this.Calculator.copy_generate(eff1);
		this.Effect2 = this.Calculator.copy_generate(eff2);
		
		this.En_xs = new double[this.NumberOfParticle * this.sysDim * this.tsda.allTime];
		this.En_xp = new double[this.NumberOfParticle * this.sysDim * this.tsda.allTime];
		this.En_xp_no_noise = new double[this.NumberOfParticle * this.sysDim * this.tsda.allTime];
		this.En_xf = new double[this.NumberOfParticle * this.sysDim * this.tsda.allTime];
		
		this.En_x0_ = new double[this.NumberOfParticle * this.sysDim * this.tsda.repSize];
		
		this.wi_x0 = new double[this.sysDim];
		this.Calculator.setvalue(this.wi_x0, Math.sqrt(var_mu));
		this.storedResidual = new double[this.NumberOfParticle * this.sysDim * this.tsda.observationalTimeNum];
		this.Vxx_inv = new double[this.tsda.observationalTimeNum][this.sysDim][this.sysDim];
		this.Vxxxx_inv = new double[this.tsda.observationalTimeNum][this.sysDim][this.sysDim];
		
		/* Count the number of updating parameters */
		this.point1 = new int[this.sysDim+1];
		this.point2 = new int[this.sysDim+1];
		
    	/* Get the number of parameters */
    	int end = 0;
    	for (int i = 0; i < this.sysDim; i++) {
    		if(end == this.Effect1.length) {
    			this.point1[i] = end;
    		} else if(this.Effect1[end][0] == i) {
				for (int c2 = end + 1; c2 < this.Effect1.length + 1; c2++) {
					if(c2 == this.Effect1.length){
						end = this.Effect1.length;
						this.point1[i + 1] = this.Effect1.length;
						break;
					}
					if(this.Effect1[c2][0] != i) {
						end = c2;
						this.point1[i + 1] = end;
						break;
					} 
				} 
			} else {
	    			this.point1[i + 1] = end;
    		}
		}
    	for (int i = 0; i < this.sysDim; i++) {
    		int num = this.point1[i + 1] - this.point1[i];
			if(num > 1) {
				num = (num * (num - 1)) / 2;
				this.point2[i + 1] = this.point2[i] + num;
			} else this.point2[i + 1] = this.point2[i];
		}
    	this.point1[this.sysDim] = this.Effect1.length;
    	this.point2[this.sysDim] = this.Effect2.length;
	}
	
	/*
	 * initialize for updating Expectation of Smoothers More Than Second order
	 */
	protected void initializeExpExtended(boolean DRUG) {
		this.Txxxx_mmmm = new double[this.Effect2.length][this.Effect2.length];
		this.Txxx_mm = new double[this.sysDim][this.Effect2.length];
		this.Txxx_mmm = new double[this.sysDim][this.Effect2.length];
		this.vec = new double[this.Effect2.length];
		if(DRUG){
			this.Txxz_mmm = new double[this.Effect2.length][this.G[0].length];
		}
		this.Cor = new double[this.sysDim][this.sysDim];
		this.lagCor = new double[this.sysDim][this.sysDim];
	}
	
	/*
	 * Kalman Filter (Prediction and Filter)
	 */
	protected void KalmanFilter(boolean DRUG, boolean INPUT) {
		for (int rep = 0; rep < this.tsda.repSize; rep++) {
			int count = 0;
			for (int j = 0; j < this.tsda.maxTime; j++) {

				/* time = {replicate(1, 2, ...) - 1} * maxTime + j = time point */
				int time = this.tsda.maxTime * rep + j;
				
				/* Prediction */
				Prediction(time, rep, j == 0, DRUG, INPUT);
				
				/* Filter */
				count = Filter(time, this.tsda.TimeArray.get(rep).get(count) == j + 1, 
						this.tsda.Validity.get(rep).get(count) || this.testPredictionAbility, rep, count);
				
				if (this.tsda.TimeArray.get(rep).size() - 1 < count) {
					count = 0;
					break;
				}
			}
		}
	}

	protected void Prediction(final int time, final int rep, boolean zero, boolean DRUG, boolean INPUT) {
		if(zero) {
			this.Calculator.multAx(this.F, this.x0_[rep], this.x_p[time]);
			if(DRUG) this.Calculator.multAddAx(this.G, this.z0_[rep], this.x_p[time]);
			if(INPUT) this.Calculator.add(this.x_p[time], this.U);
			this.Calculator.multABCt(this.F, this.v0, this.F, this.v_p[time]);			
		} else {
			this.Calculator.multAx(this.F, this.x_f[time - 1], this.x_p[time]);
			if(DRUG) this.Calculator.multAddAx(this.G, this.z[time-1], this.x_p[time]);
			if(INPUT) this.Calculator.add(this.x_p[time], this.U);
			this.Calculator.multABCt(this.F, this.v_f[time - 1], this.F, this.v_p[time]);
		}
		this.Calculator.add(this.v_p[time], this.Q);
	}	
	
	/**
	 *  override
	 */
	protected int Filter(int time, boolean observed, boolean validity, int rep_num, int count){
		System.out.println("Filtering is not over-riden!!");
		System.exit(0);
		return count;
	}
	
	protected void KalmanSmoother() {
		Calculate_J();
		Smoother();
	}

	protected void Calculate_J() {
		for (int i = 0; i < this.tsda.repSize; i++) {
			int MaxTimeOfRep = 0;
			for (int j = this.tsda.TimeArray.get(i).size() - 1; j >= 0 ; j--) {
				if(this.tsda.Validity.get(i).get(j)) {
					MaxTimeOfRep = this.tsda.TimeArray.get(i).get(j);
					break;
				}
			}
			
			for (int j = 0; j < MaxTimeOfRep; j++) {
				int time = i * this.tsda.maxTime + j;
				
				this.Calculator.symmetricInverse(this.v_p[time], this.v_pInv_ws);
				this.Calculator.changesymmetric(this.v_pInv_ws);
				
				if (j == 0)	this.Calculator.multABC(this.v0, this.Ftrans, this.v_pInv_ws, this.J0[i]);
				else this.Calculator.multABC(this.v_f[time - 1], this.Ftrans, this.v_pInv_ws, this.J[time - 1]);
			}
		}
	}

	protected void Smoother() {
		for (int i = 0; i < this.tsda.repSize; i++) {
			int MaxTimeOfRep = 0;
			for (int j = this.tsda.TimeArray.get(i).size() - 1; j >= 0 ; j--) {
				if(this.tsda.Validity.get(i).get(j)) {
					MaxTimeOfRep = this.tsda.TimeArray.get(i).get(j);
					break;
				}
			}
			for (int j = 1; j <= MaxTimeOfRep; j++) {
				//time = (i + 1) * this.tsda.maxTime - j; // 20180129 improved to address different repIDs
				int time = i * this.tsda.maxTime - j + MaxTimeOfRep;
				if (j == MaxTimeOfRep) {
					/* x_0 */ 
					this.Calculator.copy(this.sm_ws_v, this.x_s[time]);
					this.Calculator.sub(this.sm_ws_v, this.x_p[time]);
					this.Calculator.multAddAx(this.J0[i], this.sm_ws_v, this.x0_s[i]);
					/* v_0 */
					this.Calculator.copy(this.sm_ws, this.v_s[time]);
					this.Calculator.sub(this.sm_ws, this.v_p[time]);
					this.Calculator.multAddABCt(this.J0[i], this.sm_ws, this.J0[i], this.v0_s[i]);
					/* cov */
					this.Calculator.setvalue(this.covariance[time], 0);
					this.Calculator.transpose(this.J0[i], this.Jtrans);
					this.Calculator.multAB(this.v_s[time], this.Jtrans, this.covariance[time]);
				} else {
					/* x_s */ 
					this.Calculator.copy(this.sm_ws_v, this.x_s[time]);
					this.Calculator.sub(this.sm_ws_v, this.x_p[time]);
					this.Calculator.multAddAx(this.J[time - 1], this.sm_ws_v, this.x_s[time - 1]);
					/* v_s */
					this.Calculator.copy(this.sm_ws, this.v_s[time]);
					this.Calculator.sub(this.sm_ws, this.v_p[time]);
					this.Calculator.multAddABCt(this.J[time - 1], this.sm_ws, this.J[time - 1],	this.v_s[time - 1]);
					/* cov_s */
					this.Calculator.setvalue(this.covariance[time], 0);
					this.Calculator.transpose(this.J[time - 1], this.Jtrans);
					this.Calculator.multAB(this.v_s[time], this.Jtrans, this.covariance[time]);
				}
			}
		}
	}
	
	/*
	 * when NaN or OverEstimation happens, stop the process and return false
	 */
	protected void MissResult(){
		System.err.println("Wrong Estimation: NaN");
		this.convergence = true;
		this.currentLogLikelihood = (-1) * Double.POSITIVE_INFINITY;
		this.previousLogLikelihood = (-1) * Double.POSITIVE_INFINITY;
		this.Criterion = Double.POSITIVE_INFINITY;
	}
	
	protected void calculateConditionalExpectation(boolean DRUG, boolean INPUT) {
		this.usedAllHiddenTime = 0;
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
				this.usedAllHiddenTime++;
				time = this.tsda.maxTime * rep + j;
				if (this.tsda.TimeArray.get(rep).get(count) == j + 1 && this.tsda.Validity.get(rep).get(count)) {
					ytime = this.tsda.repIDSum.get(rep) + count;
					this.Calculator.multAddxyt(this.x_s[time], this.x_s[time], this.Txx_obs);
					this.Calculator.add(this.Txx_obs, this.v_s[time]);
					this.Calculator.multAddxyt(this.y[ytime], this.x_s[time], this.Tyx);
					count++;
					if (count == this.tsda.TimeArray.get(rep).size()) count = 0;
				}
				if (j == 0) {
					this.Calculator.multAddxyt(this.x_s[time], this.x_s[time], this.Txx);
					this.Calculator.add(this.Txx, this.v_s[time]);
					this.Calculator.multAddxyt(this.x_s[time], this.x0_s[rep], this.Txx_m);
					this.Calculator.add(this.Txx_m, this.covariance[time]);
					this.Calculator.multAddxyt(this.x0_s[rep], this.x0_s[rep], this.Txx_mm);
					this.Calculator.add(this.Txx_mm, this.v0_s[rep]);
					
					if(DRUG){
						this.Calculator.multAddxyt(this.x_s[time], this.z0_[rep], Txz_m);
						this.Calculator.multAddxyt(this.x0_s[rep], this.z0_[rep], Txz_mm);
					}
					if(INPUT){
						this.Calculator.add(this.sum_x, this.x_s[time]);
						this.Calculator.add(this.sum_xm, this.x0_s[rep]);
					}
				} else {
					this.Calculator.multAddxyt(this.x_s[time], this.x_s[time], this.Txx);
					this.Calculator.add(this.Txx, this.v_s[time]);
					this.Calculator.multAddxyt(this.x_s[time], this.x_s[time - 1], this.Txx_m);
					this.Calculator.add(this.Txx_m, this.covariance[time]);
					this.Calculator.multAddxyt(this.x_s[time - 1], this.x_s[time - 1], this.Txx_mm);
					this.Calculator.add(this.Txx_mm, this.v_s[time - 1]);
					
					if(DRUG){
						this.Calculator.multAddxyt(this.x_s[time], this.z[time-1], this.Txz_m);
						this.Calculator.multAddxyt(this.x_s[time-1], this.z[time-1], this.Txz_mm);
					}
					if(INPUT){
						this.Calculator.add(this.sum_x, this.x_s[time]);
						this.Calculator.add(this.sum_xm, this.x_s[time-1]);
					}
				}
			}
		}
	}
	
	/*
	 * Get Criterion(Cri : 0.0=BIC, 1.0:AIC) from current LogLikelihood
	 */
	protected double calculateBIC(boolean useF, int upMu, boolean INPUT, boolean DRUG, double upQ, boolean upD, 
			boolean R_rI, boolean upH, boolean skew, boolean kurtosis, double Cri) {
		
		/* R */
		int c = 1;
		if (!R_rI) c += this.tsda.elementNum - 1;
		
		/* Mu */
		if (upMu==1) c += this.sysDim;
		if (upMu==2) c += this.sysDim * this.tsda.repSize;
		
		/* F */
		if(useF) c += this.sm.getActiveCount(this.F);
		else c += this.sm.getActiveCount(this.A);
		
		/* H */
		if(upH) c += this.sm.getActiveCount(this.H);
		
		/* Q */
		if(upQ==0) c += this.sysDim;
		else if(upQ < 0) c++;
		
		/* G */
		if(DRUG) c += this.sm.getActiveCount(this.G);
		
		/* Input */
		if(INPUT) c += this.sysDim;
		
		/* Degradation */
		if(upD) c += this.sysDim;

		/* Skew and Kurtosis */
		if(skew) c += this.tsda.elementNum;
		if(kurtosis) c += this.tsda.elementNum;
		
		/* Criterion */
		if(Cri==1.0) return this.Criterion = -2 * this.currentLogLikelihood + (2 * this.usedObservationalTimeNum * c) / (this.usedObservationalTimeNum - c + 1); 
		else return this.Criterion = -2 * this.currentLogLikelihood + (c * Math.log(this.usedObservationalTimeNum));
	}
	
	protected double calculateBIC_Ensemble(double logLikelihood, boolean Drug, int mu0update){
		int length = this.Effect1.length + this.Effect2.length + this.D.length + this.U.length + this.Q.length + this.R.length;
		if(mu0update == 1) length += this.sysDim;
		if(mu0update == 2) length += this.sysDim * this.tsda.repSize;
		if(Drug) length += G.length * G[0].length;
		double currentBic = -2 * logLikelihood + (length * Math.log(this.usedObservationalTimeNum));
		return currentBic;
	}
	
	public double getLoglikelihood() {
		assert (this.currentLogLikelihood != 0);
		return this.currentLogLikelihood;
	}
	
	/*
	 * If H matrix is null, set logLikelihood(this.Calculator.makeMatrixFromDiag(this.Rinv), null, null)
	 */
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
					this.Calculator.multxtD(ws, this.Rinv);
					logLikelihood -= 0.5 *this.Calculator.dotProduct(ws, yx);
					
					// (iv)
					if(RiH == null) this.Calculator.multAx(this.Rinv, yx, ws2);
					else this.Calculator.multAtx(RiH, yx, ws2);
					logLikelihood += 0.5 * this.Calculator.multxtAx(this.v_f[time], ws2);
					
					/**
					 * |I +UVt| = |I + UtV|
					 *  |(H*Vp*Ht + R)| = |(H*Vp*Ht*Ri + I) * R| = |(H*Vp*Ht*Ri + I)| * |R|
					 *  	= |(Vp*Ht*Ri*H + I)| * |R|
					 */
					// LU = I + Vp*HtRiH
					this.Calculator.copy(LU, this.I);
					/* must set HtRiH*/
					this.Calculator.multAddAtB(this.v_p[time], HtRiH, LU);
					this.Calculator.LUDecomposition(LU, L_d, U_d);
					logLikelihood -= 0.5 * this.tsda.elementNum * Math.log(2 * Math.PI);
					logLikelihood -= 0.5 * this.Calculator.diagLogDeterminant(U_d);
					logLikelihood -= 0.5 * this.LogRDet;

					count++;
					if (count == this.tsda.TimeArray.get(i).size())	count = 0;
				}
			}
		}
		return logLikelihood;
	}
	
	protected void logLikelihood_Ensemble(int ytime, int time) {
		double sum = 0;
		double log2Pi = this.sysDim * Math.log(2 * Math.PI);
		for (int n = 0; n < this.NumberOfParticle; n++) {
			double l = 0;
			double[] yx = new double[this.sysDim];
			for (int i = 0; i < yx.length; i++) {
				yx[i] =  this.y[ytime][i] - this.En_xp[this.Calculator.getNumberOfEn(time, i, n, this.sysDim, this.NumberOfParticle)];
			}
			
			l -= this.Calculator.multxtDx(this.Rinv, yx);
			l -= log2Pi;
			l -= this.LogRDet;
			l *= 0.5;
			this.En_temporal_logLikelihood[this.Calculator.getNumberOfEn(ytime, 0, n, 1, this.NumberOfParticle)] = l;
			this.En_logLikelihood[n] += l;
			sum += Math.exp(l);
		}
		this.currentLogLikelihood += Math.log(sum+Double.MIN_VALUE);
	}
	
	protected void storePosterior(double[] x, int dim, ArrayList<double[]> List){	
		for (int n = 0; n < 1000; n++) {
			double[] temp = new double[dim];
			
			for (int i = 0; i < dim; i++) {
				temp[i] = x[this.Calculator.getNumberOfEn(0, i, (int)(this.sfmt.NextUnif() * this.NumberOfParticle * 0.99999), dim, this.NumberOfParticle)];
			}
			
			List.add(temp);
		}
	}
}
