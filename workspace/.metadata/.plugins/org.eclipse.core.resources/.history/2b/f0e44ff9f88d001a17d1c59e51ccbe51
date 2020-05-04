package Hasegawa.TimeSeries.Linear.VARSSM;

import java.util.ArrayList;
import java.util.Collections;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.TimeSeries.Inference;
import Hasegawa.TimeSeries.Linear.VARSSM.vsSetting;
import Hasegawa.matrix.Matrix;
import Hasegawa.stat.simpleMath;
import RandomGenerator.Sfmt;

public class vsInference extends Inference{

	public vsSetting vSet;
	public vsStorage vsSto;
	public vsParameter vPa;
	protected vsParameter tempPrePa;
	
	protected double[] variance_xs;
	protected double[] variance_z;
	
	protected double[][] WeightA;
	protected double[][] WeightG;
	
	/* Workspace */
	protected double[][] v_fRinv;
	protected double[] y_ws;
	protected double[] x_rescale_ws;
	protected double[] z_rescale_ws;
	protected double[][] A_up_ws;
	protected double[][] F_up_ws;
	protected double[][] G_up_ws;
	protected double[][] H_up_ws;
	protected double[] D_up_ws;
	protected double[] R_up_ws;
	protected double[] U_up_ws;
	protected double[][] x_up_ws;
	protected double[][] v_up_ws;
	protected double[] Q_up_ws;
	protected double[][] yx_ws;

	protected int iteration;
	public int updatingRow;
	public int[] calculationOrder;
	protected boolean continueThisCalculation = false;
	
	protected double th_EdgeNum = 0;
	public double L1_UpRate = 0.995;
	private int L1atStart = 0;
	protected int miss_update_eternal = 0;
	
	protected ArrayList<Integer>[] activeSetAList;
	protected ArrayList<Integer>[] activeSetGList;
	protected ArrayList<Integer>[] activeSetHList;
	
	protected ArrayList<Integer> Finished = new ArrayList<Integer>();
	protected double argMinA = Double.MAX_VALUE;
	protected double argMinG = Double.MAX_VALUE;
	protected double argMinH = Double.MAX_VALUE;
	protected ArrayList<Integer> activeSetA = new ArrayList<Integer>();
	protected ArrayList<Integer> activeSetG = new ArrayList<Integer>();
	protected ArrayList<Integer> activeSetH = new ArrayList<Integer>();
	protected ArrayList<Double> nextA = new ArrayList<Double>();
	protected ArrayList<Double> nextG = new ArrayList<Double>();
	protected ArrayList<Double> nextH = new ArrayList<Double>();
	
	protected double[] Delta;
	protected double[] Delta_up_ws;
	protected double[] Nu;
	
	protected double[][] Skew_mean_t;
	protected double[][] Skew_sq_mean_t;
	protected double[][] Delta_skew_mean_t;

	protected double[][] Kurtosis_mean_t;
	
	protected double[][] Kurtosis_mean_Rinv_t;
	protected double[][][] Htrans_Kurtosis_mean_Rinv_H_t;
	protected double[][][] Kurtosis_mean_Rinv_H_t;
	
	protected double[] Kurtosis_skew_sq_mean_sum;
	protected double[] Dif_y_x_skew_mean_sum;
	protected double[] Dif_y_x_abs_max_sign;
	protected double[] Dif_y_x_sq_sum;
	protected double[][] Skew_mean_x_sum;
	protected double[] Kurtosis_H_Sigma_Htrans_sum;
	
	public vsInference(final int SYD, vsSetting SET, TimeSeriesDataArray TSDA, Matrix MP, Sfmt Sf, vsStorage vStor) {
		this.sysDim = SYD;
		this.vSet = SET;
		this.tsda = TSDA;
		this.Calculator = MP;
		this.sfmt = Sf;
		this.vsSto = vStor;
	}
	
	@SuppressWarnings("unchecked")
	protected void initialize() {
		
		this.initilizeCommons(this.vSet.Degradation > 0.0, this.vSet.Drug, this.vSet.Input, !(this.sysDim == this.tsda.elementNum));
		if(this.vSet.Drug) this.initilizeDrug(this.vSet.Input);
		
		/* Set Unique Parameters */
		this.updatingRow = 0;
		this.continueThisCalculation=false;
		this.convergence = false;
		
		/* Get Variance of Drug */
		this.variance_xs = new double[this.sysDim];
		this.v_fRinv = new double[this.sysDim][this.sysDim];
		//Recently Changed from this.sysDim to this.tsda.elementNum
		this.y_ws = new double[this.tsda.elementNum];
		this.x_rescale_ws = new double[this.sysDim];
		
		/* Workspace */
		this.A_up_ws = new double[this.sysDim][this.sysDim];
		this.F_up_ws = new double[this.sysDim][this.sysDim];
		this.G_up_ws = null;
		if(this.vSet.Drug) this.G_up_ws = new double[this.sysDim][this.tsda.drugMulRepSize[0].length];
		this.H_up_ws = null;
		if(this.H != null) this.H_up_ws = new double[this.tsda.elementNum][this.sysDim];
		this.D_up_ws = null;
		if(this.vSet.Degradation > 0) this.D_up_ws = new double[this.sysDim];
		this.R_up_ws = new double[this.tsda.elementNum];
		this.U_up_ws=null;
		if(this.vSet.Input) this.U_up_ws = new double[this.sysDim];
		this.x_up_ws = new double[this.tsda.repSize][this.sysDim];
		this.v_up_ws = new double[this.sysDim][this.sysDim];
		this.Q_up_ws = new double[this.Q.length];
		
		this.Delta = new double[this.tsda.elementNum];
		this.Delta_up_ws = new double[this.tsda.elementNum];
		this.Nu = new double[this.tsda.elementNum];
		
		this.Skew_mean_t = new double[this.tsda.allTime][this.tsda.elementNum];
		this.Skew_sq_mean_t = new double[this.tsda.allTime][this.tsda.elementNum];
		this.Delta_skew_mean_t = new double[this.tsda.allTime][this.tsda.elementNum];
		
		this.Kurtosis_mean_t = new double[this.tsda.allTime][this.tsda.elementNum];

		this.Kurtosis_mean_Rinv_t = new double[this.tsda.allTime][this.tsda.elementNum];
		this.Htrans_Kurtosis_mean_Rinv_H_t = new double[this.tsda.allTime][this.sysDim][this.sysDim];
		this.Kurtosis_mean_Rinv_H_t = new double[this.tsda.allTime][this.tsda.elementNum][this.sysDim];

		this.Kurtosis_skew_sq_mean_sum = new double[this.tsda.elementNum];
		this.Dif_y_x_skew_mean_sum = new double[this.tsda.elementNum];
		this.Dif_y_x_sq_sum = new double[this.tsda.elementNum];
		this.Dif_y_x_abs_max_sign = new double[this.tsda.elementNum];
		this.Skew_mean_x_sum = new double[this.tsda.elementNum][this.sysDim];
		this.Kurtosis_H_Sigma_Htrans_sum = new double[this.tsda.elementNum];
		
		this.yx_ws = new double[this.tsda.elementNum][this.sysDim];		

		if(this.vSet.Drug){
			this.z_rescale_ws = new double[this.z[0].length];
			this.variance_z= new double[this.z[0].length];
			
			/* Mean of Drug*/
			double[] mean_z = new double[this.z[0].length];
			double[] val_z= new double[this.z[0].length];
			double[] zz = new double[this.z[0].length];
			for (int rep = 0; rep < this.tsda.repSize; rep++) {
				this.Calculator.add(mean_z, this.z0_[rep]);
				this.Calculator.multxx(this.z0_[rep], zz);
				this.Calculator.add(val_z, zz);
			}
			for (int j = 0; j < this.z.length; j++) {
				if((j + 1) % this.tsda.maxTime == 0) continue; 
				this.Calculator.add(mean_z, this.z[j]);
				this.Calculator.multxx(this.z[j], zz);
				this.Calculator.add(val_z, zz);
			}
			this.Calculator.rescale(mean_z, 1.0/(this.z.length));
			this.Calculator.rescale(val_z, 1.0/(this.z.length));
			
			/* Variance of Drug*/
			for (int i = 0; i < this.variance_z.length; i++) {
				this.variance_z[i] = val_z[i] - mean_z[i] * mean_z[i];
			}
		}

		/* Get Threshold */
		if(this.th_EdgeNum == 0){
			this.Edge();
		}
		
		/* Set Calculation Order for each row */
		this.calculationOrder = this.sm.getIthOrder(this.tsda.elementNum);
		
		/* Set Weight for Weighted LASSO */
		this.WeightA = new double[this.A.length][this.A[0].length];
		this.Calculator.setvalue(this.WeightA, 1);
		if(vSet.Drug) {
			this.WeightG = new double[this.G.length][this.G[0].length];
			this.Calculator.setvalue(this.WeightG, 1);
		}
		
		this.activeSetAList = new ArrayList[this.sysDim];
		if(this.vSet.Drug) this.activeSetGList = new ArrayList[this.sysDim];
		if(this.H != null) {
			this.activeSetHList = new ArrayList[this.tsda.elementNum];
			for (int i = 0; i < this.H.length; i++) {
				this.activeSetHList[i] = new ArrayList<Integer>();
			}
		}
		for (int i = 0; i < this.sysDim; i++) {
			this.activeSetAList[i] = new ArrayList<Integer>();
			if(this.vSet.Drug) this.activeSetGList[i] = new ArrayList<Integer>();
			if(this.H != null) this.activeSetHList[i] = new ArrayList<Integer>();
		}
	}
	
	protected void Edge(){
		this.th_EdgeNum = 1 + Math.log(this.sysDim) / Math.log(2) * 1.5;
		if(this.vSet.Drug) this.th_EdgeNum += Math.log(this.G[0].length) / Math.log(2);
		//if(this.H != null) this.th_EdgeNum = this.sysDim;
		if(this.vSet.Degradation==0) this.th_EdgeNum++;
		if(this.sm.getSumOfColLength(this.F, this.G) - 1 <= this.th_EdgeNum) {
			this.th_EdgeNum = this.sm.getSumOfColLength(this.F, this.G) - 2;
		}
	}
	
	/**
	 * As Preparation of Performing The Main Algorithm, Trying to Get Initial States
	 * 1. Get Initial Profiles Using Full Matrix 
	 * 2. Get Initial Profiles Using Sparse Matrix but Un-fixing only 1 row 
	 * @param max_loop: Max loop num. of separate update
	 * @param max_edge: Max num. of regulations
	 */
	protected void preparationPhase(int max_loop, int max_edge) {
		if(this.sysDim - 1 < max_edge) max_edge = sysDim - 1;
		
		/* 1. Get Initial Profiles Using Full Matrix */
		System.out.println("Get Initial Profiles Using Full Matrix");
		this.initialUpdate(1.0e-20, this.vSet.Degradation != 0.0 && false);
		
		/* 2. Get Initial Profiles Using Sparse Matrix but Un-fixing only 1 row */
		System.out.println("Get Initial Profiles Using Sparse Matrix but Un-fixing only 1 row");
		final long timeAtStart = System.currentTimeMillis();
		this.currentLogLikelihood = (-1) * Double.MAX_VALUE;
		for (int i = 0; (double)i / this.sysDim < max_loop	&& System.currentTimeMillis() - timeAtStart < this.vSet.timer; i++) {
			
			/* Consider Updating Order for Each Loop */
			if(i % this.sysDim == 0) {
				this.calculationOrder = this.getOrder();
			}
			
			/* Sequential Parameters Update */
			this.previousLogLikelihood = this.currentLogLikelihood;
			if(i < this.sysDim) this.vsSto.Criterion = Double.MAX_VALUE;
			this.separateUpDate(this.calculationOrder[i % this.sysDim], max_edge, i / this.sysDim);
		}
	}
	
	/* 1. Get Initial Profiles Using Full Matrix */
	protected void initialUpdate(double L1, boolean oneLoop){
		this.vPa = new vsParameter(this.sysDim, this.sfmt, this.vSet, this.tsda, this.Calculator, 0.0);
		this.tempPrePa = new vsParameter(this.sysDim, this.sfmt, this.vSet, this.tsda, this.Calculator, 0.0);
		this.initialize();
		
		this.currentLogLikelihood = -1.0 * Double.MAX_VALUE;
		this.Calculator.setvalue(this.vPa.L1, L1);
		if(this.H != null) this.Calculator.setvalue(this.vPa.L1h, L1);
		double temporalWeight = this.vSet.WeightLASSO;
		this.vSet.WeightLASSO = 0.0;
		this.convergence = false;
		while(Math.abs(this.currentLogLikelihood - this.previousLogLikelihood) > this.vSet.Condition_of_Convergence) {
			this.getParameters();
			this.KalmanFilter(this.vSet.Drug, this.vSet.Input);
			this.KalmanSmoother();
			this.previousLogLikelihood = this.currentLogLikelihood;
			if(this.H == null) this.currentLogLikelihood = this.logLikelihood(this.RinvMatrix, null, null);
			else this.currentLogLikelihood = this.logLikelihood(this.HtransRinvH, this.H, this.RinvH);
			if(this.vSet.Print_Progress) System.err.println("Initial Loop: " + this.currentLogLikelihood);
			if(this.currentLogLikelihood < this.previousLogLikelihood) this.miss_update_eternal++;
			this.initializeExp(this.vSet.Drug, this.vSet.Input);
			if(oneLoop) break;
			this.calculateConditionalExpectation(this.vSet.Drug, this.vSet.Input);
			if(this.Calculator.checkNaN(this.Txx) 
					|| Double.isNaN(this.currentLogLikelihood) 
						|| this.Calculator.checkAbsValues(this.A, 0.9999)){
				break;
			}
			this.Update(false, true, false, false, false);
			this.iteration++;
		}
		this.iteration = 0;
		this.vSet.WeightLASSO = temporalWeight;
	}
	
	/* 2. Get Initial Profiles Using Sparse Matrix but Un-fixing only 1 row */
	private void separateUpDate(int row, int maxEdge, int num){
		vsSetting clonevSET = new vsSetting();
		clonevSET.copy(this.vSet);
		clonevSET.Drug = true;
		clonevSET.Spacom = true;
		clonevSET.Separate = true;
		if(num == 0) clonevSET.WeightLASSO = 0.0;
		
		TimeSeriesDataArray cloneTSDA = new TimeSeriesDataArray(this.Calculator);
		cloneTSDA.copy(this.tsda);
		
		/* Drug (Fixing values are treated as drug) */
		int gLength = 0;
		if(this.vSet.Drug) gLength = this.tsda.drugMulRepSize[0].length;
		
		/* Targeting Nodes and Drug */
		ArrayList<Integer> Targeting = new ArrayList<Integer>();
		Targeting.add(row);
		Collections.sort(Targeting);
		ArrayList<Integer> nonTargeting = new ArrayList<Integer>();
		for (int i = 0; i < this.sysDim; i++) {
			if(Targeting.indexOf(i) == -1) nonTargeting.add(i);
		}
		int tLength = Targeting.size();
		int ntLength = nonTargeting.size();
		if(this.H == null) cloneTSDA.elementNum = tLength;
		
		/* Get Matrix (Since All elements were set by the initial process, then store 0 at first) */
		cloneTSDA.givenF = new double[tLength][tLength];
		for (int i = 0; i < cloneTSDA.givenF.length; i++) {
			for (int j = 0; j < cloneTSDA.givenF[0].length; j++) {
				cloneTSDA.givenF[i][j] = this.tsda.givenF[Targeting.get(i)][Targeting.get(j)];
			}
		}
		cloneTSDA.givenG = new double[tLength][gLength + ntLength]; 
		for (int i = 0; i < cloneTSDA.givenG.length; i++) {
			for (int j = 0; j < gLength; j++) {
				cloneTSDA.givenG[i][j] = this.tsda.givenG[Targeting.get(i)][j];
			}
			for (int j = 0; j < ntLength; j++) {
				cloneTSDA.givenG[i][gLength + j] = this.tsda.givenF[Targeting.get(i)][nonTargeting.get(j)];
			}
		}
		
		/* Set Smoothed x to Drug */
		cloneTSDA.drugMulRepSize = new double[this.x_s.length][gLength + ntLength];
		cloneTSDA.drug0_ = new double[this.tsda.repSize][gLength + ntLength];
		
		/* Set Drug Profiles */
		for (int i = 0; i < gLength; i++) {
			for (int j = 0; j < cloneTSDA.drugMulRepSize.length; j++) {
				cloneTSDA.drugMulRepSize[j][i] = this.tsda.drugMulRepSize[j][i];
			}
		}
		for (int rep = 0; rep < this.tsda.repSize; rep++) {
			for (int i = 0; i < gLength; i++) {
				cloneTSDA.drug0_[rep][i] = this.tsda.drug0_[rep][i];
			}
		}
		
		/* Set non-Targeting Profiles as Drug */
		for (int i = 0; i < ntLength; i++) {
			for (int j = 0; j < cloneTSDA.drugMulRepSize.length; j++) {
				cloneTSDA.drugMulRepSize[j][gLength + i] = this.x_s[j][nonTargeting.get(i)];
			}
		}
		for (int rep = 0; rep < this.tsda.repSize; rep++) {
			for (int i = 0; i < gLength; i++) {
				cloneTSDA.drug0_[rep][gLength + i] = this.x0_[rep][nonTargeting.get(i)];
			}
		}
		
		/* Set Observation Data, but Original Observation will be Used  When H is not null*/
		if(this.H == null){
			cloneTSDA.ObservationData = new double[this.tsda.ObservationData.length][tLength];
			for (int i = 0; i < this.tsda.ObservationData.length; i++) {
				for (int j = 0; j < tLength; j++) {
					cloneTSDA.ObservationData[i][j] = this.tsda.ObservationData[i][Targeting.get(j)];
				}
			}	
		}
		
		/* Set Parameter for Inference */
		vsStorage clonevSTO = new vsStorage();
		clonevSTO.initializeVSStorage(tLength, clonevSET, cloneTSDA, this.Calculator, 0.0, this.sfmt);
		vsInference vINF = new vsInference(tLength, clonevSET, cloneTSDA, this.Calculator, this.sfmt, clonevSTO);
		vINF.th_EdgeNum = maxEdge;
		vINF.vPa = new vsParameter(tLength, this.sfmt, clonevSET, cloneTSDA, this.Calculator, 0.0);
		vINF.tempPrePa = new vsParameter(tLength, this.sfmt, clonevSET, cloneTSDA, this.Calculator, 0.0);
		vINF.initialize();
		for (int i = 0; i < Targeting.size(); i++) {
			if(this.vSet.Degradation != 0.0) vINF.vPa.setD(this.D[Targeting.get(i)], i);
			if(this.vSet.upDateQ==0) vINF.vPa.setQ(this.Q[Targeting.get(i)], i);
			if(this.H == null) vINF.vPa.setR(this.R[Targeting.get(i)], i);
			else vINF.vPa.setR(this.R);
			if(vSet.Input) vINF.vPa.setU(this.U[Targeting.get(i)], i);
		}
			
		/* Inference */
		boolean recursive = true;
		while(recursive){
			recursive = vINF.mainRun();
			this.setNewParameter(Targeting, nonTargeting, vINF.vPa, row, tLength, gLength, ntLength);
			this.KalmanFilter(this.vSet.Drug, this.vSet.Input);
			if(this.H == null) this.currentLogLikelihood = this.logLikelihood(this.RinvMatrix, null, null);
			else this.currentLogLikelihood = this.logLikelihood(this.HtransRinvH, this.H, this.RinvH);
			this.Criterion = this.calculateBIC(false, this.vSet.Mu0_Update, this.vSet.Input, this.vSet.Drug, this.vSet.upDateQ, 
					this.vSet.Degradation > 0.0, this.vSet.R_rI <= 0.0, this.H != null, this.Delta != null, this.Nu != null, this.vSet.Criterion);
			this.currentLogLikelihood = this.penalization(this.currentLogLikelihood, this.A, this.G, this.H, this.vPa.L1, this.vPa.L1h);
			if(this.vsSto.Criterion > this.Criterion || true){
				this.storeCurrentSettings(this.vsSto);
			}
		}
		this.recallPreviousSettings(this.vsSto);
		this.setParameters();
		this.setActiveSets(row);
		this.KalmanFilter(this.vSet.Drug, this.vSet.Input);
		this.KalmanSmoother();
		this.calculateConditionalExpectation(this.vSet.Drug, this.vSet.Input);
		if(this.H == null) this.currentLogLikelihood = this.logLikelihood(this.RinvMatrix, null, null);
		else this.currentLogLikelihood = this.logLikelihood(this.HtransRinvH, this.H, this.RinvH);
		this.Criterion = this.calculateBIC(false, this.vSet.Mu0_Update, this.vSet.Input, this.vSet.Drug, this.vSet.upDateQ, 
				this.vSet.Degradation > 0.0, this.vSet.R_rI <= 0.0, this.H != null, this.Delta != null, this.Nu != null, this.vSet.Criterion);
		this.currentLogLikelihood = this.penalization(this.currentLogLikelihood, this.A, this.G, this.H, this.vPa.L1, this.vPa.L1h);
		this.storeCurrentSettings(this.vsSto);
		if(Double.isNaN(this.x0_s[0][0])) {
			this.Calculator.copy(this.x0_s, this.x0_);
			this.Calculator.copy(this.x_s, this.x_f);
		}
	}
		
	protected double penalization(double cLL, double[][] cA, double[][] cG, double[][] cH, double[] L1, double[] L1h){
		double pLL = cLL;
		for (int i = 0; i < cA.length; i++) {
			for (int j = 0; j < cA[0].length; j++) {
				pLL -= Math.abs(cA[i][j]) * L1[i] / (this.WeightA[i][j] * this.tsda.givenF[i][j]);
			}
			if(cG != null){
				for (int j = 0; j < cG[0].length; j++) {
					pLL -= Math.abs(cG[i][j]) * L1[i] / (this.WeightG[i][j] * this.tsda.givenG[i][j]);
				}
			}
		}
		for (int i = 0; i < this.tsda.elementNum; i++) {
			if(cH != null){
				for (int j = 0; j < cH[0].length; j++) {
					pLL -= Math.abs(cH[i][j]) * L1h[i];
				}
			}
		}
		return pLL;
	}
	
	private void setNewParameter(ArrayList<Integer> Targeting, ArrayList<Integer> nonTargeting, vsParameter vPa, 
								int row, int tLength, int gLength, int ntLength){
		
		/* Update parameter */
		for (int i = 0; i < vPa.getF().length; i++) {
			for (int j = 0; j < vPa.getF()[0].length; j++) {
				this.A[Targeting.get(i)][Targeting.get(j)] = vPa.getA()[i][j];
				this.F[Targeting.get(i)][Targeting.get(j)] = vPa.getF()[i][j];
			}
		}
		for (int i = 0; i < vPa.getG().length; i++) {
			for (int j = 0; j < gLength; j++) {
				this.G[Targeting.get(i)][j] = vPa.getG()[i][j];
			}
			for (int j = 0; j < ntLength; j++) {
				this.A[Targeting.get(i)][nonTargeting.get(j)] = vPa.getG()[i][j + gLength];
				this.F[Targeting.get(i)][nonTargeting.get(j)] = vPa.getG()[i][j + gLength];
			}
		}
		
		/* Update Other Parameters */
		for (int i = 0; i < tLength; i++) {
			if(this.vSet.Degradation != 0.0) this.D[Targeting.get(i)] = vPa.getD()[i]; 
			if(this.vSet.upDateQ == 0) this.Q[Targeting.get(i)] = vPa.getQ()[i];
			if(this.H == null) this.R[Targeting.get(i)] = vPa.getR()[i];
			else {
				for (int row_H = 0; row_H < this.tsda.elementNum; row_H++) {
					this.H[row_H][Targeting.get(i)] = vPa.getH()[row_H][i];
				}
			}
			if(this.vSet.Input) this.U[Targeting.get(i)] = vPa.getU()[i];
			if(this.vSet.Mu0_Update > 1) {
				for (int rep = 0; rep < this.tsda.repSize; rep++) {
					this.x0_[rep][Targeting.get(i)] = vPa.getx0()[rep][i];					
				}
			}
			//this.v0[Targeting.get(i)][Targeting.get(i)] = vPa.getv0()[i][i];
			this.vPa.L1[Targeting.get(i)] = vPa.L1[i];
		}
		if(this.H != null) {
			this.Calculator.copy(this.R, vPa.getR());
		}
	}
	
	public boolean mainRun() {
		
		/* This is used to avoid a specific problem in using HGC */
		double base_L1_UpDateRate = 0.999;
		if(this.continueThisCalculation) {
			this.continueThisCalculation = false;
		} else {
			this.previousLogLikelihood = (-1) * Double.MAX_VALUE;
			this.L1_UpRate = base_L1_UpDateRate;
		}
		
		boolean previousConvergence = false;
		while(this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]) < this.th_EdgeNum){
			
			/* Minimum L1 */
			if(this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] < 0.1) break;
			
			/* Main Process */
			this.getParameters();
			this.KalmanFilter(this.vSet.Drug, this.vSet.Input);
			this.KalmanSmoother();
			this.initializeExp(this.vSet.Drug, this.vSet.Input);
			this.calculateConditionalExpectation(this.vSet.Drug, this.vSet.Input);
			
			/* Calculate logLikelihood, used for BIC, and penalized logLikelihood */
			if(this.H == null) this.currentLogLikelihood = this.logLikelihood(this.RinvMatrix, null, null);
			else this.currentLogLikelihood = this.logLikelihood(this.HtransRinvH, this.H, this.RinvH);
			this.Criterion = this.calculateBIC(false, this.vSet.Mu0_Update, this.vSet.Input, this.vSet.Drug, this.vSet.upDateQ, 
					this.vSet.Degradation > 0.0, this.vSet.R_rI <= 0.0, this.H != null, this.Delta != null, this.Nu != null, this.vSet.Criterion);
			this.currentLogLikelihood = this.penalization(this.currentLogLikelihood, this.A, this.G, this.H, this.vPa.L1, this.vPa.L1h);
			
			/* Print Progress */
			if (this.vSet.Print_Progress) {
				System.err.println("Updating Num: " + this.updatingRow + " (Row " + this.calculationOrder[this.updatingRow] + ") - itr "+ this.iteration);
				System.err.println("logLikelihood= " + this.currentLogLikelihood);
				System.err.println("Cri= " + this.Criterion);
				System.err.println("L1= " + this.vPa.L1[this.calculationOrder[(int)this.updatingRow]]);
				System.err.println("Num. of Updating Edges= " + this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]) + " of " + this.th_EdgeNum);
				System.err.println("Num. of Total Edges = " + (this.sm.getActiveCount(this.A) + this.sm.getActiveCount(this.G)));
				System.err.println();
			}
		
			/* Check Convergence and Then Update L1 */
			if(this.currentLogLikelihood - this.previousLogLikelihood < this.vSet.Condition_of_Convergence 
					&& this.currentLogLikelihood >= this.previousLogLikelihood){
				
				/* Decrease L1_UpRate When Continuously Converged, Otherwise Set base_L1_UpDateRate */
				if(previousConvergence) this.L1_UpRate *= 0.99;
				else this.L1_UpRate = base_L1_UpDateRate;
				if(this.L1_UpRate < 0.95) this.L1_UpRate = 0.95;
				
				/* Update L1 According to L1_UpRate When Having Edges but 0.99 When Having No Edge */
				if(this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]) > 0) {
					this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] *= this.L1_UpRate;
				} else {
					this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] *= 0.99;
				}
				previousConvergence = true;
				this.convergence = true;
			} else {
				previousConvergence = false;
				this.convergence = false;
			}			
			
			/* Store Best BIC */
			if( this.Criterion < this.vsSto.Criterion &&  
					this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] > 0.05 &&
						this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]) < this.th_EdgeNum){
				this.storeCurrentSettings(this.vsSto);
			}
			
			/* At Current < Previous : errors in calculation */
			if(this.previousLogLikelihood > this.currentLogLikelihood && this.iteration > 5){
				this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] *= 0.995;
			}
			
			/* Update Parameters */	
			this.previousLogLikelihood = this.currentLogLikelihood;
			this.tempPrePa.setParameters(this.vPa, this.Calculator);
			this.Update(false, false, false, this.iteration==0, false);
			if(this.checkNaN(Double.MAX_VALUE)){
				this.MissResult();
				break;
			}
			
			/* On the HGC-Super Computer */
			this.iteration++;
			if((this.convergence && this.vSet.Separate) || (this.iteration % 30 == 29 && !this.vSet.Separate && this.vSet.Spacom)){
				this.continueThisCalculation = true;
				break;
			}
		}
		return this.continueThisCalculation;
	}

	protected boolean checkNaN(double A_th){
		/* Check NaN */
		if(this.Calculator.checkNaN(this.Txx)) {
			System.out.println("Txx");
			for (int i = 0; i < this.Txx.length; i++) {
				for (int j = 0; j < this.Txx[0].length; j++) {
					System.out.print(Txx[i][j] + " ");
				}
				System.out.println();
			}
			return true;
		}
		if( Double.isNaN(this.currentLogLikelihood)) {
			System.out.println("logLikelihood " + this.currentLogLikelihood);
			return true;			
		}
		if(this.Calculator.checkAbsValues(this.A, A_th)){
			System.out.println("A");
			for (int i = 0; i < this.A.length; i++) {
				for (int j = 0; j < this.A[0].length; j++) {
					System.out.print(this.A[i][j] + " ");
				}
				System.out.println();
			}
			return true;
		}
		if(this.Calculator.checkValues(this.Q, 1.0e5, -0.01)){
			System.out.println("Q");
			for (int i = 0; i < this.Q.length; i++) {
				System.out.print(this.Q[i] + " ");
			}
			System.out.println();
			return true;
		}
		if(this.Calculator.checkValues(this.R, 1.0e5, -0.01)){
			System.out.println("R");
			for (int i = 0; i < this.R.length; i++) {
				System.out.print(this.R[i] + " ");
			}
			System.out.println();
			return true;
		}
		return false;
	}
	
	protected void initializeUpdatingRow(boolean initial){
		this.iteration = 0;
		if(initial){
			this.vPa.L1[this.calculationOrder[this.updatingRow]] = 1.0e30;
			this.Calculator.setvalue(this.A[this.calculationOrder[this.updatingRow]], 0);
			this.Calculator.setvalue(this.F[this.calculationOrder[this.updatingRow]], 0);
			if(this.vSet.Degradation != 0.0) this.F[this.calculationOrder[this.updatingRow]][this.calculationOrder[this.updatingRow]] = 1 - this.D[this.calculationOrder[this.updatingRow]]; 
			if(vSet.Drug) this.Calculator.setvalue(this.G[this.calculationOrder[this.updatingRow]], 0);
		} else {
			this.L1atStart = 1;
			this.vPa.L1[this.calculationOrder[this.updatingRow]] *= 5;
		}
	}
	
	public void setActiveSetsAll(){
		for (int i = 0; i < this.sysDim; i++) {
			this.activeSetAList[this.calculationOrder[i]] = this.sm.getActiveRowElements(this.calculationOrder[i], this.vsSto.vPa.getA());
			if(this.vSet.Drug) this.activeSetGList[this.calculationOrder[i]] = this.sm.getActiveRowElements(this.calculationOrder[i], this.vsSto.vPa.getG());
		}
		for (int i = 0; i < this.tsda.elementNum; i++) {
			if(this.H != null) this.activeSetHList[i] = this.sm.getActiveRowElements(i, this.vsSto.vPa.getH());
		}
	}

	public void setActiveSets(int row){
		this.activeSetAList[row] = this.sm.getActiveRowElements(row, this.vsSto.vPa.getA());
		if(this.vSet.Drug) this.activeSetGList[row] = this.sm.getActiveRowElements(row, this.vsSto.vPa.getG());
	}

	public void setActiveSetsH(int row){
		if(this.H != null) this.activeSetHList[row] = this.sm.getActiveRowElements(row, this.vsSto.vPa.getH());
	}
	
	protected void getParameters() {
		this.LogRDet = this.vPa.getParameters(this.x0_, this.v0, this.A, this.F, this.Ftrans, this.H, this.G, 
				this.RinvH, this.HtransRinvH, this.D, this.Q, this.Qinv, this.R, this.Rinv, this.RinvMatrix, this.U, this.I, 
				null, null, null, null, null, null, null, null, null, null, null);
		for (int rep = 0; rep < this.tsda.repSize; rep++) {
			this.Calculator.copy(this.x0_s[rep], this.x0_[rep]); 
			this.Calculator.copy(this.v0_s[rep], this.v0); 
		}
		this.Calculator.setvalue(this.predictionError, 0);
		this.Calculator.setvalue(this.predictionAbility, 0);
	    	if(this.H != null) {
				this.Calculator.multAB(this.Rinv, this.H, this.RinvH);
				this.Calculator.multAtB(this.H, this.RinvH, this.HtransRinvH);
	    	}
	}
	
	protected int Filter(int time, boolean observed, boolean validity, int rep_num, int count){
		if(this.H != null) this.Calculator.multAx(this.H, this.x_p[time], this.y_p[time]);
		else this.Calculator.copy(this.y_p[time], this.x_p[time]);
		if (observed && validity) {
			if (this.H == null){
				this.Calculator.inversionTheorem(this.v_p[time], this.Rinv, this.v_f[time]);
				this.Calculator.changesymmetric(this.v_f[time]);
				this.Calculator.copy(this.x_f[time], this.x_p[time]);
				this.Calculator.multAB(this.v_f[time], this.Rinv, this.v_fRinv);
				this.Calculator.sub(this.y[this.tsda.repIDSum.get(rep_num) + count], this.x_p[time], this.y_ws);
				this.Calculator.multAddAx(this.v_fRinv, this.y_ws, this.x_f[time]);
			} else {
				double[][] IKH = new double[this.sysDim][this.sysDim];
				this.Calculator.symmetricInverse(this.v_p[time], IKH);
				this.Calculator.changesymmetric(IKH);
				this.Calculator.add(IKH, this.HtransRinvH);
				this.Calculator.symmetricInverse(IKH, this.v_f[time]);
				this.Calculator.changesymmetric(this.v_f[time]);

				this.Calculator.copy(this.x_f[time], this.x_p[time]);
				this.Calculator.multAddABtx(this.v_f[time], this.RinvH, this.y[this.tsda.repIDSum.get(rep_num) + count], this.x_f[time]);
				this.Calculator.multSubABx(this.v_f[time], this.HtransRinvH, this.x_p[time], this.x_f[time]);
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
				int time = this.tsda.maxTime * i + j;
				
				if(this.H == null){
					this.Calculator.copy(this.v_pInv_ws, this.Qinv);
					if(j==0) this.Calculator.inversionTheorem(this.v0, this.F, this.Qinv, this.sm_ws);
					else this.Calculator.inversionTheorem(this.v_f[time-1], this.F, this.Qinv, this.sm_ws);
					
					this.Calculator.changesymmetric(this.sm_ws);
					this.Calculator.multAB(this.Qinv, this.F,  this.sm_ws2);
					this.Calculator.multSubABCt(this.sm_ws2, this.sm_ws, this.sm_ws2, this.v_pInv_ws);
				} else {
					this.Calculator.symmetricInverse(this.v_p[time], this.v_pInv_ws);
					this.Calculator.changesymmetric(this.v_pInv_ws);
				}
				
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
				
				int time = i * this.tsda.maxTime - j + MaxTimeOfRep;
				if (j == MaxTimeOfRep) {
					//x_s0
					this.Calculator.copy(this.sm_ws_v, this.x_s[time]);
					this.Calculator.sub(this.sm_ws_v, this.x_p[time]);
					this.Calculator.multAddAx(this.J0[i], this.sm_ws_v, this.x0_s[i]);
					//v_s0
					this.Calculator.copy(this.sm_ws, this.v_s[time]);
					this.Calculator.sub(this.sm_ws, this.v_p[time]);
					this.Calculator.multAddABCt(this.J0[i], this.sm_ws, this.J0[i], this.v0_s[i]);
					//cov
					this.Calculator.setvalue(this.covariance[time], 0);
					this.Calculator.transpose(this.J0[i], this.Jtrans);
					this.Calculator.multAB(this.v_s[time], this.Jtrans, this.covariance[time]);
				} else {
					//x_s
					this.Calculator.copy(this.sm_ws_v, this.x_s[time]);
					this.Calculator.sub(this.sm_ws_v, this.x_p[time]);
					this.Calculator.multAddAx(this.J[time - 1], this.sm_ws_v, this.x_s[time - 1]);
					//v_s
					this.Calculator.copy(this.sm_ws, this.v_s[time]);
					this.Calculator.sub(this.sm_ws, this.v_p[time]);
					this.Calculator.multAddABCt(this.J[time - 1], this.sm_ws, this.J[time - 1], this.v_s[time - 1]);
					//cov
					/*if (j==1000) {
						if(this.H == null){
							//cov_T,T-1
							this.Calculator.multAB(this.v_f[time], this.Rinv, this.v_fRinv);
							this.Calculator.copy(this.sm_ws, this.v_fRinv);
							this.Calculator.SubDA(this.I, this.sm_ws);
							this.Calculator.multABC(this.sm_ws, this.F, this.v_f[time-1], this.covariance[time]);
						} else {
							double[][] IKH = new double[this.sysDim][this.sysDim];
							this.Calculator.symmetricInverse(this.v_p[time], IKH);
							this.Calculator.changesymmetric(IKH);
							this.Calculator.add(IKH, this.HtransRinvH);
							this.Calculator.symmetricInverse(IKH);
							this.Calculator.multABC(IKH, this.F, this.v_f[time-1], this.covariance[time]);
						}
					} else {*/
						this.Calculator.setvalue(this.covariance[time], 0);
						this.Calculator.transpose(this.J[time-1], this.Jtrans);
						this.Calculator.multAB(this.v_s[time], this.Jtrans, this.covariance[time]);
					//} 20180129 Improved to estimate the parameter values
				}
			}
		}
	}
	
	public void Weight(){
		if(this.vSet.WeightLASSO <= 0) return;
		this.Calculator.setvalue(this.WeightA, 1.0);
		if(vSet.Drug){
			this.Calculator.setvalue(this.WeightG, 1.0);
		}
		
		if(this.vSet.WeightLASSO==3 || this.vSet.WeightLASSO==4) {
			/* Weight using Ridge */
			for (int i = 0; i < this.sysDim; i++) {
				double[] aws = new double[this.A.length];
				this.Calculator.copy(aws, this.Txx_m[i]);
				if(this.vSet.Drug)	this.Calculator.multSubAx(this.Txz_mm, this.G[i], aws);
				if(this.vSet.Input) {
					this.Calculator.copy(this.x_rescale_ws, this.sum_xm);
					this.Calculator.rescale(this.x_rescale_ws, this.U[i]);
					this.Calculator.sub(aws, this.x_rescale_ws);
				}
				if(this.vSet.Degradation != 0.0) {
					this.Calculator.copy(this.x_rescale_ws, this.Txx_mm[i]);
					this.Calculator.rescale(this.x_rescale_ws, 1-this.D[i]);
					this.Calculator.sub(aws, this.x_rescale_ws);
				}
				
				double[][] Inv = new double[this.A.length][this.A.length];
				double[][] Gamma = new double[this.A.length][this.A.length];
				if(this.vSet.WeightLASSO==3) this.Calculator.setDiagvalue(Gamma, 1.0e-5);
				else if(this.vSet.WeightLASSO==4) this.Calculator.setDiagvalue(Gamma, 1.0e5);
				this.Calculator.add(Gamma, this.Txx_mm);
				this.Calculator.symmetricInverse(Gamma, Inv);
	    		this.Calculator.multAx(Inv, aws, this.WeightA[i]);
				for (int j = 0; j < this.A[0].length; j++) {
					this.WeightA[i][j]=1/Math.log(Math.abs(this.WeightA[i][j])+Math.E);
				}
			}
			if(this.vSet.Drug){
				for (int i = 0; i < this.sysDim; i++) {
					double[] gws = new double[this.G[0].length];
					this.Calculator.copy(gws, this.Txz_m[i]);
					this.Calculator.transpose(this.Txz_mm, this.Txz_mm_trans);
					this.Calculator.multSubAx(this.Txz_mm_trans, this.F[i], gws);
					if(this.vSet.Input) {
						this.Calculator.copy(this.z_rescale_ws, this.sum_zm);
						this.Calculator.rescale(this.z_rescale_ws, this.U[i]);
						this.Calculator.sub(gws, this.z_rescale_ws);
					}
					double[][] InvG = new double[this.G[0].length][this.G[0].length];
					double[][] Gamma = new double[this.G[0].length][this.G[0].length];
					if(this.vSet.WeightLASSO==3) this.Calculator.setDiagvalue(Gamma, 1.0e-5);
					else if(this.vSet.WeightLASSO==4) this.Calculator.setDiagvalue(Gamma, 1.0e5);
					this.Calculator.add(Gamma, this.Tzz_mm);
					this.Calculator.symmetricInverse(Gamma, InvG);
		    		this.Calculator.multAx(InvG, gws, this.WeightG[i]);
					for (int j = 0; j < this.G[0].length; j++) {
						WeightG[i][j]=1/Math.log(Math.abs(WeightG[i][j])+Math.E);	
					}
				}
			}
			return;
		}
		
		/* Mean of Expression*/
		this.variance_xs= new double[this.x_s[0].length];
		double[] mean_x = new double[this.x_s[0].length];
		double[] val_x= new double[this.x_s[0].length];
		double[] xx = new double[this.x_s[0].length];
		for (int i = 0; i < this.x0_s.length; i++) {
			this.Calculator.add(mean_x, this.x0_s[i]);
			this.Calculator.multxx(this.x0_s[i], xx);
			this.Calculator.add(val_x, xx);
		}
		this.Calculator.rescale(mean_x, (double)this.tsda.repSize);
		this.Calculator.rescale(val_x, (double)this.tsda.repSize);
		
		for (int j = 0; j < this.x_s.length; j++) {
			if((j + 1) % this.tsda.maxTime == 0) continue; 
			this.Calculator.add(mean_x, this.x_s[j]);
			this.Calculator.multxx(this.x_s[j], xx);
			this.Calculator.add(val_x, xx);
		}
		this.Calculator.rescale(mean_x, 1.0/(this.x_s.length));
		this.Calculator.rescale(val_x, 1.0/(this.x_s.length));
		
		/* Variance of Expression*/
		for (int i = 0; i < this.variance_xs.length; i++) {
			this.variance_xs[i] = val_x[i] - mean_x[i] * mean_x[i];
		}
		
		/* Update Weight */
		for (int j = 0; j <this.variance_xs.length; j++) {
			for (int i = 0; i < WeightA.length; i++) {
				/* Each Node */
				if(this.vSet.WeightLASSO==1.0) this.WeightA[i][j] = 1.0 / Math.sqrt(this.variance_xs[j]);
				if(this.vSet.WeightLASSO==2.0) this.WeightA[i][j] = Math.log((1.0 / this.variance_xs[j]) + Math.E - 1);
			}
		}
		if(!this.vSet.Drug) return;
		for (int j = 0; j < this.variance_z.length; j++) {
			for (int i = 0; i < this.WeightG.length; i++) {
				/* Each Drug*/
				if(this.vSet.WeightLASSO==1.0) this.WeightG[i][j] = 1.0/Math.sqrt(this.variance_z[j]);
				if(this.vSet.WeightLASSO==2.0) this.WeightG[i][j] = Math.log((1.0 / this.variance_z[j]) + Math.E - 1);
				if(this.vSet.WeightLASSO==5.0) this.WeightG[i][j] = 1.0/(this.variance_z[j]/this.variance_z[new simpleMath().getMaxAbsValueIndex(this.variance_z)]);
				if(this.vSet.WeightLASSO==6.0) this.WeightG[i][j] = Math.log((1.0 / this.variance_z[j]) + Math.E - 1);
				if(this.vSet.WeightLASSO==7.0) this.WeightG[i][j] = 1.0/this.variance_z[j];
			}
		}
	}
	
	protected void Update(boolean variable_cutting, boolean allUpdate, boolean update_H, boolean initialLoop, 
			boolean without_L1_fix_variables) {
		
		/*
		 * update A, F and G
		 */
		this.Calculator.setvalue(this.A_up_ws, 0);
		this.Calculator.setvalue(this.F_up_ws, 0);
		if(this.vSet.Drug) this.Calculator.setvalue(this.G_up_ws, 0);
		double[] aws = new double[this.sysDim];
		double[] gws = null;
		if(this.vSet.Drug) gws = new double[this.tsda.drugMulRepSize[0].length];
		for (int i = 0; i < this.F.length; i++) {
			/* A (= F + diag(1-d)) */
			this.Calculator.copy(aws, this.Txx_m[i]);
			if(this.vSet.Drug)	this.Calculator.multSubAx(this.Txz_mm, this.G[i], aws);
			if(this.vSet.Input) {
				this.Calculator.copy(this.x_rescale_ws, this.sum_xm);
				this.Calculator.rescale(this.x_rescale_ws, this.U[i]);
				this.Calculator.sub(aws, this.x_rescale_ws);
			}
			if(this.vSet.Degradation != 0.0) {
				this.Calculator.copy(this.x_rescale_ws, this.Txx_mm[i]);
				this.Calculator.rescale(this.x_rescale_ws, 1 - this.D[i]);
				this.Calculator.sub(aws, this.x_rescale_ws);
			}
			if(this.vSet.Degradation != 0.0) aws[i]=0;
			
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
			
			if(allUpdate){
				ArrayList<Integer> tempActiveSetA = new ArrayList<Integer>();
				for (int j = 0; j < this.F[0].length; j++) {
					if(!(this.vSet.Degradation != 0.0 && i == j)) tempActiveSetA.add(j);
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
			} else if(this.calculationOrder[(int)this.updatingRow]==i && !variable_cutting) {
				if(this.iteration == 0 && this.L1atStart == 0) {
					this.vPa.L1[i] = Math.abs(this.getSelectedValueL1(aws, gws, i, this.L1atStart) * 1.0001) * this.Qinv[i];
			    }
				
				/* Test Existing Active Set A */
				this.Finished.clear();
				ArrayList<Integer> tempActiveSetA = new ArrayList<Integer>();
				ArrayList<Integer> nextNonActiveSetA = new ArrayList<Integer>(this.activeSetAList[i]);
				this.calcArgMin(tempActiveSetA, nextNonActiveSetA, aws, i, 0, nextNonActiveSetA.size(), 0, this.Q[i], this.Qinv[i]);
				
				/* Test Active Set A +-1 */
				this.Finished.clear();
				for (int remove = -1; remove < this.sm.getActiveRowElements(i, this.A).size(); remove++) {
					tempActiveSetA = new ArrayList<Integer>(this.sm.getActiveRowElements(i, this.A));
					if(remove > -1) tempActiveSetA.remove(remove);
					this.calcArgMin(tempActiveSetA, aws, i, 0, this.Q[i], this.Qinv[i]);
					if(without_L1_fix_variables) break;
					nextNonActiveSetA = new ArrayList<Integer>();
					for (int j = 0; j < aws.length; j++) {
						if(i == j && this.vSet.Degradation != 0.0) continue;
						if(!tempActiveSetA.contains(j)) nextNonActiveSetA.add(j);
					}
					this.calcArgMin(tempActiveSetA, nextNonActiveSetA, aws, i, 0, tempActiveSetA.size() + 1, 0, this.Q[i], this.Qinv[i]);
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
					this.calcArgMin(tempActiveSetG, nextNonActiveSetG, gws, i, 0, nextNonActiveSetG.size(), 1, this.Q[i], this.Qinv[i]);
					
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
					
					this.activeSetGList[i] = new ArrayList<Integer>(this.activeSetG);
				}
			} else {
				/* A */
				/* Test Existing Active Set A */
				this.Finished.clear();
				ArrayList<Integer> tempActiveSetA = new ArrayList<Integer>();
				ArrayList<Integer> nextNonActiveSetA = new ArrayList<Integer>(this.activeSetAList[i]);
				this.calcArgMin(tempActiveSetA, nextNonActiveSetA, aws, i, 0, nextNonActiveSetA.size(), 0, this.Q[i], this.Qinv[i]);
				
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
					this.calcArgMin(tempActiveSetG, nextNonActiveSetG, gws, i, 0, nextNonActiveSetG.size(), 1, this.Q[i], this.Qinv[i]);
					
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
				for (int j = 0; j < activeSetA.size(); j++) {
					A_up_ws[i][activeSetA.get(j)] = this.nextA.get(j);
				}
			}
			
			if(this.vSet.Drug && this.activeSetG.size()!=0){
				for (int j = 0; j < activeSetG.size(); j++) {
					G_up_ws[i][activeSetG.get(j)] = this.nextG.get(j);
				}
			}
		}
		this.vPa.setA(A_up_ws);
		if(this.vSet.Drug) this.vPa.setG(G_up_ws);
		
		/*
		 * update D
		 */
		if(this.vSet.Degradation != 0.0){
			for (int i = 0; i < this.D.length; i++) {
				D_up_ws[i] = this.Txx_m[i][i]-this.Calculator.dotProduct(A_up_ws[i], this.Txx_mm[i]);
				if(this.vSet.Drug) D_up_ws[i] -= this.Calculator.dotProduct(G_up_ws[i], this.Txz_mm[i]);
				if(this.vSet.Input) D_up_ws[i] -= this.U[i]*this.sum_xm[i];
				D_up_ws[i]=1-(D_up_ws[i]/this.Txx_mm[i][i]);
				if(D_up_ws[i] < this.vSet.Degradation) D_up_ws[i] = this.vSet.Degradation;
				if(D_up_ws[i] > 0.4) D_up_ws[i] = 0.4;
			}
			this.vPa.setD(D_up_ws);
		}
		
		/* 
		 * update F
		 */
		this.Calculator.copy(F_up_ws, A_up_ws);
		if(this.vSet.Degradation != 0.0){
			for (int i = 0; i < this.F.length; i++) {
				F_up_ws[i][i]+=(1-D_up_ws[i]);
			}
		}
		this.vPa.setF(F_up_ws);
		
		// H
		double[][] Htemp = null;
		if(this.H != null){
			double[][] Txx_obsInv = new double[this.sysDim][this.sysDim];
			this.Calculator.symmetricInverse(this.Txx_obs, Txx_obsInv);
			this.Calculator.changesymmetric(Txx_obsInv);

			Htemp = new double[this.tsda.elementNum][this.sysDim];
			this.Calculator.multAB(this.Tyx, Txx_obsInv, Htemp);
			this.vPa.setH(Htemp);
		}
		
		/* 
		 * update R 
		 */
		for (int i = 0; i < this.tsda.elementNum; i++) {
			if(this.H == null) R_up_ws[i] = this.Txx_obs[i][i] - 2 * this.Tyx[i][i] + this.Tyy[i];
			else R_up_ws[i] = this.Calculator.multxtAx(Txx_obs, Htemp[i]) - 2 * this.Calculator.dotProduct(Htemp[i], this.Tyx[i]) + this.Tyy[i];	
			R_up_ws[i] /= this.usedObservationalTimeNum;
			if(R_up_ws[i] < 0) 	R_up_ws[i] = this.R[i];
			if(R_up_ws[i] > Math.abs(this.vSet.R_rI)) R_up_ws[i] = Math.abs(this.vSet.R_rI);
		}
		if (this.vSet.R_rI<=0.0) {
			this.Calculator.setvalue(R_up_ws, this.Calculator.sumofVector(R_up_ws) / (double) R_up_ws.length);
		}
		this.vPa.setR(R_up_ws);

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
			this.vPa.setx0(this.x_up_ws);
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
				if(this.Q_up_ws[i] < 0) this.Q_up_ws[i] = this.Q[i];
			}
			if (this.vSet.upDateQ < 0.0) {
				this.Calculator.setvalue(this.Q_up_ws, this.Calculator.sumofVector(this.Q_up_ws) / (double) this.Q_up_ws.length);
			}
			this.vPa.setQ(this.Q_up_ws);	
		}
	}
	
	protected void calcArgMin(ArrayList<Integer> Active, ArrayList<Integer> NonActive, double[] ws, int u_row, double parent_min, int limitLength, int A0G1H2, double coef, double coef_inv){
		
		Collections.sort(Active);
		Collections.sort(NonActive);
		
		/* update */
		double[] subWS = new double[Active.size() + 1];
		double[] nextWS = new double[Active.size() + 1];
		double[] nextWS2 = new double[Active.size() + 1];
		double[] sign = new double[Active.size() + 1];
		double[][] T_mmInvActive = new double[Active.size() + 1][Active.size() + 1];
		double[][] T_mmInvActive_ws = new double[Active.size() + 1][Active.size() + 1];
		for (int n = 0; n < NonActive.size() ; n++) {
			/* Prepare Temporal Array For Next */
			ArrayList<Integer> nextActive = new ArrayList<Integer>(Active);
			ArrayList<Integer> nextNonActive = new ArrayList<Integer>(NonActive);
			nextActive.add(NonActive.get(n));
			nextNonActive.remove(n);
			Collections.sort(nextActive);
			
			/* Skip if the Next Set was already Evaluated, or Store Next Set */
			if(Finished.contains(this.sm.getEncodedIntegerList(nextActive))) {
				continue;
			} else {
				this.Finished.add(this.sm.getEncodedIntegerList(nextActive));
			}
			
			for (int j = 0; j < subWS.length; j++) {
				subWS[j] = ws[nextActive.get(j)];
			}
			this.Calculator.copy(nextWS, subWS);
			
			if(A0G1H2 == 0) {
				this.Calculator.getSubMatrix(this.Txx_mm, nextActive, T_mmInvActive);
				this.Calculator.symmetricInverseWithWS(T_mmInvActive, T_mmInvActive_ws);
			} else if(A0G1H2 ==1){
				this.Calculator.getSubMatrix(this.Tzz_mm, nextActive, T_mmInvActive);
				this.Calculator.symmetricInverseWithWS(T_mmInvActive, T_mmInvActive_ws);
			} else if(A0G1H2 ==2){
				this.Calculator.getSubMatrix(this.Txx, nextActive, T_mmInvActive);
				this.Calculator.symmetricInverseWithWS(T_mmInvActive, T_mmInvActive_ws);
			}
			
			for (int i = 0; i < nextActive.size(); i++) {
				if(A0G1H2 == 0) {
					if(this.A[u_row][nextActive.get(i)] > 0) sign[i] = 1;
					else if(this.A[u_row][nextActive.get(i)] < 0) sign[i] = -1;
					else sign[i] = nextWS[i]/Math.abs(nextWS[i]);
				} else if(A0G1H2 == 1) {
					if(this.G[u_row][nextActive.get(i)] > 0) sign[i] = 1;
					else if(this.G[u_row][nextActive.get(i)] < 0) sign[i] = -1;
					else sign[i] = nextWS[i]/Math.abs(nextWS[i]);
				} else if(A0G1H2 == 2) {
					if(this.H[u_row][nextActive.get(i)] > 0) sign[i] = 1;
					else if(this.H[u_row][nextActive.get(i)] < 0) sign[i] = -1;
					else sign[i] = nextWS[i]/Math.abs(nextWS[i]);
				}
			}
			
			/* Search When Changing at most two sign */
			boolean possible = false;
			int pointer_1 = -1;
			while(pointer_1 < nextActive.size()){
				/* Skip if High Dimension */
				if(Active.size() > 30 && pointer_1 > -1) break;
				
				/* Change Sign */
				if(pointer_1 > -1) sign[this.sm.getNthMinAbsValueIndex(subWS, pointer_1 + 1)] *= -1;
				
				int pointer_2 = -1;
				while(pointer_2 < nextActive.size()){
					/* At first, check original sign */
					if(pointer_2 <= pointer_1 && pointer_1 != -1) {
						pointer_2++;
						continue;
					}
					possible = true;
					
					/* Change Sign */
					if(pointer_2 > -1) sign[this.sm.getNthMinAbsValueIndex(subWS, pointer_2 + 1)] *= -1;
					this.Calculator.copy(nextWS, subWS);
					for (int j = 0; j < nextWS.length; j++) {
						if(A0G1H2 == 0) nextWS[j] -= sign[j] * this.vPa.L1[u_row] * coef 
								/ (Math.abs(this.WeightA[u_row][nextActive.get(j)] * Math.abs(this.tsda.givenF[u_row][nextActive.get(j)])));
						else if(A0G1H2 == 1) nextWS[j] -= sign[j] * this.vPa.L1[u_row] * coef 
								/ (Math.abs(this.WeightG[u_row][nextActive.get(j)] * Math.abs(this.tsda.givenG[u_row][nextActive.get(j)])));
						else if(A0G1H2 == 2) nextWS[j] -= sign[j] * this.vPa.L1h[u_row] *  this.R[u_row];
					}
					
					/* Calculate Times Inverse */
					this.Calculator.multAx(T_mmInvActive, nextWS, nextWS2);
					
					/* Check Sign */
					for (int j = 0; j < nextWS2.length; j++) {
						if(nextWS2[j] * sign[j] < 0 ) possible = false;
						if(Math.abs(nextWS2[j]) > 0.5 && this.vSet.Degradation != 0.0 && this.H == null) possible = false;
					}
					
					/* Fine if sign is matched or regulation term is very small */
					if(A0G1H2 < 2){
						if(possible | this.vPa.L1[u_row] <= 1.0e-20) break;
					} else if(A0G1H2 ==2){
						if(possible | this.vPa.L1h[u_row] <= 1.0e-20) break;
					}
					
					/* Change Back Sign */
					if(pointer_2 > -1) sign[this.sm.getNthMinAbsValueIndex(subWS, pointer_2 + 1)] *= -1;
					pointer_2++;
				}
				if(A0G1H2 < 2){
					if(possible | this.vPa.L1[u_row] <= 1.0e-20) break;
				} else if(A0G1H2 ==2){
					if(possible | this.vPa.L1h[u_row] <= 1.0e-20) break;
				}
				if(pointer_1 > -1) sign[this.sm.getNthMinAbsValueIndex(subWS, pointer_1 + 1)] *= -1;				
				pointer_1++;
			}
			if(A0G1H2 < 2){
				if(!(possible | this.vPa.L1[u_row] <= 1.0e-20)) continue;
			} else if(A0G1H2 ==2){
				if(!(possible | this.vPa.L1h[u_row] <= 1.0e-20)) continue;
			}
			
			/* Calculate Arg.Min */
			double min = 0;
			if(A0G1H2 == 0) {
				this.Calculator.getSubMatrix(this.Txx_mm, nextActive, T_mmInvActive);
				min = this.Calculator.multxtAx(T_mmInvActive, nextWS2)
						- (2 * this.Calculator.dotProduct(subWS, nextWS2));
			} else if(A0G1H2 == 1) { 
				this.Calculator.getSubMatrix(this.Tzz_mm, nextActive, T_mmInvActive);
				min = this.Calculator.multxtAx(T_mmInvActive, nextWS2) 
						- (2 * this.Calculator.dotProduct(subWS, nextWS2));
			} else if(A0G1H2 == 2) { 
				this.Calculator.getSubMatrix(this.Txx, nextActive, T_mmInvActive);
				min = this.Calculator.multxtAx(T_mmInvActive, nextWS2) 
						- (2 * this.Calculator.dotProduct(subWS, nextWS2));
			}
			
			if(A0G1H2 < 2){
				min *=  coef_inv;
			} else if(A0G1H2 == 2){
				min *=  coef_inv;				
			}
			for (int j = 0; j < nextWS2.length; j++) {
				if(A0G1H2 == 0) min += 2 * this.vPa.L1[u_row] * Math.abs(nextWS2[j])
						 / (Math.abs(this.WeightA[u_row][nextActive.get(j)] * Math.abs(this.tsda.givenF[u_row][nextActive.get(j)])));
				 else if(A0G1H2 == 1) min += 2 * this.vPa.L1[u_row] * Math.abs(nextWS2[j]) 
						 / (Math.abs(this.WeightG[u_row][nextActive.get(j)] * Math.abs(this.tsda.givenG[u_row][nextActive.get(j)])));
				 else if(A0G1H2 == 2) min += 2 * this.vPa.L1h[u_row] * Math.abs(nextWS2[j]);
			}
			
			/* Skip the case that (Parent < Next), (0 < Next) or (Root_size=1 < Next)  */
			if(min > 0 || min > parent_min) {
				continue;
			}
			if(A0G1H2 == 0 && argMinA > min) {
				this.argMinA = min;
				this.activeSetA = new ArrayList<Integer>(nextActive);
				this.nextA.clear();
				for (int i = 0; i < nextActive.size(); i++) {
					this.nextA.add(nextWS2[i]);
				}
				Collections.sort(this.activeSetA);
			} else if(A0G1H2 == 1 && argMinG > min) {
				this.argMinG = min;
				this.activeSetG = new ArrayList<Integer>(nextActive);
				this.nextG.clear();
				for (int i = 0; i < nextActive.size(); i++) {
					this.nextG.add(nextWS2[i]);
				}
				Collections.sort(this.activeSetG);
			} else if(A0G1H2 == 2 && argMinH > min) {
				this.argMinH = min;
				this.activeSetH = new ArrayList<Integer>(nextActive);
				this.nextH.clear();
				for (int i = 0; i < nextActive.size(); i++) {
					this.nextH.add(nextWS2[i]);
				}
				Collections.sort(this.activeSetH);
			}
			
			if(nextActive.size() < limitLength && nextActive.size() < this.th_EdgeNum) {
				this.calcArgMin(nextActive, nextNonActive, ws, u_row, min, limitLength, A0G1H2, coef, coef_inv);
			}
		}
	}
	
	/**
	 * Calculate argMin of A and Z
	 * @param Active
	 * @param ws
	 * @param u_row
	 * @param A_TRUE
	 */
	protected void calcArgMin(ArrayList<Integer> Active, double[] ws, int u_row, int A0G1H2, double coef, double coef_inv){
		Collections.sort(Active);
		double[] subWS = new double[Active.size()];
		double[] nextWS = new double[Active.size()];
		double[] nextWS2 = new double[Active.size()];
		double[] sign = new double[Active.size()];
		double[][] T_mmInvActive = new double[Active.size()][Active.size()];
		double[][] T_mmInvActive_ws = new double[Active.size()][Active.size()];
		for (int j = 0; j < subWS.length; j++) {
			subWS[j] = ws[Active.get(j)];
		}
		this.Calculator.copy(nextWS, subWS);
		
		if(A0G1H2 ==0) {
			this.Calculator.getSubMatrix(this.Txx_mm, Active, T_mmInvActive);
			this.Calculator.symmetricInverseWithWS(T_mmInvActive, T_mmInvActive_ws);
		} else if(A0G1H2 ==1){
			this.Calculator.getSubMatrix(this.Tzz_mm, Active, T_mmInvActive);
			this.Calculator.symmetricInverseWithWS(T_mmInvActive, T_mmInvActive_ws);
		} else if(A0G1H2 ==2){
			this.Calculator.getSubMatrix(this.Txx_obs, Active, T_mmInvActive);
			this.Calculator.symmetricInverseWithWS(T_mmInvActive, T_mmInvActive_ws);
		}

		/* Set Sign */
		for (int i = 0; i < Active.size(); i++) {
			if(A0G1H2 == 0){
				if(this.A[u_row][Active.get(i)] > 0) sign[i] = 1;
				else if(this.A[u_row][Active.get(i)] < 0) sign[i] = -1;
				else sign[i] = nextWS[i]/Math.abs(nextWS[i]);
			} else if(A0G1H2 == 1){
				if(this.G[u_row][Active.get(i)] > 0) sign[i] = 1;
				else if(this.G[u_row][Active.get(i)] < 0) sign[i] = -1;
				else sign[i] = nextWS[i]/Math.abs(nextWS[i]);
			} else if(A0G1H2 == 2){
				if(this.H[u_row][Active.get(i)] > 0) sign[i] = 1;
				else if(this.H[u_row][Active.get(i)] < 0) sign[i] = -1;
				else sign[i] = nextWS[i]/Math.abs(nextWS[i]);
			}
		}
		
		/* Search When Changing at most two sign */
		boolean possible = false;
		int pointer_1 = -1;
		while(pointer_1 < Active.size()){
			/* Skip if High Dimension */
			if(Active.size() > 30 && pointer_1 > -1) break;
			
			/* Change Sign */
			if(pointer_1 > -1) sign[this.sm.getNthMinAbsValueIndex(subWS, pointer_1 + 1)] *= -1;
			
			int pointer_2 = -1;
			while(pointer_2 < Active.size()){
				/* At first, check original sign */
				if(pointer_2 <= pointer_1 && pointer_1 != -1) {
					pointer_2++;
					continue;
				}
				possible = true;
				
				/* Change Sign */
				if(pointer_2 > -1) sign[this.sm.getNthMinAbsValueIndex(subWS, pointer_2 + 1)] *= -1;
				this.Calculator.copy(nextWS, subWS);
				for (int j = 0; j < nextWS.length; j++) {
					if(A0G1H2 == 0) nextWS[j] -= sign[j] * this.vPa.L1[u_row] * coef
							/ (Math.abs(this.WeightA[u_row][Active.get(j)] * Math.abs(this.tsda.givenF[u_row][Active.get(j)])));
					else if(A0G1H2 == 1) nextWS[j] -= sign[j] * this.vPa.L1[u_row] * coef 
							/ (Math.abs(this.WeightG[u_row][Active.get(j)] * Math.abs(this.tsda.givenG[u_row][Active.get(j)])));
					else if(A0G1H2 == 2) nextWS[j] -= sign[j] * this.vPa.L1h[u_row] * coef;
				}
				
				/* Calculate Times Inverse */
				this.Calculator.multAx(T_mmInvActive, nextWS, nextWS2);

				/* Check Sign */
				for (int j = 0; j < nextWS2.length; j++) {
					if(nextWS2[j] * sign[j] < 0 ) possible = false;
					if(Math.abs(nextWS2[j]) > 0.5 && this.vSet.Degradation != 0.0 && this.H == null) possible = false;
				}
				
				/* Fine if sign is matched or regulation term is very small */
				if(A0G1H2 < 2){
					if(possible | this.vPa.L1[u_row] <= 1.0e-20) break;
				} else if(A0G1H2 ==2){
					if(possible | this.vPa.L1h[u_row] <= 1.0e-20) break;
				}
				
				/* Change Back Sign */
				if(pointer_2 > -1) sign[this.sm.getNthMinAbsValueIndex(subWS, pointer_2 + 1)] *= -1;
				pointer_2++;
			}
			/* Fine if sign is matched or regulation term is very small */
			if(A0G1H2 < 2){
				if(possible | this.vPa.L1[u_row] <= 1.0e-20) break;
			} else if(A0G1H2 ==2){
				if(possible | this.vPa.L1h[u_row] <= 1.0e-20) break;
			}
			
			/* Change Back Sign */
			if(pointer_1 > -1) sign[this.sm.getNthMinAbsValueIndex(subWS, pointer_1 + 1)] *= -1;				
			pointer_1++;
		}
		/* Fine if sign is matched or regulation term is very small */		
		if(A0G1H2 < 2){
			if(!(possible | this.vPa.L1[u_row] <= 1.0e-20)) return;
		} else if(A0G1H2 ==2){
			if(!(possible | this.vPa.L1h[u_row] <= 1.0e-20)) return;
		}
		
		/* Calculate Arg.Min */
		double min = 0;
		if(A0G1H2 == 0) {
			this.Calculator.getSubMatrix(this.Txx_mm, Active, T_mmInvActive);
			min = this.Calculator.multxtAx(T_mmInvActive, nextWS2)
					- (2 * this.Calculator.dotProduct(subWS, nextWS2));
		} else if(A0G1H2 == 1){ 
			this.Calculator.getSubMatrix(this.Tzz_mm, Active, T_mmInvActive);
			min = this.Calculator.multxtAx(T_mmInvActive, nextWS2) 
					- (2 * this.Calculator.dotProduct(subWS, nextWS2));
		} else if(A0G1H2 == 2) { 
			this.Calculator.getSubMatrix(this.Txx_obs, Active, T_mmInvActive);
			min = this.Calculator.multxtAx(T_mmInvActive, nextWS2) 
					- (2 * this.Calculator.dotProduct(subWS, nextWS2));
		}
		
		if(A0G1H2 < 2){
			min *=  coef_inv;
		} else if(A0G1H2 == 2){
			min *= coef_inv;				
		}
		
		for (int j = 0; j < nextWS2.length; j++) {
			 if(A0G1H2 == 0) min += 2 * this.vPa.L1[u_row] * Math.abs(nextWS2[j])
					 / (Math.abs(this.WeightA[u_row][Active.get(j)] * Math.abs(this.tsda.givenF[u_row][Active.get(j)])));
			 else if(A0G1H2 == 1) min += 2 * this.vPa.L1[u_row] * Math.abs(nextWS2[j]) 
					 / (Math.abs(this.WeightG[u_row][Active.get(j)] * Math.abs(this.tsda.givenG[u_row][Active.get(j)])));
			 else if(A0G1H2 == 2) min += 2 * this.vPa.L1h[u_row] * Math.abs(nextWS2[j]);
		}
		
		/* Skip the case that (Parent < Next), (0 < Next) or (Root_size=1 < Next)  */
		if(A0G1H2 == 0 && argMinA > min) {
			this.argMinA = min;
			this.activeSetA = new ArrayList<Integer>(Active);
			this.nextA.clear();
			for (int i = 0; i < Active.size(); i++) {
				this.nextA.add(nextWS2[i]);
			}
			Collections.sort(this.activeSetA);
		} else if(A0G1H2 == 1 && argMinG > min) {
			this.argMinG = min;
			this.activeSetG = new ArrayList<Integer>(Active);
			this.nextG.clear();
			for (int i = 0; i < Active.size(); i++) {
				this.nextG.add(nextWS2[i]);
			}
			Collections.sort(this.activeSetG);
		} else if(A0G1H2 == 2 && argMinH > min) {
			this.argMinH = min;
			this.activeSetH = new ArrayList<Integer>(Active);
			this.nextH.clear();
			for (int i = 0; i < Active.size(); i++) {
				this.nextH.add(nextWS2[i]);
			}
			Collections.sort(this.activeSetH);
		}
	}
	
	/*
	 * Selecting L1 value having "num" elements 
	 */
	protected double getSelectedValueL1(double[] aws, double[] gws, int row, int num){
		int gLength = 0;
		if(gws!=null) gLength = gws.length;
		double[] ws = new double[aws.length+gLength];
		for (int i = 0; i < aws.length; i++) {
			ws[i]=Math.abs(aws[i]) * Math.abs(this.WeightA[row][i]) * Math.abs(this.tsda.givenF[row][i]);
		}
		for (int i = 0; i < gLength; i++) {
			ws[aws.length+i]=Math.abs(gws[i]) * Math.abs(this.WeightG[row][i]) * Math.abs(this.tsda.givenG[row][i]);
		}
		this.sm.sort(ws, 0, ws.length-1);
		return ws[num];
	}
	
	/*
	 * Store current parameters to vSto
	 */
	protected void storeCurrentSettings(vsStorage vSto){
		/* Parameters */
		vSto.vPa.setParameters(this.x0_, this.v0, this.A, this.F, this.H, this.G, this.D, this.Q, 
							  this.R, this.U, this.I, null, null, null, null, null, null, null);
		
		/* Profiles */
		vSto.setExpectations(this.Txx_obs, this.Txx, this.Tyx, this.Txx_m, this.Txx_mm, this.Txz_mm, 
							 this.Txz_m, this.sum_x, this.sum_xm, this.sum_zm, null, null, null, null, 
							 null, null, null, null, null, null, null, null, this.Calculator);
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
				this.RinvH, this.HtransRinvH, this.D, this.Q, this.Qinv, this.R, this.Rinv, this.RinvMatrix, this.U, this.I, 
				null, null, null, null, null, null, null, null, null, null, null);
		
		/* Profiles */
		vSto.getExpectations(this.Txx_obs, this.Txx, this.Tyx, this.Txx_m, this.Txx_mm, this.Txz_mm, this.Txz_m, 
							 this.sum_x, this.sum_xm, this.sum_zm, null, null, null, null, null, null, 
							 null, null, null, null, null, null, this.Calculator);
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
	
	public int[] getOrder(){
		final int separate = 20;
		int forEach = 5;
		if(forEach > this.sysDim) forEach = this.sysDim;
		
		double[][] temp_x_s =this.Calculator.copy_generate(this.x_s);
		this.Calculator.normalizeCol(temp_x_s);
		double[][] variance = new double[separate][this.sysDim];
		double[] baseLine = new double[this.sysDim];
		for (int r = 0; r < this.tsda.repSize; r++) {
			for (int t = tsda.maxTime - 2; t < tsda.maxTime; t++) {
				this.Calculator.add(baseLine, temp_x_s[t + (this.tsda.maxTime * r)]);
			}
		}
		this.Calculator.rescale(baseLine, 1.0/(tsda.repSize * 2));
		
		for (int i = 0; i < this.sysDim; i++) {
			for (int r = 0; r < this.tsda.repSize; r++) {
				for (int t = 0; t < tsda.maxTime; t++) {
					variance[(int) Math.floor(t/((double)tsda.maxTime/separate))][i] 
						+= (temp_x_s[t + (this.tsda.maxTime * r)][i] - baseLine[i]) 
							* (temp_x_s[t + (this.tsda.maxTime * r)][i] - baseLine[i]);
				}
			}
		}
		double[][] part_ord = new double[separate][this.sysDim];
		for (int s = 0; s < separate; s++) {
			for (int i = 0; i < this.sysDim; i++) {
				part_ord[s][i]=i;
			}
			this.sm.sort(variance[s], part_ord[s], 0, this.sysDim - 1);
		}
		
		ArrayList<Integer> order = new ArrayList<Integer>();
		for (int s = 0; s < separate; s++) {
			int[] random = this.sm.getRandomOrder(0, forEach-1, this.sfmt);
			for (int i = 0; i < forEach; i++) {
				if(!order.contains((int)part_ord[s][random[i]]))order.add((int)part_ord[s][random[i]]);
			}
		}
		for (int i = 0; i < this.sysDim; i++) {
			if(!order.contains(i))order.add(i);
		}
		return this.Calculator.ArrayListToArrayInt(order);
	}
	
	public void setParameters(){
		if(vSet.R_rI <= 0) this.Calculator.setvalue(this.R, this.Calculator.sumofVector(this.R)/this.R.length);
		if(vSet.upDateQ < 0) this.Calculator.setvalue(this.Q, this.Calculator.sumofVector(this.Q)/this.Q.length);
		this.vPa.setParameters(this.x0_, this.v0, this.A, this.F, this.H, this.G, this.D, this.Q, this.R, this.U, this.I, 
				null, null, null, null, null, null, null);
	}
	
	public void variableSelection(){
		variableSelection(this.updatingRow);
	}
	
	public void variableSelection(int up) {
		this.updatingRow = up;
		
		while (true) {
			/* Recall Previous Setting */
			this.recallPreviousSettings(this.vsSto);
			this.setParameters();
			
			/* Get Pruning Row */
			int previous_length = this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]);
			int row = this.calculationOrder[this.updatingRow];
			if(this.sm.getActiveRowCount(this.A, this.G, row) < 5) {
				break;
			}
			
			/* Call */
			if(this.vSet.Print_Progress){
				System.out.println("Current Row = " + row + ", Updating Row = " + this.updatingRow);
				for (int i = 0; i < this.activeSetAList[row].size(); i++) {
					System.out.print(this.activeSetAList[row].get(i) + " ");
				}
				System.out.println();
			}
			
			/* Store Min_Criterion Setting */
			vsStorage vSTO_minCriterion = new vsStorage();
			vSTO_minCriterion.initializeVSStorage(this.sysDim, this.vSet, this.tsda, this.Calculator, 0.0, this.sfmt);
			this.storeCurrentSettings(vSTO_minCriterion);
			
			/* Pruning A */
			for (int j = -1; j < this.activeSetAList[row].size(); j++) {
				int candidate = -1;
				if(j != -1){
					candidate = this.sm.getNthMinAbsValueIndex(this.vPa.getA()[row], this.activeSetAList[row].size() - j, 0);
					this.vPa.getA()[row][candidate] = 0;
					this.vPa.getF()[row][candidate] = 0;
					this.activeSetAList[row].remove(this.activeSetAList[row].indexOf(candidate));
				}
				this.previousLogLikelihood = (-1) * Double.MAX_VALUE;
				this.miss_update_eternal = 0;
				this.iteration = 0;
				while(true){
					this.getParameters();
					this.KalmanFilter(this.vSet.Drug, this.vSet.Input);
					this.KalmanSmoother();
					this.initializeExp(this.vSet.Drug, this.vSet.Input);
					this.calculateConditionalExpectation(this.vSet.Drug, this.vSet.Input);
					if(this.H == null) this.currentLogLikelihood = this.logLikelihood(this.RinvMatrix, null, null);
					else this.currentLogLikelihood = this.logLikelihood(this.HtransRinvH, this.H, this.RinvH);
					
					this.Criterion = this.calculateBIC(false, this.vSet.Mu0_Update, this.vSet.Input, this.vSet.Drug, this.vSet.upDateQ, 
							this.vSet.Degradation > 0.0, this.vSet.R_rI <= 0.0, this.H != null, this.Delta != null, this.Nu != null, this.vSet.Criterion);
					this.currentLogLikelihood = this.penalization(this.currentLogLikelihood, this.A, this.G, this.H, this.vPa.L1, this.vPa.L1h);
					if( this.Criterion < vSTO_minCriterion.Criterion && !this.Calculator.checkAbsValues(this.A, 0.999999)){
						this.storeCurrentSettings(vSTO_minCriterion);
					}
					
					/* Print Progress */
					if (this.vSet.Print_Progress) {
						System.err.println("Updating Num: " + this.updatingRow + " (Row " + this.calculationOrder[this.updatingRow] + ") - itr "+ this.iteration);
						System.err.println("logLikelihood= " + this.currentLogLikelihood);
						System.err.println("Cri= " + this.Criterion);
						System.err.println("L1= " + this.vPa.L1[this.calculationOrder[(int)this.updatingRow]]);
						System.err.println("Num. of Updating Edges= " + this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]) + " of " + this.th_EdgeNum);
						System.err.println("Num. of Total Edges = " + (this.sm.getActiveCount(this.A) + this.sm.getActiveCount(this.G)));
						System.err.println();
					}
				
					/* Check Convergence */
					if(this.currentLogLikelihood - this.previousLogLikelihood < this.vSet.Condition_of_Convergence 
							&& this.currentLogLikelihood > this.previousLogLikelihood
								&& this.iteration > 1){
						break;
					} else {
						this.convergence = false;
					}
					
					/* At Current < Previous : errors in calculation */
					if(this.previousLogLikelihood > this.currentLogLikelihood && 
					   this.previousLogLikelihood - this.currentLogLikelihood > 1.0e-5 && 
					   this.iteration > 5){
						miss_update_eternal++;
						if(miss_update_eternal > 400) break;
					}
					this.previousLogLikelihood = this.currentLogLikelihood;
					
					/* Check whether the update should be the next step */
					this.Update(true, false, false, false, false);
					if(this.checkNaN(Double.MAX_VALUE)){
						this.MissResult();
						break;
					}
					this.iteration++;
				}
				this.recallPreviousSettings(this.vsSto);
				this.setParameters();
				this.setActiveSets(row);			
			}
			
			/* Pruning G*/
			if(this.vSet.Drug){
				/* Call */
				if(this.vSet.Print_Progress){
					for (int i = 0; i < this.activeSetGList[row].size(); i++) {
						System.out.print(this.activeSetGList[row].get(i) + " ");
					}
					System.out.println();
				}
				for (int j = 0; j < this.activeSetGList[row].size(); j++) {
					int candidate = this.sm.getNthMinAbsValueIndex(this.vPa.getG()[(int)row], this.activeSetGList[row].size() - j, 0);
					this.vPa.getG()[(int)row][candidate] = 0;
					this.activeSetGList[row].remove(this.activeSetGList[row].indexOf(candidate));
					this.previousLogLikelihood = (-1) * Double.MAX_VALUE;
					int rev = 0;
					this.iteration = 0;
					while(true){
						this.getParameters();
						this.KalmanFilter(this.vSet.Drug, this.vSet.Input);
						this.KalmanSmoother();
						this.initializeExp(this.vSet.Drug, this.vSet.Input);
						this.calculateConditionalExpectation(this.vSet.Drug, this.vSet.Input);
						if(this.H == null) this.currentLogLikelihood = this.logLikelihood(this.RinvMatrix, null, null);
						else this.currentLogLikelihood = this.logLikelihood(this.HtransRinvH, this.H, this.RinvH);
						this.Criterion = this.calculateBIC(false, this.vSet.Mu0_Update, this.vSet.Input, this.vSet.Drug, this.vSet.upDateQ, 
								this.vSet.Degradation > 0.0, this.vSet.R_rI <= 0.0, this.H != null, this.Delta != null, this.Nu != null, this.vSet.Criterion);
						this.currentLogLikelihood = this.penalization(this.currentLogLikelihood, this.A, this.G, this.H, this.vPa.L1, this.vPa.L1h);
						/* v1 */
						if(this.Criterion < vSTO_minCriterion.Criterion && !this.Calculator.checkAbsValues(this.A, 0.999999)){
							this.storeCurrentSettings(vSTO_minCriterion);
						}
						
						/* Print Progress */
						if (this.vSet.Print_Progress) {
							System.err.println("Updating Num: " + this.updatingRow + " (Row " + this.calculationOrder[this.updatingRow] + ") - itr "+ this.iteration);
							System.err.println("logLikelihood= " + this.currentLogLikelihood);
							System.err.println("Cri= " + this.Criterion);
							System.err.println("L1= " + this.vPa.L1[this.calculationOrder[(int)this.updatingRow]]);
							System.err.println("Num. of Updating Edges= " + this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]) + " of " + this.th_EdgeNum);
							System.err.println("Num. of Total Edges = " + (this.sm.getActiveCount(this.A) + this.sm.getActiveCount(this.G)));
							System.err.println();
						}
					
						/* Check Convergence */
						if(this.currentLogLikelihood - this.previousLogLikelihood < this.vSet.Condition_of_Convergence 
								&& this.currentLogLikelihood > this.previousLogLikelihood
									&& this.iteration > 1){
							break;
						} else {
							this.convergence = false;
						}
						
						/* At Current < Previous : errors in calculation */
						if(this.previousLogLikelihood > this.currentLogLikelihood && this.previousLogLikelihood - this.currentLogLikelihood > 1.0e-5 && this.iteration > 5){
							rev++;
							if(rev > 10) break;
						}
						this.previousLogLikelihood = this.currentLogLikelihood;
						
						/* Check whether the update should be the next step */
						this.Update(true, false, false, false, false);
						if(this.checkNaN(Double.MAX_VALUE)){
							this.MissResult();
							break;
						}
						this.iteration++;
					}
					this.recallPreviousSettings(this.vsSto);				
					this.setParameters();
					this.setActiveSets(row);
				}
			}
			
			/* Get Optimal Set*/
			this.recallPreviousSettings(vSTO_minCriterion);				
			this.setParameters();
			this.storeCurrentSettings(this.vsSto);
			this.setActiveSets(row);
			this.getParameters();
			
			if(this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]) == previous_length) break;
		}
	}
	
	public void variableSelectionH(){
		variableSelectionH(this.updatingRow);
	}
	
	public void variableSelectionH(int up) {
		this.updatingRow = up;
		while (true) {
			/* Recall Previous Setting */
			this.recallPreviousSettings(this.vsSto);				
			this.setParameters();
			
			/* Get Pruning Row */
			int previous_length = this.sm.getActiveRowCount(this.H, null, this.updatingRow);
			int row = this.updatingRow;
			if(this.sm.getActiveRowCount(this.H, null, row) < 2) {
				break;
			}
			
			/* Call */
			if(this.vSet.Print_Progress){
				System.out.println("Current Row = " + row + ", Updating Row = " + this.updatingRow);
				for (int i = 0; i < this.activeSetHList[row].size(); i++) {
					System.out.print(this.activeSetHList[row].get(i) + " ");
				}
				System.out.println();
			}
			
			/* Store Min_Criterion Setting */
			vsStorage vSTO_minCriterion = new vsStorage();
			vSTO_minCriterion.initializeVSStorage(this.sysDim, this.vSet, this.tsda, this.Calculator, 0.0, this.sfmt);
			this.storeCurrentSettings(vSTO_minCriterion);
			
			/* Pruning H */
			for (int j1 = -1; j1 < this.activeSetHList[row].size(); j1++) {
				for (int j2 = -1; j2 < Math.min(this.activeSetHList[row].size() - 1, 2); j2++) {
					if(j1 != -1){
						if(j1 <= j2) continue;
						int candidate = -1;
						candidate = this.sm.getNthMinAbsValueIndex(this.vPa.getH()[row], this.activeSetHList[row].size() - j1, 0);
						this.vPa.getH()[row][candidate] = 0;
						this.activeSetHList[row].remove(this.activeSetHList[row].indexOf(candidate));
						if(j2 != -1) {
							candidate = this.sm.getNthMinAbsValueIndex(this.vPa.getH()[row], this.activeSetHList[row].size() - j2, 0);
							this.vPa.getH()[row][candidate] = 0;
							this.activeSetHList[row].remove(this.activeSetHList[row].indexOf(candidate));							
						}
					}
					this.previousLogLikelihood = (-1) * Double.MAX_VALUE;
					this.miss_update_eternal = 0;
					this.iteration = 0;
					while(true){
						this.getParameters();
						this.KalmanFilter(this.vSet.Drug, this.vSet.Input);
						this.KalmanSmoother();
						this.initializeExp(this.vSet.Drug, this.vSet.Input);
						this.calculateConditionalExpectation(this.vSet.Drug, this.vSet.Input);
						if(this.H == null) this.currentLogLikelihood = this.logLikelihood(this.RinvMatrix, null, null);
						else this.currentLogLikelihood = this.logLikelihood(this.HtransRinvH, this.H, this.RinvH);
						
						this.Criterion = this.calculateBIC(false, this.vSet.Mu0_Update, this.vSet.Input, this.vSet.Drug, this.vSet.upDateQ, 
								this.vSet.Degradation > 0.0, this.vSet.R_rI <= 0.0, this.H != null, this.Delta != null, this.Nu != null, this.vSet.Criterion);
						this.currentLogLikelihood = this.penalization(this.currentLogLikelihood, this.A, this.G, this.H, this.vPa.L1, this.vPa.L1h);
						if( this.Criterion < vSTO_minCriterion.Criterion && !this.Calculator.checkAbsValues(this.A, 0.999999)){
							this.storeCurrentSettings(vSTO_minCriterion);
						}
						
						/* Print Progress */
						if (this.vSet.Print_Progress) {
							System.err.println("Updating Num: " + this.updatingRow + " (Row " + this.updatingRow + ") - itr "+ this.iteration);
							System.err.println("logLikelihood= " + this.currentLogLikelihood);
							System.err.println("Cri= " + this.Criterion);
							System.err.println("L1= " + this.vPa.L1h[(int)this.updatingRow]);
							System.err.println("Num. of Updating Edges= " + this.sm.getActiveRowCount(this.H, null, this.updatingRow) + " of " + this.th_EdgeNum);
							System.err.println("Num. of Total Edges = " + this.sm.getActiveCount(this.H));
							System.err.println();
						}
						
						/* Check Convergence */
						if(this.currentLogLikelihood - this.previousLogLikelihood < this.vSet.Condition_of_Convergence 
								&& this.currentLogLikelihood > this.previousLogLikelihood
									&& this.iteration > 1){
							break;
						} else {
							this.convergence = false;
						}
						
						/* At Current < Previous : errors in calculation */
						if(this.previousLogLikelihood > this.currentLogLikelihood && this.previousLogLikelihood - this.currentLogLikelihood > 1.0e-5 && this.iteration > 5){
							miss_update_eternal++;
							if(miss_update_eternal > 400) break;
						}
						this.previousLogLikelihood = this.currentLogLikelihood;
						
						/* Check whether the update should be the next step */
						this.Update(true, false, true, false, false);
						if(this.checkNaN(Double.MAX_VALUE)){
							this.MissResult();
							break;
						}
						this.iteration++;
					}
					this.recallPreviousSettings(this.vsSto);
					this.setParameters();
					this.setActiveSetsH(row);			
				}
			}
			
			/* Get Optimal Set*/
			this.recallPreviousSettings(vSTO_minCriterion);
			this.setParameters();
			this.storeCurrentSettings(this.vsSto);
			this.setActiveSetsH(row);
			this.getParameters();
			
			if(this.sm.getActiveRowCount(this.H, null, this.updatingRow) == previous_length) break;
		}
	}
	
	public void override_vPa(vsParameter vPa_to) {
		this.vPa = vPa_to;
	}
	
	public void checkPredictionAbility() {
		this.getParameters();
		this.testPredictionAbility = true;
		this.KalmanFilter(this.vSet.Drug, this.vSet.Input);
		this.currentLogLikelihood = this.logLikelihood(this.RinvMatrix, null, null);
		this.Calculator.copy(this.vsSto.predictionAbility, this.predictionAbility);
	}

}