package Hasegawa.TimeSeries.Linear.SSM;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.TimeSeries.Inference;
import Hasegawa.matrix.Matrix;
import Hasegawa.stat.simpleMath;
import RandomGenerator.Sfmt;

public class ssmInference extends Inference{

	protected ssmSetting sSet;
	protected ssmStorage sSto;
	protected ssmParameter sPa;
	
	protected double[][] ssmD;
	protected double[][] Phi;
	
	protected Boolean monotone;
	
	public ssmInference(final int syd, ssmSetting SET, TimeSeriesDataArray TSDA, Matrix MP, Sfmt Sf, ssmStorage STO) {
		this.sysDim = syd;	
		this.sSet = SET;
		this.tsda = TSDA;
		this.Calculator = MP;
		this.sfmt = Sf;
		this.sSto = STO;
		this.sm = new simpleMath();
	}
	
	public void PreparePermutation(double[][] ph, double[][] pg) {
		this.sSto.perPhi = new double[this.tsda.elementNum][this.tsda.elementNum];
		this.sSto.result = new double[this.tsda.elementNum][this.tsda.elementNum];
		this.Calculator.copy(this.sSto.perPhi, ph);
		if(this.sSet.Drug){
			this.sSto.perDrug = new double[this.tsda.elementNum][this.tsda.drugMulRepSize[0].length];
			this.sSto.resultDrug = new double[this.tsda.elementNum][this.tsda.drugMulRepSize[0].length];	
			this.Calculator.copy(this.sSto.perDrug, pg);
		}
		this.tsda.ReserveData();
	}

	
	protected void generateNextParameter(){
		this.sPa = new ssmParameter(this.sysDim, this.sfmt, this.sSet, this.tsda, this.Calculator);
	}

	protected void initialize() {
		
		this.initilizeCommons(false, this.sSet.Drug, false, true);
		if(this.sSet.Drug) this.initilizeDrug(this.sSet.Input);
		
		/*
		 * Unique Parameters. 
		 */
		this.monotone = true;
		this.RinvH = new double[this.tsda.elementNum][this.sysDim];
		this.HtransRinvH = new double[this.sysDim][this.sysDim];	
	}

	protected void getParameters(ssmParameter indicator) {
		this.x0_ = indicator.getx0();
		this.v0 = indicator.getv0();
		for (int rep = 0; rep < this.tsda.repSize; rep++) {
			this.Calculator.copy(this.x0_s[rep], this.x0_[rep]);
			this.Calculator.copy(this.v0_s[rep], this.v0);
		}
		this.F = indicator.getF();		
		this.H = indicator.getH();
		if(this.sSet.Drug){
			this.G = indicator.getG();
		}
		this.Q = indicator.getQ();
		this.R = indicator.getR();
		this.Calculator.copy(this.Ftrans, this.F);
		this.Calculator.transpose(this.Ftrans);
		this.Calculator.copy(this.Rinv, this.R);
		this.Calculator.DiagInvese(this.Rinv);
		this.Calculator.multAB(this.Rinv, this.H, this.RinvH);
		this.Calculator.multAtB(this.H, this.RinvH, this.HtransRinvH);
		this.I = indicator.getI();
		this.LogRDet = this.Calculator.diagLogDeterminant(this.R);
	}
	
	public void run() {
		
		if(this.sSet.Permutation) this.tsda.MakeBootstrapSample(this.sfmt);
		
		this.initialize();
		this.generateNextParameter();
		
		for (int i = 0; i < this.sSet.maxLoop; i++) {
			
			this.getParameters(this.sPa);
			this.KalmanFilter(this.sSet.Drug, false);
			this.KalmanSmoother();
			this.initializeExp(this.sSet.Drug, false);
			this.calculateConditionalExpectation(this.sSet.Drug, false);
			this.currentLogLikelihood = this.logLikelihood(this.HtransRinvH, this.H, this.RinvH);
			
			if (this.currentLogLikelihood < this.previousLogLikelihood) {
				this.monotone = false;
			}
			if(this.currentLogLikelihood < this.previousLogLikelihood){
				System.out.println("ssssss");
			}
			if (Math.abs(this.currentLogLikelihood - this.previousLogLikelihood) < this.sSet.Condition_of_Convergence) {
				this.convergence = true;
				break;
			}
			
			//Check Convergence and whether H becomes 0 or not. 
			if(this.Calculator.checkNaN(this.F) || 
					this.Calculator.checkAbsValues(this.F, 1.0e2, 1.0e-8) || 
					this.Calculator.checkNaN(this.H) || 
					this.Calculator.checkAbsValues(this.H, 1.0e2, 1.0e-8)){
				this.MissResult();
				break;
			}
				
			this.previousLogLikelihood = this.currentLogLikelihood;
			
			if (this.sSet.Print_Progress) {
				System.err.println("iter=" + i);
				System.err.println("logLikelihood=" + this.currentLogLikelihood);
			}
			this.Update();
		}
		this.calculateBIC(true, this.sSet.Mu0_Update, false, this.sSet.Drug, 1.0, false, 
				this.sSet.R_rI<=0.0, true, false, false, 0.0);
		System.out.println("SystemDimension=" + this.sysDim);
		System.out.println("LogLikelihood=" + this.currentLogLikelihood);
		System.out.println("BIC=" + this.Criterion);
	}
	
	protected int Filter(int time, boolean observed, boolean validity, int rep_num, int count){
		
		if(this.H != null) this.Calculator.multAx(this.H, this.x_p[time], this.y_p[time]);
		else this.Calculator.copy(this.y_p[time], this.x_p[time]);
		if (observed) {
			double[][] IKH = new double[this.sysDim][this.sysDim];
			this.Calculator.symmetricInverse(this.v_p[time], IKH);
			this.Calculator.changesymmetric(IKH);
			this.Calculator.add(IKH, this.HtransRinvH);
			this.Calculator.symmetricInverse(IKH, this.v_f[time]);
			this.Calculator.changesymmetric(this.v_f[time]);

			this.Calculator.copy(this.x_f[time], this.x_p[time]);
			this.Calculator.multAddABtx(this.v_f[time], this.RinvH, this.y[this.tsda.repIDSum.get(rep_num) + count], this.x_f[time]);
			this.Calculator.multSubABx(this.v_f[time], this.HtransRinvH, this.x_p[time], this.x_f[time]);
			
			count++;
			
		} else {
			this.Calculator.copy(this.x_f[time], this.x_p[time]);
			this.Calculator.copy(this.v_f[time], this.v_p[time]);
		}
		
		this.Calculator.copy(this.x_s[time], this.x_f[time]);
		this.Calculator.copy(this.v_s[time], this.v_f[time]);
		
		return count;
	}

	protected void Update() {
		//Mu
		if (this.sSet.Mu0_Update > 0) {
			double[] x_same = new double[this.sysDim];
			double[][] x_sep = new double[this.tsda.repSize][this.sysDim];
			double[][] v = new double[this.sysDim][this.sysDim];

			for (int rep = 0; rep < this.tsda.repSize; rep++) {
				this.Calculator.add(x_same, this.x0_s[rep]);
				this.Calculator.add(x_sep[rep], this.x0_s[rep]);
				this.Calculator.add(v, this.v0_s[rep]);
			}
			this.Calculator.rescale(x_same, (1.0/this.tsda.repSize));
			this.Calculator.rescale(v, (1.0 / this.tsda.repSize));
			if(this.sSet.Mu0_Update == 1){
				for (int i = 0; i < x_sep.length; i++) {
					for (int j = 0; j < x_sep[0].length; j++) {
						x_sep[i][j] = x_same[j];
					}
				}
				
			}
			this.sPa.setx0(x_sep);
			//this.sPa.setv0(v);
		}

		double[][] Txx_obsInv = new double[this.sysDim][this.sysDim];
		double[][] Txx_mmInv = new double[this.sysDim][this.sysDim];;
		this.Calculator.symmetricInverse(this.Txx_obs, Txx_obsInv);
		this.Calculator.changesymmetric(Txx_obsInv);
		this.Calculator.symmetricInverse(this.Txx_mm, Txx_mmInv);
		this.Calculator.changesymmetric(Txx_mmInv);
		
		// H
		double[][] Htemp = new double[this.tsda.elementNum][this.sysDim];
		this.Calculator.multAB(this.Tyx, Txx_obsInv, Htemp);
		this.sPa.setH(Htemp);
		
		//R
		double[] Rtemp = new double[this.tsda.elementNum];
		for (int i = 0; i < this.tsda.elementNum; i++) {
			Rtemp[i] = this.Calculator.multxtAx(Txx_obs, Htemp[i]) - 2 * this.Calculator.dotProduct(Htemp[i], this.Tyx[i]) + this.Tyy[i];
			Rtemp[i] /= this.tsda.observationalTimeNum;
		}
		if(this.sSet.R_rI <= 0) this.Calculator.setvalue(Rtemp, this.Calculator.sumofVector(Rtemp) / (double) Rtemp.length);
		this.sPa.setR(Rtemp);
		
		//F
		double[][] Ftemp = new double[this.sysDim][this.sysDim];
		double[][] Fws = new double[this.sysDim][this.sysDim];
		this.Calculator.copy(Fws, this.Txx_m);
		if(this.sSet.Drug)	this.Calculator.multSubABt(this.G, this.Txz_mm, Fws);
		this.Calculator.multAB(Fws, Txx_mmInv, Ftemp);
		this.sPa.setF(Ftemp);
		
		//G
		if(this.sSet.Drug){
			double[][] Gtemp = new double[this.sysDim][this.tsda.drugMulRepSize[0].length];
			double[][] Gws = new double[this.sysDim][this.tsda.drugMulRepSize[0].length];
			this.Calculator.copy(Gws, this.Txz_m);
			this.Calculator.multSubAB(this.F, this.Txz_mm, Gws);			
			this.Calculator.multAB(Gws, this.Tzz_mmInv, Gtemp);
			this.sPa.setG(Gtemp);
		}
		if (this.sSet.R_rI <= 0.0) {
			double value = 0;
			for (int i = 0; i < Rtemp.length; i++)
				value += Rtemp[i];
			this.Calculator.setvalue(Rtemp, value / (double) Rtemp.length);
			this.sPa.setR(Rtemp);
		}
	}

	public boolean storeBestParameters(boolean isLast){
		if (this.sSto.Criterion > this.Criterion || this.sSet.Permutation || isLast) {
			this.sSto.logLikelihood = this.currentLogLikelihood;
			this.sSto.Criterion = this.Criterion;
			this.sSto.Monotone = this.monotone;
			this.sSto.Convergence = this.convergence;
			this.sSto.sPa.storeParameters(this.sPa);
			if(isLast) return false;
		}
		if (this.sSet.Permutation){
			this.run_CSSM();			
			if(this.Calculator.checkAbsValues(this.Phi, 1, 0)) return true;
			this.ComparePhiMatrix(this.sSto.result, this.sSto.perPhi, this.sPa, this.sSto.perDrug, this.sSto.resultDrug);
		}
		this.sSto.SuccessCount++;
		return false;
	}
	
	public void run_CSSM(){
		System.out.println("------CSSM Start------");
		System.out.println("LogLikelihood(before)=" + this.sSto.logLikelihood);
		
		this.getParameters(this.sSto.sPa);
		this.CSSM();
		
		this.getParameters(this.sPa);
		this.KalmanFilter(this.sSet.Drug, false);
		this.KalmanSmoother();
		
		this.currentLogLikelihood = this.logLikelihood(this.HtransRinvH, this.H, this.RinvH);
		this.calculateBIC(true, this.sSet.Mu0_Update, false, this.sSet.Drug, 1.0, false, 
				this.sSet.R_rI <= 0.0, true, false, false, 0.0);

		this.sSto.setProfiles(this.x0_s, this.x_p, this.x_f, this.x_s, this.y_p);
		
		System.out.println("LogLikelihood(after)=" + this.currentLogLikelihood);
		System.out.println("------CSSM End------");	
	}
	
	protected void CSSM() {
		
		double[] rootR = new double[this.tsda.elementNum];
		for (int i = 0; i < rootR.length; i++) rootR[i] = Math.sqrt(R[i]);
		double[] rootRinv = new double[this.tsda.elementNum];
		for (int i = 0; i < rootRinv.length; i++) rootRinv[i] = Math.sqrt(Rinv[i]);
		
		double[][] rootRinvH = new double[this.tsda.elementNum][this.sysDim];
		this.Calculator.multAB(rootRinv, this.H, rootRinvH);
		
		double[][] u = new double[this.tsda.elementNum][this.sysDim];
		double[] s;
		if (this.tsda.elementNum > this.sysDim) s = new double[this.sysDim];
		else s = new double[this.tsda.elementNum];
		double[][] v = new double[this.sysDim][this.sysDim];
		this.Calculator.SVDecomposition(rootRinvH, u, v, s);
		double[][] Dec = new double[this.sysDim][this.sysDim];
		double[] eigen = new double[this.sysDim];
		for (int i = 0; i < Dec.length; i++) {
			Dec[i][i] = s[i];
			eigen[i] = s[i] * s[i];
		}
	
		double[][] tempH = new double[this.tsda.elementNum][this.sysDim];
		this.Calculator.multAB(u, Dec, tempH);
		this.Calculator.multAB(rootR, tempH);
		this.sPa.setH(tempH);
		
		double[][] tempF = new double[this.sysDim][this.sysDim];
		this.Calculator.transpose(v);
		this.Calculator.multAddABCt(v, F, v, tempF);
		this.sPa.setF(tempF);
		
		if(this.sSet.Drug){
			double[][]tempG = new double[this.sysDim][this.tsda.drugMulRepSize[0].length];
			this.Calculator.multAB(v, this.G, tempG);
			this.sPa.setG(tempG);
		}
		for (int rep = 0; rep < this.tsda.repSize; rep++) {
			this.Calculator.multAx(v, this.x0_[rep], this.x0_s[rep]);
		}
		this.sPa.setx0(this.x0_s);			
	
		double[] eigeninv = new double[this.sysDim];
		this.Calculator.copy(eigeninv, eigen);
		this.Calculator.DiagInvese(eigeninv);
			
		this.ssmD = new double[this.sysDim][this.tsda.elementNum];
		for (int i = 0; i < this.sysDim; i++){
			for (int j = 0; j < this.tsda.elementNum; j++){
				this.ssmD[i][j] = tempH[j][i] * eigeninv[i] * rootRinv[j];
			}
		}
		this.sPa.setssmD(this.ssmD);
		
		Phi = new double[this.tsda.elementNum][this.tsda.elementNum];
		double[][] Dtrans = new double[this.tsda.elementNum][this.sysDim];
		this.Calculator.multAB(eigen, tempF);
		this.Calculator.transpose(this.ssmD, Dtrans);
		this.Calculator.multABC(Dtrans, tempF, this.ssmD, this.Phi);
		this.sPa.setPh(this.Phi);
			
		if(this.sSet.Drug){
			double[][] rootRinvHG = new double[this.tsda.elementNum][G[0].length];
			this.Calculator.multAB(rootRinvH, this.G, rootRinvHG);
			this.sPa.setrootRinvHG(rootRinvHG);			
		}
	}
	
	public void ComparePhiMatrix(double[][] resultPhi, double[][] originalPhi, ssmParameter ip, double[][] originalG, double[][] resultG) {		
		double[][] phi = new double[this.tsda.elementNum][this.tsda.elementNum];
		phi = ip.getPhi();
		for (int i = 0; i < originalPhi.length; i++)
			for (int j = 0; j < originalPhi[0].length; j++)
				if (originalPhi[i][j] > phi[i][j])
					resultPhi[i][j]++;
		if(originalG!=null){
			double[][] rG=new double[ip.getrootRinvHG().length][ip.getrootRinvHG()[0].length];
			rG = ip.getrootRinvHG();
			for (int i = 0; i < originalG.length; i++)
				for (int j = 0; j < originalG[0].length; j++)
					if (originalG[i][j] > rG[i][j])
						resultG[i][j]++;
		}
	}
}
