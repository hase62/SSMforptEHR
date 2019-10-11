package Hasegawa.TimeSeries.Linear.VARSSMEx;

import java.util.ArrayList;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.TimeSeries.Linear.VARSSM.vsSetting;
import Hasegawa.TimeSeries.Linear.VARSSM.vsStorage;
import Hasegawa.TimeSeries.Linear.VARSSM.vsInference;
import Hasegawa.matrix.Matrix;
import RandomGenerator.Sfmt;

public class vseInference extends vsInference{

	protected boolean delta_changed = true;

	public vseInference(final int SYD, vsSetting SET, TimeSeriesDataArray TSDA, Matrix MP, Sfmt Sf, vsStorage vStor) {
		super(SYD, SET, TSDA, MP, Sf, vStor);
	}
	
	public void preparationPhase(int max_loop, int max_edge) {
		System.out.println("Get Initial Profiles Using Full Matrix");
		this.initialUpdate(1.0e-20, this.vSet.Degradation != 0.0 && false);
		this.Calculator.setvalue(this.vPa.L1, 1.0e-20);
		if(this.H != null) this.Calculator.setvalue(this.vPa.L1h, 1.0e-20);
	}
	
	public boolean mainRun() {
		/* This is used to avoid a specific problem in using HGC */
		double base_L1_UpDateRate = 1.01;
		boolean previousConvergence = true;
		if(this.continueThisCalculation) {
			this.continueThisCalculation = false;
		} else {
			this.previousLogLikelihood = (-1) * Double.MAX_VALUE;
			this.L1_UpRate = base_L1_UpDateRate;
			this.miss_update_eternal = 0;
			previousConvergence = false;
		}
		
		//int miss_update = 0;
		while(this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]) != 0){
			 if(this.iteration > 20000) break;
			 if(this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]) < 2 && this.iteration > 10000) break;
			 
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
			double currentLogLikelihood_nonReg = this.currentLogLikelihood;
			this.currentLogLikelihood = this.penalization(this.currentLogLikelihood, this.A, this.G, this.H, this.vPa.L1, this.vPa.L1h);
			
			/* Print Progress */
			int num_of_updating_edge = this.sm.getActiveRowCount(this.A, this.G, this.calculationOrder[this.updatingRow]);
			if (this.vSet.Print_Progress) {
				System.err.println("Updating Num: " + this.updatingRow + " (Row " + this.calculationOrder[this.updatingRow] + ") - itr "+ this.iteration);
				System.err.println("Miss-update: " + this.miss_update_eternal);
				System.err.println("logLikelihood= " + this.currentLogLikelihood);
				System.err.println("(non-reg logLikelihood= " + currentLogLikelihood_nonReg + ")");
				System.err.println("Cri= " + this.Criterion);
				System.err.println("L1= " + this.vPa.L1[this.calculationOrder[(int)this.updatingRow]]);
				System.err.println("Num. of Updating Edges= " + num_of_updating_edge + " of " + this.th_EdgeNum);
				System.err.println("Num. of Total Edges = " + (this.sm.getActiveCount(this.A) + this.sm.getActiveCount(this.G)));
				System.err.println();
			}
		
			/* At Current < Previous : errors in calculation */
			if(this.previousLogLikelihood > this.currentLogLikelihood && !previousConvergence && iteration > 5){
				if(this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] < 10) {
					this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] += 0.05;
				} else {
					this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] *= L1_UpRate;
				}
				this.miss_update_eternal++;
				previousConvergence = true;
			} else if(this.currentLogLikelihood - this.previousLogLikelihood < this.vSet.Condition_of_Convergence 
					&& this.currentLogLikelihood >= this.previousLogLikelihood){
				
				/* Decrease L1_UpRate When Continuously Converged, Otherwise Set base_L1_UpDateRate */
				boolean increase_L1 = previousConvergence || this.currentLogLikelihood - this.previousLogLikelihood < this.vSet.Condition_of_Convergence * 0.1;
				if(increase_L1 && this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] > 10) this.L1_UpRate *= 1.05;
				else this.L1_UpRate = base_L1_UpDateRate;
				if(this.L1_UpRate > 1.4) this.L1_UpRate = 1.4;
				
				/* Update L1 According to L1_UpRate When Having Edges but 0.99 When Having No Edge */
				if(this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] < 10) {
					this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] += 0.05;
				} else {
					this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] *= this.L1_UpRate;
				}
				previousConvergence = true;
				this.convergence = true;
			} else {
				previousConvergence = false;
				this.convergence = false;
			}
			
			/* Store Best BIC */
			int num_of_updating_edge_best = this.sm.getActiveRowCount(this.vsSto.vPa.getA(), this.vsSto.vPa.getG(), this.calculationOrder[this.updatingRow]);
			int full_edge_num = this.A.length;
			if(this.G != null) full_edge_num += this.G[0].length;
			if(this.iteration > 10 && 
			   (
			    (this.Criterion < this.vsSto.Criterion &&  
				 this.vPa.L1[this.calculationOrder[(int)this.updatingRow]] > 0.05 &&
				 (this.vsSto.Criterion > (Double.MAX_VALUE * 0.1) |  
				 num_of_updating_edge < this.th_EdgeNum |
				 num_of_updating_edge <= num_of_updating_edge_best)) 
			   | 
			    (num_of_updating_edge_best > full_edge_num - 2 && 
				 this.miss_update_eternal < 5)
			   )
			  ){
				this.storeCurrentSettings(this.vsSto);
			}
			
			/* Update Parameters */	
			this.previousLogLikelihood = this.currentLogLikelihood;
			this.tempPrePa.setParameters(this.vPa, this.Calculator);

			this.Update(false, false, false, this.iteration == 0, false);
			if(this.checkNaN(Double.MAX_VALUE) || miss_update_eternal >= 300){
				this.MissResult();
				break;
			}
			
			/* On the HGC-Super Computer */
			this.iteration++;
			if(this.iteration % 200 == 99 && this.vSet.Spacom){
				this.continueThisCalculation = true;
				break;
			}
		}
		return this.continueThisCalculation;
	}
	
	public void initializeUpdatingRow(boolean AG){
		this.iteration = 0;
		if(AG) {
			this.th_EdgeNum = 1 + Math.log(this.sysDim) / Math.log(2) * 1.5;
			if(this.vSet.Drug) this.th_EdgeNum += Math.log(this.G[0].length) / Math.log(2);
			if(this.vSet.Degradation==0) this.th_EdgeNum++;
			if(this.sm.getSumOfColLength(this.F, this.G) - 1 <= this.th_EdgeNum) {
				this.th_EdgeNum = this.sm.getSumOfColLength(this.F, this.G) - 2;
			}
		}
		else this.th_EdgeNum = this.sysDim;

		ArrayList<Integer> list = new ArrayList<Integer>();
		if(AG){
			this.vPa.L1[this.calculationOrder[this.updatingRow]] = 1.0e-20;
			for (int i = 0; i < this.A[this.calculationOrder[this.updatingRow]].length; i++) {
				if(i == this.calculationOrder[this.updatingRow] && this.vSet.Degradation > 0) continue;
				list.add(i);
			}
			this.activeSetAList[this.calculationOrder[this.updatingRow]] = list;
			if(this.vSet.Drug) {
				list = new ArrayList<Integer>();
				for (int i = 0; i < this.G[this.calculationOrder[this.updatingRow]].length; i++) {
					list.add(i);
				}
				this.activeSetGList[this.calculationOrder[this.updatingRow]] = list;
			}
		} else if(!AG){
			this.vPa.L1h[this.calculationOrder[this.updatingRow]] = 1.0e-20;
			for (int i = 0; i < this.H[this.calculationOrder[this.updatingRow]].length; i++) {
				list.add(i);
			}
			this.activeSetHList[this.calculationOrder[this.updatingRow]] = list;
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
		if(this.miss_update_eternal < 300) {
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
		if(this.miss_update_eternal < 300) {
			this.vPa.setF(F_up_ws);
		}
		
		/*
		 * update H
		 */
		H_up_ws = null;
		if((this.H != null && !update_H)) {
			double[][] Txx_obsInv = new double[this.sysDim][this.sysDim];
			this.Calculator.symmetricInverse(this.Txx_obs, Txx_obsInv);
			this.Calculator.changesymmetric(Txx_obsInv);
			
			H_up_ws = new double[this.tsda.elementNum][this.sysDim];
			this.Calculator.multAB(this.Tyx, Txx_obsInv, H_up_ws);
			if(this.miss_update_eternal < 300) {
				this.vPa.setH(H_up_ws);
			}
		} else if(this.H != null && update_H) {
			H_up_ws = new double[this.tsda.elementNum][this.sysDim];
			this.Calculator.setvalue(this.H_up_ws, 0);
			double[] hws = new double[this.sysDim];
			for (int i = 0; i < this.H.length; i++) {
				this.Calculator.copy(hws, this.Tyx[i]);
				
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
			if(this.miss_update_eternal < 300) {
				this.vPa.setH(H_up_ws);
			}
		}
		
		/* 
		 * update R 
		 */
		for (int i = 0; i < this.tsda.elementNum; i++) {
			if(this.H == null) R_up_ws[i] = this.Txx_obs[i][i] - 2 * this.Tyx[i][i] + this.Tyy[i];
			else R_up_ws[i] = this.Calculator.multxtAx(Txx_obs, H_up_ws[i]) - 2 * this.Calculator.dotProduct(H_up_ws[i], this.Tyx[i]) + this.Tyy[i];	
			R_up_ws[i] /= this.usedObservationalTimeNum;
			if(R_up_ws[i] < 0) 	R_up_ws[i] = this.R[i];
			if(R_up_ws[i] > Math.abs(this.vSet.R_rI)) R_up_ws[i] = Math.abs(this.vSet.R_rI);
		}
		if (this.vSet.R_rI<=0.0) {
			this.Calculator.setvalue(R_up_ws, this.Calculator.sumofVector(R_up_ws) / (double) R_up_ws.length);
		}
		if(this.miss_update_eternal < 200) {
			this.vPa.setR(R_up_ws);
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
			if(this.miss_update_eternal < 200) {
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
			if(this.miss_update_eternal < 30) {
				this.vPa.setQ(this.Q_up_ws);	
			}
		}
	}
	
	public boolean hRun() {
		/* This is used to avoid a specific problem in using HGC */
		double base_L1_UpDateRate = 1.01;
		boolean previousConvergence = true;
		if(this.continueThisCalculation) {
			this.continueThisCalculation = false;
		} else {
			this.previousLogLikelihood = (-1) * Double.MAX_VALUE;
			this.L1_UpRate = base_L1_UpDateRate;
			this.miss_update_eternal = 0;
			previousConvergence = false;
		}
		
		//int miss_update = 0;
		while(this.sm.getActiveRowCount(this.H, null, this.calculationOrder[this.updatingRow]) != 0  && this.iteration < 20000){
			 if(this.iteration > 20000) break;
			 if(this.sm.getActiveRowCount(this.H, null, this.calculationOrder[this.updatingRow]) < 3 && this.iteration > 10000) break;

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
			double currentLogLikelihood_nonReg = this.currentLogLikelihood;
			this.currentLogLikelihood = this.penalization(this.currentLogLikelihood, this.A, this.G, this.H, this.vPa.L1, this.vPa.L1h);
			
			/* Print Progress */
			int num_of_updating_edge = this.sm.getActiveRowCount(this.H, null, this.calculationOrder[this.updatingRow]);
			if (this.vSet.Print_Progress) {
				System.err.println("Updating Num: " + this.updatingRow + " (Row " + this.calculationOrder[this.updatingRow] + ") - itr "+ this.iteration);
				System.err.println("Miss-update: " + this.miss_update_eternal);
				System.err.println("logLikelihood= " + this.currentLogLikelihood);
				System.err.println("(non-reg logLikelihood= " + currentLogLikelihood_nonReg + ")");
				System.err.println("Cri= " + this.Criterion);
				System.err.println("L1h= " + this.vPa.L1h[this.calculationOrder[(int)this.updatingRow]]);
				System.err.println("Num. of Total Edges (A, G) = " + (this.sm.getActiveCount(this.A) + this.sm.getActiveCount(this.G)));
				System.err.println("Num. of Updating Edges (H) = " + num_of_updating_edge);
				System.err.println("Num. of Total Edges (H) = " + (this.sm.getActiveCount(this.H)));
				System.err.println();
			}
		
			/* At Current < Previous : errors in calculation */
			if(this.previousLogLikelihood > this.currentLogLikelihood && !previousConvergence && iteration > 5){
				if(this.vPa.L1h[this.calculationOrder[(int)this.updatingRow]] < 10) {
					this.vPa.L1h[this.calculationOrder[(int)this.updatingRow]] += 0.05;
				} else {
					this.vPa.L1h[this.calculationOrder[(int)this.updatingRow]] *= this.L1_UpRate;
				}
				//miss_update ++;
				this.miss_update_eternal++;
				previousConvergence = true;
			} else if(this.currentLogLikelihood - this.previousLogLikelihood < this.vSet.Condition_of_Convergence 
					&& this.currentLogLikelihood >= this.previousLogLikelihood){
				
				/* Decrease L1_UpRate When Continuously Converged, Otherwise Set base_L1_UpDateRate */
				boolean increase_L1 = previousConvergence || this.currentLogLikelihood - this.previousLogLikelihood < this.vSet.Condition_of_Convergence * 0.1;
				if(increase_L1 && this.vPa.L1h[this.calculationOrder[(int)this.updatingRow]] > 10) this.L1_UpRate *= 1.05;
				else this.L1_UpRate = base_L1_UpDateRate;
				if(this.L1_UpRate > 1.4) this.L1_UpRate = 1.4;
				
				/* Update L1 According to L1_UpRate When Having Edges but 0.99 When Having No Edge */
				if(this.vPa.L1h[this.calculationOrder[(int)this.updatingRow]] < 10) {
					this.vPa.L1h[this.calculationOrder[(int)this.updatingRow]] += 0.5;
				} else {
					this.vPa.L1h[this.calculationOrder[(int)this.updatingRow]] *= this.L1_UpRate;
				}
				//miss_update = 0;
				previousConvergence = true;
				this.convergence = true;
			} else {
				previousConvergence = false;
				this.convergence = false;
			}			
			
			/* Store Best BIC */
			if( this.Criterion < this.vsSto.Criterion && this.iteration > 0){
				this.storeCurrentSettings(this.vsSto);
			}
			
			/* Update Parameters */	
			this.previousLogLikelihood = this.currentLogLikelihood;
			this.tempPrePa.setParameters(this.vPa, this.Calculator);
			
			this.Update(false, false, true, this.iteration == 0, false);
			if(this.checkNaN(Double.MAX_VALUE)  || miss_update_eternal >= 300){
				this.MissResult();
				break;
			}
			
			/* On the HGC-Super Computer */
			this.iteration++;
			if(this.iteration % 200 == 99 && this.vSet.Spacom){
				this.continueThisCalculation = true;
				break;
			}
		}
		return this.continueThisCalculation;
	}
	
	public void without_L1_Run() {
		/* This is used to avoid a specific problem in using HGC */
		this.setParameters();
		this.currentLogLikelihood = -1.0 * Double.MAX_VALUE;
		this.Calculator.setvalue(this.vPa.L1, 0);
		if(this.H != null) this.Calculator.setvalue(this.vPa.L1h, 0);
		this.vSet.WeightLASSO = 0.0;
		this.convergence = false;
		while(Math.abs(this.currentLogLikelihood - this.previousLogLikelihood) > this.vSet.Condition_of_Convergence) {
			/* Main Process */
			this.getParameters();
			this.KalmanFilter(this.vSet.Drug, this.vSet.Input);
			this.KalmanSmoother();
			this.initializeExp(this.vSet.Drug, this.vSet.Input);
			this.calculateConditionalExpectation(this.vSet.Drug, this.vSet.Input);
			
			/* Calculate logLikelihood, used for BIC, and penalized logLikelihood */
			this.previousLogLikelihood = this.currentLogLikelihood;
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
				System.err.println("L1h= " + this.vPa.L1h[this.calculationOrder[(int)this.updatingRow]]);
				System.err.println("Num. of Total Edges (A, G) = " + (this.sm.getActiveCount(this.A) + this.sm.getActiveCount(this.G)));
				System.err.println("Num. of Total Edges (H)= " + (this.sm.getActiveCount(this.H)));
				System.err.println();
			}
			/* Store Best BIC */
			if( this.Criterion < this.vsSto.Criterion && this.iteration > 0){
				this.storeCurrentSettings(this.vsSto);
			}
			
			/* Update Parameters */	
			this.tempPrePa.setParameters(this.vPa, this.Calculator);
			this.Update(false, false, true, this.iteration == 0, true);
			if(this.checkNaN(Double.MAX_VALUE)){
				this.MissResult();
				break;
			}
		}
	}

	public double getRatioXZ() {
		return(this.Calculator.getDiagSum(this.Txx_mm) / this.Tzz_mm[0][0]);
	}
}