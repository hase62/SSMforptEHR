package Hasegawa.TimeSeries.Linear.VARSSMEx;

import java.io.File;
import java.io.IOException;
import java.util.Locale;import Hasegawa.IO.Reader;
import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.TimeSeries.Linear.VARSSM.vsSetting;
import Hasegawa.TimeSeries.Linear.VARSSM.vsStorage;
import Hasegawa.matrix.Matrix;
import RandomGenerator.Sfmt;

public class mainVARSSMEx {
	public static void main(String[] args) throws IOException {
		
		System.out.println("Start!!");
		final int argsLength = 6;
		
		if (args.length < argsLength) {
			System.err.println("java -jar VARSSM.jar <root> <data> <Structure> <id> <varssm.set> <output> <drug> <drug.z0> <drug structure>");
			/* Users/takanorihasegawa/Dropbox/Program/java/input/SSM 
			 * /Users/takanorihasegawa/Dropbox/Program/java/input/SSM Input/WNT5A-Stable90-obsData_norm_10nodes_Noise0.1-0.3ID0-10.txt NO 0 Set/varssm_simpletest.set NO NO NO NO
			 */
			System.exit(0);
		}

		final long timeAtStart = System.currentTimeMillis();
		Locale.setDefault(Locale.ENGLISH);
		
		Matrix Calculator = new Matrix();
		/* Read Time-series */
		final TimeSeriesDataArray TimeSeriesData = new TimeSeriesDataArray(Calculator);
		
		/* Read Setting */
		TimeSeriesData.Set(args[0] + "/" + args[1], "\t");
		final int ID = Integer.parseInt(args[3]);
		final vsSetting Setting = new vsSetting();
		Setting.set_vs(Setting.settingReader(args[0] + "/" + args[4], "\t"));
		
		/* Random Seed */
		Sfmt Sf = new Sfmt(ID);
		
		/* Set System Dimension */
		int systemDimension = Setting.systemDimension;
		if(systemDimension > TimeSeriesData.elementNum) {
			systemDimension = TimeSeriesData.elementNum;
		}
		
		/* Set Weight Matrix */
		TimeSeriesData.givenF = new double[systemDimension][systemDimension];
		if(Setting.GivenRegulation && systemDimension == TimeSeriesData.elementNum) {
			TimeSeriesData.givenF = Reader.ReadMatrixDouble(args[0] + "/" + args[2], "\t");
		} else {
			Calculator.setvalue(TimeSeriesData.givenF, 1.0);	
		}
		TimeSeriesData.givenG = null;
		if(Setting.Drug) {
			TimeSeriesData.readDrugProfiles(args[0]+ "/" +args[6], args[0] + "/" + args[7], "\t");
			if(Setting.GivenRegulation) {
				try {
					TimeSeriesData.givenG = Reader.ReadMatrixDouble(args[0] + "/" + args[8], "\t");
				} catch (Exception e) {
					System.err.println("Please prepare <drug structure>");
					System.exit(0);
				}
			} else {
				TimeSeriesData.givenG = new double[systemDimension][TimeSeriesData.drug0_[0].length];
				Calculator.setvalue(TimeSeriesData.givenG, 1.0);		
			}
		}
		
		/* Set Steady State of The Time-series as 0 */
		if(Setting.LongSteadyState)	TimeSeriesData.changeDataForSteadyState();

		/* Set Simulation Interval */
		if(Setting.shortenInterval!=1.0) {
			TimeSeriesData.shortenSimulationInterval((int)Setting.shortenInterval);
		}
		
		/* Out Put Information */
		System.out.println("Data=" + args[0] + "/" + args[1]);
		System.out.println("Max Update Count=" + Setting.maxLoop);
		System.out.println("Maximum Calculation Time=" + TimeSeriesData.maxTime);
		System.out.println("System Dimension=" + systemDimension);
		System.out.println("The Amount of Elements=" + TimeSeriesData.elementNum);
		System.out.println("Maximum Time = " + TimeSeriesData.maxTime);
		System.out.println("Replicate Size = " + TimeSeriesData.repSize);
		System.out.println("Amount of Time = " + TimeSeriesData.observationalTimeNum);
		System.out.println("Random Seed ID = " + ID);
		System.out.println("Criterion(0:BIC, 1:AIC)=" + Setting.Criterion);
		System.out.println("Simulation Interval=*" + 1/Setting.shortenInterval);
		System.out.println("-----Program is Executed-----");
		
		/* Main Process */
		vsStorage vSTO = new vsStorage();
		vSTO.initializeVSStorage(systemDimension, Setting, TimeSeriesData, Calculator, 0.0, Sf);
		vseInference vINF = new vseInference(systemDimension, Setting, TimeSeriesData, Calculator, Sf, vSTO);
		
		/* Set Initial Parameter Values */
		System.out.println("Initial Loop 1 - Fixing Upate ...... ");
		vINF.preparationPhase(0, 0);

		/* Adjust Scale for Z */
		double z_multi = 1;		
		if(Setting.Drug){
			double[][] givenF_tmp = Calculator.copy_generate(TimeSeriesData.givenF);
			Calculator.setvalue(TimeSeriesData.givenF, 1.0);
			double[][] givenG_tmp = null;
			if(Setting.Drug) {
				givenG_tmp = Calculator.copy_generate(TimeSeriesData.givenG);
				Calculator.setvalue(TimeSeriesData.givenG, 1.0);
			}
			z_multi = Math.sqrt(vINF.getRatioXZ()) * 1.0;
			System.out.println("Initial Loop 2 - Fixing Upate ...... " + z_multi);
			Sf = new Sfmt(ID);
			TimeSeriesData.drugModifyVariance(z_multi);		
			vSTO = new vsStorage();
			vSTO.initializeVSStorage(systemDimension, Setting, TimeSeriesData, Calculator, 0.0, Sf);
			vINF = new vseInference(systemDimension, Setting, TimeSeriesData, Calculator, Sf, vSTO);
			
			vINF.preparationPhase(0, 0);
			Calculator.copy(TimeSeriesData.givenF, givenF_tmp);
			if(Setting.Drug) {
				Calculator.copy(TimeSeriesData.givenG, givenG_tmp);				
			}
		}
		System.out.println("Initial Loop 2 - Normalization Done ...... " + Math.sqrt(vINF.getRatioXZ()));	
		vINF.vsSto.Criterion = Double.MAX_VALUE;
		
		/* Temp
		File mkdir_ = new File(args[0] + "/" + args[5] + "/UKF/Cri" +Setting.Criterion +"ID" + ID);
		if (!mkdir_.exists());
		String pass2_ = mkdir_.getPath();
		Calculator.write(pass2_ + "/z_multi" + "-" + ID + ".txt", z_multi);
		mkdir_ = new File(args[0] + "/" + args[5] + "/UKF/Cri" +Setting.Criterion +"ID" + ID + "/WithoutL1");
		pass2_ = mkdir_.getPath();
		Calculator.write(pass2_ + "/z_multi" + "-" + ID + ".txt", z_multi);		
		*/
		
		/* Prediction Error */
		if(Setting.testPredictionAbility != 0){
			for (int r = 0; r < TimeSeriesData.Validity.size(); r++) {
				TimeSeriesData.Validity.get(r).set((int) Setting.testPredictionAbility, false);
			}
		}
		for (int l = 0; l < vINF.vSet.maxLoop; l++) {
			if(System.currentTimeMillis() - timeAtStart > (long)(Setting.timer * 1.0 / 2.0)) l = (int)vINF.vSet.maxLoop - 1;
			/* Set Calculation Order, i.e., Calculate Genes Varied among Early Time-steps First */
			vINF.calculationOrder = vINF.getOrder();
			/* For Each Gene */
			for (int i = 0; i < vINF.calculationOrder.length; i++) {
				if(System.currentTimeMillis() - timeAtStart > (long)(Setting.timer * 1.0 / 2.0) && l != (int)vINF.vSet.maxLoop - 1) break;

				System.out.println("Loop:" + l + " Row:" + i + " Start "+ (System.currentTimeMillis() - timeAtStart) / (1000 * 60));
				vINF.Weight();
				vINF.updatingRow = i;
				vINF.initializeUpdatingRow(true);
				vINF.setParameters();
				while(true) {
					if(!vINF.mainRun()) break;
				}
				vINF.recallPreviousSettings(vINF.vsSto);
				vINF.setParameters();
				vINF.setActiveSetsAll();
				System.out.println("Criterion " + vINF.vsSto.Criterion 
								   + ", Log-likelihood" + vINF.vsSto.logLikelihood 
								   + ", Prediction Error " + Calculator.sumofVector(vINF.vsSto.predictionAbility));
				
				/* Variable Selection */
				if(l > vINF.vSet.maxLoop - 3){
					if(System.currentTimeMillis() - timeAtStart > (long)(Setting.timer * 1.0 / 2.0) && l != (int)vINF.vSet.maxLoop - 1) break;
					vINF.variableSelection();
					vINF.setActiveSetsAll();
					System.out.println("After Checking All Combination:\nCriterion " + vINF.vsSto.Criterion + ", Prediction Error " + Calculator.sumofVector(vINF.vsSto.predictionAbility));
				}
			}
		}

		vINF.calculationOrder = new int[TimeSeriesData.elementNum];
		for (int i = 0; i < TimeSeriesData.elementNum; i++) {
			vINF.calculationOrder[i] = i;	
		}
		int starting_row = 0;
		if(vINF.vSet.maxLoop > 1) starting_row = 1;
		if(vINF.vSet.maxLoop > 4) starting_row = 2;
		for (int l = starting_row; l < vINF.vSet.maxLoop; l++) {
			if(vINF.vPa.getH() == null) break;
			if(System.currentTimeMillis() - timeAtStart > (long)(Setting.timer * 3.0 / 4.0)) l = (int)vINF.vSet.maxLoop - 1;

			/* For Each Gene */
			for (int i = 0; i < TimeSeriesData.elementNum; i++) {
				if(System.currentTimeMillis() - timeAtStart > (long)(Setting.timer * 3.0 / 4.0) && l != (int)vINF.vSet.maxLoop - 1) break;

				System.out.println("Loop:" + l + " Row:" + i + " Start "+ (System.currentTimeMillis() - timeAtStart) / (1000 * 60));
				vINF.Weight();
				vINF.updatingRow = i;
				vINF.initializeUpdatingRow(false);
				vINF.setParameters();
				while(true) {
					if(!vINF.hRun()) break;
				}
				vINF.recallPreviousSettings(vINF.vsSto);
				vINF.setParameters();
				vINF.setActiveSetsAll();
				System.out.println("Criterion " + vINF.vsSto.Criterion 
								   + ", Log-likelihood" + vINF.vsSto.logLikelihood 
								   + ", Prediction Error " + Calculator.sumofVector(vINF.vsSto.predictionAbility));
				
				/* Variable Selection */
				if(l > vINF.vSet.maxLoop - 2){
					if(System.currentTimeMillis() - timeAtStart > (long)(Setting.timer * 3.0 / 4.0) && l != (int)vINF.vSet.maxLoop - 1) break;
					vINF.variableSelectionH();
					vINF.setActiveSetsAll();
					System.out.println("After Checking All Combination:\nCriterion " + vINF.vsSto.Criterion + ", Prediction Error " + Calculator.sumofVector(vINF.vsSto.predictionAbility));
				}
			}
		}
		if(Setting.testPredictionAbility != 0){
			vINF.checkPredictionAbility();
		}
				
		/* Write Results */
		File mkdir = new File(args[0] + "/" + args[5] + "/UKF/Cri" +Setting.Criterion +"ID" + ID);
		if (mkdir.mkdir());
		else mkdir.mkdirs();
		String pass2 = mkdir.getPath();
		Calculator.write(pass2 + "/Loglikelihood"  + "-" + ID+ ".txt", vSTO.logLikelihood);
		Calculator.write(pass2 + "/Criterion" + "-" + ID + ".txt", vSTO.Criterion);
		Calculator.write(pass2 + "/PredictionError" + "-" + ID + ".txt", Calculator.sumofVector(vSTO.predictionError));
		
		double calcTime = (System.currentTimeMillis() - timeAtStart) / 60 / 1000;
		if(Setting.testPredictionAbility == 0){
			
			Calculator.write(pass2 + "/F" + "-" + ID + ".txt", vSTO.vPa.getF());
			if(systemDimension != TimeSeriesData.elementNum) Calculator.write(pass2 + "/H" + "-" + ID + ".txt", vSTO.vPa.getH());
			Calculator.write(pass2 + "/Q" + "-" + ID + ".txt", vSTO.vPa.getQ());
			if(Setting.Degradation > 0.0) Calculator.write(pass2 + "/D" + "-" + ID + ".txt",vSTO.vPa.getD());
			Calculator.write(pass2 + "/R" + "-" + ID + ".txt", vSTO.vPa.getR());
			if(Setting.Input) Calculator.write(pass2 + "/U" + "-" + ID + ".txt",vSTO.vPa.getU());
			Calculator.write(pass2 + "/x0" + "-" + ID + ".txt",vSTO.vPa.getx0());
			Calculator.write(pass2 + "/xp" + "-" + ID + ".txt",vSTO.getx_p());
			Calculator.write(pass2 + "/xf" + "-" + ID + ".txt",vSTO.getx_f());
			Calculator.write(pass2 + "/xs" + "-" + ID + ".txt",vSTO.getx_s());
			Calculator.write(pass2 + "/yp" + "-" + ID + ".txt",vSTO.gety_p());
			Calculator.write(pass2 + "/Time" + "-" + ID + ".txt", calcTime);			
			Calculator.write(pass2 + "/v0" + "-" + ID + ".txt",vSTO.vPa.getv0());
			Calculator.write(pass2 + "/z_multi" + "-" + ID + ".txt", z_multi);
			if(Setting.Drug) Calculator.write(pass2 + "/G" + "-" + ID + ".txt",vSTO.vPa.getG());
		} else {
			Calculator.write(pass2 + "/PredictionAbility" + "-" + ID + ".txt", Calculator.sumofVector(vSTO.predictionAbility));
		}
		
		System.out.println("Best MarginalLogLikelihood = " + vSTO.logLikelihood);
		System.out.println("Best Criteria = " + vSTO.Criterion + "\n");
		System.out.println("");
		System.out.println("Program End " + calcTime);

		/* Result Without L1 */
		vINF.without_L1_Run();
		
		/* Write Results */
		mkdir = new File(args[0] + "/" + args[5] + "/UKF/Cri" +Setting.Criterion +"ID" + ID + "/WithoutL1");
		if (mkdir.mkdir());
		else mkdir.mkdirs();
		pass2 = mkdir.getPath();
		Calculator.write(pass2 + "/Loglikelihood"  + "-" + ID+ ".txt", vSTO.logLikelihood);
		Calculator.write(pass2 + "/Criterion" + "-" + ID + ".txt", vSTO.Criterion);
		Calculator.write(pass2 + "/PredictionError" + "-" + ID + ".txt", Calculator.sumofVector(vSTO.predictionError));
		
		calcTime = (System.currentTimeMillis() - timeAtStart) / 60 / 1000;
		if(Setting.testPredictionAbility == 0){
			
			Calculator.write(pass2 + "/F" + "-" + ID + ".txt", vSTO.vPa.getF());
			if(systemDimension != TimeSeriesData.elementNum) Calculator.write(pass2 + "/H" + "-" + ID + ".txt", vSTO.vPa.getH());
			Calculator.write(pass2 + "/Q" + "-" + ID + ".txt", vSTO.vPa.getQ());
			if(Setting.Degradation > 0.0) Calculator.write(pass2 + "/D" + "-" + ID + ".txt",vSTO.vPa.getD());
			Calculator.write(pass2 + "/R" + "-" + ID + ".txt", vSTO.vPa.getR());
			if(Setting.Input) Calculator.write(pass2 + "/U" + "-" + ID + ".txt",vSTO.vPa.getU());
			Calculator.write(pass2 + "/x0" + "-" + ID + ".txt",vSTO.vPa.getx0());
			Calculator.write(pass2 + "/xp" + "-" + ID + ".txt",vSTO.getx_p());
			Calculator.write(pass2 + "/xf" + "-" + ID + ".txt",vSTO.getx_f());
			Calculator.write(pass2 + "/xs" + "-" + ID + ".txt",vSTO.getx_s());
			Calculator.write(pass2 + "/yp" + "-" + ID + ".txt",vSTO.gety_p());
			Calculator.write(pass2 + "/Time" + "-" + ID + ".txt", calcTime);			
			Calculator.write(pass2 + "/v0" + "-" + ID + ".txt",vSTO.vPa.getv0());
			Calculator.write(pass2 + "/z_multi" + "-" + ID + ".txt", z_multi);
			if(Setting.Drug) Calculator.write(pass2 + "/G" + "-" + ID + ".txt",vSTO.vPa.getG());
		} else {
			Calculator.write(pass2 + "/PredictionAbility" + "-" + ID + ".txt", Calculator.sumofVector(vSTO.predictionAbility));
		}
		
		System.out.println("Best MarginalLogLikelihood = " + vSTO.logLikelihood);
		System.out.println("Best Criteria = " + vSTO.Criterion + "\n");
		System.out.println("");
		System.out.println("Program End " + calcTime);
	}
}
