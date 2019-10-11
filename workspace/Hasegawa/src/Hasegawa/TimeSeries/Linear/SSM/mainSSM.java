package Hasegawa.TimeSeries.Linear.SSM;

import java.io.File;
import java.io.IOException;
import java.util.Locale;

import Hasegawa.IO.Reader;
import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.IO.Writer;
import Hasegawa.matrix.sMatrix;
import Hasegawa.matrix.Matrix;
import RandomGenerator.Sfmt;

public class mainSSM {
	public static void main(String[] args) throws IOException {

		System.out.println("Start!!");
		final int argsLength = 7;

		if (args.length < argsLength) {
			System.err.println("java -jar mainSSM.jar <data> <systemDimension> <repeatnumber> <id> <Drug.txt> <ssm.set> <output>");
			System.err.println("ex) /Users/takaorihasegawa_mpro/Dropbox/Program/java/input/SSM/Input/hirosaki/TSDA.SCALE.tsda_full.MALE.30.60.txt 10 3 1 Z.zt_full.MALE.30.60.SCALE.txt.selected.txt ./../../Set/ssm.set Output");
			System.err.println("if you do not use drug-profile, please fill in <Drug.txt> as /NoDrug ");
			System.err.println("if you want to execute Permutation Test, ");
			System.err.println("java -jar mainSSM.jar <data> <systemDimension> <repeatnumber> <id> <Drug.txt> <ssm.set> <output> <Phi> (<PhiG>)~for DRUG");
			System.exit(0);
		}

		final long timeAtStart = System.currentTimeMillis();
		Locale.setDefault(Locale.ENGLISH);

		final File data = new File(args[0]);
		final int systemDimension = Integer.parseInt(args[1]);
		final int repeatNumber = Integer.parseInt(args[2]);
		final int ID = Integer.parseInt(args[3]);

		ssmSetting Setting = new ssmSetting();
		String pass = data.getParent();
		Setting.set_ssm(Setting.settingReader(pass + "/" + args[5], "\t"));
		
		if (Setting.Permutation && args.length < argsLength + 1 ) {
			System.err.println("java -jar mainSSM.jar <data> <systemDimension> <repeatnumber> <id> <Drug.txt> <ssm.set> <output> <Phi> (<PhiG>)~for DRUG");
			System.exit(0);
		}

		Sfmt sfmt = new Sfmt(ID);
		TimeSeriesDataArray TimeSeriesData = new TimeSeriesDataArray(new Matrix());
		TimeSeriesData.Set(data.getPath(), "\t");
		if(Setting.Drug) {
			TimeSeriesData.readDrugProfiles(pass + "/" + args[4], pass+"/Z0.txt", "\t");
		}
		
        /*
         * Permutation Test
         */
		double[][] perPhi = null;
		double[][] perG = null;		
		if(Setting.Permutation){
	    	perPhi = Reader.ReadMatrixDouble(pass + "/" + args[7], "\t");
	    	if(Setting.Drug) {
	    		perG = Reader.ReadMatrixDouble(pass + "/" + args[8], "\t");
	    	}
		}
		
		System.out.println("Amount of data=" + TimeSeriesData.elementNum);
		System.out.println("Maximum Time=" + TimeSeriesData.maxTime);
		System.out.println("Replicate Size=" + TimeSeriesData.repSize);
		System.out.println("Amount of Time=" + TimeSeriesData.observationalTimeNum);
		System.out.println("Random Seed=" + ID);
		System.out.println("Permutation=" + Setting.Permutation);
		System.out.println("System Dimension=" + systemDimension);

		/*
		 * Thread Process
		 */
		ssmInference[] sINF = new ssmInference[(int) Setting.Thread];
		ssmStorage[] sSTO = new ssmStorage[(int) Setting.Thread];
		ssmSetting[] sSET = new ssmSetting[(int) Setting.Thread];
		TimeSeriesDataArray[] TSDA = new TimeSeriesDataArray[(int) Setting.Thread];
		Matrix[] CALCULATOR = new Matrix[(int) Setting.Thread];
		ssmMultiProcess[] sMP = new ssmMultiProcess[(int) Setting.Thread];
		Sfmt[] SFMT = new Sfmt[(int) Setting.Thread];
		for (int i = 0; i < (int) Setting.Thread; i++) {
			sSTO[i] = new ssmStorage();
			sSET[i] = new ssmSetting();
			sSET[i].copy(Setting);
			SFMT[i] = new Sfmt((int) (sfmt.NextUnif() * 1.0e4));
			CALCULATOR[i] = new Matrix();
			TSDA[i] = new TimeSeriesDataArray(CALCULATOR[i]);
			TSDA[i].copy(TimeSeriesData);
			if (Setting.SameDimensionAmongThreads){
				sSTO[i].initializeSSMStorage(systemDimension, sSET[i], TSDA[i], CALCULATOR[i], SFMT[i]);
				sINF[i] = new ssmInference(systemDimension, sSET[i], TSDA[i], CALCULATOR[i], SFMT[i], sSTO[i]);
			} else {
				sSTO[i].initializeSSMStorage(systemDimension + i, sSET[i], TSDA[i], CALCULATOR[i], SFMT[i]);
				sINF[i] = new ssmInference(systemDimension + i, sSET[i], TSDA[i], CALCULATOR[i], SFMT[i], sSTO[i]);
			}
			if(Setting.Permutation) sINF[i].PreparePermutation(perPhi, perG);

			sMP[i] = new ssmMultiProcess(repeatNumber, sINF[i], sSET[i].timer);
			sMP[i].start();
		}

		try {
			for (int i = 0; i < sMP.length; i++){
				sMP[i].join();
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		if(Setting.Permutation){
	        double success = 0;
	        double result[][] = new double[TimeSeriesData.elementNum][TimeSeriesData.elementNum];
	        for (int i = 0; i < sMP.length; i++){
	        	success += (double) sSTO[i].SuccessCount;
				sMatrix.sMP.add(result, sSTO[i].result);
	        }
			sMatrix.sMP.write(pass + "/" + args[6] + "/Success-"+systemDimension+"-"+ID+".txt", success);
			sMatrix.sMP.write(pass + "/" + args[6] + "/Permutation-"+systemDimension+"-"+ID+".txt", result);
			if(Setting.Drug){
		        result = new double[TimeSeriesData.elementNum][TimeSeriesData.drugMulRepSize[0].length];
		        for (int i = 0; i < sMP.length; i++) sMatrix.sMP.add(result, sSTO[i].resultDrug);{
		        	sMatrix.sMP.write(pass+"/PermutationDrug-"+systemDimension+"-"+ID+".txt", result);
		        }
			}
		} else {
			int optimum = 0;
			for (int i = 0; i < sSTO.length; i++){
				if (sSTO[i].Criterion < sSTO[optimum].Criterion) optimum = i;
			}

			if (Setting.SameDimensionAmongThreads) {
				for (int i = 0; i < Setting.Thread; i++) {
					if(sSTO[i].SuccessCount==0) {
						System.out.println("NO SUCCESS: " + i);
						continue;
					}
					File mkdir = new File(pass + "/" + args[6] + "/SD" + systemDimension + "Repeat" + repeatNumber + "ID" + ID + "GShift-"+Setting.Drug+"-"+i);
					if (mkdir.mkdir());
					else mkdir.mkdirs();
					String pass2 = mkdir.getPath();
					sMatrix.sMP.write(pass2 + "/Loglikelihood-" + systemDimension + "-" + ID+ ".txt", sSTO[i].logLikelihood);
					sMatrix.sMP.write(pass2 + "/BIC-"+ systemDimension + "-" + ID + ".txt", sSTO[i].Criterion);
					sMatrix.sMP.write(pass2 + "/Success-"+ systemDimension + "-" + ID + ".txt",(double) sSTO[i].SuccessCount);
					sMatrix.sMP.write(pass2 + "/F-"+ systemDimension + "-" + ID + ".txt", sSTO[i].sPa.getF());
					sMatrix.sMP.write(pass2 + "/H-"+ systemDimension + "-" + ID + ".txt", sSTO[i].sPa.getH());
					sMatrix.sMP.write(pass2 + "/R-"+ systemDimension + "-" + ID + ".txt", sSTO[i].sPa.getR());
					sMatrix.sMP.write(pass2 + "/x_p-"+ systemDimension + "-" + ID + ".txt", sSTO[i].getx_p());
					sMatrix.sMP.write(pass2 + "/x_f-"+ systemDimension + "-" + ID + ".txt", sSTO[i].getx_f());
					sMatrix.sMP.write(pass2 + "/x_s-"+ systemDimension + "-" + ID + ".txt", sSTO[i].getx_s());
					sMatrix.sMP.write(pass2 + "/y_p-"+ systemDimension + "-" + ID + ".txt", sSTO[i].gety_p());
					sMatrix.sMP.write(pass2 + "/ssmD-"+ systemDimension + "-" + ID + ".txt", sSTO[i].sPa.getssmD());
					sMatrix.sMP.write(pass2 + "/x0-"+ systemDimension + "-" + ID + ".txt", sSTO[i].sPa.getx0());
					sMatrix.sMP.write(pass2 + "/v0-"+ systemDimension + "-" + ID + ".txt", sSTO[i].sPa.getv0());
					sMatrix.sMP.write(pass2 + "/Phi-"+ systemDimension + "-" + ID + ".txt", sSTO[i].sPa.getPhi());
					if(Setting.Drug){
						Writer.write(pass2 + "/G-"+ systemDimension + "-" + ID + ".txt", sSTO[i].sPa.getG());
						sMatrix.sMP.write(pass2 + "/rootRinvHG-"+ systemDimension + "-" + ID + ".txt", sSTO[i].sPa.getrootRinvHG());
					}
				}
			} else {
				File mkdir = new File(pass + "/" + args[6] + "/SD" + systemDimension + "Repeat" + repeatNumber + "ID" + ID + "GShift-"+Setting.Drug);
				if (mkdir.mkdir());
				else mkdir.mkdirs();
				String pass2 = mkdir.getPath();
				String[] s = new String[sSTO.length];
				for (int i = 0; i < sSTO.length; i++) {
					int c = systemDimension + i;
					s[i] = "System:" + c + "\t" + sSTO[i].Criterion + "\t" + sSTO[i].logLikelihood;
				}
				Writer.write(pass2 + "/BICs-" + systemDimension + "-" + ID + ".txt", s);
			}
			System.out.println("BestLogLikelihood=" + sSTO[optimum].logLikelihood);
			System.out.println("BestBIC=" + sSTO[optimum].Criterion);
			System.out.println("Monotone=" + sSTO[optimum].Monotone);
			System.out.println("Convergence=" + sSTO[optimum].Convergence);
		}

		System.out.println("");
		System.out.println("program end "+ (System.currentTimeMillis() - timeAtStart) / 60 / 1000);
	}
}