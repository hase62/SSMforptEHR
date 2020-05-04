package Hasegawa.IO;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.io.File;

import Hasegawa.matrix.Matrix;
import Hasegawa.stat.simpleMath;
import RandomGenerator.Sfmt;

public class TimeSeriesDataArray {
	
	private Matrix Calculator;

	// TimeArray:(ex)TimeArray[1]={1,2,4,8,16,32},[2]={1,2,4,8,16,32},[3]={1,2,4,8}
	public ArrayList<ArrayList<Integer>> TimeArray = new ArrayList<ArrayList<Integer>>();
	protected ArrayList<ArrayList<Integer>> TimeArray_reserve = new ArrayList<ArrayList<Integer>>();
	
	// Validity:(ex)Validity[1]={true,true,false,...},[2]={true,true,false,...},[3]={false,false,true...}
	public ArrayList<ArrayList<Boolean>> Validity = new ArrayList<ArrayList<Boolean>>();

	// repID:(ex) if replicate is {1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3}, repID is
	// {6,6,4}
	public ArrayList<Integer> repID = new ArrayList<Integer>();

	// repIDSum:{0,6,12}
	public ArrayList<Integer> repIDSum = new ArrayList<Integer>();

	// repSize:3
	public int repSize = 0;

	// ObservationData[time][gene]
	public double[][] ObservationData;
	
	// ObservationData[time][gene]
	protected double[][] ObservationData_reserve;
	
	// observationalTimeNum:The number of time series
	public int observationalTimeNum = 0;

	// maxTime:Maximum time (ex) time(1,2,4,1,2,3,1,3),maxTime=4
	public int maxTime = 0;

	// allTime:maxTime * repID.size();
	public int allTime = 0;

	// elementNum:The amount of genes
	public int elementNum = 0;
	
	// The List of Gene Name
	public ArrayList<String> NameList = new ArrayList<String>();
	
	// The List of Gene Name.
	public String[] Name;
	
	// Z0
	public double[][] drug0_;
	
	// Z_n
	//public double[][] drugRaw;
	
	// Z_n
	public double[][] drugMulRepSize;
	
	// Given Regulation for F
	public double[][] givenF;
	
	// Given Regulation for G
	public double[][] givenG;
	
	public TimeSeriesDataArray(Matrix Calculator){
		this.Calculator = Calculator;
	}

	public void copy(TimeSeriesDataArray TSDA){
		
		this.TimeArray = new ArrayList<ArrayList<Integer>>(TSDA.TimeArray);
		if(TSDA.TimeArray_reserve != null) {
			this.TimeArray_reserve = new ArrayList<ArrayList<Integer>>(TSDA.TimeArray_reserve);
		}
		this.Validity = new ArrayList<ArrayList<Boolean>>(TSDA.Validity);	
		
		this.repID = new ArrayList<Integer>(TSDA.repID);
		this.repIDSum = new ArrayList<Integer>(TSDA.repIDSum);
		this.repSize = TSDA.repSize;

		this.ObservationData = new double[TSDA.observationalTimeNum][TSDA.elementNum];
		this.Calculator.copy(this.ObservationData, TSDA.ObservationData);
		if(TSDA.ObservationData_reserve!=null) {
			this.ObservationData_reserve = new double[TSDA.observationalTimeNum][TSDA.elementNum];
			this.Calculator.copy(this.ObservationData_reserve, TSDA.ObservationData_reserve);
		}

		this.observationalTimeNum = TSDA.observationalTimeNum;
		this.maxTime = TSDA.maxTime;
		this.allTime = TSDA.allTime;
		this.elementNum = TSDA.elementNum;
		
		this.NameList = new ArrayList<String>(TSDA.NameList);
		this.Name = (String[])this.NameList.toArray(new String[0]);
		
		if(TSDA.drugMulRepSize!=null){
			this.drugMulRepSize = new double[TSDA.drugMulRepSize.length][TSDA.drugMulRepSize[0].length];
			this.drug0_ = new double[TSDA.repSize][TSDA.drug0_[0].length];
			this.Calculator.copy(this.drugMulRepSize, TSDA.drugMulRepSize);
			this.Calculator.copy(this.drug0_, TSDA.drug0_);
		}
		
		if(TSDA.givenF!=null) {
			this.givenF = this.Calculator.copy_generate(TSDA.givenF);
		}
		if(TSDA.givenG!=null){
			this.givenG = this.Calculator.copy_generate(TSDA.givenG);
		}
	}
	
	/**
	 * Main routine of TimeSeriesDataArray
	 * Instead of using this method, you can use copy or makeSelectiongSample
	 */
	public void Set(String str, String tor) throws IOException {
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		ArrayList<ArrayList<Double>> ObservationData_list = new ArrayList<ArrayList<Double>>();
		String temp;
		while ((temp = fRead.readLine()) != null) {
			StringTokenizer st = new StringTokenizer(temp, tor);
			String st_temp = st.nextToken();
			
			/* Rep */
			if (st_temp.toLowerCase().equals("\"repid\"") 
					|| st_temp.toLowerCase().equals("@repid") || st_temp.toLowerCase().equals("rep_id")
						||st_temp.toLowerCase().equals("repid")||st_temp.toLowerCase().equals("rep")) {
				int count = 0;
				int pre = 1;
				repIDSum.add(0);
				while (st.hasMoreTokens() == true) {
					String nt = st.nextToken();					
					if ((int)Double.parseDouble(nt) == pre) count++;
					else {
						this.repID.add(count);
						if(pre!=1) repIDSum.add(repIDSum.get(repIDSum.size()-1)+repID.get(repID.size()-2));
						count = 1;
						pre = (int)Double.parseDouble(nt);
					}
				}
				this.repID.add(count);
				if(pre!=1) repIDSum.add(repIDSum.get(repIDSum.size()-1)+repID.get(repID.size()-2));
				this.repSize = this.repID.size();	
				continue;
			}
			
			/* Time */
			if (st_temp.toLowerCase().equals("\"time\"") 
					|| st_temp.toLowerCase().equals("@time") 
						|| st_temp.toLowerCase().equals("time")) {
				ArrayList<Integer> EachRepTime = new ArrayList<Integer>();
				ArrayList<Boolean> val = new ArrayList<Boolean>();
				int Time;
				String s =st.nextToken();
				Time = (int)Double.parseDouble(s);
				EachRepTime.add(Time);
				val.add(true);
				this.maxTime=0;
				this.observationalTimeNum=1;
				int preTime = Time;
				while(st.hasMoreTokens()) {					
						Time = (int)Double.parseDouble(st.nextToken());
						this.observationalTimeNum++;
						if(preTime>=Time) {
							this.TimeArray.add(EachRepTime);
							this.Validity.add(val);
							EachRepTime = new ArrayList<Integer>();
							val = new ArrayList<Boolean>();
						}
						EachRepTime.add(Time);
						val.add(true);
						preTime=Time;
						if (this.maxTime < Time) this.maxTime = Time;
				}
				this.TimeArray.add(EachRepTime);		
				this.Validity.add(val);
				continue;
			}
			
			this.NameList.add(st_temp);
			ArrayList<Double> TimeRow = new ArrayList<Double>();
			while(st.hasMoreTokens())
				TimeRow.add(Double.parseDouble(st.nextToken()));
			ObservationData_list.add(TimeRow);
			this.elementNum++;
		}
		this.Name = (String[])this.NameList.toArray(new String[0]);
		fInput.close();

		this.ObservationData = new double[observationalTimeNum][elementNum];
		for (int i = 0; i < this.observationalTimeNum; i++){
			for (int j = 0; j < this.elementNum; j++){
				this.ObservationData[i][j] = ObservationData_list.get(j).get(i);
			}
		}
		
		this.allTime = this.maxTime * this.repID.size();
	}
	
	/**
	 * Shorten TimeSeries (s)
	 * ex)1, 2, 3, 5, 8 -> 1, 3, 5, 9, 15 (s=2)
	 */
	public void shortenSimulationInterval(int s){
		
		ArrayList<ArrayList<Integer>> tempTimeArray = new ArrayList<ArrayList<Integer>>();
		for (int j = 0; j < this.repSize; j++) {
			ArrayList<Integer> temp = new ArrayList<Integer>();
			for (Integer integer : this.TimeArray.get(j)) {
				temp.add((integer*s)-(s-1));
			}
			tempTimeArray.add(temp);
		}
		this.TimeArray = new ArrayList<ArrayList<Integer>>(tempTimeArray);
		
		int previous_maxTime = this.maxTime;
		this.maxTime = (this.maxTime * s) - (s - 1);
		this.allTime = this.maxTime * repSize;
		
		if(this.drugMulRepSize!=null){
			double[][] tempDrug = new double[this.maxTime * this.repSize][this.drugMulRepSize[0].length];
			for (int rep = 0; rep < this.repSize; rep++) {
				for (int i = 0; i < previous_maxTime; i++) {
					int tPreviousStart = rep * previous_maxTime + i;
					int tPreviousEnd =  rep * previous_maxTime + i + 1;
					
					int tCurrentStart = rep * maxTime + i * s;
					for (int j = 0; j < tempDrug[i].length; j++) {
						if(i == previous_maxTime - 1){
							tempDrug[tCurrentStart][j] = this.drugMulRepSize[tPreviousStart][j];
							continue;
						}
						for (int n = 0; n < s; n++) {
							//if(i != this.drugMulRepSize.length - 1 && this.drugMulRepSize[i][j] == 1 && this.drugMulRepSize[i + 1][j] == 0) tempDrug[(i * s) + n][j] = 1;
							//else 
							tempDrug[tCurrentStart + n][j] = this.drugMulRepSize[tPreviousStart][j] - 
									(this.drugMulRepSize[tPreviousStart][j] - this.drugMulRepSize[tPreviousEnd][j]) * (double) n / s;
						}
					}
				}
			}
			this.drugMulRepSize = this.Calculator.copy_generate(tempDrug);
		}
	}
	
	/**
	 * Must MakeReserve before MakeBootstrapSample
	 * This method store observed data to Observation_reserve, and make Bootstrap samples to ObsevationData
	 */
	public void MakeBootstrapSample(Sfmt sf){
		if(this.ObservationData_reserve == null) this.ReserveData();
		ArrayList<Integer> sw = new simpleMath().getSwap(this.TimeArray.get(0).size(), sf);
		this.TimeArray = new ArrayList<ArrayList<Integer>>();
		int counter = 0;
		for (int i = 0; i < this.TimeArray_reserve.size(); i++) {
			ArrayList<Integer> ta = new ArrayList<Integer>();
			for (int j = 0; j < sw.size(); j++) {
				int n = this.TimeArray_reserve.get(0).get(sw.get(j));
				if(this.TimeArray_reserve.get(i).contains(n)) {
					this.Calculator.copy(this.ObservationData[counter], 
						this.ObservationData_reserve[this.repIDSum.get(i) 
						    + this.TimeArray_reserve.get(i).indexOf(n)]);
					ta.add(this.TimeArray_reserve.get(0).get(j));
					counter++;
				}
			}
			this.TimeArray.add(ta);
		}
	}
	
	public void ReserveData(){
		this.ObservationData_reserve = new double[observationalTimeNum][elementNum];
		this.Calculator.copy(this.ObservationData_reserve, this.ObservationData);
		this.TimeArray_reserve = new ArrayList<ArrayList<Integer>>(this.TimeArray);
	}
	
	public void nPredictionErrorDataSet(Sfmt sf, int n){
		ArrayList<Integer> Test = new ArrayList<Integer>();
		for (int i = 0; i < this.TimeArray.get(0).size(); i++) {
			//if(i % (int)(this.TimeArray.get(0).size() / n) 
					//== (this.TimeArray.get(0).size() - 1) % (int)(this.TimeArray.get(0).size() / n)) 
			if(i >= this.TimeArray.get(0).size() - n)	Test.add(this.TimeArray.get(0).get(i));
		}
		
		for (int i = 0; i < this.Validity.size(); i++) {
			for (int j = 0; j < this.Validity.get(0).size(); j++) {
				if(Test.contains(this.TimeArray.get(i).get(j))) this.Validity.get(i).set(j, false); 
			}
		}
	}
	
	public void setValidityTrue(){
		for (int i = 0; i < this.Validity.size(); i++) {
			for (int j = 0; j < this.Validity.get(0).size(); j++) {
				this.Validity.get(i).set(j, true); 
			}
		}
	}
	
	/**
	 * Read Drug Profile
	 */
	public void readDrugProfiles(String pass1, String pass2, String split) throws NumberFormatException, IOException{
		double[][] Read = Reader.ReadMatrixDouble(pass1, split);
		double[][] drugRaw = new double[Read[0].length][Read.length];
		this.Calculator.transpose(Read, drugRaw);
		File file = new File(pass2);
		if (file.exists()) {
			double[][] drug0Raw = Reader.ReadMatrixDouble(pass2, split);
			this.drug0_ = new double[this.repSize][drug0Raw.length];
			for (int r = 0; r < this.repSize; r++) {
				for (int i = 0; i < drug0Raw.length; i++) {
					if(drug0Raw[0].length < this.repSize) {
						this.drug0_[r][i] = drug0Raw[i][0];
					} else{
						this.drug0_[r][i] = drug0Raw[i][r];	
					}
				}
			}
		} else {
			System.out.println("Z0 is replaced by Z[1], ... ");
			this.drug0_ = new double[this.repSize][drugRaw[0].length];
			for (int r = 0; r < this.repSize; r++) {
				for (int i = 0; i < drugRaw[0].length; i++) {
					this.drug0_[r][i] = drugRaw[r * this.maxTime][i];
				}
			}
			for (int r = 0; r < this.repSize; r++) {
				for (int j = 0; j < this.maxTime - 1; j++) {
					for (int i = 0; i < drugRaw[0].length; i++) {
						drugRaw[r * this.maxTime + j][i] = drugRaw[r * this.maxTime + j + 1][i];
					}
				}
			}
		}
		this.drugMulExe(drugRaw);	
	}
	
	/**
	 * Make Drug Profile having multiple replicate
	 */
	protected void drugMulExe(double[][] drugRaw){
		if(drugRaw.length < this.repSize * this.maxTime) {
			this.drugMulRepSize = new double[this.maxTime * this.repSize][drugRaw[0].length];
			for (int r = 0; r < this.repSize; r++) {
				for (int i = 0; i < this.maxTime; i++) {
					for (int j = 0; j < drugRaw[0].length; j++) {
						this.drugMulRepSize[r * this.maxTime + i][j] = drugRaw[i][j];
					}
				}
			}
		} else {
			this.drugMulRepSize = this.Calculator.copy_generate(drugRaw);
		}
	}
	
	public void drugModifyVariance(double r){
		for (int i = 0; i < this.drugMulRepSize.length; i++) {
			for (int j = 0; j < this.drugMulRepSize[0].length; j++) {
				this.drugMulRepSize[i][j] = this.drugMulRepSize[i][j] * r;
			}
		}
		for (int i = 0; i < this.drug0_.length; i++) {
			for (int j = 0; j < this.drug0_[0].length; j++) {
				this.drug0_[i][j] = this.drug0_[i][j] * r;
			}
		}
	}
	
	public void addGaussianNoiseToObservation(double var, Sfmt s){
		for (int i = 0; i < this.ObservationData.length; i++) {
			for (int j = 0; j < this.ObservationData[0].length; j++) {
				this.ObservationData[i][j] += Math.sqrt(var) * s.NextNormal();
			}
		}
	}
	
	/*
	 * Only Linear Effect
	 * ex) y = 1.5a - 0.5b + 1
	 */
	public double[][] toRegressionInputData(Boolean DRUG, Boolean Input){
		int timePoints = this.observationalTimeNum - this.repSize;
		int option = 0;
		if(Input) option++;
		if(DRUG) option += this.drugMulRepSize[0].length;
		double[][] X = new double[timePoints][this.elementNum + option];
		int counter = 0;
		
		for (int i = 0; i < this.TimeArray.size(); i++) {
			for (int j = 0; j + 1 < this.TimeArray.get(i).size(); j++) {
				for (int j2 = 0; j2 < this.elementNum; j2++) {
					X[counter][j2] = this.ObservationData[this.repIDSum.get(i) + j][j2];
				}
				if(DRUG)  {
					for (int d = 0; d < this.drugMulRepSize[0].length; d++) {
						X[counter][this.elementNum + d] = this.drugMulRepSize[this.TimeArray.get(i).get(j) - 1][d];
					}
				}
				if(Input) X[counter][X[0].length - 1] = 1;
				counter++;
			}
		}
		return X;
	}
	
	/*
	 * Only Linear Effect
	 * ex) y = 1.5a - 0.5b + 1
	 */
	public double[][] toRegressionInputDataIndicatingColum(ArrayList<Integer> regulator, Boolean DRUG, Boolean Input){
		int timePoints = this.observationalTimeNum - this.repSize;
		int option = 0;
		if(Input) option++;
		if(DRUG) option += this.drugMulRepSize[0].length;
		double[][] X = new double[timePoints][regulator.size() + option];
		int counter = 0;
		
		for (int i = 0; i < this.TimeArray.size(); i++) {
			for (int j = 0; j + 1 < this.TimeArray.get(i).size(); j++) {
				for (int j2 = 0; j2 < regulator.size(); j2++) {
					X[counter][j2] = this.ObservationData[this.repIDSum.get(i) + j][regulator.get(j2)];
				}
				if(DRUG)  {
					for (int d = 0; d < this.drugMulRepSize[0].length; d++) {
						X[counter][regulator.size() + d] = this.drugMulRepSize[this.TimeArray.get(i).get(j) - 1][d];
					}
				}
				if(Input) X[counter][X[0].length - 1] = 1;
				counter++;
			}
		}
		return X;
	}
	
	/*
	 * Only Linear Effect
	 * ex) y = 1.5a - 0.5b + 1
	 */
	public double[][] toRegressionInputDataIndicatingColumWithAdditionalValues(ArrayList<Integer> regulator, Boolean DRUG, Boolean Input, double[] Self){
		int timePoints = this.observationalTimeNum - this.repSize;
		int option = 0;
		if(Input) option++;
		if(DRUG) option += this.drugMulRepSize[0].length;
		double[][] X = new double[timePoints][regulator.size() + option + 1];
		int counter = 0;
		
		for (int i = 0; i < this.TimeArray.size(); i++) {
			for (int j = 0; j + 1 < this.TimeArray.get(i).size(); j++) {
				for (int j2 = 0; j2 < regulator.size(); j2++) {
					X[counter][j2] = this.ObservationData[this.repIDSum.get(i) + j][regulator.get(j2)];
				}
				if(DRUG)  {
					for (int d = 0; d < this.drugMulRepSize[0].length; d++) {
						X[counter][regulator.size() + d] = this.drugMulRepSize[this.TimeArray.get(i).get(j) - 1][d];
					}
				}
				if(Input) X[counter][X[0].length - 2] = 1;
				X[counter][X[0].length - 1] = Self[this.repIDSum.get(i) + j];
				counter++;
			}
		}
		return X;
	}
	
	/*
	 * Include Cross Effect
	 * ex) y = 1.5a - 0.5b -0.15ab + 1
	 */

	public double[][] toRegressionInputData2(Boolean DRUG, Boolean Input){
		int timePoints = this.observationalTimeNum - this.repSize;
		int option = 0;
		int col = this.elementNum + (this.elementNum * (this.elementNum - 1) / 2);
		if(Input) option++;
		if(DRUG) option += this.drugMulRepSize[0].length;
		double[][] X = new double[timePoints][col + option];
		int counter = 0;
		
		for (int i = 0; i < this.TimeArray.size(); i++) {
			for (int j = 0; j + 1 < this.TimeArray.get(i).size(); j++) {
				for (int j2 = 0; j2 < this.elementNum; j2++) {
					X[counter][j2] = this.ObservationData[this.repIDSum.get(i) + j][j2];
				}
				int cross = this.elementNum;
				for (int c1 = 0; c1 < this.elementNum; c1++) {
					for (int c2 = c1+1; c2 < this.elementNum; c2++) {
						X[counter][cross] = this.ObservationData[this.repIDSum.get(i) + j][c1] 
								* this.ObservationData[this.repIDSum.get(i) + j][c2];
						cross++;
					}
				}
				if(DRUG)  {
					for (int d = 0; d < this.drugMulRepSize[0].length; d++) {
						X[counter][col + d] = this.drugMulRepSize[this.TimeArray.get(i).get(j) - 1][d];
					}
				}
				if(Input) X[counter][X[0].length - 1] = 1;
				counter++;
			}
		}
		return X;
	}
	
	/*
	 * Include Cross Effect
	 * ex) y = 1.5a - 0.5b -0.15ab + 1
	 */
	public double[][] toRegressionInputData2IndicatingColum(ArrayList<Integer> regulator, Boolean DRUG, Boolean Input){
		double[][] Obs = this.Calculator.getPartOfColVectors(this.ObservationData, regulator);
		int elementNum_ = Obs[0].length;
		int timePoints = this.observationalTimeNum - this.repSize;
		int option = 0;
		int col = elementNum_ + (elementNum_ * (elementNum_ - 1) / 2);
		if(Input) option++;
		if(DRUG) option += this.drugMulRepSize[0].length;
		double[][] X = new double[timePoints][col + option];
		int counter = 0;
		
		for (int i = 0; i < this.TimeArray.size(); i++) {
			for (int j = 0; j + 1 < this.TimeArray.get(i).size(); j++) {
				for (int j2 = 0; j2 < elementNum_; j2++) {
					X[counter][j2] = Obs[this.repIDSum.get(i) + j][j2];
				}
				int cross = elementNum_;
				for (int c1 = 0; c1 < elementNum_; c1++) {
					for (int c2 = c1+1; c2 < elementNum_; c2++) {
						X[counter][cross] = Obs[this.repIDSum.get(i) + j][c1] * Obs[this.repIDSum.get(i) + j][c2];
						cross++;
					}
				}
				if(DRUG)  {
					for (int d = 0; d < this.drugMulRepSize[0].length; d++) {
						X[counter][col + d] = this.drugMulRepSize[this.TimeArray.get(i).get(j) - 1][d];
					}
				}
				if(Input) X[counter][X[0].length - 1] = 1;
				counter++;
			}
		}
		return X;
	}
	
	/*
	 * Include Cross Effect
	 * ex) y = 1.5a - 0.5b -0.15ab + 1
	 */
	public double[][] toRegressionInputData2IndicatingColumWithAdditionalValues(ArrayList<Integer> regulator, Boolean DRUG, Boolean Input, double[] Self){
		double[][] Obs = this.Calculator.getPartOfColVectors(this.ObservationData, regulator);
		int elementNum_ = Obs[0].length;
		int timePoints = this.observationalTimeNum - this.repSize;
		int option = 0;
		int col = elementNum_ + (elementNum_ * (elementNum_ - 1) / 2);
		if(Input) option++;
		if(DRUG) option += this.drugMulRepSize[0].length;
		double[][] X = new double[timePoints][col + option + 1];
		int counter = 0;
		
		for (int i = 0; i < this.TimeArray.size(); i++) {
			for (int j = 0; j < this.TimeArray.get(i).size() - 1 ; j++) {
				for (int j2 = 0; j2 < elementNum_; j2++) {
					X[counter][j2] = Obs[this.repIDSum.get(i) + j][j2];
				}
				int cross = elementNum_;
				for (int c1 = 0; c1 < elementNum_; c1++) {
					for (int c2 = c1+1; c2 < elementNum_; c2++) {
						X[counter][cross] = Obs[this.repIDSum.get(i) + j][c1] * Obs[this.repIDSum.get(i) + j][c2];
						cross++;
					}
				}
				if(DRUG)  {
					for (int d = 0; d < this.drugMulRepSize[0].length; d++) {
						X[counter][col + d] = this.drugMulRepSize[this.TimeArray.get(i).get(j) - 1][d];
					}
				}
				if(Input) X[counter][X[0].length - 2] = 1;
				X[counter][X[0].length - 1] = Self[this.repIDSum.get(i) + j];
				counter++;
			}
		}
		return X;
	}
	
	public double[] toRegressionOutputData(int node){
		int timePoints = this.observationalTimeNum - this.repSize;
		double[] y = new double[timePoints];
		int counter = 0;
		
		for (int i = 0; i < this.TimeArray.size(); i++) {
			for (int j = 1; j < this.TimeArray.get(i).size(); j++) {
				y[counter] = this.ObservationData[this.repIDSum.get(i) + j][node];
				counter++;
			}
		}
		return y;
	}
	
	public double[] toRegressionOutputDataOnlyFluctuation(int node, double rate){
		int timePoints = this.observationalTimeNum - this.repSize;
		double[] y = new double[timePoints];
		int counter = 0;
		
		for (int i = 0; i < this.TimeArray.size(); i++) {
			for (int j = 1; j < this.TimeArray.get(i).size(); j++) {
				y[counter] = this.ObservationData[this.repIDSum.get(i) + j][node] - (rate * this.ObservationData[this.repIDSum.get(i) + j - 1][node]);
				counter++;
			}
		}
		return y;
	}

	public void changeDataForSteadyState() {
		int r = 0;
		double[] ss = new double[this.elementNum];
		for (int rep = 0; rep < this.repSize; rep++) {
			r += this.TimeArray.get(rep).size();
			for (int i = 0; i < ss.length; i++) {
				ss[i] += this.ObservationData[r-1][i] / this.repSize;
			}
		}
		for (int t = 0; t < this.ObservationData.length; t++) {
			for (int i = 0; i < this.ObservationData[0].length; i++) {
				this.ObservationData[t][i] -= ss[i];
			}
		}
	}
	
	public double[] calc_kurtosis(double nu_kurtosis) {
		double[] mean = new double[this.elementNum];
		double[] kurtosis = new double[this.elementNum];
		for (int t = 0; t < this.ObservationData.length; t++) {
			for (int i = 0; i < this.elementNum; i++) {
				mean[i] += this.ObservationData[t][i] / (double)this.ObservationData.length;
			}
		}
		for (int t = 0; t < this.ObservationData.length; t++) {
			for (int i = 0; i < this.elementNum; i++) {
				kurtosis[i] += Math.pow(this.ObservationData[t][i] - mean[i], 4) / (double)this.ObservationData.length;
			}
		}
		for (int i = 0; i < kurtosis.length; i++) {
			if(kurtosis[i]< 4) {
				kurtosis[i] = nu_kurtosis;
			} else {
				kurtosis[i] = nu_kurtosis / Math.sqrt(kurtosis[i] - 3);				
				if(kurtosis[i] < 2) kurtosis[i] = 2;
			}
		}
		return kurtosis;
	}
}