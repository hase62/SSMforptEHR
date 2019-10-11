package Hasegawa.TimeSeries.Linear.SSM;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.TimeSeries.Storage;
import Hasegawa.matrix.Matrix;
import RandomGenerator.Sfmt;

public class ssmStorage extends Storage{
		
	public ssmParameter sPa;
	public double[][] perPhi;
	public double[][] result;
	public double[][] perDrug;
    public double[][] resultDrug;
    public boolean Monotone;
    public boolean Convergence;
    public int SuccessCount = 0;
	 
	public void initializeSSMStorage(int sysDim, ssmSetting set, TimeSeriesDataArray tsda, Matrix Calculator, Sfmt Sf){
		this.Criterion = Double.MAX_VALUE;
		this.logLikelihood = (-1)*Double.MAX_VALUE;
		this.SuccessCount = 0;
		
		this.sPa = new ssmParameter(sysDim, Sf, set, tsda, Calculator);
		this.Initialize(sysDim, tsda, set);
	}
    
}
