package Hasegawa.TimeSeries.Linear.VARSSM;

import Hasegawa.IO.TimeSeriesDataArray;
import Hasegawa.TimeSeries.Storage;
import Hasegawa.matrix.Matrix;
import RandomGenerator.Sfmt;

public class vsStorage extends Storage{
	
	public vsParameter vPa;
	double nonPenalizedLogLikelihood;
	
	public void initializeVSStorage(int sysDim, vsSetting set, TimeSeriesDataArray tsda, Matrix Calculator, double initial, Sfmt Sf){
		this.vPa = new vsParameter(sysDim, Sf, set, tsda, Calculator, initial);
		this.Initialize(sysDim, tsda, set);
	}
}
