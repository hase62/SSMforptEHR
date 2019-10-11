package Hasegawa.TimeSeries;

import Hasegawa.IO.SettingReader;

public class Setting extends SettingReader{
    
	public double Mean_of_Mu0 = 0;
	public double Var_of_Mu0 = 1;
	public int Mu0_Update = 0;
	public boolean Print_Progress = false;
	public long timer = Integer.MAX_VALUE;
			
	public int Thread = 1;
	public double R_rI = 1;
	public double Criterion=0;
	public double shortenInterval = 1;
	public boolean Drug = false;
	
	public boolean Input = false;
	public double Condition_of_Convergence = 1.0E-3;
	public double maxLoop = 300000;
	public double upDateQ=0.0;
	public int systemDimension = 1000;
	
	public double Degradation=0.0;
	public boolean GivenRegulation = false;
	public boolean LongSteadyState = false;
}