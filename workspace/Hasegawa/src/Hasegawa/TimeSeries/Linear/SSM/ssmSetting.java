package Hasegawa.TimeSeries.Linear.SSM;

import java.util.ArrayList;
import Hasegawa.TimeSeries.Setting;

public class ssmSetting extends Setting{
	
	public double Mean_of_H = 0;
	public double SD_of_H = 1;
	public double Lower_bound_of_F = -1.5;
	public double UPPer_bound_of_F = 1.5;
	
	public boolean Absolute_Convergence = false;
	public boolean Permutation = false;
	public boolean SameDimensionAmongThreads = false;
	
	public void set_ssm(final ArrayList<Double> s){	
    	this.Mean_of_H = s.get(0);
    	this.SD_of_H = s.get(1);
    	this.Lower_bound_of_F = s.get(2);
    	this.UPPer_bound_of_F = s.get(3);
    	this.Mean_of_Mu0 = s.get(4);
    	this.Var_of_Mu0 = s.get(5);
    	this.Condition_of_Convergence = s.get(6);
    	this.maxLoop = s.get(7);
    	this.R_rI = s.get(8);
    	this.Mu0_Update = (int)(double)s.get(9);
    	if(s.get(10)==0.0) this.Print_Progress =true;
    	if(s.get(11)==0.0) this.Absolute_Convergence=true;
    	if(s.get(12)==0.0) this.Permutation = true;
		this.timer =  (long)(s.get(13) * 60 * 1000);
		if(this.timer==0) this.timer = Long.MAX_VALUE;
		this.Thread = (int)(double)s.get(14);
    	if(s.get(15)==0.0) this.SameDimensionAmongThreads =true; 
    	if(s.get(16)==0.0) this.Drug = true;
    }
    
	public void copy(ssmSetting SET){
    	
		this.Mean_of_H = SET.Mean_of_H;
    	this.SD_of_H = SET.SD_of_H;
    	this.Lower_bound_of_F = SET.Lower_bound_of_F;
    	this.UPPer_bound_of_F = SET.UPPer_bound_of_F;
    	this.Mean_of_Mu0 = SET.Mean_of_Mu0;
    	
    	this.Var_of_Mu0 = SET.Var_of_Mu0;
    	this.Condition_of_Convergence = SET.Condition_of_Convergence;
    	this.maxLoop = SET.maxLoop;
    	this.R_rI = SET.R_rI;
    	this.Mu0_Update = SET.Mu0_Update;

    	this.Print_Progress = SET.Print_Progress;
    	this.Absolute_Convergence = SET.Absolute_Convergence;
    	this.Permutation = SET.Permutation;
    	this.timer = SET.timer;
    	
    	this.Thread = SET.Thread;
    	this.SameDimensionAmongThreads = SET.SameDimensionAmongThreads; 
		this.Drug = SET.Drug;
	}
}