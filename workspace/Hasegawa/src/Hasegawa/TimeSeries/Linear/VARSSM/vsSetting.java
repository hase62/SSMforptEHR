package Hasegawa.TimeSeries.Linear.VARSSM;

import java.util.ArrayList;

import Hasegawa.TimeSeries.Setting;

public class vsSetting extends Setting{

	public double WeightLASSO = 0;
	public boolean Spacom = true;
	public boolean Separate=true;
	public double testPredictionAbility = 0;
	
	public void set_vs(final ArrayList<Double> s){
		this.Mean_of_Mu0 = s.get(0);
		this.Var_of_Mu0 = s.get(1);
		this.Condition_of_Convergence = s.get(2);
		this.maxLoop = s.get(3);
		this.Mu0_Update = (int)(double)s.get(4);
		
		if(s.get(5)==0.0) this.Print_Progress =true;
		this.timer =  (long)(s.get(6) * 60 * 1000);
		if(this.timer==0) this.timer = Long.MAX_VALUE;
		//this.Thread = (int)(double)s.get(7);
		if(s.get(7)==0.0) this.Drug = true;
		this.R_rI = s.get(8);
		if(s.get(9)==0.0) this.Input = true;
		
		this.Degradation = s.get(10);
		this.Criterion = s.get(11);
		this.shortenInterval =s.get(12);
		this.upDateQ = s.get(13);
		if(s.get(14)==0.0) this.GivenRegulation = true;
		
		this.testPredictionAbility = s.get(15);
		if(s.get(16)==0.0) this.LongSteadyState = true;
		this.systemDimension = (int)(double) s.get(17);
	}
	
	public void copy(vsSetting SET){
    	this.Mean_of_Mu0 = SET.Mean_of_Mu0;
    	this.Var_of_Mu0 = SET.Var_of_Mu0;
    	this.Condition_of_Convergence = SET.Condition_of_Convergence;
    	this.maxLoop = SET.maxLoop;
    	this.R_rI = SET.R_rI;
    	this.Mu0_Update = SET.Mu0_Update;
    	this.Print_Progress = SET.Print_Progress;
    	this.timer = SET.timer;
    	//this.Thread = SET.Thread;
		this.Drug = SET.Drug;
		this.Input = SET.Input;
		this.Degradation = SET.Degradation;
		this.Criterion = SET.Criterion;
		//this.WeightLASSO = SET.WeightLASSO;
		//this.Spacom = SET.Spacom;
		//this.appropriateOrder =SET.appropriateOrder;
		this.shortenInterval = SET.shortenInterval;
		this.upDateQ = SET.upDateQ;
		//this.Separate = SET.Separate;
		this.GivenRegulation = SET.GivenRegulation;
		this.testPredictionAbility = SET.testPredictionAbility;
		this.LongSteadyState = SET.LongSteadyState;
	}
}
