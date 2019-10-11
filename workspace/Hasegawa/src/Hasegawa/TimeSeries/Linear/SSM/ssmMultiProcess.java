package Hasegawa.TimeSeries.Linear.SSM;

public class ssmMultiProcess extends Thread {

	protected int repeatNum;
	protected ssmInference sInf;
	protected double timer = 0;
	
	public ssmMultiProcess(int repeat, ssmInference inference, double time) {
		this.repeatNum = repeat;
		this.sInf = inference;
		this.timer = time;
	}
	
	public void run() {
		final long timeAtStart = System.currentTimeMillis();
		for (int i = 0; i < this.repeatNum && (System.currentTimeMillis() - timeAtStart < this.timer);) {
			this.sInf.run();
			if (!Double.isInfinite(this.sInf.getLoglikelihood()) && !Double.isNaN(this.sInf.getLoglikelihood())) {
				if(this.sInf.storeBestParameters(false)) continue;
				i++;
				System.out.println("Count: "+i);
			}
		}
		this.sInf.run_CSSM();
		this.sInf.storeBestParameters(true);
		System.out.println("End");
	}
}