package fem;


public class TimeFunction {
	
	public TimeFunction(int id1, double amp,double per,double phase1){
		id=id1;
		amplitude=amp;
		period=per;
		phase=phase1;
	}

	public int id;
	public double amplitude,period,phase;
	
	public double getValue(double t){
		
		double val= amplitude*Math.cos(2*Math.PI/period*t+phase*Math.PI/180);
		
		return val;
	}

}