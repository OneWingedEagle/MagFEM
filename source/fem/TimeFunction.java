package fem;

import math.util;


public class TimeFunction {
	
	public TimeFunction(int id1, double amp,double per,double phase1){
		id=id1;
		amplitude=amp;
		period=per;
		phase=phase1;
		type=0;
	}

	public TimeFunction(int id1, double a0,double per,double a1, double a2,double a3){
		id=id1;
		Cdc=a0;
		period=per;
		Ccos=a1;
		Csin=a2;
		Texp=a3;
		type=1;
	}

	
	public int id,type;
	public double amplitude,period,phase;
	public double Cdc,Ccos,Csin,Cexp,Texp;
	
	public double getValue(double t){
		
		double val=0;
		
		if(type==0)
			val=amplitude*Math.cos(2*Math.PI/period*t+phase*Math.PI/180);
		else if(type==1){
			double wt=2*Math.PI/period*t;
			val=Cdc+Math.exp(Texp*t)*(Ccos*Math.cos(wt)+Csin*Math.sin(wt));
		}

		return val;
	}

}