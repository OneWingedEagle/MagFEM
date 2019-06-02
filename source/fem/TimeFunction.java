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

	public TimeFunction(int id1, double a0,double a1,double per,double a2, double a3,double a4){
		id=id1;
		Cdc=a0;
		Cramp=a1;
		period=per;
		Ccos=a2;
		Csin=a3;
		Texp=a4;
		type=1;
	}

	
	public int id,type;
	public double amplitude,period,phase;
	public double Cdc,Cramp,Ccos,Csin,Cexp,Texp;
	
	public double getValue(double t){
		
		double val=0;
		
		if(type==0)
			val=amplitude*Math.cos(2*Math.PI/period*t+phase*Math.PI/180);
		else if(type==1){
			double wt=2*Math.PI/period*t;
			val=Cdc+Cramp*t+Math.exp(Texp*t)*(Ccos*Math.cos(wt)+Csin*Math.sin(wt));
		}

		return val;
	}

}