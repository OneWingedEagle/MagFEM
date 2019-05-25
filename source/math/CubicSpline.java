package math;

import static java.lang.Math.pow;
public class CubicSpline {

	CubicSpline(){
		
	}
	
public static void main(String[] args) throws Exception{
		
	CubicSpline csp=new CubicSpline ();
	Vect x=new Vect();
	x=x.linspace(1, 30,30);
	double[][] xy=new double[x.length][2];
	for(int i=0;i<x.length;i++)
	{
		xy[i][0]=x.el[i];
		xy[i][1]=.5*pow(x.el[i],3)+2*pow(x.el[i],2)-x.el[i]+4;
		System.out.println(xy[i][1]);
		
	}
	//util.show(xy);
	double[][] Ci=csp.getCoefs(xy);
	util.show(Ci);
	Vect vx=new Vect();
	vx=vx.linspace(7, 8, 10);
	Vect vy=new Vect(vx.length);
	for(int i=0;i<vx.length;i++)
	{
		
		vy.el[i]=Ci[7][0]*pow(vx.el[i]-7,3)+Ci[7][1]*pow(vx.el[i]-7,2)+Ci[7][2]*(vx.el[i]-7)+Ci[7][3];

		
	}
	vx.show();
	vy.show();
		
	}

	public double[][] getCoefs(double[][] xy){
		int I=xy.length-1;
		double[][] coefs=new double[I][4];
		int L=I-1;
		double[] h=new double[I];
		for(int i=0;i<I;i++)
			h[i]=xy[i+1][0]-xy[i][0];

		Vect b=new Vect(L);
		for(int i=0;i<b.length;i++)
			b.el[i]=6*(xy[i][1]+xy[i+2][1]-2*xy[i+1][1])/Math.pow(h[i],2);

		SpMat A=new SpMat(L);
			A.row[0]=new SpVect(L,2);
			A.row[0].el[0]=4;
			A.row[0].index[0]=0;
			A.row[0].el[1]=1;
			A.row[0].index[1]=1;
			for(int i=1;i<L-1;i++){
				A.row[i]=new SpVect(L,3);
				A.row[i].el[0]=1;
				A.row[i].el[1]=4;
				A.row[i].el[2]=1;
				A.row[i].index[0]=i-1;
				A.row[i].index[1]=i;
				A.row[i].index[2]=i+1;
			}
			A.row[L-1]=new SpVect(L,2);
			A.row[L-1].el[0]=1;
			A.row[L-1].index[0]=L-2;
			A.row[L-1].el[1]=4;
			A.row[L-1].index[1]=L-1;
			
			double error=1e-6;
			int N=1000;
			
			SpMatSolver solver=new SpMatSolver();
			Vect M1=solver.CG(A, b, error, N,new Vect(L));
			double[] M=new double[xy.length];
			M[0]=M1.el[0]; M[xy.length-1]=M1.el[L-1];
			for(int i=0;i<L;i++){
				M[i+1]=M1.el[i];
			}

			for(int i=0;i<I;i++){
				coefs[i][0]=(M[i+1]-M[i])/(6*h[i]);
				coefs[i][1]=M[i]/2;
				coefs[i][2]=(xy[i+1][1]-xy[i][1])/h[i]-h[i]*(M[i+1]+2*M[i])/6;
				coefs[i][3]=xy[i][1];
			}
			
			
		return coefs;
	}
}
