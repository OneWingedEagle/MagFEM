package math;
import java.util.Arrays;
import static java.lang.Math.*;
public class Tajimi {

	
public Tajimi(){}
	
public static void main(String[] args){
	
	double D=50;
	
	double F=20e3;
	
	Vect x=new Vect().linspace(0, 10, 100);
	
	Vect y=new Vect(x.length);
	for(int i=0;i<y.length;i++)
		y.el[i]=Math.sqrt(3)/2*x.el[i];
	

	//util.plot(x,y);

}
	public static void main2(String[] args) throws Exception{
		

		int K=400;
		double[] f=new double[K+1];
		double[] fi=new double[K];
		double[] fb=new double[K];
		double[] SK=new double[K];
		double[] A=new double[K];
		double[] tt=new double[K];
		
		f[K]=7;
		
		double S0=94.76*4*Math.PI;
		double fg=2.92;
		double kesaig=0.34;
		double kes2=pow(kesaig,2);
		double T=40;
		int N=400;
		int nn=300;
		
		double[] aa=new double[K];
		double[][] ff=new double[K][nn];
	
		//monte carlo simulation %%%%%%%%%%%%%
		double dt=T/N;
		double fs =1/dt;


		for(int jj=0;jj<nn;jj++){
		
			for(int k=1;k<K;k++)	
			    f[k]= f[0] + (f[K]- f[0])*Math.random();
			Arrays.sort(f);
			
		for(int k=0;k<K;k++){	
		    fi[k]=2*Math.PI*Math.random();
		    fb[k]=Math.sqrt(f[k]*f[k+1]);
		    double xx=pow(fb[k]/fg,2);
		    SK[k]=S0*(1+4*kes2*xx)/(pow(1-pow(fb[k]/fg,2),2)+4*kes2*xx);
		    A[k]=sqrt(2*SK[k]*(f[k+1]-f[k]));
		}
		
		for(int j=0;j<N;j++){
		    double x=0;
		   double t=j*dt;
		   tt[j]=t;
			for(int k=0;k<K;k++)
		    x=x+A[k]*sin(2*PI*fb[k]*t+fi[k]);
  
			aa[0]=0;
		    aa[j]=x;
	
		}
		   
		
		//util.hshow(aa);
		for(int k=0;k<K;k++)
		ff[k][jj]=aa[k];


		}
		
		util.plot(tt,aa);
		
		//util.show(ff);
		//xlswrite('ag.xls',ff,'b2:ko4000')




	}
}
