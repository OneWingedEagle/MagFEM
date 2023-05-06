package math;
import java.util.Arrays;
import static java.lang.Math.*;
public class TajimiNonst {

	
public TajimiNonst(){}
	
	public static void main(String[] args) throws Exception{
		

		int K=40;
		double[] f=new double[K+1];
		double[] fi=new double[K];
		double[] fb=new double[K];
		double[] SK=new double[K];
		double[] A=new double[K];
		
		
		f[K]=7;
		
		double S0=94.76*4*Math.PI;
		double fg=2.92;
		double kesaig=0.34;
		double kes2=pow(kesaig,2);
		double T=40;
		
		int N=400;
		int nn=30;
		
		double a1=0.085;
		double b1=0.17;
		a1=0;
		b1=1;
		double c=pow(a1/b1,a1/(b1-a1))-pow(a1/b1,b1/(b1-a1));
		double[] aa=new double[N];
		double[] a=new double[N];
		double[][] ff=new double[K][nn];
		double[] tt=new double[N];
	
		//monte carlo simulation %%%%%%%%%%%%%
		double dt=T/N;
		double fs =1/dt;


		for(int jj=0;jj<nn;jj++){
		
			for(int k=1;k<K;k++)	
			    f[k]= f[0] + (f[K]- f[0])*Math.random();
			Arrays.sort(f);
			
		for(int k=0;k<K;k++){	
		    fi[k]=2*PI*Math.random();
		    fb[k]=Math.sqrt(f[k]*f[k+1]);
		    double xx=pow(fb[k]/fg,2);
		    SK[k]=S0*(1+4*kes2*xx)/(pow(1-pow(fb[k]/fg,2),2)+4*kes2*xx);
		    A[k]=sqrt(2*SK[k]*(f[k+1]-f[k]));
		}
		
		for(int j=0;j<N;j++){
		    double x=0;
		   double t=j*dt;
		   tt[j]=t;
		   double ck=(exp(-a1*t)-exp(-b1*t))/c;


			for(int k=0;k<K;k++){
		    x=x+A[k]*sin(2*PI*fb[k]*t+fi[k]);

				}
					    
		    x*=ck;
		    
		    aa[0]=0;
		    aa[j]=x;
		    a[j]=x;
	
		}
		   
		
		//util.hshow(aa);
		for(int k=0;k<K;k++)
		ff[k][jj]=aa[k];


		}
		
		util.plot(tt,a);
		
		//util.show(ff);
		//xlswrite('ag.xls',ff,'b2:ko4000')




	}
}
