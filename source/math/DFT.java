package math;
import  static java.lang.Math.*;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import io.Loader;
public class DFT {
	
  public static void main(String argv[]) {
	  String arrFile="arr7.txt";
	  int N=1800;
	  int J=7;

	  double[][] arr=new Loader().loadArrays(N,J,arrFile);
	  
	  double[][] a=new Mat(arr).transp().el;
	  
	  Complex[][] Y=new Complex[J][N];
	  for(int j=0;j<J;j++){

	 Y[j]=dft(a[j]);
	  }

	 // N=1;
	  

	  for(int k=0;k<N;k++){
		  for(int j=0;j<J;j++){
			  System.out.print(Y[j][k].norm()/N+"\t");
		  }
		  util.pr("");
	  }
	 // util.pr("");
/*	  util.pr(Y.length);
util.pr(s);*/


}
  
  public static void main33(String argv[]) {
	  

	  double[] a=new Loader().loadArray();
	//  double[] b=new double[a.length/2];
/*	  for(int k=0;k<b.length;k++){
			b[k]=a[b.length+k];
			  }*/
	 Complex[] Y=dft(a);

	 int N=Y.length;

	// N=1;
	 double s=0;
	  for(int k=0;k<Y.length;k++){
	 util.pr(Y[k].norm()/N);
	 if( k<100)
		 //s+=Math.pow(b[k],2);
	 s+=Math.pow(Y[k].norm(),2);
	  }

	 // util.pr("");
/*	  util.pr(Y.length);
util.pr(s);*/


}
  
  public static Complex[] dft(double[] x){
	  int N=x.length;
	  Complex[] X=new Complex[N];
	  
	  for(int k=0;k<N;k++){
		  X[k]=new Complex();
		  for(int n=0;n<N;n++){
			  double c=cos(2*PI*n*k/N);
			  double s=-sin(2*PI*n*k/N);
			  Complex z=new Complex(c,s);
			  X[k]=X[k].add(z.times(x[n]));
		  }
		  
	  }
	  
	  return X;
	  
  }
  
  public static Complex[] dft(Complex[] x){
	  int N=x.length;
	  Complex[] X=new Complex[N];
	  
	  for(int k=0;k<N;k++){
		  X[k]=new Complex();
		  for(int n=0;n<N;n++){
			  double c=cos(2*PI*n*k/N);
			  double s=-sin(2*PI*n*k/N);
			  Complex z=new Complex(c,s);
			  X[k]=X[k].add(z.times(x[n]));
		  }
		  
	  }
	  
	  return X;
	  
  }
  
  public static Complex[] idft(Complex[] X){
	  int N=X.length;
	  double rN=1.0/N;
	  Complex[] x=new Complex[N];
	  
	  for(int n=0;n<N;n++){
		 x[n]=new Complex();
		  for(int k=0;k<N;k++){
			  double c=cos(2*PI*n*k/N);
			  double s=sin(2*PI*n*k/N);
			  Complex z=new Complex(c,s);
			  x[n]=x[n].add(z.times(X[k]));
		  }
		  x[n]= x[n].times(rN);
		  
	  }
	  
	  return x;
	  
  }
  
	
	public static double[][] fftr(double[] x)
	{
	
		double pi=4*atan(1);
		int N,n,k;

		N=x.length;
		
		double X[][]= new double[N][2];

		if(N==1)
		{
			
			X[0][0]=x[0];
			X[0][1]=0;
		}
			else {
		
			double[] xe=new double[N/2];	
			double[] xo=new double[N/2];
			double[][] Xe= new double[N/2][2];
			double[][] Xo= new double[N/2][2];
			
			for(n=0;n<N/2;n++)
			{
				xe[n]=x[2*n];
				xo[n]=x[2*n+1];

			}
	
			Xe=fftr(xe);

			Xo=fftr(xo);

		    for(k=0;k<N/2;k++){
		    X[k][0]=Xe[k][0]+cos(2*pi*k/N)*Xo[k][0]+sin(2*pi*k/N)*Xo[k][1];    
		    X[k][1]=Xe[k][1]+cos(2*pi*k/N)*Xo[k][1]-sin(2*pi*k/N)*Xo[k][0];
		    X[k+N/2][0]=Xe[k][0]-cos(2*pi*k/N)*Xo[k][0]-sin(2*pi*k/N)*Xo[k][1];    
		    X[k+N/2][1]=Xe[k][1]-cos(2*pi*k/N)*Xo[k][1]+sin(2*pi*k/N)*Xo[k][0];	
		    	
		    }

		    }
			


		return X;

	}
	
// Class FFTC: FFT of a complex vector
	
	public static double[][] fftc(double[][] x)
	{
	
		double pi=4*atan(1);
		int N,n,k;

		N=x.length;
	//	System.out.println("****  "+x.length);
		double X[][]= new double[N][2];
		
		
		
		if(N==1)
		{
			
			X=x;
		}
			else {
		
			double[][] xe=new double[N/2][2];	
			double[][] xo=new double[N/2][2];
			double[][] Xe= new double[N/2][2];
			double[][] Xo= new double[N/2][2];
			
			for(n=0;n<N/2;n++)
			{
				xe[n]=x[2*n];
				xo[n]=x[2*n+1];

			}
	
			Xe=fftc(xe);

			Xo=fftc(xo);

		    for(k=0;k<N/2;k++){
		    X[k][0]=Xe[k][0]+cos(2*pi*k/N)*Xo[k][0]+sin(2*pi*k/N)*Xo[k][1];    
		    X[k][1]=Xe[k][1]+cos(2*pi*k/N)*Xo[k][1]-sin(2*pi*k/N)*Xo[k][0];
		    X[k+N/2][0]=Xe[k][0]-cos(2*pi*k/N)*Xo[k][0]-sin(2*pi*k/N)*Xo[k][1];    
		    X[k+N/2][1]=Xe[k][1]-cos(2*pi*k/N)*Xo[k][1]+sin(2*pi*k/N)*Xo[k][0];	
		    	
		    }

		    }
			


		return X;

	}
	
	public static double[] ifft(double[][] X, int L)
	{
		double[][] nx=new double[L][2];
		double[] x=new double[L];
	    double Ld=L; 
			nx=nifft(X);
	for(int i=0;i<L;i++)
		x[i]=nx[i][0]/Ld;
	
	return x;	
	}
	
//Class  CIFFT:  Inverse FFT with a complex output vector
	
	public static double[][] cifft(double[][] X, int L)
	{
		double[][] nx=new double[L][2];
		double[][] x=new double[L][2];
		nx=nifft(X);
		double Ld=L;
	for(int i=0;i<L;i++){
		x[i][0]=nx[i][0]/Ld;
		x[i][1]=nx[i][1]/Ld;
	}
	return x;	
	}
	
// Class NIFT: Auxiliary inverse FFT with a complex output
	
	public static double[][] nifft(double[][] X)
	{
	
		double pi=4*atan(1);
		int N,n;

		N=X.length;

		double x[][]= new double[N][2];

		
		if(N==1)
		{
			
			x=X;
		}
			else {
				
			double[][] Xe= new double[N/2][2];
			double[][] Xo= new double[N/2][2];
			double[][] xe=new double[N/2][2];	
			double[][] xo=new double[N/2][2];
			
			
			for(n=0;n<N/2;n++)
			{
				Xe[n]=X[2*n];
				Xo[n]=X[2*n+1];

			}
	
			xe=nifft(Xe);

			xo=nifft(Xo);

		    for(n=0;n<N/2;n++){
		    	x[n][0]=xe[n][0]+cos(2*pi*n/N)*xo[n][0]-sin(2*pi*n/N)*xo[n][1];    
			    x[n][1]=xe[n][1]+cos(2*pi*n/N)*xo[n][1]+sin(2*pi*n/N)*xo[n][0];
			    x[n+N/2][0]=xe[n][0]-cos(2*pi*n/N)*xo[n][0]+sin(2*pi*n/N)*xo[n][1];    
			    x[n+N/2][1]=xe[n][1]-cos(2*pi*n/N)*xo[n][1]-sin(2*pi*n/N)*xo[n][0];
		    	
		    }

		    }
			


		return x;

	}
	
	
	 public static Complex[] fft(double [] x){
		  int N=x.length;
		  
		  double[] x2=new double[2*N];
		  for(int k=0;k<N;k++)
			  x2[2*k]=x[k];
		  
		  DoubleFFT_1D X1= new  DoubleFFT_1D(N);
		  X1.complexForward(x2);
		
		  Complex[] X=new Complex[N];
		  
		  for(int k=0;k<N;k++){
			  X[k]=new Complex(x2[2*k],x2[2*k+1]);
}
		  
		  return X;
		  
	  }

	 public static Complex[] ifft(Complex [] X){
		  int N=X.length;
		  
		  double[] X1=new double[2*N];
		  for(int k=0;k<N;k++){
			  X1[2*k]=X[k].re;
			  X1[2*k+1]=X[k].im;
		  }
		  

		  DoubleFFT_1D fft1= new  DoubleFFT_1D(N);
		  fft1.complexInverse(X1,true);
		
		  Complex[] x=new Complex[N];
		  
		  for(int k=0;k<N;k++){
			  x[k]=new Complex(X1[2*k],X1[2*k+1]);
}
		  
		  return x;
		  
	  }
	 
		public static Complex[][] FFT2D(double[][] x){
			
			int I=x.length;
			int J=x[0].length;
			
			double[][] X=new double[I][2*J];
			for(int i=0;i<I;i++)
				for(int j=0;j<J;j++){
					X[i][2*j]=x[i][j];
				
			}

			
			DoubleFFT_2D ft2=new DoubleFFT_2D(I,J);
			ft2.complexForward(X);
			
			Complex[][] Y=new Complex[I][J];
			
			for(int i=0;i<I;i++)
				for(int j=0;j<J;j++){
					Y[i][j]=new Complex(X[i][2*j],X[i][2*j+1]);
				
			}
			return Y;
		}
		
		public static Complex[][] iFFT2D(Complex[][] x){
			
			int I=x.length;
			int J=x[0].length;
			
			double[][] X=new double[I][2*J];
			for(int i=0;i<I;i++)
				for(int j=0;j<J;j++){
					X[i][2*j]=x[i][j].re;
					X[i][2*j+1]=x[i][j].im;
				
			}

			DoubleFFT_2D ft2=new DoubleFFT_2D(I,J);
			ft2.complexInverse(X,true);
			
			
			Complex[][] Y=new Complex[I][J];
			
			for(int i=0;i<I;i++)
				for(int j=0;j<J;j++){
					Y[i][j]=new Complex(X[i][2*j],X[i][2*j+1]);
				
			}
			return Y;
		}
}
