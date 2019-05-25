package emsolution;

import java.util.Random;

import static java.lang.Math.*;
import math.Mat;
import math.Vect;
import math.util;

public class Preisach1D {

	public Random r;
	public int M,kk,dim;
	public long seed;
	public double cfm,cfw,mean,width,Hs,Bs,Bseff,DB;
	public boolean[] on;
	public double[][] a;
	public double[] K;
	public double  epsdB=1e-2;
	public double  epsdBdH=1e-5;


	public Preisach1D(){}


	public Preisach1D(int M, double mean,double width,double Hmax,double Bs,double Bseff, long seed){

		this.M=M;
		this.mean=mean;
		this.width=width;
		this.seed=seed;
		this.Hs=Hmax;
		this.Bs=Bs;
		this.Bseff=Bseff;
		this.DB=Bs/M;
		
		
		dim=1;

		 double phideg=10;
		 
		 double phirad=phideg*PI/180;
		 
			cfm=0;
			cfw=0;
			 
		r=new Random(3564656);

		K=new double[M];
		a=new double[M][2];
		on=new boolean[M];



		for(int j=0;j<M;j++){

			K[j]=1;

			double am=mean*r.nextGaussian()*(1+cfm*abs(sin(phirad)));
			double d=width*abs(r.nextGaussian())*(1+cfw*abs(sin(phirad)));
			a[j][0]=am-d/2;

			a[j][1]=am+d/2;


		}


	}
	
	public void updateAni(double c1, double c2){


				cfm=c1;
				cfw=c2;
				 
			r=new Random(3564656);

			for(int j=0;j<M;j++){



				double am=mean*r.nextGaussian()*(1+cfm);
				double d=width*abs(r.nextGaussian())*(1+cfw);
				a[j][0]=am-d/2;

				a[j][1]=am+d/2;


			}


		}
	
	
	public Preisach1D deepCopy(){
		
		Random r=new Random();
		long ss=this.seed;
		Preisach1D pr=new Preisach1D(this.M, this.mean,this.width,this.Hs,this.Bs,this.Bseff,ss);
		
		return pr;
	}

	public static void main2(String[] args)
	{

		



	}


	public Mat demagnetize(){
		return demagnetize(Hs,20);
	}
	public Mat demagnetize(double Hm){
		return demagnetize(Hm,20);
	}

	public Mat demagnetize(int nCycles){
		return demagnetize(Hs,nCycles);
	}

	
	
	public Mat demagnetize(double Hm,int nCycles){

		int L=200;
		double t=0;
		double dt=1./L;

		//double cc=2.0/nCycles;

		Vect H=new Vect(nCycles*L);

		int  ix=0;
		for(int i=0;i<nCycles;i++)
			for(int j=0;j<L;j++){

					double x=(1-(t+dt)/nCycles)*Math.sin(2*Math.PI*t-PI/2);
					//double x=Math.exp(-cc*t)*Math.sin(2*Math.PI*t);
				t+=dt;
				H.el[ix]=x*Hm;
				ix++;
			}


		Mat BH=this.getCurve(H);

		return BH;

	}

	
public Mat magnetizeUpTo(double Bpeak,int L){
		

		Vect H=new Vect().linspace(0, Hs, L);
		Vect B=new Vect(L);

		int ix=0;
		for(int i=0;i<L;i++){
			if(i==0){
				B.el[i]=this.getRes();
				H.el[i]=H.el[0];

				continue;
			}

			double dB=0;
			for(int j=0;j<M;j++)
			{

				if(H.el[i]>H.el[i-1] && H.el[i]>a[j][1]){
					if(!on[j]){
						dB+=this.DB*K[j];
						on[j]=true;
					}
				}
				else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0] ){

					if(on[j]){
						dB-=this.DB*K[j];
						on[j]=false;
					}
				}

			}

			B.el[i]=B.el[i-1]+2*dB;
			
			ix++;
			if(B.el[i]>Bpeak || Math.abs(B.el[i]-Bpeak)<=epsdB) {
				B.el[i]=Bpeak;
				break;
			}
			

			
		}


		Mat BH1=new Mat(ix,2);
		for(int i=0;i<ix;i++){
			BH1.el[i][0]=H.el[i];
			BH1.el[i][1]=B.el[i];
		}
		
		 Mat BH=this.distill(BH1);

		return BH;

	
	}
	
public Mat magnetizeDownTo(double Bpeak,int L){
		

		Vect H=new Vect().linspace(0, -Hs, L);
		Vect B=new Vect(L);

		int ix=0;
		for(int i=0;i<L;i++){
			if(i==0){
				B.el[i]=this.getRes();
				H.el[i]=H.el[0];

				continue;
			}

			double dB=0;
			for(int j=0;j<M;j++)
			{

				if(H.el[i]>H.el[i-1] && H.el[i]>a[j][1]){
					if(!on[j]){
						dB+=this.DB*K[j];
						on[j]=true;
					}
				}
				else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0] ){

					if(on[j]){
						dB-=this.DB*K[j];
						on[j]=false;
					}
				}

			}

			B.el[i]=B.el[i-1]+2*dB;
			
			ix++;
			if(B.el[i]<Bpeak || Math.abs(B.el[i]-Bpeak)<=epsdB) {
				B.el[i]=Bpeak;
				break;
			}
			
			
		}


		Mat BH1=new Mat(ix+1,2);
		for(int i=0;i<=ix;i++){
			BH1.el[i][0]=H.el[i];
			BH1.el[i][1]=B.el[i];
		}
		 Mat BH=this.distill(BH1);

		return BH;

	
	}


	

/*
public Mat getCurveX(Vect H){
	
	Mat H1=new Mat(H.length,2);
	H1.setCol(H, 0);
	
	Mat BH=getLocus(H1);
	
	Mat BH1=new Mat(H.length,2);
	
	BH1.setCol(H, 0);
	BH1.setCol(BH.getColVect(2), 1);
	
	return BH1;
	
}*/

	public Mat getCurve(Vect H){



		int L=H.length;
		Vect B=new Vect(L);
		
		


		for(int i=0;i<L;i++){
			if(i==0){
				B.el[i]=this.getRes();

				continue;
			}

			double dB=0;
			for(int j=0;j<M;j++)
			{

				if(H.el[i]>H.el[i-1] && H.el[i]>a[j][1]){
					if(!on[j]){
						dB+=this.DB*K[j];
						on[j]=true;
					}
				}
				else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0] ){

					if(on[j]){
						dB-=this.DB*K[j];
						on[j]=false;
					}
				}

			}

			B.el[i]=B.el[i-1]+2*dB;
		}


		Mat BH=new Mat(L,2);
		for(int i=0;i<L;i++){
			BH.el[i][0]=H.el[i];
			BH.el[i][1]=B.el[i];
		}

		return BH;

	}
	

	

	public double getRes(){
		double Br=0;

		for(int j=0;j<M;j++){

			if(on[j]) 
				Br+=this.DB;
			else	
				Br+=-this.DB;


		}

		return Br;
	}
	

	public  Mat initial(int L){
		this.demagnetize();

		return 	magnetizeUpTo(this.Bseff,L);
	}


	public  Mat symMajorDesc(int L){
		return symDesc(Bseff,L);
	}


	public  Mat symMajorAsc(int L){
		
		return symAsc(-Bseff,L);
	}
	
	public  Mat symMajorFull(int L){
		
		return symFull(Bseff,L);
	}
	
	public  Mat symDesc(double Bpeak,int L){

		this.demagnetize();

		Mat BH1=this.magnetizeUpTo(Bpeak,L);

		int L1=BH1.nRow;
		double H1=BH1.el[L1-1][0];
		
	
	
		Vect seqH=new Vect().linspace(H1, -Hs, L);
	

		Mat BH2=this.getCurve(seqH);
		
		Mat BH3=this.distill(BH2);
		
		int ix=0;
	
		for(int i=0;i<BH3.nRow;i++){
			if(BH3.el[i][1]>=-Bpeak)
				ix++;
		}

		
		Mat BH4=new Mat(ix,2);
		
		for(int i=0;i<ix;i++){
			BH4.el[i]=BH3.el[i];
	
		}


		return BH4;
	}


	public  Mat symAsc(double Bpeak,int L){

		this.demagnetize();

		Mat BH1=this.magnetizeDownTo(-Bpeak,L);
		
		int L1=BH1.nRow;
		double H1=BH1.el[L1-1][0];
		
		Vect seqH=new Vect().linspace(H1, Hs, L);
	

		Mat BH2=this.getCurve(seqH);
		
		Mat BH3=this.distill(BH2);
		
		int ix=0;
	
		for(int i=0;i<BH3.nRow;i++){
			if(BH3.el[i][1]<=Bpeak)
				ix++;
		}

		
		Mat BH4=new Mat(ix,2);
		
		for(int i=0;i<ix;i++){
			BH4.el[i]=BH3.el[i];
	
		}
		

		return BH4;

	}
	
	public  Mat symFull(double Bpeak,int L){
		
		Mat BH1=this.symAsc(Bpeak,L);


		Mat BH2=this.symDesc(Bpeak,L);
		int L1=BH1.nRow;
		int L2=BH2.nRow;
		Mat BH3=new Mat(L1+L2+1,2);
		for(int i=0;i<BH1.nRow;i++){
			BH3.el[i]=BH1.el[i];
		}
		for(int i=0;i<BH2.nRow;i++){
			BH3.el[i+L1]=BH2.el[i];
		}

		BH3.el[L1+L2]=BH1.el[0];
		
		return BH3;
		
	}
	
	public  Mat revDesc(double Bpeak,int L){
		
	this.demagnetize();
	this.magnetizeDownTo(-Bseff,L);
	
	Vect H=new Vect().linspace(-Hs, Hs, L);
	
	Vect B=new Vect(L);

	int ix=0;
	for(int i=0;i<L;i++){
		if(i==0){
			B.el[i]=this.getRes();
			H.el[i]=H.el[0];

			continue;
		}

		double dB=0;
		for(int j=0;j<M;j++)
		{

			if(H.el[i]>H.el[i-1] && H.el[i]>a[j][1]){
				if(!on[j]){
					dB+=this.DB*K[j];
					on[j]=true;
				}
			}
			else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0] ){

				if(on[j]){
					dB-=this.DB*K[j];
					on[j]=false;
				}
			}

		}

		B.el[i]=B.el[i-1]+2*dB;
		
		ix++;
		if(B.el[i]>Bpeak || Math.abs(B.el[i]-Bpeak)<=epsdB) {
			B.el[i]=Bpeak;
			break;
		}
		
		
	}

	double H1=H.el[ix-1];

	H=new Vect().linspace(H1, -Hs, L);

	Mat BH1=this.getCurve(H);

	 Mat BH=this.distill(BH1);

	
	return BH;
	
	}
	
	public  Mat revAsc(double Bpeak,int L){
		
		this.demagnetize();
		this.magnetizeUpTo(Bseff,L);
		
		Vect H=new Vect().linspace(Hs, -Hs, L);
	
		Vect B=new Vect(L);

		int ix=0;
		for(int i=0;i<L;i++){
			if(i==0){
				B.el[i]=this.getRes();
				H.el[i]=H.el[0];

				continue;
			}

			double dB=0;
			for(int j=0;j<M;j++)
			{

				if(H.el[i]>H.el[i-1] && H.el[i]>a[j][1]){
					if(!on[j]){
						dB+=this.DB*K[j];
						on[j]=true;
					}
				}
				else if(H.el[i]<H.el[i-1] && H.el[i]<a[j][0] ){

					if(on[j]){
						dB-=this.DB*K[j];
						on[j]=false;
					}
				}

			}

			B.el[i]=B.el[i-1]+2*dB;
			
			ix++;
			if(B.el[i]<Bpeak || Math.abs(B.el[i]-Bpeak)<=epsdB) {
				B.el[i]=Bpeak;
				break;
			}
			
			
		}
		


		double H1=H.el[ix-1];

		H=new Vect().linspace(H1, Hs, L);

		Mat BH1=this.getCurve(H);

		 Mat BH=this.distill(BH1);

		return BH;
		
		}

	public Mat distill(Mat BH1){

		int Leff=1;

		int L=BH1.nRow;
		boolean[] skip=new boolean[L];
		for(int i=1;i<BH1.nRow;i++){


			if( Math.abs(BH1.el[i][1]-BH1.el[i-1][1])/Math.abs(BH1.el[i][0]-BH1.el[i-1][0])<epsdBdH){
				skip[i]=true;
			
				continue;
			}


			Leff++;
		}

		Mat BH=new Mat(Leff,2);

		int ix=0;
		for(int i=0;i<BH1.nRow;i++){

			if(!skip[i]){
				BH.el[ix][0]=BH1.el[i][0];
				BH.el[ix][1]=BH1.el[i][1];
				ix++;
				
			}


		}




		return BH;
	}

}
