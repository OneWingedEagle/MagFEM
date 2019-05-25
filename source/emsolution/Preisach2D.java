package emsolution;

import java.util.Random;

import static java.lang.Math.*;
import math.Mat;
import math.Vect;
import math.util;

public class Preisach2D {

	public Random r;
	public int M,nphi,dim,kRotated;
	public long seed;
	public double cfm,cfw,mean,width,Hs,BsM,Bseff,DB2D,dphiRad;
	public boolean[][] on;
	public double[] phi;
	public double[][][] a;
	public double[][] K;
	public double  epsdB=1e-2;
	public double  epsdBdH=1e-5;


	public Preisach2D(){}


	public Preisach2D(int M, double mean,double width,double Hmax,double BsM,double Bseff, long seed){

		this.M=M;
		this.mean=mean;
		this.width=width;
		this.seed=seed;
		this.Hs=Hmax;
		this.BsM=BsM;
		this.Bseff=Bseff;

		dim=2;

		cfm=4;
		cfw=4;

		int nh=9;

		nphi=2*nh+1;

		if(nphi==1)  this.dphiRad=0;
		else
			this.dphiRad=Math.PI/(nphi-1);


		double sum=1;
		for(int i=1;i<nphi;i++){
			sum+=sin(i*dphiRad);
		}


		sum/=nphi;



		DB2D=BsM/M/nphi/sum;;


		r=new Random(3564656);

		K=new double[M][nphi];
		a=new double[M][2][nphi];
		on=new boolean[M][nphi];
		phi=new double[nphi];

		double dphiDeg=180.0/(nphi-1);

		for(int k=0;k<nphi;k++){

			phi[k]=k*dphiDeg;

			double phirad=k*this.dphiRad;

			for(int j=0;j<M;j++){

				K[j][k]=1-.0*sin(2*phirad);

				double am=mean*r.nextGaussian()*(1+cfm*abs(sin(phirad)));
				//double am=2*mean*(.5-r.nextDouble())*(1+cfm*abs(sin(phirad)));
				double d=width*abs(r.nextGaussian())*(1+cfw*abs(sin(phirad)));

				a[j][0][k]=am-d/2;

				a[j][1][k]=am+d/2;




			}
		}


	}

	public Preisach2D deepCopy(){

		//Random r=new Random();
		long ss=this.seed;
		Preisach2D pr=new Preisach2D(this.M, this.mean,this.width,this.Hs,this.BsM,this.Bseff,ss);
		pr.dphiRad=this.dphiRad;

		return pr;
	}

	public static void main2(String[] args)
	{

	}

	public  Mat initial(double Bpeak,int L){

		this.demagnetize();
		Mat BH1=this.magnetizeUpTo(Bpeak,L);


		return BH1;
	}

	public Mat getCurveAlt(Vect H){

		int L=H.length;
		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);


		for(int k=0;k<nphi;k++){


			int kr=(k+kRotated)%nphi;


			Vect Halt=H.times(cos(k*dphiRad));


			Vect er=new Vect(cos(k*dphiRad),sin(k*dphiRad));

			for(int i=0;i<L;i++){

				if(i==0){

					B[i][k]=B[i][k].add(this.getRes(kr));
					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(Halt.el[i]>Halt.el[i-1] && Halt.el[i]>a[j][1][kr]){
						if(!on[j][kr]){
							dB=dB.add(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=true;
						}
					}
					else if(Halt.el[i]<=Halt.el[i-1] && Halt.el[i]<a[j][0][kr] ){

						if(on[j][kr]){
							dB=dB.sub(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=false;
						}
					}

				}

				B[i][k]=B[i-1][k].add(dB.times(2));
			}
		}

		Vect[] Bsum=new Vect[L];

		for(int i=0;i<L;i++){
			Bsum[i]=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bsum[i]=Bsum[i].add(B[i][j]);

		}


		Mat BH=new Mat(L,3);
		for(int i=0;i<L;i++){
			BH.el[i][0]=H.el[i];

			BH.el[i][1]=Bsum[i].el[0];
			BH.el[i][2]=Bsum[i].el[1];
		}


		return BH;


	}



	public Mat getLocusBRotation(double Hm,int Lc,int Nc){

		int L=Nc*Lc;
		Mat Hp=new Mat(L,2);
		for(int j=0;j<L;j++){
			Hp.el[j][0]=Hm*Math.cos(4*j*Math.PI/L);
			Hp.el[j][1]=Hm*Math.sin(4*j*Math.PI/L);
		}

		return getLocus(Hp);



	}

	public Mat getLocus(Mat H){

		int L=H.nRow;
		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);

		Vect Hr=new Vect(L);


		for(int k=0;k<nphi;k++){

			
			int kr=(k+kRotated)%nphi;
	
			Vect er=new Vect(cos(k*dphiRad),sin(k*dphiRad));

			for(int i=0;i<L;i++)
				Hr.el[i]=new Vect(H.el[i]).dot(er);


			for(int i=0;i<L;i++){


				if(i==0){

					B[i][k]=B[i][k].add(this.getRes(kr));

					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(Hr.el[i]>Hr.el[i-1] && Hr.el[i]>a[j][1][kr]){
						if(!on[j][kr]){
							dB=dB.add(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=true;
						}
					}
					else if(Hr.el[i]<=Hr.el[i-1] && Hr.el[i]<a[j][0][kr] ){

						if(on[j][kr]){
							dB=dB.sub(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=false;
						}
					}

				}

				B[i][k]=B[i-1][k].add(dB.times(2));
			}
		}

		Vect[] Bsum=new Vect[L];

		for(int i=0;i<L;i++){
			Bsum[i]=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bsum[i]=Bsum[i].add(B[i][j]);
		}


		Mat BH=new Mat(L,4);
		for(int i=0;i<L;i++){
			BH.el[i][0]=H.el[i][0];
			BH.el[i][1]=H.el[i][1];
			BH.el[i][2]=Bsum[i].el[0];
			BH.el[i][3]=Bsum[i].el[1];
		}



		return BH;

	}
	

	
	public Mat getLocusb(Mat H){ // some problem

		int L=H.nRow;
		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);

		double Hr=0,Hrp;


		for(int i=0;i<L;i++){

			Hrp=Hr;
			
			for(int k=0;k<nphi;k++){


				int kr=(k+kRotated)%nphi;

				Vect er=new Vect(cos(k*dphiRad),sin(k*dphiRad));

				Hr=new Vect(H.el[i]).dot(er);

				if(i==0){

					B[i][k]=B[i][k].add(this.getRes(kr));

					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(i>0 && Hr>Hrp && Hr>a[j][1][kr]){
						if(!on[j][kr]){
							dB=dB.add(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=true;
						}
					}
					else if(i>0 && Hr<=Hrp && Hr<a[j][0][kr] ){

						if(on[j][kr]){
							dB=dB.sub(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=false;
						}
					}

				}


				B[i][k]=B[i-1][k].add(dB.times(2));
			}
		}

		Vect[] Bsum=new Vect[L];

		for(int i=0;i<L;i++){
			Bsum[i]=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bsum[i]=Bsum[i].add(B[i][j]);
		}


		Mat BH=new Mat(L,4);
		for(int i=0;i<L;i++){
			BH.el[i][0]=H.el[i][0];
			BH.el[i][1]=H.el[i][1];
			BH.el[i][2]=Bsum[i].el[0];
			BH.el[i][3]=Bsum[i].el[1];
		}



		return BH;

	}

	public Mat getLoopBinputBisec(double Br){


		this.demagnetize();



		int iang=0;
		double phiRad=iang*Math.PI/18;;

		int steps = 36;
		Mat  X= new Mat(steps, 2);


		double yxRatio = 1;
		double a11 = Math.cos(phiRad);
		double a21 = Math.sin(phiRad);
		double a12 = -a21*yxRatio;
		double a22 = a11*yxRatio;


		double Xm = 10;

		for (int i = 0; i < steps; ++i){
			double wt = 4 * Math.PI*i / steps;
			double xx = Xm* Math.cos(wt);
			double yy = Xm* Math.sin(wt);

			double xxt = a11*xx + a12*yy;
			double yyt = a21*xx + a22*yy;

			X.el[i][0] = xxt;
			X.el[i][1] = yyt;

		}

	//	Mat Y=X.times(Br/Xm);

		int P=5;


		int L=steps*P;
		Vect[] H=new Vect[L];

		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);


		double Hrp,Hr=0;

		double Bm=0;
		double eps=.1;


		int ix=-1;

		double Bm1=.01,Bm2=Bseff;
		double magn1=10,magn2=Hs, magn=0;

		magn1=10;magn2=1000;

		boolean reached;


		for(int i=0;i<steps;i++){
			
			Vect eH=new Vect(X.el[i]).normalized();
			
			int px=0;
			while(px<P){

				ix++;
				px++;

				reached=false;

				magn=(magn1+magn2)/2;
				
				H[ix]=eH.times(magn);

			Hrp=Hr;
			
			for(int k=0;k<nphi;k++){


				int kr=(k+kRotated)%nphi;

				Vect er=new Vect(cos(k*dphiRad),sin(k*dphiRad));

				Hr=H[ix].dot(er);

				if(i==0){

					B[ix][k]=B[ix][k].add(this.getRes(kr));

					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(i>0 && Hr>Hrp && Hr>a[j][1][kr]){
						if(!on[j][kr]){
							dB=dB.add(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=true;
						}
					}
					else if(i>0 && Hr<=Hrp && Hr<a[j][0][kr] ){

						if(on[j][kr]){
							dB=dB.sub(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=false;
						}
					}

				}


				B[ix][k]=B[ix-1][k].add(dB.times(2));
			}
			
			Vect Bt=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bt=Bt.add(B[ix][j]);
			Bm=Bt.norm();



			if((Bm-Br)*(Bm1-Br)<0)
			{

				magn2=magn1;
				Bm2=Bm1;
				magn1=magn;
				Bm1=Bm;

			}
			else if((Bm-Br)*(Bm2-Br)<0){

				magn1=magn;
				Bm1=Bm;


			}
			else if(Bm<Br){

				magn1=magn;
				Bm1=Bm;
				magn2*=2;
				Bm2*=2;
			}
			else if(Bm>Br){

				magn1=magn;
				Bm1=Bm;
				magn2/=2;
				Bm2/=2;
			}

			util.pr(magn1+" "+magn2+ "  magn "+magn+"     <><>      "+Bm1+" "+Bm2+" Bm "+Bm);

			if(Math.abs(Bm-Br)<eps) reached=true;
				
		}
		}

		Vect[] Bsum=new Vect[L];

		for(int i=0;i<L;i++){
			Bsum[i]=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bsum[i]=Bsum[i].add(B[i][j]);
		}


		Mat BH=new Mat(L,4);
		for(int i=0;i<L;i++){
			BH.el[i][0]=H[i].el[0];
			BH.el[i][1]=H[i].el[1];
			BH.el[i][2]=Bsum[i].el[0];
			BH.el[i][3]=Bsum[i].el[1];
		}


		return BH;

	}


	public Mat getLoopBinput(double Br,double cf){


		this.demagnetize();
		
		Mat BH1=this.magnetizeUpTo(Br, 1000);
		this.demagnetize();
		this.getCurveAlt(BH1.getColVect(0)); // this is to reset the model to the given Bpeak;

		double  H0=1;

		int iang=0;
		double phiRad=iang*Math.PI/18;;

		int steps = 360*2;

		Mat  X= new Mat(steps, 2);



		double yxRatio = 1;
		double a11 = Math.cos(phiRad);
		double a21 = Math.sin(phiRad);
		double a12 = -a21*yxRatio;
		double a22 = a11*yxRatio;


		double Xm = H0;

		for (int i = 0; i < steps; ++i){
			double wt = 4 * Math.PI*i / steps;
			double xx = Xm* Math.cos(wt);
			double yy = Xm* Math.sin(wt);

			double xxt = a11*xx + a12*yy;
			double yyt = a21*xx + a22*yy;

			X.el[i][0] = xxt;
			X.el[i][1] = yyt;

		}


		int P=100;


		int L=steps*P;
		Vect[] H=new Vect[L];

		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);

		Vect[] Hf=new Vect[steps];
		Vect[] Bf=new Vect[steps];
		Vect Bt=new Vect(2);

		double Hrp,Hr=0;

		double Bm=Br;
		double eps=.1;

		int ix=-1;

		boolean reached;

		double cr=1;

		for(int i=0;i<steps;i++){
			//cr=0;
			Vect eH=new Vect(X.el[i]).normalized();
			
			int px=0;
			while(px<P){

				ix++;
				px++;

				reached=false;
				
				if(Bm<Br)
				cr*=1+cf*(Br-Bm);
				else
					cr/=1+cf*(Bm-Br);
				
/*				if(Bm<Br)
				cr+=cf*(Br-Bm);
				else
					cr-=cf*(Bm-Br);	
			
				if(cr<=0) cr=1;*/
				
				
				H[ix]=new Vect(X.el[i]).add(eH.times(cr));
				if(ix==0)H[ix].hshow();

				if(ix==0) Hrp=Hr-1;
				else 	Hrp=Hr;
			
			for(int k=0;k<nphi;k++){


				int kr=(k+kRotated)%nphi;

				Vect er=new Vect(cos(k*dphiRad),sin(k*dphiRad));

				Hr=H[ix].dot(er);

				if(i==0){

					B[ix][k]=B[ix][k].add(this.getRes(kr));

					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(i>0 && Hr>Hrp && Hr>a[j][1][kr]){
						if(!on[j][kr]){
							dB=dB.add(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=true;
						}
					}
					else if(i>0 && Hr<=Hrp && Hr<a[j][0][kr] ){

						if(on[j][kr]){
							dB=dB.sub(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=false;
						}
					}

				}


				B[ix][k]=B[ix-1][k].add(dB.times(2));
			}
			
			Bt=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bt=Bt.add(B[ix][j]);
			Bm=Bt.norm();
			
			if(px>1 && Math.abs(Br-Bm)<eps){
				break;
			}

		}
			
			util.pr(px+"  "+Bm);

			Bf[i]=Bt.deepCopy();
			
			Hf[i]=new Vect(X.el[i]).add(eH.times(cr));
		}

		Mat R=new Mat();
		Vect v=new Vect(steps/2);
		Mat BH=new Mat(steps/2,4);
		for(int i=0;i<BH.nRow;i++){
			int j=i+steps/2;
			
		v.el[i]=1*util.getAng(Hf[j])-1*util.getAng(Bf[j]);
		
		double diff=Math.acos(Math.cos(v.el[i]))/Math.PI*180;
		

				
		//util.pr((int)(180/Math.PI*Math.floor(util.getAng(Bf[j]))));
		
		if(i==0)
			 R=util.rotMat2D(v.el[i]);
		
		v.el[i]=diff;

			Vect Hrotated=R.mul(Hf[j]);
			Vect Brotated=R.mul(Bf[j]);
			
			BH.el[i][0]=Hrotated.el[0];
			BH.el[i][1]=Hrotated.el[1];
			
			BH.el[i][2]=Brotated.el[0];
			BH.el[i][3]=Brotated.el[1];
			
			
/*			BH.el[i][0]=Hf[j].el[0];
			BH.el[i][1]=Hf[j].el[1];
			BH.el[i][2]=Bf[j].el[0];
			BH.el[i][3]=Bf[j].el[1];*/
			
	
			
		}
		//v.show();

	//	util.plot(v);

		return BH;

	}




	public Mat demagnetize(){
		return demagnetize(Hs,100,20);
	}
	public Mat demagnetize(double Hm){
		return demagnetize(Hm,100,20);
	}

	public Mat demagnetize(int nCycles){
		return demagnetize(Hs,100,nCycles);
	}



	public Mat demagnetize(double Hm,int L,int nCycles){



		for(int k=0;k<this.nphi;k++)
			for(int j=0;j<this.M;j++)
				on[j][k]=false;

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

		Mat Hp=new Mat(nCycles*L,2);

		ix=0;
		t=0;
		for(int i=0;i<nCycles;i++)
			for(int j=0;j<L;j++){

				double x=(1-(t+dt)/nCycles)*Math.sin(2*Math.PI*t-PI/2);
				double y=(1-(t+dt)/nCycles)*Math.cos(2*Math.PI*t-PI/2);
				//double x=Math.exp(-cc*t)*Math.sin(2*Math.PI*t);
				t+=dt;
				Hp.el[ix][0]=x*Hm;
				Hp.el[ix][1]=y*Hm;
				ix++;
			}


		Mat BH=this.getLocus(Hp);

		return BH;

	}


	public Mat magnetizeUpTo(double Bpeak,int L){

		Vect H=new Vect().linspace(0, Hs, L);


		Vect[][] B=new Vect[L][nphi];
		for(int i=0;i<L;i++)
			for(int k=0;k<nphi;k++)
				B[i][k]=new Vect(2);

		Vect Hr=new Vect(L);


		for(int k=0;k<nphi;k++){

			int kr=(k+kRotated)%nphi;

			Vect er=new Vect(cos(k*dphiRad),sin(k*dphiRad));


			Hr=H.times(cos(k*dphiRad));

			for(int i=0;i<L;i++){


				if(i==0){
					//B[i][k]=er.times(this.getRes());
					B[i][k]=B[i][k].add(this.getRes(kr));

					continue;
				}

				Vect dB=new Vect(2);
				for(int j=0;j<M;j++)
				{

					if(Hr.el[i]>Hr.el[i-1] && Hr.el[i]>a[j][1][kr]){
						if(!on[j][kr]){
							dB=dB.add(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=true;
						}
					}
					else if(Hr.el[i]<=Hr.el[i-1] && Hr.el[i]<a[j][0][kr] ){

						if(on[j][kr]){
							dB=dB.sub(er.times(this.DB2D*K[j][kr]));
							on[j][kr]=false;
						}
					}

				}

				B[i][k]=B[i-1][k].add(dB.times(2));


			}
		}


		Vect[] Bsum=new Vect[L];

		for(int i=0;i<L;i++){
			Bsum[i]=new Vect(2);
			for(int j=0;j<nphi;j++)
				Bsum[i]=Bsum[i].add(B[i][j]);
		}

		int ix=0;
		for(int i=0;i<L;i++){
			if(Bsum[i].el[0]>Bpeak) break;
			ix++;
		}
		
		
		
	

		Mat BH1=new Mat(ix,3);
		for(int i=0;i<ix;i++){
			BH1.el[i][0]=H.el[i];
			BH1.el[i][1]=Bsum[i].el[0];
			BH1.el[i][2]=Bsum[i].el[1];
		}

		
		Mat BH=this.distill(BH1);


		return BH;

	}

	public Vect getRes(){

		Vect Br=new Vect(2);


		for(int k=0;k<nphi;k++){


			Vect er=new Vect(cos(k*dphiRad),sin(k*dphiRad));

			for(int j=0;j<M;j++)
			{

				if(on[j][k])
					Br=Br.add(er.times(this.DB2D*K[j][k]));
				else
					Br=Br.add(er.times(-this.DB2D*K[j][k]));

			}

		}


		return Br;
	}

	public Vect getRes(int k){

		Vect Br=new Vect(2);

		Vect er=new Vect(cos(k*dphiRad),sin(k*dphiRad));

		for(int j=0;j<M;j++)
		{

			if(on[j][k])
				Br=Br.add(er.times(this.DB2D*K[j][k]));
			else
				Br=Br.add(er.times(-this.DB2D*K[j][k]));

		}




		return Br;
	}


	public  Mat symMajorDesc(int L){
		return symDesc(Bseff,L);
	}

	public  Mat symDesc(double Bpeak,int L){

		Mat[] m=new Mat[3];

		this.demagnetize();

		Mat BH1=this.magnetizeUpTo(Bpeak,L);

		int L1=BH1.nRow;
		double H1=BH1.el[L1-1][0];

		this.demagnetize();

		Vect seqH=new Vect().linspace(0, H1, L).aug(new Vect().linspace(H1, -Hs, L));
		BH1=this.getCurveAlt(seqH);

		//util.plot(BH1.getColVect(0),BH1.getColVect(1));

		int jx=0;

		while(jx<BH1.nRow && BH1.el[jx+1][0]>=BH1.el[jx][0]){jx++;}





		Mat BH2=new Mat(BH1.nRow-jx,3);


		for(int i=0;i<BH2.nRow;i++)
			BH2.el[i]=BH1.el[i+jx];



		//util.plot(BH2.getColVect(0),BH2.getColVect(1));


		Mat BH3=this.distill(BH2);

		int ix=0;

		while(ix<BH3.nRow && BH3.el[ix][1]>=-Bpeak){ix++;}

		Mat BH4=new Mat(ix,3);

		for(int i=0;i<ix;i++){
			BH4.el[i][0]=BH3.el[i][0];
			BH4.el[i][1]=BH3.el[i][1];
			BH4.el[i][2]=BH3.el[i][2];
		}

		//util.plot(BH4.getColVect(0),BH4.getColVect(1));	 

		return BH4;
	}



	public Mat distill(Mat BH1){

		int Leff=1;
		boolean col3=(BH1.nCol==3);

		int L=BH1.nRow;
		boolean[] skip=new boolean[L];
		for(int i=1;i<BH1.nRow;i++){

			if( Math.abs(BH1.el[i][1]-BH1.el[i-1][1])/Math.abs(BH1.el[i][0]-BH1.el[i-1][0])<epsdB)
				if( Math.abs(BH1.el[i][1]-BH1.el[i-1][1])/Math.abs(BH1.el[i][0]-BH1.el[i-1][0])<epsdBdH)
				{
					skip[i]=true;

					continue;
				}


			Leff++;
		}

		int ncol=2;
		if(col3) ncol=3;
		Mat BH=new Mat(Leff,ncol);

		int ix=0;
		for(int i=0;i<BH1.nRow;i++){

			if(!skip[i]){
				BH.el[ix][0]=BH1.el[i][0];
				BH.el[ix][1]=BH1.el[i][1];
				if(col3)
					BH.el[ix][2]=BH1.el[i][2];
				ix++;

			}


		}




		return BH;
	}



}
