package math;

import static java.lang.Math.*;
import io.Loader;
import io.Writer;

import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Random;

import javax.swing.JFileChooser;
import javax.swing.JFrame;

import materialData.CurrentWaveForm;
import triangulation.ConvexHull;

import org.apache.commons.math3.special.Erf;
import org.math.plot.Plot2DPanel;

import fem.Element;
import fem.Model;
import fem.Node;


public class util {
	static DecimalFormat df=new DecimalFormat("0.00000E00");
	static String regex="[ : ,=\\t]+";

	static String logFile=null;
	public util(){}

	public static void main(String[] args) {
		int N=10;
		int P=3,Q=4;
		int M=P*Q;
		
		double [][] data=new double[N][M];
		for(int i=0;i<N;i++){
			int jx=0;
			for(int j=0;j<P;j++)
				for(int k=0;k<Q;k++){
					data[i][jx]=(int)(100*Math.exp(-abs(i-jx)*1./N));
					jx++;
				}
		}

		
		util.pr(N+"   "+M);
		for(int i=0;i<M;i++){
			util.ph("A10"+i+"\t");
		}
		util.pr("");
		for(int i=0;i<N;i++){
			for(int j=0;j<M;j++){

			util.ph(data[i][j]+"\t");
				}
			util.ph("B10"+i+"\t");
			util.pr("");
			
		}


	}


	public static void ContaminTransport() {
		// 1 D 

		boolean steady=false;
		
		double x0=0, xn=10000;
		int N=200;


		double dx=(xn-x0)/N;
		int nT1=10;
		int nTsub=10;
		if(steady){
			 nT1=1;
			 nTsub=1;
		}
		
		int nT=nT1*nTsub;

		double t0=0;
		double t1=nT1*86400;
		double dt=(t1-t0)/nT;
		

		
		int nunx=0;
		

		double src_duration=1e6*dt;
		
		double[] knownValues=new double[N];
		int[] un_kn_map=new int[N];
		int[] unIndex=new int[N];
		
		for(int i=0;i<N;i++){
			un_kn_map[i]=-1;
			unIndex[i]=-1;
		}

		for(int i=0;i<N;i++){

			if(i>=0 && i<1){
				knownValues[i]=1;
			}else
				if(i>=N-1 && i<-N){
					knownValues[i]=.5;
				}else{
				unIndex[i]=nunx;
				un_kn_map[nunx]=i;
				nunx++;
			}
			
		}
		
		Vect qC0=new Vect(N);
		
		for(int i=0;i<N;i++){

		//	if(i>=95 && i<=105)
		///	qC0.el[100]=-1.6e-4;
					
		}
		//qC0.el[N-1]=-1.2e-3;
		
		//qC0.el[50]=-1.0e-4;
		
		double V=4e-3;
		double D=2e-1;
		double Cr =V*dt/dx; // Cre should be less than 1.
		double Pe =V*dx/D; // Pe should be less than 2.
		util.pr("Cr = "+Cr+ "     *****  Cre = V*dt/dx  should be <= 1.");
		util.pr("Pe = "+Pe+ "     *****  Pe = V*dx/D  should be <= 2.");
	//	util.pr(Pe);
		double C0=1;
		
			double dx2=dx*dx;
		double at=1/dt;
		double am=-(V/(2*dx)+D/dx2);
		double a=2*D/dx2;
		double an=(V/(2*dx)-D/dx2);

		Mat H=new Mat(N, N);
		
		for(int i=0;i<N;i++){

			if(i>0)
				H.el[i][i-1]=am;
			
			H.el[i][i]=a;
			
		if(i<N-1) H.el[i][i+1]=an;

			if(nT>1)
				H.el[i][i]+=at;

		}
		

		int nun=nunx;
	//	H.show();

		Mat H2=new Mat(nun, nun);
	//	util.show(unIndex);	
//util.hshow(H.size());
//util.hshow(H2.size());

		for(int k=0;k<N;k++){
			if(unIndex[k]>=0) {
			Vect v=new Vect(nun);

			for(int j=0;j<N;j++){
				if(unIndex[j]>=0){
				v.el[unIndex[j]]=H.el[k][j];
				H.el[k][j]=0;
				}
			}
			//v.hshow();
			H2.setRow(v, unIndex[k]);
			}
			}
//H2.show();
//util.hshow(un_kn_map);
//util.hshow(unIndex);
		Vect Cx=new Vect(N);
		for(int i=0;i<N;i++){
			if(unIndex[i]<0)
				Cx.el[i]=knownValues[i];
		}
	
//		Cx.show();
//H.show();
		Vect bx=H.mul(Cx);
		

	

		
//	bx.show();
		Vect b0=new Vect(nun);
		Vect qC=new Vect(nun);
		for(int i=0;i<N;i++){
			if(unIndex[i]>=0){
				b0.el[unIndex[i]]=-bx.el[i];
				qC.el[unIndex[i]]=-qC0.el[i];
			}
		}

		Mat M=H2;
		
		Vect[] C=new Vect[nT1];
		
		
		C[0]= new Vect(N);
		for(int k=0;k<N;k++)
			C[0].el[k]=Cx.el[k];
		
		Vect sol=new Vect(nun);
			
		Vect b=new Vect(nun);
		Mat Minv=M.inv();

		int m=0;
		for(int m1=0;m1<nT1;m1++)
			for(int m2=0;m2<nTsub;m2++){
			m++;
			double t=m*dt;
			double f=1;//exp(-1*t*.000002);//+.9*cos(t*m/N*2*4*PI);
			//if(t>5*86400) f=0;
		//	if(t>src_duration)
			//	b=sol.times(at);
		//	else
				b=b0.add(qC.times(f)).add(sol.times(1*at));	
			
			sol=Minv.mul(b);

	//		sol.hshow();

			if(m2==0 && (steady ||m1>0)){
			 C[m1]= new Vect(N);
			// if(t<=src_duration)

				for(int k=0;k<N;k++){
					if(unIndex[k]>=0)
						C[m1].el[k]=sol.el[unIndex[k]];
					else
						C[m1].el[k]=knownValues[k];
				}
			 
		//	 C[m1].hshow();
			}
			}
//
		//analytical

		
		Mat y=new Mat(nT1,N);
	
		for(int m1=0;m1<nT1;m1++){
			double t=m1*nTsub*dt;
			if(t==0) t=1e-6;
		for(int k=0;k<N;k++){
			double x=k*dx;
			double tt=C0*.5*(Erf.erfc((x-V*t)/(2*sqrt(D*t)))+Math.exp(V*x/D)*Erf.erfc((x+V*t)/(2*sqrt(D*t))));;
	
			y.el[m1][k]=tt;
		}
		}
		
	Mat AN=new Mat(N,1+nT1);
	Mat NUM=new Mat(N,nT1+1);
	
	for(int k=0;k<N;k++){
		AN.el[k][0]=k*dx;
		NUM.el[k][0]=k*dx;
		for(int m1=0;m1<nT1;m1++){
			AN.el[k][m1+1]=y.el[m1][k];
		}
		
		for(int m1=0;m1<nT1;m1++){
			NUM.el[k][m1+1]=C[m1].el[k];
		}
	//	X.el[k][1+nT]=y.el[k];
	}
	//X.show();
	

	
	String[] titles1=new String[AN.nCol];
	for(int k=0;k<titles1.length;k++)
		titles1[k]="AN, t= "+(k)*nTsub*dt;
	String[] titles2=new String[AN.nCol];
	for(int k=0;k<titles2.length;k++)
		titles2[k]="t= "+(k)+" days";
	
	//titles[X.nCol-2]="Analytical steady ";
	
	//util.plotBunch(X.el,titles);
	
	//util.plotBunch(AN.el,titles1);
	//NUM.show();

	util.plotBunch(NUM.el,titles2);
	

	}
	
	public static void gmeshFunc() 		{
		double height=10;
		//DFTdecopler();

		int nP=4;
		Vect[] corners=new Vect[nP];
		corners[0]=new Vect(-1,-1,0);
		corners[1]=new Vect(1,-1,0);
		corners[2]=new Vect(1,1,0);
		corners[3]=new Vect(-1,1,0);
		for(int k=0;k<nP;k++){

			util.pr("Point("+(k+1)+") = {"+corners[k].el[0]+","+corners[k].el[1]+","+corners[k].el[2]+"};");
		}
		int nL=4;
		int[][] lines=new int[nP][2];
		lines[0][0]=1;	lines[0][1]=2;
		lines[1][0]=2;	lines[1][1]=3;
		lines[2][0]=3;	lines[2][1]=4;
		lines[3][0]=4;	lines[3][1]=1;
		for(int k=0;k<nL;k++){

			util.pr("Line("+(k+1)+") = {"+lines[k][0]+","+lines[k][1]+"};");
		}


		util.pr("Curve Loop(1) = {4, 1, 2, 3};");						
		util.pr("Plane Surface(1) = {1};");						


		int nCircles1=6;

		double r=.1;
		double r1=.21;


		Vect[] cents1=new Vect[nCircles1];

		for(int j=0;j<nCircles1;j++){
			cents1[j]=new Vect(r1*cos(j*PI/3+PI/6),r1*sin(j*PI/3+PI/6),0);
		}




		int nCircles2=6;
		double r2=.65;

		Vect[] cents2=new Vect[nCircles2];

		for(int j=0;j<nCircles2;j++){
			cents2[j]=new Vect(r2*cos(j*PI/3+PI/6),r2*sin(j*PI/3+PI/6),0);
		}

		int nCircs=nCircles1*nCircles1;

		Vect[] cents=new Vect[nCircs];


		for(int i=0;i<nCircles2;i++){
			Mat R=util.rotMat2D(i*PI/3+PI/6);
			for(int j=0;j<nCircles1;j++){

				cents[i*nCircles1+j]=R.mul(cents1[j].add(new Vect(r2,0,0)).v2()).v3();
			}
		}

		//util.pr("SetFactory(\"OpenCASCADE\");");
		for(int k=0;k<nCircs;k++){

			util.pr("Point("+(1000+k*1000+1)+") = {"+(cents[k].el[0]+r)+","+cents[k].el[1]+","+cents[k].el[2]+"};");
			util.pr("Point("+(1000+k*1000+2)+") = {"+cents[k].el[0]+","+(cents[k].el[1]+r)+","+cents[k].el[2]+"};");
			util.pr("Point("+(1000+k*1000+3)+") = {"+(cents[k].el[0]-r)+","+cents[k].el[1]+","+cents[k].el[2]+"};");
			util.pr("Point("+(1000+k*1000+4)+") = {"+cents[k].el[0]+","+(cents[k].el[1]-r)+","+cents[k].el[2]+"};");
			util.pr("Point("+(1000+k*1000+5)+") = {"+cents[k].el[0]+","+cents[k].el[1]+","+cents[k].el[2]+"};");

			util.pr("Circle("+(4*k+nL+1)+") = {"+(1000+k*1000+1)+","+(1000+k*1000+5)+","+(1000+k*1000+2)+"};");
			util.pr("Circle("+(4*k+nL+1+1)+") = {"+(1000+k*1000+2)+","+(1000+k*1000+5)+","+(1000+k*1000+3)+"};");
			util.pr("Circle("+(4*k+nL+1+2)+") = {"+(1000+k*1000+3)+","+(1000+k*1000+5)+","+(1000+k*1000+4)+"};");
			util.pr("Circle("+(4*k+nL+1+3)+") = {"+(1000+k*1000+4)+","+(1000+k*1000+5)+","+(1000+k*1000+1)+"};");

			//util.pr("Circle("+(k+nL+1)+") = {"+cents[k].el[0]+","+cents[k].el[1]+","+cents[k].el[2]+","+r+"0, 2*Pi};");

			util.pr("Curve Loop("+(k+nL+1)+") = {"+(4*k+nL+1)+","+(4*k+nL+2)+","+(4*k+nL+3)+","+(4*k+nL+4)+"};");

			//util.pr("Transfinite Curve {"+(k+nL+1)+"} = 20 Using Progression 1");	

			util.pr("Plane Surface("+(k+nL+1)+") = {"+(k+nL+1)+"};");			
		}



		for(int j=0;j<nCircles2;j++){
			util.pr("Extrude {{0,0,"+height+"},{0,0,1},{"+cents2[j].el[0]+","+cents2[j].el[1]+","+cents2[j].el[2]+"}, 2*Pi}{");	
			util.pr("Surface {"+(nL+j*nCircles2+1)+"}; Surface {"+(nL+j*nCircles2+2)+"}; Surface {"+(nL+j*nCircles2+3+"}; Surface {"+(nL+j*nCircles2+4)+"}; Surface {"+(nL+j*nCircles1+5)+"}; Surface {"+(nL+j*nCircles2+6))+"}; Layers{20}; Recombine;}");

		}

		util.pr("Extrude {{0,0,"+height+"},{0,0,1},{0,0,0}, 0}{");	
		util.pr("Surface {1}; Layers{20}; Recombine;}");

		for(int i=0;i<800;i++){
			util.pr("Transfinite Curve {"+(i+1)+"} = 20 Using Progression 1");	

		}
	}

	public static void DFTdecopler() {

		long startTime = System.currentTimeMillis();

		int N=2;


		Vect x=new Vect(N);

		int M=2*N-1;

		Mat Rm=new Mat(N,N);

		for(int i=0;i<N;i++){
			//if(i==0)
			x.el[i]=i+1./(i+2);//100./(1+abs(i-.5*N)/N);
			//R.el[i]=(1+pow(i-N/2.,2)/(N*N));
			for(int j=0;j<N;j++){
				Rm.el[i][j]=1+1*Math.abs(i-j);//1./(1.+abs(i-j));
			}
		}
		//	b=Rm.inv().mul(b);;
		x.hshow();
		//	Rm=Rm.inv();
		//	Rm.show();
		//	Vect I1=Rm.inv().mul(b);
		Vect y1=Rm.mul(x);
		//Rm=Rm.inv().deepCopy();
		//y1.hshow();


		Vect y2=new Vect(2*M);
		for(int j=0;j<N;j++){
			y2.el[j]=y1.el[j];

		}
		//y2.el[0]=y2.el[M-1];
		y2.hshow();

		Vect x2=new Vect(M);
		for(int j=0;j<N;j++){
			x2.el[j]=x.el[j];

		}

		Vect Rj=new Vect(M);

		for(int i=0;i<N;i++){

			for(int j=0;j<N;j++){
				Rj.el[i-j+N-1]=Rm.el[i][j];
			}
		}



		Complex[] Rk=DFT.dft(Rj.el);
		//Complex[] x2k=DFT.dft(x2.el);
		Complex[] y2k=DFT.dft(y2.el);
		Complex[] x2k=new Complex[M];

		/*					util.pr("Rk --------------");

					for(int j=0;j<Rk.length;j++){
						util.pr(String.format("%10.5f",Rk[j].norm()));
						}
					util.pr("x2k --------------");
					for(int j=0;j<x2k.length;j++){
					//	util.pr(String.format("%10.5f",x2k[j].norm()));
						util.pr(String.format("%10.5f",x2k[j].re)+String.format("%10.5f",x2k[j].re));

						}
		 */

		//util.pr("y2k --------------");
		for(int j=0;j<M;j++){
			//	y2k[j]=x2k[j].times(Rk[j]);
			//I2k[j]=I1k[j].times(Rk[j]);
		}

		for(int j=0;j<M;j++){
			//util.pr(String.format("%10.5f",y2k[j].norm()));
			//util.pr(String.format("%10.5f",Xk.re)+", "+String.format("%10.5f",Xk.im));
			//df.format(A[i][j])+"\t"
		}

		for(int j=0;j<M;j++){
			//	I2k[j]=I2k[j].times(Rk[j].inv());
			//I2k[j]=I1k[j].times(Rk[j]);
		}

		Complex[] I2c=DFT.idft(y2k);

		Vect I2=new Vect(M);


		for(int j=0;j<M;j++){
			I2.el[j]=I2c[j].re;
		}

		I2.hshow();
		Vect I3=new Vect(N);
		for(int j=0;j<N;j++){
			I3.el[j]=I2c[j+N-1].re;

		}
		//I3.el[0]=I2c[2*N-2].re;

		I3.hshow();

		//	util.pr("y2k --------------");
		for(int j=0;j<M;j++){

			x2k[j]=y2k[j].times(Rk[j].inv());
			//I2k[j]=I1k[j].times(Rk[j]);
		}

		util.pr("x2k --------------");
		for(int j=0;j<x2k.length;j++){
			util.pr(String.format("%10.5f",x2k[j].re)+String.format("%10.5f",x2k[j].re));
		}

		Complex[] x3c=DFT.idft(x2k);

		Vect x3=new Vect(M);


		for(int j=0;j<M;j++){
			x3.el[j]=x3c[j].re;
		}

		x3.hshow();

	}

	public static void main3(String[] args) throws Exception{

	}
	public static void main4(String[] args) {

		Mat M=new Mat();

		M.symrand(2);;
		for(int i=0;i<2;i++)
			for(int j=0;j<2;j++)			
				M.el[i][j]=.5-Math.random();
		M.show();

		//Vect x=M.eigQR(1e-6);
		//Mat Q=M.eigVectQR(1e-6);

		//	x.hshow();
		//Q.show();
		/*				M.el[0][0]=5;
				M.el[0][1]=2;
				M.el[1][0]=-1+M.el[0][1];
				M.el[1][1]=4;*/

		Vect v=new Vect().rand(2, 0, 1).sub(new Vect(.5,.5));

		//v=Q.getColVect(1);
		//v.el[0]+=.005;
		Vect v0=v.deepCopy();
		Vect vpr=v.deepCopy();

		Model model=new Model("C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\matrix pole\\quad.txt");
		Mat R=util.rotMat2D(v, new Vect(0,1));
		for(int j=1;j<=model.numberOfNodes;j++){
			Vect z=model.node[j].getCoord();

			Vect zr=R.mul(z);
			model.node[j].setCoord(zr);
			model.node[j].setDeformable(true);

		}
		Model model0=model.deepCopy();

		model.writeMesh("C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\matrix pole\\bun.txt");

		int N=50;
		for(int i=0;i<N;i++){
			vpr=v;
			v=M.mul(v);
			v.timesVoid(1./v.norm());
			//	v.hshow();
			Mat R2=util.rotMat2D(v, vpr);

			for(int j=1;j<=model.numberOfNodes;j++){
				Vect z0=model0.node[j].getCoord();
				Vect z=model.node[j].getCoord();
				Vect zr=R2.mul(z);


				Vect u=zr.sub(z0);

				model.node[j].setU(u.times(1e-7));		
				model.node[j].setCoord(zr);		

			}
			//util.pr(model.node[4].getCoord().sub(model.node[1].getCoord()).norm());	
			//model.writeMesh("C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\matrix pole\\bun"+i+".txt");
			model.writeNodalField("C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\matrix pole\\disp"+i+".txt", -1);
		}


	}

	private void invertYparam(){
		int M=25;
		int N=M*4;

		String file="D:\\Works\\2018\\2-Port 3D Trans\\cln\\1.8reg\\imp";

		Loader loder=new Loader();
		double[][] xx=loder.loadArrays(N,3,file);

		Mat[] Yr=new Mat[M];
		Mat[] Ym=new Mat[M];
		//Mat[] Zr=new Mat[M];
		//Mat[] Zm=new Mat[M];
		for(int i=0;i<M;i++){
			Yr[i]=new Mat(2,2);
			Ym[i]=new Mat(2,2);
			Yr[i].el[0][0]=xx[i][1];
			Ym[i].el[0][0]=xx[i][2];
			Yr[i].el[0][1]=xx[M+i][1];
			Ym[i].el[0][1]=xx[M+i][2];
			Yr[i].el[1][1]=xx[2*M+i][1];
			Ym[i].el[1][1]=xx[2*M+i][2];
			Yr[i].el[1][0]=xx[3*M+i][1];
			Ym[i].el[1][0]=xx[3*M+i][2];

			Complex Y11=new Complex(Yr[i].el[0][0],Ym[i].el[0][0]);
			Complex Y12=new Complex(Yr[i].el[0][1],Ym[i].el[0][1]);
			Complex Y21=new Complex(Yr[i].el[1][0],Ym[i].el[1][0]);
			Complex Y22=new Complex(Yr[i].el[1][1],Ym[i].el[1][1]);

			Complex detY=Y11.times(Y22).sub(Y12.times(Y21));
			Complex detYinv=detY.inv();

			Complex Z11=Y22.times(detYinv);
			Complex Z12=Y21.times(detYinv).times(-1);
			Complex Z21=Y12.times(detYinv).times(-1);
			Complex Z22=Y11.times(detYinv);

			util.pr(xx[i][0]+"\t"+Z11.re+"\t"+Z11.im);



		}

	}

	public static void main2(String[] args) throws Exception{
		int n=5;
		int[] xx=new int[n];
		for(int i=1;i<9;i++)
			for(int j=1;j<9;j++)
				for(int k=1;k<9;k++)
					for(int p=1;p<9;p++)
						for(int q=1;q<9;q++)
							//for(int r=1;r<9;r++)
						{
							xx[0]=i;
							xx[1]=j;
							xx[2]=k;
							xx[3]=p;
							xx[4]=q;


							int num=xx[0];
							for(int i1=1;i1<n;i1++){
								num*=10;
								num+=xx[i1];
							}


							int num2=0;
							//	for(int i1=0;i1<3;i1++)
							//for(int j1=0;j1<3;j1++)
							//	for(int k1=0;k1<2;k1++){

							num2=xx[0];
							num2+=(int)Math.pow(xx[1],xx[2]);
							//num2+=(int)Math.pow(xx[2],xx[3]);
							num2*=xx[3];


							int num3=0;
							num3+=(int)Math.pow(xx[0],xx[1]);
							num3+=(int)Math.pow(xx[2],xx[3]);
							//num3+=xx[4];

							num3=xx[0];
							num3+=(int)Math.pow(xx[1]*xx[2],1);
							num3+=(int)Math.pow(xx[3],xx[4]);
							//num3*=xx[4];
							//num3+=xx[3]*xx[4];
							//num3+=xx[4];
							if(num3==num) 
								//util.pr(num);
								util.pr(num+" = "+xx[0]+"+"+xx[1]+"*"+xx[2]+" + "+xx[3]+"^"+xx[4]);

							if(num2==num ) 
								//util.pr(num);
								util.pr(num+" = ("+xx[0]+"+"+xx[1]+"^"+xx[2]+" )* "+xx[3]);
							//	}

						}


		//double h=1;
		//double app=20*PI/180;
		//double b=2*tan(app/2)*h;

		double b=1;
		double rad80=PI*80./180;

		double x1=0;
		double x2=b;


		int nn=100;

		Vect ang_x_vals=new Vect(nn);

		Vect ang_C1_vals=new Vect(nn);

		double dx=b/2/100;

		double AD=sin(rad80)*b/sin(rad80/2);

		Vect Dcoord=new Vect(x1+AD*cos(PI/3),AD*sin(PI/3));
		Vect Ccoord=new Vect(x2,0);

		Vect AC=new Vect(b,0);

		for(int i=0;i<-nn;i++){

			double x=x1+i*dx;

			double y=x*tan(rad80);

			Vect Ecoord=new Vect(x,y);

			Vect ED=Dcoord.sub(Ecoord);
			Vect EC=Ccoord.sub(Ecoord);



			double x_angle=acos(ED.dot(EC)/(ED.norm()*EC.norm()));


			double C1_angle=acos(EC.dot(AC)/(AC.norm()*EC.norm()));


			ang_x_vals.el[i]=x_angle*180/PI;
			//ang_x_vals.el[i]=tan(x_angle);

			ang_C1_vals.el[i]=C1_angle*180/PI;;


			//util.pr(ang_x_vals.el[i]);

		}

		//ang_C1_vals.hshow();
		//ang_x_vals.hshow();

		//plot(ang_x_vals);
		//plot(ang_C1_vals,ang_x_vals);

		double s=0;

		int kx=0;
		int ky=10;
		int kz=10;

		Vect M=new Vect(1,0,0);
		Vect H=new Vect(3);

		double d=1;
		double dz=1;
		double vol=d*d*dz;

		for(int k=-kz;k<=kz;k++)
			for(int i=-kx;i<=kx;i++)
				for(int j=-ky;j<=ky;j++){
					Vect r=new Vect(i*d,j*d,k*dz);	
					double rn=r.norm();
					//	if(rn>2.) continue;
					if(j==0) continue;

					Vect m=M.times(vol);
					double mdotr=m.dot(r);

					double r2=r.dot(r);
					double r3=rn*r2;
					double r5=r3*r2;

					H=H.add(r.times(mdotr*3/r5).add(m.times(-1/r3)));

					s+=(3*i*d*i*d/r5-1./r3);

				}

		H=H.times(1./(4*PI));

		s/=4*PI;
		//H=H.add(new Vect(2./3,0,0));

		H.hshow();

		util.pr(s);


	}


	public static void mainForCoilMotionTest(String[] args) throws Exception{


		CurrentWaveForm inductance=new CurrentWaveForm("emf//inductance2.txt");
		util.plot(inductance.TI);
		/*		double T=.01;
		double wsin=2*PI/T;
		double Tsw=1./7200;
		int N=1800;
		double[] y=new double[N];
		double dt=T/N;
		for(int i=0;i<y.length;i++){
			double t=i*dt;
			y[i]=3*Math.sqrt(2)*Math.cos(wsin*t+4*PI/3)*(1+.1*util.triangWave(Tsw,t)+.1*util.triangWave(Tsw/2,t));
			util.pr(t+"\t"+y[i]);
		}*/

		int N=121;
		double[] current=new double[N];
		double dt=5e-3;
		double L0=.4;
		L0=1.83e-3;
		double C=.196;
		C=-5.2e-4;
		double R=5.0;
		R=1.e-1;
		double V=1e0;
		double T=2e-1;
		double omega=2*Math.PI/T;
		current[0]=8;
		for(int i=1;i<current.length;i++){
			double t=i*dt;
			//double L=(L0+C*sin(omega*t));
			double L=inductance.getI(t);
			//double L=(L0+10*C*t);
			//double dLdt=omega*C*cos(omega*t);
			double dLdt=(inductance.getI(t)-inductance.getI(t-dt))/dt;

			double denum=L/dt+(R+dLdt);
			//	V=100*sin(omega*t);
			double num=V+L*current[i-1]/dt;
			current[i]=num/denum;
		}
		double[] current2=new double[100];
		for(int i=0;i<current2.length;i++){
			//current2[i]=current[N-1-current2.length+i];
		}
		plot(current);

		util.show(current);

	}


	public static double max(double[] x){
		double max=x[0];
		for(int i=1;i<x.length;i++)
			if(x[i]>max)
				max=x[i];
		return max;
	}
	public static double max(final double a, final double b){
		double result =a;;
		if(b>a) result =b;
		return result;
	}


	public static int max(int[] x){
		int max=x[0];
		for(int i=1;i<x.length;i++)
			if(x[i]>max)
				max=x[i];
		return max;
	}	
	public static int indmax(double[] x){
		int indmax=0;
		double max=x[0];
		for(int i=1;i<x.length;i++)
			if(x[i]>max)
				indmax=i;
		return indmax;
	}
	public static int indpiv(Mat x,int j){

		int indmax=j;
		double max=abs(x.el[0][j]);
		for(int i=j;i<x.nRow;i++)
			if(abs(x.el[i][j])>max){
				max=abs(x.el[i][j]);
				indmax=i;
			}
		return indmax;
	}

	public static double[] linspace(double a, double b,int N){
		double[] v=new double[N];
		double d=(b-a)/(N-1);
		for(int i=0;i<N;i++)
			v[i]=a+i*d;
		return v;
	}

	public static double[][] copy(double[][] a){
		int I=a.length;
		int J=a[0].length;
		double[][] a1=new double[I][J];
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				a1[i][j]=a[i][j];			

		return a1;
	}

	public static double[][] cubicSpl(double[][] xy){
		int I=xy.length-1;
		double[][] coefs=new double[I][4];
		int L=I-1;
		double[] h=new double[I];
		for(int i=0;i<I;i++)
			h[i]=xy[i+1][0]-xy[i][0];

		Vect b=new Vect(L);
		for(int i=0;i<b.length;i++)
			b.el[i]=6*((xy[i+2][1]-xy[i+1][1])/h[i+1]-(xy[i+1][1]-xy[i][1])/h[i]);

		Mat A=new Mat(L,L);

		A.el[0][0]=2*(h[0]+h[1]);
		A.el[0][1]=h[1];

		for(int i=1;i<L-1;i++){

			A.el[i][i-1]=h[i];
			A.el[i][i]=2*(h[i]+h[i+1]);
			A.el[i][i+1]=h[i+1];

		}

		A.el[L-1][L-2]=h[L-1];
		A.el[L-1][L-1]=2*(h[L-2]+h[L-1]);

		Vect M1=gaussel(A, b);
		double[] M=new double[xy.length];
		M[0]=0; M[xy.length-1]=0;
		for(int i=1;i<M.length-1;i++){
			M[i]=M1.el[i-1];
		}

		for(int i=0;i<coefs.length;i++){
			coefs[i][0]=(M[i+1]-M[i])/(6*h[i]);
			coefs[i][1]=M[i]/2;
			coefs[i][2]=(xy[i+1][1]-xy[i][1])/h[i]-h[i]*(M[i+1]+2*M[i])/6;
			coefs[i][3]=xy[i][1];
		}


		return coefs;
	}

	public static Vect gaussel(Mat A, Vect b){
		int[] dim=A.size();
		if(dim[0]!=dim[1]) throw new IllegalArgumentException("Matrix is not square");
		int I=dim[0];
		Mat Ab=new Mat();
		Ab=A.aug(b);
		Ab.low0();
		Vect x=new Vect(I);
		x=solveup(Ab);
		return x;

	}

	public static Mat rotEuler(Vect rotAx,double alpha)
	{
		double e1,e2,e3,e4;

		e1=rotAx.el[0]*sin(alpha/2);
		e2=rotAx.el[1]*sin(alpha/2);
		e3=rotAx.el[2]*sin(alpha/2);
		e4=cos(alpha/2);

		Mat M=new Mat(3,3);
		M.el[0][0]=pow(e1,2)-pow(e2,2)-pow(e3,2)+pow(e4,2);
		M.el[0][1]=2*(e1*e2-e3*e4);
		M.el[0][2]=2*(e1*e3+e2*e4);
		M.el[1][0]=2*(e1*e2+e3*e4);
		M.el[1][1]=-pow(e1,2)+pow(e2,2)-pow(e3,2)+pow(e4,2);
		M.el[1][2]=2*(e2*e3-e1*e4);
		M.el[2][0]=2*(e1*e3-e2*e4);
		M.el[2][1]=2*(e2*e3+e1*e4);
		M.el[2][2]=-pow(e1,2)-pow(e2,2)+pow(e3,2)+pow(e4,2);
		return M;
	}


	public static Mat rotMat(Vect newAx,Vect oldAx){

		if(newAx.length==2) return rotMat2D(newAx,oldAx);

		Mat M=new Mat(3,3);

		double newAxn=newAx.norm();
		if(newAxn==0){M.eye(); return M;}
		double oldAxn=oldAx.norm();
		double alpha,cos;
		Vect rotAx=oldAx.cross(newAx);
		if(rotAx.norm()==0){

			M.eye();
			return M;

		}

		rotAx.normalize();


		cos=oldAx.dot(newAx)/(newAxn*oldAxn);
		if(cos>=1)
			alpha=0;
		else if(cos<=-1)alpha=PI;
		else
			alpha=acos(cos);

		return rotEuler(rotAx,alpha);
	}

	public static Mat rotMat(Vect newAx,Vect oldAx, Vect rotAx){

		if(newAx.length==2) return rotMat2D(newAx,oldAx);

		Mat M=new Mat(3,3);

		double newAxn=newAx.norm();
		double oldAxn=oldAx.norm();
		double rotAxn=rotAx.norm();
		if(newAxn==0|| oldAxn==0|| rotAxn==0){M.eye(); return M;}

		double alpha=0;
		double  cos=oldAx.dot(newAx)/(newAxn*oldAxn);
		if(cos>=1)
			alpha=0;
		else if(cos<=-1)alpha=PI;
		else
			alpha=acos(cos);


		return rotEuler(rotAx,alpha);
	}


	public static Mat rotMat2D(Vect newAx,Vect oldAx){

		double ang1=getAng(oldAx);
		double ang2=getAng(newAx);

		return rotMat2D(ang2-ang1);
	}

	public static Mat rotMat2D(double rad){


		Mat M=new Mat(2,2);
		M.el[0][0]=cos(rad);
		M.el[1][1]=	M.el[0][0];
		M.el[0][1]=-sin(rad);
		M.el[1][0]=-M.el[0][1];

		return M;
	}

	public static Mat tensorize(Vect v){
		int dim=(v.length+3)/3;
		Mat S=new Mat(dim,dim);
		for(int i=0;i<dim;i++)
			S.el[i][i]=v.el[i];
		if(dim==2) {
			S.el[0][1]=v.el[2];
			S.el[1][0]=v.el[2];
		}
		else {
			S.el[0][1]=v.el[3];
			S.el[1][0]=v.el[3];
			S.el[1][2]=v.el[4];
			S.el[2][1]=v.el[4];
			S.el[0][2]=v.el[5];
			S.el[2][0]=v.el[5];
		}
		return S;

	}

	public static Vect vectorize(Mat S){
		int dim=S.nCol;
		int L=3*(dim-1);
		Vect v=new Vect(L);
		for(int i=0;i<dim;i++)
			v.el[i]=S.el[i][i];

		if(dim==2) {
			v.el[2]=S.el[0][1];
		}
		else {
			v.el[3]=S.el[0][1];
			v.el[4]=S.el[1][2];
			v.el[5]=S.el[0][2];

		}
		return v;

	}

	public static Mat rotMatrix(Mat Q1,Mat Q2){
		Mat T=new Mat(3,3);
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				T.el[i][j]=Q1.getColVect(i).dot(Q2.getColVect(j));
		return T;
	}

	public static Vect solveup(Mat Ab){
		int I=Ab.nRow;
		int J=Ab.nCol;
		if(I!=J-1) throw new IllegalArgumentException("Matrix is not square");
		Vect x=new Vect(I);
		x.el[I-1]=Ab.el[I-1][J-1]/Ab.el[I-1][I-1];

		for(int i=I-2;i>=0;i--){
			double s=0;
			for(int j=i+1;j<I;j++)
				s=s+Ab.el[i][j]*x.el[j];
			x.el[i]=(Ab.el[i][J-1]-s)/Ab.el[i][i];
		}

		return x;
	}

	public static double getAng(Vect v){
		double ang=0;
		if(v.norm()==0) return ang;
		else if(v.el[0]>=0 && v.el[1]>=0) ang=atan(abs(v.el[1]/v.el[0]));
		else if(v.el[0]<=0 && v.el[1]>=0) ang=PI-atan(abs(v.el[1]/v.el[0]));
		else if(v.el[0]>=0 && v.el[1]<=0) ang=2*PI-atan(abs(v.el[1]/v.el[0]));

		else ang=atan(abs(v.el[1]/v.el[0]))+Math.PI;

		return ang;
	}

	public  static File getJFile(int mode){
		JFileChooser fileChooser = new JFileChooser();
		File theDirectory = new File(System.getProperty("user.dir"));
		fileChooser.setCurrentDirectory(theDirectory);
		int returnValue;
		if(mode==0)
			returnValue = fileChooser.showOpenDialog(null);
		else
			returnValue = fileChooser.showSaveDialog(null);

		if (returnValue == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			System.out.println(selectedFile.getPath());


			return selectedFile;

		}

		return null;
	}



	public  static String getFile(int mode){
		String filePath="";
		FileDialog fd;
		Frame f=new Frame();
		if(mode==0)
			fd= new FileDialog(f,"Select  file",FileDialog.LOAD);
		else
			fd= new FileDialog(f,"Select  file",FileDialog.SAVE);
		fd.setVisible(true);
		fd.toFront();
		String Folder=fd.getDirectory();
		String File = fd.getFile();
		if(Folder!=null && File!=null)
		{

			filePath=Folder+"\\"+File;

		}
		f.dispose();
		fd.dispose();

		return filePath;
	}

	public  static String getFile(){

		return getFile(0);
	}


	public static void shuffle(int[] ar){
		Random rnd = new Random();
		for (int i = ar.length - 1; i > 0; i--)
		{
			int index = rnd.nextInt(i + 1);
			// Simple swap
			int a = ar[index];
			ar[index] = ar[i];
			ar[i] = a;
		}
	}

	public static int[] sortind(int[] a){
		int[] ind=new int[a.length];
		int[][] v=new int[a.length][2];
		for(int i=0;i<a.length;i++){
			v[i][0]=a[i];
			v[i][1]=i;
		}
		int[] temp=new int[2];
		for(int i=0;i<a.length-1;i++){
			for(int j=0;j<a.length-i-1;j++)
				if(v[j+1][0]<v[j][0]){
					temp=v[j];    
					v[j]=v[j+1];
					v[j+1]=temp;

				}

		}
		for(int i=0;i<a.length;i++)
			ind[i]=v[i][1];
		return ind;
	}

	public static double[][][] grid(double[] x, double[] y){
		double[][][] grid= new double[2][y.length][x.length];
		for(int i=0;i<y.length;i++)
			for(int j=0;j<x.length;j++){
				grid[0][i][j]=x[j];
				grid[1][i][j]=y[i];}
		return grid;
	}
	public static void show(double[][] A){
		for(int i=0;i<A.length;i++){
			for(int j=0;j<A[0].length;j++)

				System.out.print(df.format(A[i][j])+"\t");
			System.out.println();
		}
		System.out.println();
	}
	public static void show(double[] v){
		for(int i=0;i<v.length;i++)
			System.out.println(df.format(v[i])+"\t");
		System.out.println();
	}

	public static void hshow(double[] v){
		for(int i=0;i<v.length;i++)
			System.out.print(df.format(v[i])+"\t");
		System.out.println();
	}

	public static double[] times(double[] v,double a){
		double[] y=new double[v.length];
		for(int i=0;i<v.length;i++)
			y[i]=a*v[i];
		return y;
	}

	public static double[][] times(double[][] M,double a){
		double[][] y=new double[M.length][M[0].length];
		for(int i=0;i<M.length;i++)
			for(int j=0;j<M[0].length;j++)
				y[i][j]=M[i][j]*a;
		return y;
	}


	public static double saw(double t)
	{

		double s=sin(t);
		return s;
		/*	double c=cos(t);

		double y=0;

		if(s>0 && c>0)
			y_
		int k=(int)(t);
		double rem=t-k;


			double y=0;
			if(rem<=.25)
			y=4*rem;
			else if(rem<=.75)
				y=2-4*rem;
			else 
				y=-4+4*rem;



		return y;*/
	}



	public static double triangWave(double cycle, double t)
	{


		double k=floor(t/cycle);
		double rem=t/cycle-k;


		double y=0;
		if(rem<=.25)
			y=4*rem;
		else if(rem<=.75)
			y=2-4*rem;
		else 
			y=-4+4*rem;



		return y;
	}

	public static void show(int[][] A){
		for(int i=0;i<A.length;i++){
			for(int j=0;j<A[0].length;j++)
				System.out.format("%d\t",A[i][j]);
			System.out.println();
		}
		System.out.println();
	}
	public static void show(int[] v){
		for(int i=0;i<v.length;i++)
			System.out.format("%d\n",v[i]);
		System.out.println();
	}
	public static void hshow(int[] v){
		for(int i=0;i<v.length;i++)
			System.out.format("%d\t",v[i]);
		System.out.println();
	}

	public static void show(byte[] v){
		for(int i=0;i<v.length;i++)
			System.out.format("%d\t",v[i]);
		System.out.println();
	}

	public static void show(boolean[] v){
		for(int i=0;i<v.length;i++)
			System.out.format("%s\t",v[i]);
		System.out.println();
	}

	public static void show(boolean[][] A){
		for(int i=0;i<A.length;i++){
			for(int j=0;j<A[0].length;j++)
				System.out.format("%s\t",A[i][j]);
			System.out.println();
		}
	}


	public static void show(String[] s){
		for(int i=0;i<s.length;i++){
			System.out.format("%s\n",s[i]);
		}
	}
	public static void hshow(String[] s){
		for(int i=0;i<s.length;i++){
			System.out.format("%s  ",s[i]);
		}
		System.out.println();
	}

	public static void pr(double a){

		System.out.println(a);

	}


	public static void pr(String a){

		System.out.println(a);
		if(logFile!=null){
			try{
				PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(logFile,true)));
				pw.println(a);
				pw.close();
			} catch(IOException e){System.err.println("IOException: " + e.getMessage());}
		}

	}
	public static void pr(int a){

		System.out.println(a);

	}
	public static void pr(boolean b){

		System.out.println(b);

	}


	public static void ph(double a){

		System.out.print(a);

	}
	public static void ph(int a){

		System.out.print(a);

	}

	public static void ph(String a){

		System.out.print(a);

	}


	public static void plot(Vect y){
		double[] x=new double[y.length];
		for(int i=0;i<x.length;i++)
			x[i]=i;
		plot(x,y.el);
	}

	public static void plot(double[] y){
		double[] x=new double[y.length];
		for(int i=0;i<x.length;i++)
			x[i]=i;
		plot(x,y);
	}
	
	public static void plot(String dataName, double[] y,String panelName){

		plot(dataName,y,Color.black,panelName);
	}
	
	public static void plot(String dataName,double[] y,Color color ,String panelName){
		util.pr(y.length);
		double[] x=new double[y.length];
		for(int i=0;i<x.length;i++)
			x[i]=i;
		plot(dataName,color,x,y,panelName);
	}


	public static void plot(Vect x, Vect y){
		plot(x.el,y.el);
	}

	public static void plot(double[] x, double[] y){

		plot("y=f(x)",Color.black,x,y);
	}

	public static void plot(Mat M){

		plot("y=f(x)",Color.black,M.el);
	}


	public static void plot(double[][] XY){

		plot("y=f(x)",Color.black,XY);
	}

	public static void plot(String name, Color c,double[] x, double[] y){


		plot(name,c,x,y,"Plot Panel");


	}
	
	public static void plot(String name, Color c,double[] x, double[] y, String panelName){

		double[][] A=new double[x.length][2];
		for(int i=0;i<x.length;i++){
			A[i][0]=x[i];
			A[i][1]=y[i];

		}

		plot(name,c,A,panelName);


	}



	public static void deleteDir(File dir) {
		if (dir.isDirectory()) {

			String[] children = dir.list();
			for (int i=0; i<children.length; i++) {
				deleteDir(new File(dir, children[i]));

			}
		}
		else
			dir.delete();



	}

	public static void plot(String name, Color c,double[][] XY){
		 plot(name, c,XY,"Plot Panel");
	}

	public static void plot(String name, Color c,double[][] XY,String panelName){

		Plot2DPanel plot = new Plot2DPanel();

		plot.setFont( new Font("Times New Roman", 1, 13));
		plot.addLinePlot(name, c, XY);
		//plot.setFont( new Font("Times New Roman", 1, 120));
		//	 util.pr(plot.getFont().toString());
		JFrame frame = new JFrame(panelName);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setSize(500,400);
		frame.setContentPane(plot);
		frame.setVisible(true);


	}

	public static void plotBunch(double[][] data){
		//	DecimalFormat df=new DecimalFormat("00.0");

		String[] titles = new String[data[0].length-1];

		for(int j=0;j<titles.length;j++)
			titles[j]="data "+j;

			
		plotBunch(data,titles);
		

	}
	
	public static void plotBunch(double[][] data, String[] titles){
		//	DecimalFormat df=new DecimalFormat("00.0");

		Plot2DPanel plot = new Plot2DPanel();
		double[] x=new double[data.length];
		double[] y=new double[data.length];

		for(int i=0;i<data.length;i++)
			x[i]=data[i][0];

		for(int j=0;j<data[0].length-1;j++){
			for(int i=0;i<x.length;i++)
				y[i]=data[i][j+1];
			plot.addLinePlot(titles[j], x, y);

		}

		plot.setAxisLabel(0,"x");
		plot.setAxisLabel(1,"y");
		plot.addLegend("EAST");

		JFrame frame = new JFrame("plot panel");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setSize(1000,800);
		frame.setContentPane(plot);
		frame.setVisible(true);

	}

	public static void plotBunch(double[][]... data){

		int n=data.length;
		String[] name=new String[n];
	
		for(int j=0;j<n;j++)
			name[j]="curve "+j;

		plotBunch(name,data);
	}
	

	public static void plotBunch(String[] name, double[][]... data){
		//	DecimalFormat df=new DecimalFormat("00.0");

		Plot2DPanel plot = new Plot2DPanel();
		double[] x,y;

		for(int j=0;j<data.length;j++){

			x=new double[data[j].length];
			y=new double[data[j].length];
			for(int i=0;i<x.length;i++){
				x[i]=data[j][i][0];
				y[i]=data[j][i][1];
			}


			plot.addLinePlot(name[j], x, y);

		}

		plot.setAxisLabel(0,"x");
		plot.setAxisLabel(1,"y");
		plot.addLegend("EAST");

		JFrame frame = new JFrame("plot panel");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setSize(500,400);
		frame.setContentPane(plot);
		frame.setVisible(true);

	}

	public static void plotBunch(String[] name, Mat[] data){
		plotBunch(name, data, data.length);

	}

	public static void plotBunch(Mat[] data){
		plotBunch(data, data.length);

	}

	public static void plotBunch(String[] name,Mat[] data, int n){
		//	DecimalFormat df=new DecimalFormat("00.0");

		Plot2DPanel plot = new Plot2DPanel();
		double[] x,y;

		for(int j=0;j<n;j++){

			x=new double[data[j].nRow];
			y=new double[data[j].nRow];
			for(int i=0;i<x.length;i++){
				x[i]=data[j].el[i][0];
				y[i]=data[j].el[i][1];
			}


			plot.addLinePlot(name[j], x, y);

		}

		plot.setAxisLabel(0,"x");
		plot.setAxisLabel(1,"y");
		plot.addLegend("EAST");

		JFrame frame = new JFrame("plot panel");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setSize(500,400);
		frame.setContentPane(plot);
		frame.setVisible(true);

	}

	public static void plotBunch(String[] name,Mat[] data, int n1, int n2){
		//	DecimalFormat df=new DecimalFormat("00.0");

		Plot2DPanel plot = new Plot2DPanel();
		double[] x,y;

		for(int j=n1;j<=n2;j++){

			x=new double[data[j].nRow];
			y=new double[data[j].nRow];
			for(int i=0;i<x.length;i++){
				x[i]=data[j].el[i][0];
				y[i]=data[j].el[i][1];
			}


			plot.addLinePlot(name[j-n1], x, y);

		}

		plot.setAxisLabel(0,"x");
		plot.setAxisLabel(1,"y");
		plot.addLegend("EAST");

		JFrame frame = new JFrame("plot panel");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setSize(500,400);
		frame.setContentPane(plot);
		frame.setVisible(true);

	}


	public static void plotBunch(Mat[] data, int n){


		String[] name=new String[n];
		for(int j=0;j<n;j++)
			name[j]="curve "+j;

		plotBunch(name,data,n);

	}

	public static void plotBunch(Mat[] data, int n1, int n2){

		int n=n2-n1+1;

		String[] name=new String[n];
		for(int j=0;j<n;j++)
			name[j]="curve "+j;

		plotBunch(name,data,n1,n2);

	}

	public static void plot(SpMat A){

		Plot2DPanel plot = new Plot2DPanel();
		int N=A.nRow;
		for(int i=0;i<A.nRow;i++){
			// int ir=N-1-i;
			int ir=i;
			int L=A.row[ir].nzLength;

			Vect x=new Vect(L);
			Vect y=new Vect(L);
			for(int j=0;j<L;j++){
				x.el[j]=(double)A.row[ir].index[j];
				y.el[j]=i;


			}

			//if(i==N-1)
			plot.addScatterPlot("",Color.red, x.el, y.el);


		}

		JFrame frame = new JFrame("a plot panel");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setSize(800,800);
		frame.setContentPane(plot);
		frame.setVisible(true);
	}

	public static Vect Atiken(Vect v2,Vect v1, Vect v){
		Vect Av=new Vect(v.length);
		for(int i=0;i<Av.length;i++)
			Av.el[i]=(v2.el[i]*v.el[i]-v1.el[i]*v1.el[i])/(v2.el[i]-2*v1.el[i]+v.el[i]);

		return Av;
	}
	public  static double  Atiken(double x2,double x1, double x){
		double Ax;

		Ax=(x2*x-x1*x1)/(x2-2*x1+x);

		return Ax;
	}


	public static void quickSort(double[] x){
		Sort.quick(x);

	}


	public static void quickSort(double[] x, int[] ind){
		Sort.quick(x,ind);

	}
	public static int search(int[] A,int ic,int a){
		int m=-1;
		for(int i=0;i<ic+1;i++){
			if(A[i]==a){
				m=i;
				break;
			}
		}
		return m;
	}

	public static int search(int[] A,int a){
		int m=-1;
		for(int i=0;i<A.length;i++){
			if(A[i]==a){
				m=i;
				break;
			}
		}
		return m;
	}

	static public double J1(double x) {

		double ax;
		double y;
		double ans1, ans2;

		if ( (ax = Math.abs(x)) < 8.0) {
			y = x * x;
			ans1 = x * (72362614232.0 + y * ( -7895059235.0 + y * (242396853.1
					+ y * ( -2972611.439 + y * (15704.48260 + y * ( -30.16036606))))));
			ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
					+ y * (99447.43394 + y * (376.9991397 + y * 1.0))));
			return ans1 / ans2;
		} else {
			double z = 8.0 / ax;
			double xx = ax - 2.356194491;
			y = z * z;

			ans1 = 1.0 + y * (0.183105e-2 + y * ( -0.3516396496e-4
					+
					y * (0.2457520174e-5 + y * ( -0.240337019e-6))));
			ans2 = 0.04687499995 + y * ( -0.2002690873e-3
					+ y * (0.8449199096e-5 + y * ( -0.88228987e-6
							+ y * 0.105787412e-6)));
			double ans = Math.sqrt(0.636619772 / ax) *
					(Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);
			if (x < 0.0) ans = -ans;
			return ans;
		}
	}

	public static int fact(int k){
		int f=1;

		for(int i=1;i<=k;i++)
			f*=i;

		return f;
	}

	public static String first(String line){

		String[] sp=line.split(regex);
		int b=0;
		while(b<sp.length-1 &&sp[b].equals("")){b++;}

		return sp[b];
	}



	public static String dropLeadingSpaces(String line){

		int L=line.length();

		int ix=0;

		while(line.charAt(ix++)==' ') {};

		String line2=String.copyValueOf(line.toCharArray(), ix-1, L-ix+1);

		return line2;
	}


	public static void wait(int ms){	try {
		Thread.sleep(ms);
	} catch (InterruptedException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	}


	public static Vect[] getHull(Vect[] p){
		return ConvexHull.getConvexHull(p);
	}


	public static SpMat loadSpMat(	String file){
		
		return loadSpMat(file,4,false);


	}
	

	public static SpMat loadSpMat(	String file,int skip, boolean zeroBase){
		
		String regex1="[ : ,()=\\t]+";
		
		SpMat Ms=null;

		Mat M=null;

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			int Nmax=1000000;
			int[][] rc=new int[Nmax][2];
			double[] vals=new double[Nmax];

			for(int k=0;k<skip;k++)
			line=br.readLine();

			int ix=-1;
			for(int i=0;i<2000000;i++){
				line=br.readLine();
				if(line==null || line.equals("")) break;
				sp=line.split(regex1);
				int ib=0;
				if(sp[0].equals("")) ib++;

				int row=Integer.parseInt(sp[ib]);
				int col=Integer.parseInt(sp[ib+1]);
				if(!zeroBase){
					row--;
					col--;
				}
				double val=Double.parseDouble(sp[ib+2]);
				ix++;
				rc[ix][0]=row;
				rc[ix][1]=col;
				vals[ix]=val;

				//util.pr(row+" , "+col+" , "+val);
			}

			int nRows=rc[ix][1]+1;

			 M=new Mat(nRows,nRows);
			for(int i=0;i<=ix;i++){

				int row=rc[i][0];

				int col=rc[i][1];

				double val=vals[i];
				if(abs(val)>1e-6)
				M.el[row][col]=val;

			}
			//M.show();


			 Ms=new SpMat(nRows);

			for(int i=0;i<nRows;i++){
				Vect v=M.getColVect(i);
				int size=0;
				for(int j=0;j<v.length;j++)
					if(v.el[j]!=0) size++;

				Ms.row[i]=new SpVect(nRows,size);
				ix=0;
				for(int j=0;j<v.length;j++){
					if(v.el[j]!=0) {
						Ms.row[i].el[ix]=v.el[j];	
						Ms.row[i].index[ix]=j;	
						ix++;
					}

				}
		
			}



			br.close();
			fr.close();
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading model file.");
		}


		 Ms=new SpMat(M);
		return Ms;

	}
	
public static Vect loadSpVect(String file,int L){
		
		String regex1="[ : ,()=\\t]+";
		
		Vect v=null;

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			int Nmax=1000000;
			int[][] rc=new int[Nmax][2];
			double[] vals=new double[Nmax];

			line=br.readLine();
			line=br.readLine();
			line=br.readLine();
	
			line=br.readLine();

			int ix=-1;
			for(int i=0;i<1000000;i++){
				line=br.readLine();
				if(line==null || line.equals("")) break;
				sp=line.split(regex1);
				int ib=0;
				if(sp[0].equals("")) ib++;

				int row=Integer.parseInt(sp[ib])-1;
				int col=Integer.parseInt(sp[ib+1])-1;
				double val=Double.parseDouble(sp[ib+2]);
				ix++;
				rc[ix][0]=row;
				rc[ix][1]=col;
				vals[ix]=val;

				//util.pr(row+" , "+col+" , "+val);
			}


			 v=new Vect(L);
			for(int i=0;i<=ix;i++){

				int row=rc[i][0];

				double val=vals[i];
			//	if(val!=0) util.pr(row+" "+val);
				v.el[row]=val;

			}

			br.close();
			fr.close();
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading model file.");
		}



		return v;

	}

	public static void setLogFile(String path){
		
		logFile=path;
	}
	
	public static String getLogFile(){
		
		return logFile;
	}
}
