package emsolution;

import io.Loader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Random;

import math.Mat;
import math.Vect;
import math.util;

public class PlayModel2D {

	public Mat shapeFunc;
	public Preisach2D ps;
	//public int nSym, nDesc,nAsc, nInit, nMajor,nTotCurves;
	int[] nInitial,nMajor,nSymLoop,nDescending,nAscending,nTotCurves,nAni,Lani;
	public boolean anisot;
	public Mat[] BHraw,BH;
	public double[] zk,pk;
	public int dim=2,nHyst;



	public PlayModel2D()
	{	
		int kr=0;
		double BsM=1.8;
		double Bseff=1.7;
		double Hs=430;

		int Mp=1000;
		double Hs0=100;
		double mean=.2*Hs0;
		double width=.2*Hs0;

		ps=new Preisach2D(Mp,mean,width,Hs,BsM,Bseff,3564656);
		ps.kRotated=kr;


	}

	public static void main2(String[] args)
	{

		PlayModel2D pm=new PlayModel2D();

		//String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		String file="C:\\Works\\HVID\\hys_dataHy";

		pm.createData(file);
		//pm.createDataMath(file);
		//pm.rotation();

		//pm.getHRotationalBinput();
	//	pm.getBRotationalHinput();
		/*		try {
			Thread.sleep(200);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/

		//pm.loadData(file);


	}
	

public void createData(String file){


	int nSet=18;
	int nInit=1;
	int nMajor=1;
	int nSymLoops=14;
	int nDescending=0;
	int nAscending=0;
	int nAni=18;
	int Lani=18;

	double dang=180.0/nSet;

	int nTot=nInit+nMajor+nSymLoops+nDescending+nAscending;

	double Bs=ps.Bseff;
	double Hs=ps.Hs;


	int LL=1000;

	double Bseff=ps.Bseff;

	Mat[][] BHs=new Mat[nSet][nTot];
	Mat[][] BHsr=new Mat[nSet][nTot];
	Mat[]BHsAv=new Mat[nTot];
	Mat[] BHsrAv=new Mat[nTot];
	Mat[] BHani=new Mat[nAni];
	Vect Bani=new Vect();

	Vect Bdiv=new Vect().linspace(Bs, 0, nSymLoops+2);


	for(int ia=0;ia<nSet;ia++){

		ps.kRotated=ia;

		int ix=0;

		for(int i=0;i<nInit;i++){


			Mat BHtemp=ps.initial(Bs,LL);

			int L=50;
			Vect B0=new Vect().linspace(0,.1,4);
			Vect B1=new Vect().linspace(.125,Bseff,L-4);
			Vect B=B0.aug(B1);

			BHs[ia][ix]=new Mat(L,3);

			BHs[ia][ix]=new Mat(L,3);
			for(int j=0;j<L;j++){
				BHs[ia][ix].el[j][0]=getH(BHtemp,B.el[j]);
				BHs[ia][ix].el[j][1]=B.el[j];

			}

			//	util.plot(BHs[ix].getColVect(0),BHs[ix].getColVect(1));

			BHsr[ia][ix]=this.getHrBr(BHs[ia][ix]);



			ix++;

		}




		for(int i=0;i<nMajor;i++){

			Mat BHtemp=ps.symMajorDesc(LL);

			int L=100;
			Vect B=new Vect().linspace(Bdiv.el[0],-Bdiv.el[0],L);

			BHs[ia][ix]=new Mat(L,3);
			for(int j=0;j<L;j++){
				BHs[ia][ix].el[j][0]=getH(BHtemp,B.el[j]);
				BHs[ia][ix].el[j][1]=B.el[j];

			}


			BHsr[ia][ix]=this.getHBij(BHs[ia][ix],0);


			ix++;

		}

		for(int i=0;i<nSymLoops;i++){

			//double Bp=Bseff-(i+1)*Bseff/(nSymLoops+1);	

			double Bp=Bdiv.el[i+1];

			Mat BHtemp=ps.symDesc(Bp,LL);

			int L=100;
			Vect B=new Vect().linspace(Bp,-Bp,L);
			BHs[ia][ix]=new Mat(L,3);
			for(int j=0;j<L;j++){
				BHs[ia][ix].el[j][0]=getH(BHtemp,B.el[j]);
				BHs[ia][ix].el[j][1]=B.el[j];

			}


			BHsr[ia][ix]=this.getHBij(BHs[ia][ix],0);



			ix++;



		}


		for(int i=0;i<nDescending;i++){

		}

		for(int i=0;i<nAscending;i++){

		}

	}


	if(nTot>0){
		util.plotBunch(BHsr[0]);
		//	util.plotBunch(BHsr[1]);
	}


	for(int i=0;i<nTot;i++){
		Mat M=BHs[0][i].deepCopy();
		for(int j=1;j<nSet;j++)
			M=M.add(BHs[j][i]);
		BHsAv[i]=M.times(1.0/nSet);

		BHsrAv[i]=this.getHBij(BHsAv[i],0);
	}

	//	util.plotBunch(BHsrAv);


	if(nAni>0){

		Bani=new Vect().linspace(0, 1.7, Lani);

		Mat[] BHaniT=new Mat[nAni];


		for(int i=0;i<nAni;i++){

			double phiRad=i*dang*Math.PI/180;
			//phiRad=0*Math.PI/2;

			Vect Hr=new Vect().linspace(0,1.8*Hs,LL);

			/*		int Lp=2000;
			Mat Hp=new Mat(Lp,2);
			for(int j=0;j<Lp;j++){
				Hp.el[j][0]=.8*Hs*Math.cos(4*j*Math.PI/Lp)*Math.cos(phiRad);
				Hp.el[j][1]=.8*Hs*Math.sin(4*j*Math.PI/Lp)*Math.cos(phiRad);
			}*/


			Mat H=new Mat(LL,2);
			for(int j=0;j<LL;j++){
				H.el[j][0]=Hr.el[j]*Math.cos(phiRad);
				H.el[j][1]=Hr.el[j]*Math.sin(phiRad);
			}


			for(int k=0;k<ps.M;k++)
				for(int j=0;j<ps.nphi;j++){
					ps.on[k][j]=false;
				}

			ps.demagnetize();

			/*	util.pr(M.nCol);

			util.plot(getHBij(M,0).el);
			 */

			//	ps.getRes().hshow();

			//Mat BH2=ps.getLocus(Hp);
			Mat BH2=ps.getLocus(H);

			Mat BH3=this.getHrBr(BH2);

			//util.plot(getHBij(BH2,3).el);
			//	util.plot(BH2.getColVect(0),BH2.getColVect(2));


			BHaniT[i]=BH3.deepCopy();

			BHaniT[i]=new Mat(Lani,2);

			BHani[i]=new Mat(Lani,2);


			for(int j=0;j<Lani;j++){
				double Ht=this.getH(BH3, Bani.el[j]);
				BHaniT[i].el[j][0]=Ht;
				BHaniT[i].el[j][1]=Bani.el[j];
				BHani[i].el[j][0]=Ht*Math.cos(phiRad);
				BHani[i].el[j][1]=Ht*Math.sin(phiRad);

			}

		}


		//	Bani.show();




		util.plotBunch(BHaniT);
	}
	//util.plotBunch(BHs);

	boolean write =true;
	if(write){
	
		
			String file1="C:\\Works\\HVID\\hys_dataHXAv";
		this.writeHystDataAv(BHsAv, BHani, file1);
			String file2="C:\\Works\\HVID\\hys_dataHX";
		this.writeHystData(BHs, file2);
		
	}



}

public void writeHystData(Mat[][] BHs, String file)
{
	 writeHystData(BHs, file,1);
}

public void writeHystData(Mat[][] BHs, String file,int nInit){
	
	int nSet=BHs.length;

	
	Mat[] BHani=new Mat[1];
	int nAni=0;

	int nMajor=1;
	int nTot=BHs[0].length;
	int nSymLoops=nTot-nMajor-nInit;
	int nDescending=0;
	int nAscending=0;
	


	DecimalFormat dfB=new DecimalFormat("#.000");
	DecimalFormat dfH=new DecimalFormat("#.0");


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		
			
	
			double Hsefft=1000;
			double Bsefft=1;
			if(nInit==1)
				Bsefft=BHs[0][0].el[BHs[0][0].nRow-1][1];
			else
				Bsefft=BHs[0][0].el[0][1];

			
			int ioVariable=1;
			int identificationMethod=1;
			
			pwBun.println(ioVariable+"\t"+identificationMethod+"\t"+nSet);
			pwBun.println("*Bs*Hs*");
			pwBun.print(Bsefft+"\t");
			for(int ia=0;ia<nSet;ia++){

				if(nInit==1)
					Hsefft=BHs[ia][0].el[BHs[ia][0].nRow-1][0];
				else
					Hsefft=BHs[0][0].el[0][0];

				pwBun.print(dfH.format(Hsefft)+"\t");
			}
			pwBun.println();

			pwBun.println("* ---- Number of initial, major, symmetric, descending and asending currves, respectively:");
			pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);
			pwBun.println("* angles:");
			for(int ia=0;ia<nSet;ia++){
				pwBun.print(dfH.format(ia*10)+"\t");
			}
		
			pwBun.println();
			String str="";
			for(int i=0;i<nTot;i++){
				
				if(i==nInit-1){
					str="* ----Initial curve";
					//pwBun.println(str);
				}
				else if(i==nInit+nMajor-1) 
					{
					str="* ----Major loop";
					pwBun.println(str);
					}
				else if(i==nInit+nMajor) {
					str="* ----Symmetrical loops";
					pwBun.println(str);
				}
				else if(i==nInit+nMajor+nSymLoops) {
					str="* ----Descending loops";
					pwBun.println(str);
				}
				else if(i==nInit+nMajor+nSymLoops) {
					str="* ----Ascending loops";
					pwBun.println(str);
				}

				pwBun.println(BHs[0][i].nRow);
			
				for(int j=0;j<BHs[0][i].nRow;j++){
					for(int ia=0;ia<nSet;ia++){
						pwBun.print(dfH.format(BHs[ia][i].el[j][0])+"\t");
					}
					pwBun.println(dfB.format(BHs[0][i].el[j][1]));
				}
			}

			pwBun.println("* ----- Rotational hysteresis loss");
			pwBun.println("* Number of data *");
			pwBun.println("0");
			pwBun.println("* B * loss");
			pwBun.println("* ----- Anisotropy");
			pwBun.println("* Number B magnitude * Number of angles *");
			pwBun.println(0+"\t"+0); 		//	pwBun.println(Lani+"\t"+nAni);
			pwBun.println("* B * H ﾂｷﾂｷﾂｷﾂｷﾂｷ * easy axis of magnetization*");
			pwBun.println("* B * H ﾂｷﾂｷﾂｷﾂｷﾂｷ * hard axis of magnetization*");
			


			util.pr("Simulated angle-dependent hysteresis data was written to "+file+".");

			pwBun.close();
		}
		catch(IOException e){}
		

}

public void writeHystDataXY(Mat[][] BHsx,Mat[][] BHsy, String file,int nInit){
	
	int nSet=BHsx.length;

	
	Mat[] BHani=new Mat[1];
	int nAni=0;

	int nMajor=1;
	int nTot=BHsx[0].length;
	int nSymLoops=nTot-nMajor-nInit;
	int nDescending=0;
	int nAscending=0;
	


	DecimalFormat dfB=new DecimalFormat("#.000");
	DecimalFormat dfH=new DecimalFormat("#.0");


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		
			
			double Hsefft=1000;
			double Bsefft=1;
			if(nInit==1)
				Bsefft=BHsx[0][0].el[BHsx[0][0].nRow-1][1];
			else
				Bsefft=BHsx[0][0].el[0][1];

	
			pwBun.println(1+"\t"+1+"\t"+nSet);
			pwBun.println("*Bs*Hs*");
			pwBun.print(Bsefft+"\t");
			for(int ia=0;ia<nSet;ia++){

				if(nInit==1)
					Hsefft=BHsx[ia][0].el[BHsx[ia][0].nRow-1][0];
				else
					Hsefft=BHsx[0][0].el[0][0];

				pwBun.print(dfH.format(Hsefft)+"\t");
			}
			pwBun.println();

			pwBun.println("*xxx *");
			pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);
			pwBun.println("*angles *");
			for(int ia=0;ia<nSet;ia++){
				pwBun.print(dfH.format(ia*15)+"\t");
			}
			pwBun.println();
			for(int i=0;i<nTot;i++){
				pwBun.println("*xxx");
				pwBun.println(BHsx[0][i].nRow);
				for(int j=0;j<BHsx[0][i].nRow;j++){
					for(int ia=0;ia<nSet;ia++){
						pwBun.print(dfH.format(BHsx[ia][i].el[j][0])+"\t");
						pwBun.print(dfH.format(BHsy[ia][i].el[j][0])+"\t");
		
					}
					pwBun.println(dfB.format(BHsx[0][i].el[j][1]));
				}
			}

			pwBun.println("*----- xxx");
			pwBun.println("*xx");
			pwBun.println("0");
			pwBun.println("*xx");
			pwBun.println("*xx");
			pwBun.println("*xx *");
			pwBun.println(0+"\t"+0); 		//	pwBun.println(Lani+"\t"+nAni);
			pwBun.println("*xx*");
			


			util.pr("Simulated angle-dependent hysteresis data was written to "+file+".");

			pwBun.close();
		}
		catch(IOException e){}
		

}


public void createDataMath(String file){


	double c1=3;
	double c2=.01;

	double c3=.1;
	
	double c4=.02;
	
	
	int nSet=18;
	int nInit=1;
	int nMajor=1;
	int nSymLoops=14;
	int nDescending=0;
	int nAscending=0;
	int nAni=18;
	int Lani=50;

	double dang=180.0/nSet;

	int nTot=nInit+nMajor+nSymLoops+nDescending+nAscending;

	double Bs=ps.Bseff;
	double Hs=ps.Hs;


	int LL=1000;

	double Bseff=ps.Bseff;

	Mat[][] BHs=new Mat[nSet][nTot];
	Mat[][] BHsr=new Mat[nSet][nTot];
	Mat[]BHsAv=new Mat[nTot];
	Mat[] BHsrAv=new Mat[nTot];
	Mat[] BHani=new Mat[nAni];
	Vect Bani=new Vect();

	Vect Bdiv=new Vect().linspace(Bs, 0, nSymLoops+2);


	for(int ia=0;ia<nSet;ia++){

		ps.kRotated=ia;
		
		

		int ix=0;

		for(int i=0;i<nInit;i++){


			int L=50;
			Vect B=new Vect().linspace(0,Bseff,L);

			BHs[ia][ix]=new Mat(L,3);

	
			for(int j=0;j<L;j++){
				BHs[ia][ix].el[j][0]=Hs*Math.sinh(B.el[j]*(c1/(1+c2*ia)))*B.el[j];
				BHs[ia][ix].el[j][1]=B.el[j];

			}
			
			

			//	util.plot(BHs[ix].getColVect(0),BHs[ix].getColVect(1));

			BHsr[ia][ix]=BHs[ia][ix];



			ix++;

		}




		for(int i=0;i<nMajor;i++){


			int L=100;
			Vect B=new Vect().linspace(Bdiv.el[0],-Bdiv.el[0],L);

			BHs[ia][ix]=new Mat(L,3);
			for(int j=0;j<L;j++){
			//	BHs[ia][ix].el[j][0]=Hs*Math.sinh(B.el[j]*(c1/(1+c2*ia)))-c3/(1+c4*ia)*Math.pow(Bdiv.el[0]-B.el[j],2)*Math.pow(B.el[j]+Bdiv.el[0],2);
				BHs[ia][ix].el[j][0]=Hs*Math.sinh(B.el[j]*(c1/(1+c2*ia)));//-c3/(1+c4*ia)*Math.abs((Bdiv.el[0]-B.el[j])*(B.el[j]+Bdiv.el[0]));
				
				BHs[ia][ix].el[j][1]=B.el[j]+c3/(1+c4*ia)*Math.pow(Bdiv.el[0]-B.el[j],2)*Math.pow(B.el[j]+Bdiv.el[0],2);;

			}


			BHsr[ia][ix]=BHs[ia][ix];


			ix++;

		}

		for(int i=0;i<nSymLoops;i++){

			//double Bp=Bseff-(i+1)*Bseff/(nSymLoops+1);	

			double Bp=Bdiv.el[i+1];


			int L=100;
			Vect B=new Vect().linspace(Bp,-Bp,L);
			BHs[ia][ix]=new Mat(L,3);
			for(int j=0;j<L;j++){
				//BHs[ia][ix].el[j][0]=Hs*Math.sinh(B.el[j]*(c1/(1+c2*ia)))-c3/(1+c4*ia)*Math.pow(Bp-B.el[j],2)*Math.pow(B.el[j]+Bp,2);
				
				BHs[ia][ix].el[j][0]=Hs*Math.sinh(B.el[j]*(c1/(1+c2*ia)));

				BHs[ia][ix].el[j][1]=B.el[j]+c3/(1+c4*ia)*Math.pow(Bp-B.el[j],2)*Math.pow(B.el[j]+Bp,2);

			}


			BHsr[ia][ix]=BHs[ia][ix];



			ix++;



		}


		for(int i=0;i<nDescending;i++){

		}

		for(int i=0;i<nAscending;i++){

		}

	}


	if(nTot>0){
		util.plotBunch(BHsr[0]);
		//util.plotBunch(BHsr[9]);
	}


	for(int i=0;i<nTot;i++){
		Mat M=BHs[0][i].deepCopy();
		for(int j=1;j<nSet;j++)
			M=M.add(BHs[j][i]);
		BHsAv[i]=M.times(1.0/nSet);

		BHsrAv[i]=this.getHBij(BHsAv[i],0);
	}

	//	util.plotBunch(BHsrAv);


	if(nAni>0){

		Bani=new Vect().linspace(0, 1.7, Lani);

		Mat[] BHaniT=new Mat[nAni];


		for(int i=0;i<nAni;i++){

			BHaniT[i]=new Mat(Lani,2);
			double phiRad=i*dang*Math.PI/180;
			//phiRad=0*Math.PI/2;

		
			for(int j=0;j<Lani;j++){
			BHaniT[i].el[j][0]=Hs*Math.sinh(Bani.el[j]*(c1/(1+c2*i)));
			BHaniT[i].el[j][1]=Bani.el[j];
			}
			

			BHani[i]=new Mat(Lani,4); 
			
			for(int j=0;j<Lani;j++){
				BHani[i].el[j][0]=BHaniT[i].el[j][0]*Math.cos(phiRad);
				BHani[i].el[j][1]=BHaniT[i].el[j][0]*Math.sin(phiRad);
				BHani[i].el[j][2]=BHaniT[i].el[j][1]*Math.cos(phiRad);
				BHani[i].el[j][3]=BHaniT[i].el[j][1]*Math.sin(phiRad);
			}


		}


		util.plotBunch(BHaniT);
	}

	boolean write =true;
	if(write){
		
			String file1="C:\\Works\\HVID\\hys_dataHXAv";
		this.writeHystDataAv(BHsAv, BHani, file1);
			String file2="C:\\Works\\HVID\\hys_dataHX";
		this.writeHystData(BHs, file2);
		
	}


/*	if(write){


		double Hseff=BHs[0][0].el[BHs[0][0].nRow-1][0];

		DecimalFormat dfB=new DecimalFormat("#.00");
		DecimalFormat dfH=new DecimalFormat("#.0");


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		

			for(int ia=0;ia<nSet;ia++){


				pwBun.println(1+"\t"+1+"\t"+nSet+"\t"+ps.phi[ia]);
				pwBun.println("*Bs*Hs*");
				pwBun.println(Bseff+"\t"+Hseff);

				pwBun.println("* �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｾ域腸�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽW�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ[�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ[�ｿｽ�ｽｿ�ｽｽv�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ�ｾ守ｧｰ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ[�ｿｽ�ｽｿ�ｽｽv�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ~�ｿｽ�ｽｿ�ｽｽ�ｾ域腸�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ郢ｽ�ｽｸ�ｿｽ�ｽｿ�ｽｽ�ｾ域腸�ｽｿ�ｽｽ *");
				pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);

				for(int i=0;i<nTot;i++){
					pwBun.println("*xxx");
					pwBun.println(BHs[ia][i].nRow);
					for(int j=0;j<BHs[ia][i].nRow;j++)
						pwBun.println(BHs[ia][i].el[j][0]+"\t"+BHs[ia][i].el[j][1]);
				}

				pwBun.println("* ----- �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ]�ｿｽ�ｽｿ�ｽｽq�ｿｽ�ｽｿ�ｽｽX�ｿｽ�ｽｿ�ｽｽe�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽV�ｿｽ�ｽｿ�ｽｽX�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ");
				pwBun.println("* B�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *");
				pwBun.println("0");
				pwBun.println("* B * �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ");
				pwBun.println("* ----- �ｿｽ�ｽｿ�ｽｽ�ｾ吝桁�ｽｿ�ｽｽ");
				pwBun.println("* B�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽp�ｿｽ�ｽｿ�ｽｽx�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *");
				pwBun.println(0+"\t"+0); 		//	pwBun.println(Lani+"\t"+nAni);
				pwBun.println("* B * H �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *�ｿｽ�ｽｿ�ｽｽ@�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽe�ｿｽ�ｽｿ�ｽｽ�ｾ戊ｶ｣�ｽｿ�ｽｽ");

				for(int i=0;i<0*Lani;i++){
					pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
					for(int j=0;j<nAni;j++){
						pwBun.print(dfH.format(BHani[j].el[i][0])+"\t");
					}
					pwBun.println();
				}
				pwBun.println("* B * H �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *�ｿｽ�ｽｿ�ｽｽ@�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ霆ｸ");
				for(int i=0;i<0*Lani;i++){
					pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
					for(int j=0;j<nAni;j++){
						pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
					}
					pwBun.println();
				}
	
			}


			util.pr("Simulated angle-dependent hysteresis data was written to "+file+".");

			pwBun.close();
		}
		catch(IOException e){}



		String fileAv="C:\\Works\\HVID\\hys_dataHAvy";


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fileAv)));		


			pwBun.println(1+"\t"+1+"\t"+1+"\t"+0);
			pwBun.println("*Bs*Hs*");
			pwBun.println(Bseff+"\t"+Hseff);

			pwBun.println("* �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｾ域腸�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽW�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ[�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ[�ｿｽ�ｽｿ�ｽｽv�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ�ｾ守ｧｰ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ[�ｿｽ�ｽｿ�ｽｽv�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ~�ｿｽ�ｽｿ�ｽｽ�ｾ域腸�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ郢ｽ�ｽｸ�ｿｽ�ｽｿ�ｽｽ�ｾ域腸�ｽｿ�ｽｽ *");
			pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);

			for(int i=0;i<nTot;i++){
				pwBun.println("*xxx");
				pwBun.println(BHsAv[i].nRow);
				for(int j=0;j<BHsAv[i].nRow;j++)
					pwBun.println(BHsAv[i].el[j][0]+"\t"+BHsAv[i].el[j][1]);
			}

			pwBun.println("* ----- �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ]�ｿｽ�ｽｿ�ｽｽq�ｿｽ�ｽｿ�ｽｽX�ｿｽ�ｽｿ�ｽｽe�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽV�ｿｽ�ｽｿ�ｽｽX�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ");
			pwBun.println("* B�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *");
			pwBun.println("0");
			pwBun.println("* B * �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ");
			pwBun.println("* ----- �ｿｽ�ｽｿ�ｽｽ�ｾ吝桁�ｽｿ�ｽｽ");
			pwBun.println("* B�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽp�ｿｽ�ｽｿ�ｽｽx�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *");
			pwBun.println(Lani+"\t"+nAni);
			pwBun.println("* B * H �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *�ｿｽ�ｽｿ�ｽｽ@�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽe�ｿｽ�ｽｿ�ｽｽ�ｾ戊ｶ｣�ｽｿ�ｽｽ");

			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(Bani.el[i])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][0])+"\t");
				}
				pwBun.println();
			}
			pwBun.println("* B * H �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *�ｿｽ�ｽｿ�ｽｽ@�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ霆ｸ");
			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(Bani.el[i])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
				}
				pwBun.println();
			}
						pwBun.println();
				pwBun.println("End of hysteresis data set "+ia);
				pwBun.println();



			util.pr("Simulated angle-averaged hysteresis data was written to "+fileAv+".");

			pwBun.close();
		}
		catch(IOException e){}



	}
*/


}

	public void getBRotationalHinput(){

		int iang=0;
		double phiRad=iang*Math.PI/18;;

		int steps = 360;
		Mat  X= new Mat(steps, 2);


		double yxRatio = 0;
		double a11 = Math.cos(phiRad);
		double a21 = Math.sin(phiRad);
		double a12 = -a21*yxRatio;
		double a22 = a11*yxRatio;


		double Xm = 75;

		for (int i = 0; i < steps; ++i){
			double wt = 4 * Math.PI*i / steps;
		
			double xx = Xm* Math.cos(wt);
			double yy = Xm* Math.sin(wt);

			double xxt = a11*xx + a12*yy;
			double yyt = a21*xx + a22*yy;

			X.el[i][0] = xxt;
			X.el[i][1] = yyt;
	
		}
		
		//util.plot(X);

		//X.show();


		iang=0;
		
		//Mat[] MM=new Mat[ps.nphi];
		
		Mat BH=ps.getLocus(X);
		
		Mat MM=this.getHBij(BH, 0);
		Mat NN=this.getHBij(BH, 1);
		//util.plot(BH.getColVect(2));
		
		int L=BH.nRow/2;
		Mat BHr=new Mat(L,2);
		
	//	HH.show();
		
			Vect Hr=new Vect(L);
		Vect Br=new Vect(L);
		
	
			
	double ang=Math.atan(BH.el[1][1]/BH.el[1][0]);
	//util.pr(ang/Math.PI*180);
	//util.pr(ang/Math.PI*180);
	
	Vect er=new Vect(Math.cos(ang),Math.sin(ang));
	for(int j=0;j<Hr.length;j++){
		int jx=j+L;
/*		Hr.el[j]=new Vect(BH.el[jx][0],BH.el[jx][1]).dot(er);
		Br.el[j]=new Vect(BH.el[jx][2],BH.el[jx][3]).dot(er);*/
		Hr.el[j]=new Vect(BH.el[jx][0],BH.el[jx][1]).dot(er);
		Br.el[j]=new Vect(BH.el[jx][2],BH.el[jx][3]).dot(er);
		//util.pr(Hr.el[i]+"\t"+Br.el[i]);x
	}

/*	Hr=BH.getColVect(0);
	Br=BH.getColVect(2);*/

	BHr.setCol(Hr,0);
	BHr.setCol(Br,1);

	util.plot(BHr);
	
	double loss=0;
	for(int i=0;i<MM.nRow-1;i++){
		Vect dB=new Vect(NN.el[i+1]).sub(new Vect(NN.el[i]));
		loss+=new Vect(MM.el[i]).dot(dB);
	}
	util.pr("Alternating loss per cycle: "+loss);


	//	Mat HH=new Mat(new Loader().loadArrays(360,2,"C:\\Works\\HVID\\hpath.txt"));	
		//MM[i]=this.getHBij(BH,0);

	//	Mat BH=ps.getLocus(HH);

		//BH.show();

	//	util.plot(BH.getColVect(2),BH.getColVect(3));

	//	util.plot(MM);

	}



	public void getHRotationalBinput(){


		Mat BH=ps.getLoopBinput(1,.5);

	//BH.show();
	Mat MM=this.getHBij(BH,0);
	Mat NN=this.getHBij(BH,1);
	
	for(int i=0;i<0;i++){
	MM.el[i*10][0]*=2;
	MM.el[i*10][1]*=2;
	NN.el[i*10][0]*=2;
	NN.el[i*10][1]*=2;
	}
	//MM.show();
	util.plot(MM);
	
	double loss=0;
	for(int i=0;i<MM.nRow-1;i++){
		Vect dB=new Vect(NN.el[i+1]).sub(new Vect(NN.el[i]));
		loss+=new Vect(MM.el[i]).dot(dB);
	}
	util.pr("Rotational loss per cycle: "+loss);
	util.plot(NN);
	MM.show();

}




public void writeHystDataAv(Mat[] BHs,Mat[] BHani, String file){
	
	int nSet=1;
	
	double Bseff= BHs[0].el[BHs[0].nRow-1][1];
	
	int nAni=BHani.length;
	int Lani=BHani[0].nRow;
	int nInit=1;
	int nMajor=1;
	int nTot=BHs.length;
	int nSymLoops=nTot-2;
	int nDescending=0;
	int nAscending=0;
	
	double Hseff=BHs[0].el[BHs[0].nRow-1][0];

	DecimalFormat dfB=new DecimalFormat("#.00");
	DecimalFormat dfH=new DecimalFormat("#.0");


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		


			pwBun.println(1+"\t"+1+"\t"+nSet+"\t"+0);
			pwBun.println("*Bs*Hs*");
			pwBun.println(Bseff+"\t"+Hseff);

			pwBun.println("* �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｾ域腸�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽW�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ[�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ[�ｿｽ�ｽｿ�ｽｽv�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ�ｾ守ｧｰ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ[�ｿｽ�ｽｿ�ｽｽv�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ~�ｿｽ�ｽｿ�ｽｽ�ｾ域腸�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽ郢ｽ�ｽｸ�ｿｽ�ｽｿ�ｽｽ�ｾ域腸�ｽｿ�ｽｽ *");
			pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);

			for(int i=0;i<nTot;i++){
				pwBun.println("*xxx");
				pwBun.println(BHs[i].nRow);
				for(int j=0;j<BHs[i].nRow;j++)
					pwBun.println(BHs[i].el[j][0]+"\t"+BHs[i].el[j][1]);
			}

			pwBun.println("* ----- �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ]�ｿｽ�ｽｿ�ｽｽq�ｿｽ�ｽｿ�ｽｽX�ｿｽ�ｽｿ�ｽｽe�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽV�ｿｽ�ｽｿ�ｽｽX�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ");
			pwBun.println("* B�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *");
			pwBun.println("0");
			pwBun.println("* B * �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ");
			pwBun.println("* ----- �ｿｽ�ｽｿ�ｽｽ�ｾ吝桁�ｽｿ�ｽｽ");
			pwBun.println("* B�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ * �ｿｽ�ｽｿ�ｽｽp�ｿｽ�ｽｿ�ｽｽx�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *");
			pwBun.println(Lani+"\t"+nAni);
			pwBun.println("* B * H �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *�ｿｽ�ｽｿ�ｽｽ@�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽe�ｿｽ�ｽｿ�ｽｽ�ｾ戊ｶ｣�ｽｿ�ｽｽ");
			
			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
				}
				pwBun.println();
			}
			pwBun.println("* B * H �ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ *�ｿｽ�ｽｿ�ｽｽ@�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ�ｽｿ�ｽｽ�ｿｽ霆ｸ");
			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
				}
				pwBun.println();
			}
			BHani[0].show();

			util.pr("Simulated angle-dependent hysteresis data was written to "+file+".");

			pwBun.close();
		}
		catch(IOException e){}
		

}




public void rotation(){}


public void loadData(String file){}



public  double getH(Mat BH,double B){
	return getH(BH,B,0);
}

public  double getH(Mat BH,double B, int mode){

	int ih, ib;
	if(mode==0){
		ih=0;ib=1;
	}
	else{
		ih=0;ib=2;
	}

	boolean ascending=true;
	
	

	int im=BH.nRow-1;

	if (BH.el[im][ih]<BH.el[0][ih]) ascending=false;
	 
	int i1=0;
	int i2=BH.nRow-1;

	if(ascending){

		if(B<=BH.el[i1][ib]){
			double cc=(BH.el[i1+1][ih]-BH.el[i1][ih])/(BH.el[i1+1][ib]-BH.el[i1][ib]);
			return BH.el[i1][ih];//+cc*(B-BH.el[i1][ib]);
			//return BH.el[i1][ih];
		}
		if(B>=BH.el[i2][ib]){
			//return BH.el[i2][ih];
			double cc=(BH.el[i2][ih]-BH.el[i2-1][ih])/(BH.el[i2][ib]-BH.el[i2-1][ib]);
			return BH.el[i2][ih];//+cc*(B-BH.el[i2][ib]);
			//return BH.el[i][ih]+(B-BH.el[i][ib])*(BH.el[i+1][ih]-BH.el[i][ih])/(BH.el[i+1][ib]-BH.el[i][ib]);
		}
		/*i=BH.nRow-1;
		if(B>=BH.el[i][ib])
			return BH.el[i][ih]+(B-BH.el[i][ib])*(BH.el[i][ih]-BH.el[i-1][ih])/(BH.el[i][ib]-BH.el[i-1][ib]);
		 */
	}
	else{

	
		if(B<=BH.el[i2][ib]){
			double cc=(BH.el[i2][ih]-BH.el[i2-1][ih])/(BH.el[i2][ib]-BH.el[i2-1][ib]);
		
			return BH.el[i2][ih];//+cc*(B-BH.el[i2][ib]);
			//return BH.el[i1][ih];
		}
		if(B>=BH.el[i1][ib]){

			
			//return BH.el[i2][ih];
			double cc=(BH.el[i1+1][ih]-BH.el[i1][ih])/(BH.el[i1+1][ib]-BH.el[i1][ib]);


			return BH.el[i1][ih];//+cc*(B-BH.el[i1][ib]);
		}

	}

	int j=0;

	if(ascending)
		while(BH.el[j+1][ib]<B){j++;}
	else
		while(BH.el[j+1][ib]>B){j++;}



	double cc=(BH.el[j+1][ih]-BH.el[j][ih])/(BH.el[j+1][ib]-BH.el[j][ib]);

	double H=BH.el[j][ih]+(B-BH.el[j][ib])*cc;

	return H;

}

public int getj(double[] array, double p){

	int j=0;
	if(array.length==1) return j;
	while(array[j+1]<p){j++;}
	return j;
}

public double shapeFuncAt(double z, double p){

	double f=0;
	int jz=-1;
	while(this.zk[jz+1]>z){jz++;}

	int jp=-1;
	while(this.pk[jp+1]>z){jp++;}


	double fQ11=shapeFunc(jz,jp);
	double fQ12=shapeFunc(jz+1,jp);
	double fQ21=shapeFunc(jz,jp+1);
	double fQ22=shapeFunc(jz+1,jp+1);

	f+=fQ11*(this.pk[jp+1]-p)*(this.zk[jz+1]-z);

	f+=fQ21*(p-this.pk[jp])*(this.zk[jz+1]-z);

	f+=fQ12*(this.pk[jp+1]-p)*(z-this.zk[jz]);

	f+=fQ22*(p-this.pk[jp])*(z-this.zk[jz]);

	f*=1.0/(this.pk[jp+1]-this.pk[jp])*(this.zk[jz+1]-this.zk[jz]);


	return f;

}

private double shapeFunc(int jz, int jp) {
	// TODO Auto-generated method stub
	return 0;
}

public void simulateData(){


	//ps.magnetize(10);

	int L=100;
	int M=this.nHyst;
	Vect B=new Vect(2*L);
	Vect H=new Vect(2*L);
	Vect[] pk=new Vect[M];
	Vect[] sk=new Vect[M];
	double pkp=0;
	double skp=0;
	double Bp=0;
	double Xs=1.5;
	double pkps=0;

	//	this=new double[M];
	//double[] ita=new double[M];

	for(int k=0;k<M;k++){
		pk[k]=new Vect(2*L+1);
		//zeta[k]=1*k;
	}


	for(int i=-L;i<L;i++){
		int n=i+L;

		double x=Math.exp(-.01*n)*Math.sin(2*Math.PI*n*1.0/L);
		//double x=n*.2;
		//x=this.pk[i];

		B.el[n]=x;
		//double pkp=-1000;
		for(int k=0;k<this.nHyst;k++){

			if(n>0) pkp=pk[k].el[n-1];

			pk[k].el[n]=hystron1D(x,this.zk[k],pkp);
			//	pk[k].el[n]=pk(x,this.zk[k],pkp);

			H.el[n]+=this.shapeFuncAt(k,pk[k].el[n]);
			//H.el[n]+=this.func(pk[k].el[n]);

		}
	}


	util.plot(B,H);





}

/*public void distillBHData(){

		BH=new Mat[this.nTotCurves];

		int nDB=2*(this.nSym+nMajor)+1;

		Vect divBinit=new Vect().linspace(0, Bs, (this.nSym+nMajor)+1);
		Vect divB=new Vect().linspace(Bs, -Bs,nDB);



		int L0=divB.length;


		int L=0;

		int ix=0;

		for(int i=0;i<this.nInit;i++){

			L=divBinit.length;

			BH[ix]=new Mat(L,2);

			for(int j=0;j<L;j++){

				BH[ix].el[j]=this.gethbAsc(ix,divBinit.el[j]);



			}

			ix++;
		}



		for(int i=0;i<this.nMajor;i++){

			L=divB.length;

			BH[ix]=new Mat(L,2);

			for(int j=0;j<L;j++){


				BH[ix].el[j]=this.gethbDesc(ix,divB.el[j]);
			}

			ix++;
		}


		for(int i=0;i<this.nSym;i++){

			L=divB.length-2*i;

			BH[ix]=new Mat(L,2);

			for(int j=i;j<L0-i;j++){

				BH[ix].el[j-i]=this.gethbDesc(ix,divB.el[j]);

			}

			ix++;
		}




	}*/

public double[] gethbDesc(int i,double B ){
	double[] hb=new double[2];
	int L=this.BHraw[i].nRow;



	if(B>=this.BHraw[i].el[0][1]){
		hb=BHraw[i].el[0]; 


	}

	else if(B<=this.BHraw[i].el[L-1][1]){
		hb=BHraw[i].el[L-1]; 

	}
	else{	
		int j=-1;
		while(this.BHraw[i].el[j+1][1]>B){j++;}
		double m1=(B-this.BHraw[i].el[j][1])/(this.BHraw[i].el[j+1][1]-this.BHraw[i].el[j][1]);
		double dH=this.BHraw[i].el[j+1][0]-this.BHraw[i].el[j][0];

		double H=m1*dH+this.BHraw[i].el[j][0];
		hb[0]=H; 
		hb[1]=B; 


	}

	return hb;
}

public Mat getHBij(Mat BH,int k){

	Mat BHr=new Mat(BH.nRow,2);

	if(BH.nCol==4){

		if(k==0)
			for(int i=0;i<BHr.nRow;i++){
				BHr.el[i][0]=BH.el[i][0];

				BHr.el[i][1]=BH.el[i][1];
			
			}
		else if(k==1)
			for(int i=0;i<BHr.nRow;i++){
				
				BHr.el[i][0]=BH.el[i][2];

				BHr.el[i][1]=BH.el[i][3];

			}
		else if(k==2)
			for(int i=0;i<BHr.nRow;i++){
				BHr.el[i][0]=BH.el[i][0];

				BHr.el[i][1]=BH.el[i][2];
			}
		else if(k==3)
			for(int i=0;i<BHr.nRow;i++){
				BHr.el[i][0]=BH.el[i][0];

				BHr.el[i][1]=BH.el[i][3];
			}
		else {
			throw new IllegalArgumentException("Matrix has 4 columns. No more than 4 combination possible.");
		}
	}

	else if(BH.nCol==3){

		if(k==0)
			for(int i=0;i<BHr.nRow;i++){
				BHr.el[i][0]=BH.el[i][0];

				BHr.el[i][1]=BH.el[i][1];
			}
		else if(k==1)
			for(int i=0;i<BHr.nRow;i++){
				BHr.el[i][0]=BH.el[i][0];

				BHr.el[i][1]=BH.el[i][2];
			}
		else {
			throw new IllegalArgumentException("Matrix has 3 columns. No more than 2 combination possible.");
		}
	}

	return BHr;

}

public Mat getHrBr(Mat BH){

	Mat BHr=new Mat(BH.nRow,2);

	for(int i=0;i<BHr.nRow;i++)
	{
		if(BH.nCol==4)
		{
			BHr.el[i][0]=new Vect(BH.el[i][0],BH.el[i][1]).norm();

			BHr.el[i][1]=new Vect(BH.el[i][2],BH.el[i][3]).norm();
		}
		else if(BH.nCol==3){

			BHr.el[i][0]=BH.el[i][0];

			BHr.el[i][1]=new Vect(BH.el[i][1],BH.el[i][2]).norm();
		}
	}


	return BHr;

}

public double[] gethbAsc(int i,double B ){
	double[] hb=new double[2];
	int L=this.BHraw[i].nRow;

	if(B<=this.BHraw[i].el[0][1]){
		hb=BHraw[i].el[0]; 

	}

	else if(B>=this.BHraw[i].el[L-1][1]){
		hb=BHraw[i].el[L-1]; 
	}
	else{	
		int j=-1;
		while(this.BHraw[i].el[j+1][1]<B){j++;}
		double m1=(B-this.BHraw[i].el[j][1])/(this.BHraw[i].el[j+1][1]-this.BHraw[i].el[j][1]);
		double dH=this.BHraw[i].el[j+1][0]-this.BHraw[i].el[j][0];

		double H=m1*dH+this.BHraw[i].el[j][0];
		hb[0]=H; 
		hb[1]=B; 


	}

	return hb;
}


public static double func(double x)
{
	util.pr(x);
	double y=Math.cbrt(x);

	return y;

}

public static double func2(double x)
{

	double y=Math.pow(x, 3);


	return y;


}

public static double pk(double H,double zeta,double pkp)
{

	double y=0;

	y=Math.max(Math.min(pkp,H+zeta),H-zeta);

	return y;

}

public static double hystron1D(double B,double zeta,double pp)
{

	double y=0;

	if(zeta==0) y=B;
	else
		y=B-(B-pp)/Math.max(Math.abs(B-pp)/zeta,1);

	return y;

}

public static double sk(double skp,double B,double B0,double ita)
{

	double y=0;

	y=Math.max(Math.min(B-B0+skp,ita),-ita);

	return y;


}
}
