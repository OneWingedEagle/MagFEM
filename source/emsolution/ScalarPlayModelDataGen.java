package emsolution;

import io.Loader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Random;

import static java.lang.Math.*;
import math.Mat;
import math.Vect;
import math.util;

public class ScalarPlayModelDataGen {

	public double[] zk;


	public ScalarPlayModelDataGen()
	{	}

	public static void main(String[] args)
	{

		HysDataGraph pg=new HysDataGraph();
		
		pg.loadHysData();

		double Bm=1;
		Mat Hanixy=new Mat(pg.nAni[0],2);
		for(int i=0;i<pg.nAni[0];i++){
			double ang=i*PI/18;
			int j=0;
			while(pg.BHAni[0][i].el[j+1][0]<Bm+.01){j++;}
			
			Hanixy.el[i][0]=pg.BHAni[i][0].el[j][1];
			Hanixy.el[i][1]=pg.BHAni[i][0].el[j][2];
			//Vect er=new Vect(cos(ang),sin(ang));
		}
		//Hanixy.show();
		
	//	util.plot(Hanixy);
		
		ScalarPlayModelDataGen pm=new ScalarPlayModelDataGen();
		

		String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
	//	String file="C:\\Works\\HVID\\hys_dataH";
		
		//String file2="C:\\Works\\PlayModel\\hysRotation";
		String file2="C:\\Works\\PlayModel\\hysRotNew";
		Loader ld=new Loader();

		
		int nloop=15;
		
		Mat BHij=new Mat(ld.loadArrays(1024,4*nloop+1,file2));
		//BHij.show();
		
		Mat[] BB=new Mat[nloop];
		Mat[] HH=new Mat[nloop];
		
		for(int i=0;i<BB.length;i++){
			BB[i]=new Mat(BHij.nRow,2);
			HH[i]=new Mat(BHij.nRow,2);
			
/*			BB[i].setCol(BHij.getColVect(4*i+1), 0);
			BB[i].setCol(BHij.getColVect(4*i+3), 1);
			
			HH[i].setCol(BHij.getColVect(4*i+2), 0);
			HH[i].setCol(BHij.getColVect(4*i+4), 1);*/
			
			BB[i].setCol(BHij.getColVect(4*i), 0);
			BB[i].setCol(BHij.getColVect(4*i+1), 1);
			
			HH[i].setCol(BHij.getColVect(4*i+2), 0);
			HH[i].setCol(BHij.getColVect(4*i+3), 1);
		}
	
		//HH[8].show();
	
		util.plot(HH[12]);
		
		HysOutputPlot hop=new HysOutputPlot();
		hop.loadData();
		Mat[] M2=new Mat[3];
		//M2[0]=hop.HH[0];
		

		M2[1]=HH[9];
		
		
		M2[2]=BB[9];
		Vect v=new Vect(M2[1].nRow);
		for(int i=0;i<M2[1].nRow;i++){
			double ang1=util.getAng(new Vect(M2[1].el[i]));
			double ang2=util.getAng(new Vect(M2[2].el[i]));
			v.el[i]=acos(cos(ang1-ang2))*180/PI;
		}

		util.plot(v);
	
		//M2[2]=Hanixy;
		
		int L=230;
/*		for(int i=0;i<18;i++){
			HH[9].el[i*30+L][0]*=2;
			HH[9].el[i*30+L][1]*=2;
		}*/
		
		//util.plot(HH[9]);

		
		Mat[] BHani2x=new Mat[17];
		Mat[] BHani2y=new Mat[17];
		for(int j=0;j<BHani2x.length;j++)
		{
			BHani2x[j]=new Mat(18,2);
			BHani2y[j]=new Mat(18,2);
		}

		//util.pr(BB[0].nRow);
				
		for(int i=0;i<17;i++)
			for(int j=0;j<17;j++)
			{
				BHani2x[i].el[j+1][0]=new Vect(BB[j].el[i*30+L]).norm();
				BHani2x[i].el[j+1][1]=HH[j].el[i*30+L][0];
				
				BHani2y[i].el[j+1][0]=BHani2x[i].el[j][0];
				BHani2y[i].el[j+1][1]=HH[j].el[i*30+L][1];
			}
		
	//	BHani2x[0].show();
		Mat Ba1=new Mat(18,19);
		Mat Ba2=new Mat(18,19);
		
		Ba1.setCol(BHani2x[0].getColVect(0), 0);
		Ba2.setCol(BHani2x[0].getColVect(0), 0);
		for(int i=1;i<18;i++){
			Ba1.setCol(BHani2x[i-1].getColVect(1), i);
			
			Ba2.setCol(BHani2y[i-1].getColVect(1), i);
		}
		
		
	Ba1.setCol(BHani2x[0].getColVect(1).times(-1), 18);	
	Ba2.setCol(BHani2y[0].getColVect(1).times(-1), 18);	
	
	//Ba2.show();
	
	util.plotBunch(M2,2);
	
/*Mat BHani2x=new Mat(18,18);
Mat BHani2y=new Mat(18,18);
for(int j=0;j<17;j++){
	BHani2x.el[j+1][0]=new Vect(BB[j].el[0]).norm();
	BHani2y.el[j+1][0]=	BHani2x.el[j+1][0];
}

		
for(int i=0;i<17;i++)
	for(int j=0;j<17;j++)
	{
		BHani2x.el[i+1][j+1]=HH[i].el[j*10][0];
		BHani2y.el[i+1][j+1]=HH[i].el[j*10][1];
	}
	*/
//BHani2x.show();
	
		//pm.createData(file);


		//pm.loadData(file);

		//pm.simulateData();
		/*		for(int i=0;i<pm.BH.length;i++)
			pm.BH[i].transp().show();*/
		//	pm.createData();
		/*	
		Mat[] BH1=new Mat[2*pm.BH.length];

		for(int i=0;i<pm.BH.length;i++){
			BH1[i]=pm.BHraw[i];
			BH1[i+pm.BH.length]=pm.BH[i];

		}

		util.plotBunch(BH1);*/
	}




	public void createData(String file){


		double Bs=1.8;
		double Hs=1600;
		
		
		int L=1000;

		

		Vect B=new Vect(L);
		Vect H=new Vect(L);
		
		PlayModelShapeGraph pgs=new PlayModelShapeGraph();
		pgs.loadShapeFunc();
		
		Mat[] BHs=pgs.BH;
		

		Mat[] M2=new Mat[BHs.length];
		for(int i=0;i<BHs.length;i++){
			M2[i]=new Mat(BHs[0].nRow,2);
			M2[i].setCol(BHs[0].getColVect(0), 0);
			for(int j=0;j<M2[i].nRow;j++){
				double y=this.getFuncOfB(BHs[i], M2[i].el[j][0]);
				M2[i].el[j][1]=y;
			}
		}
		

		//util.plotBunch(M2);
		
		
		Vect zkmp=M2[0].getColVect(0);
		
		int im=0;
		
		while(zkmp.el[im]<0){im++;}
		
		int K=zkmp.length-im;
		
		Vect zk=new Vect(K);
		for(int i=0;i<K;i++)
			zk.el[i]=BHs[0].el[i+im][0];
		


		zk.show();
			
		L=zk.length;
	//	BHs[2].show();
		
	/*	for(int i=0;i<BHs.length;i++){
			util.pr(BHs[i].nRow);
		}*/
		
		for(int i=0;i<L;i++){
			B.el[i]=1.5*sin(i*4*Math.PI/L);
			
			for(int k=0;k<1;k++){
		
			if(i>0){
				//H.el[i]+=getFuncOfB(M2[k],this.hysteron1D(B.el[i], zk.el[k], H.el[i-1]));
				H.el[i]=M2[k].el[i][1];//,this.hysteron1D(B.el[i], zk.el[k], H.el[i-1]));
		
			}
		}
		}

		util.plot(H,B);



		
		boolean write=false;

if(write=true){



	DecimalFormat dfB=new DecimalFormat("#.00");
	DecimalFormat dfH=new DecimalFormat("#.0");


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		

			pwBun.println(1+"\t"+1);
			pwBun.println("*Bs*Hs*");
			pwBun.println(Bs+"\t"+Hs);

			pwBun.println("* 初磁化曲線数 * メジャーループ数 * 対称ループ数 * 下降曲線数 * 上昇曲線数 *");
/*			pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);

			for(int i=0;i<nTot;i++){
				pwBun.println("*xxx");
				pwBun.println(BHs[i].nRow);
				for(int j=0;j<BHs[i].nRow;j++)
					pwBun.println(BHs[i].el[j][0]+"\t"+BHs[i].el[j][1]);
			}*/

			pwBun.println("* ----- 回転ヒステリシス損");
			pwBun.println("* B数 *");
			pwBun.println("0");
			pwBun.println("* B * 損失");
			pwBun.println("* ----- 異方性");
			pwBun.println("* B数 * 角度数 *");
			pwBun.println(0+"\t"+0);
			pwBun.println("* B * H ･････ *　磁化容易軸");
			
				

			util.pr("Simulated hysteresis data was written to "+file+".");

			pwBun.close();
		}
		catch(IOException e){}


	}



	}
	


	public  double func(double x)
	{
		double y=Math.cbrt(x);

		return y;

	}

	public  double funcOfB(double zk,double p)
	{

	
	
		double y=zk*Math.pow(p,2)*p;


		return y;


	}
	
	public  double getFuncOfB(Mat BH, double p)
	{


		int ix=0;
		int iy=1;
		
		int i1=0;
		int i2=BH.nRow-1;
		


		if(p<=BH.el[i1][ix])
			return BH.el[i1][iy];
		if(p>=BH.el[i2][ix]){
			return BH.el[i2][iy];
	
		}


		int j=0;


			while(BH.el[j+1][ix]<p){j++;}

		double cc=(BH.el[j+1][iy]-BH.el[j][iy])/(BH.el[j+1][ix]-BH.el[j][ix]);

		double H=BH.el[j][iy]+(p-BH.el[j][ix])*cc;

		
		return H;
				
	}

	public  double pk(double H,double zeta,double pkp)
	{

		double y=0;

		y=Math.max(Math.min(pkp,H+zeta),H-zeta);

		return y;

	}



	public  double sk(double skp,double B,double B0,double ita)
	{

		double y=0;

		y=Math.max(Math.min(B-B0+skp,ita),-ita);

		return y;


	}
	
	public  double hysteron1D(double B,double zeta,double pp)
	{

		double y=0;

		if(zeta==0) y=B;
		else
			y=B-(B-pp)/Math.max(Math.abs(B-pp)/zeta,1);

		return y;

	}
}
