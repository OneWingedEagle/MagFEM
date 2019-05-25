package materialData;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Scanner;

import math.Vect;
import math.util;



public class BHCurve {
	
		public double[][] BH;
		private double[][] gradBH;
		public int length;
		private double eps=4e-6;
		private double mu0=4*Math.PI*1e-7;

	
	public BHCurve(String file) throws Exception{
		double[][] BH1=new double[20000][2];
		int i;
		 try{
		 Scanner scr=new Scanner(new FileReader(file));
		 while(scr.hasNext()){
	
		 while(!scr.next().startsWith("begin")){}
		 int j=0;
		 String s=scr.next();
		 while(!s.contains("end")){
		 BH1[j][0]=Double.parseDouble(s);
		 s=scr.next();
		 BH1[j++][1]=Double.parseDouble(s);
		 s=scr.next();
		 util.pr(BH1[j-1][0]+"\t "+BH1[j-1][1]);

		 }
		this.length=j;
		 }
			
		 this.BH=new double[this.length][2];
		 for(int k=0;k<this.length;k++)
		 this.BH[k]=BH1[k];
		 
		 scr.close();
		 }
		 catch(IOException fnf){
			 throw new Exception(fnf);
		 }
		 setGradBH();
			
			

			}
	
	public BHCurve(double[][] BH){
		this.length=BH.length;
		this.BH=new double[this.length][2];
		for(int i=0;i<this.length;i++)
			for(int j=0;j<2;j++)
			this.BH[i][j]=BH[i][j];
		setGradBH();
	}
			
			public void setGradBH(){
				 this.gradBH=new double[this.BH.length][2];
				 for(int k=0;k<this.gradBH.length-1;k++){
					 this.gradBH[k][0]=(this.BH[k+1][1]-this.BH[k][1]);
					 this.gradBH[k][1]=(this.BH[k+1][0]-this.BH[k][0])/this.gradBH[k][0];
				 }
				 this.gradBH[this.gradBH.length-1][0]=1E6;
				 this.gradBH[this.gradBH.length-1][1]=this.gradBH[this.gradBH.length-2][1];
				
			}
			
		
			
	
	public static void main2(String[] args) throws Exception{
		
/*		BHCurve BH=new BHCurve("H350");
		
		Curve cv=new Curve(BH,600,600);
		cv.show(true);*/

BHCurve BH2=new BHCurve("ttt");
		
		Curve cv2=new Curve(BH2,600,600);
		cv2.show(true);
		Vect v=new Vect(1030);
		
		for(int i=0;i<v.length;i++)
			util.pr(BH2.getH(i*1.8/1030));
	
		/*String bh = System.getProperty("user.dir") + "\\bh.xls";
		BH2.writexls(bh);
	*/
	}
	
	public  double getH(double B){

		if(B>=this.BH[this.length-1][1])
			//return this.BH[this.length-1][0]+(B-this.BH[this.length-1][1])*this.gradBH[this.length-1][1];
			return this.BH[this.length-1][0]+(B-this.BH[this.length-1][1])/mu0;
		
		int j=0;

		while(this.BH[j+1][1]<B){j++;}
		return (this.BH[j][0]*(this.BH[j+1][1]-B)+this.BH[j+1][0]*(B-this.BH[j][1]))/this.gradBH[j][0];
				
	}
	
	public  double getB(double H){
		
		if(H>=this.BH[this.length-1][0])
			return this.BH[this.length-1][1]+(H-this.BH[this.length-1][0])/this.gradBH[this.length-1][0];
		
		int j=0;
		while(this.BH[j+1][0]<H){j++;}
		return (this.BH[j][1]*(this.BH[j+1][0]-H)+this.BH[j+1][1]*(H-this.BH[j][0]))/(this.BH[j+1][0]-this.BH[j][0]);
				
	}
	
	public  double getdHdB(double B){

		if(B>=this.BH[this.length-1][1])
			return this.gradBH[this.length-1][1];

		int j=0;
		while(this.BH[j+1][1]<B){j++;}
		if(j>=this.length-1)
			return this.gradBH[this.length-2][1];

		double interp=((B-this.BH[j][1])*this.gradBH[j+1][1]+(this.BH[j+1][1]-B)*this.gradBH[j][1])/this.gradBH[j][0];
		return interp;
				
	}
	
	public  double getNu(double B){
	
		if(B==0) B=this.eps;
		
		return getH(B)/B;
	
				
	}
	


	
	
	public  double getNuVar(double B){
		if(B==0) B=this.eps;
		return (getdHdB(B)-getNu(B))/(B*B);
	
				
	}
	


	public  double[][] getBH(){
		double[][] BH1=new double[BH.length][BH[0].length];
		for(int i=0;i<this.length;i++)
		for(int j=0;j<2;j++)
			BH1[i][j]=this.BH[i][j];
		return BH1;
	
				
	}
	
	public  double getPemB(double B1,double B2){
		double pem=0;
		boolean partial=false;
		int i;
		for(i=0; i<length-1;i++){
			if(BH[i][1]<B1) continue;
			if(BH[i+1][1]>B2) { partial=true; break;}
		
			pem+=(BH[i+1][0]+BH[i][0])*(BH[i+1][1]-BH[i][1])/2;
			
		}

		if(partial){
		pem+=(getH(B2)+BH[i][0])*(B2-BH[i][1])/2;
		}
		else if(i==length-1)
		{
		pem+=(getH(B2)+BH[length-1][0])*(B2-BH[length-1][1])/2;
		}
	
		return pem;

}
	
	public  double getPem(double H1,double H2){
		double pem=0;
		boolean partial=false;
		int i;
		for(i=0; i<length-1;i++){
			if(BH[i][0]<H1) continue;
			if(BH[i+1][0]>H2) { partial=true; break;}
		
			pem+=(BH[i+1][1]+BH[i][1])*(BH[i+1][0]-BH[i][0])/2;
			
		}

		if(partial){
		pem+=(getB(H2)+BH[i][1])*(H2-BH[i][0])/2;
		}
		else if(i==length-1)
		{
		pem+=(getB(H2)+BH[length-1][1])*(H2-BH[length-1][0])/2;
		}
	
		return pem;

}
	
	public  double getPem(double H){
		
	
		return getPem(0,H);

}
	
/*	public  double getPem(double H1,double H2){
		double pem=0;
		
		int i;
		for(i=0; i<length-1;i++){
			if(BH[i][0]<H1) continue;
			if(BH[i+1][0]>H2) break;
		
			pem+=(BH[i+1][1]+BH[i][1])*(BH[i+1][0]-BH[i][0])/2;
			
		}
		if(H2>BH[i][0])
		pem+=(getB(H2)+BH[i][1])*(H2-BH[i][0])/2;
	
		return pem;

}*/
	
	public int geti(double B){
		
		int i=0;
		if(this.BH.length==1) return i;
		while(this.BH[1][i+1]<B){i++;}
		return i;
	}
	
	public void showCurve(){
		Curve cv=new Curve(this,800,600);
		cv.show(true);

	}

	
}
