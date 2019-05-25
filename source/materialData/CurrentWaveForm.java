package materialData;

import java.io.FileReader;
import java.util.Scanner;

import io.Writer;
import math.Mat;
import math.Vect;
import math.util;



public class CurrentWaveForm {
	
		public double[][] TI;
		public int length;
		boolean periodic=true;

	
	public CurrentWaveForm(String f) throws Exception{
		double[][] TI1=new double[20000][2];
		
		String file = System.getProperty("user.dir") + "\\"+f;

		int i;
		// try{
		 Scanner scr=new Scanner(new FileReader(file));
		 while(scr.hasNext()){
	
		 while(!scr.next().equals("begin")){}
		 int j=0;
		 String s=scr.next();
		 while(!s.equals("end")){
		 TI1[j][0]=Double.parseDouble(s);
		 s=scr.next();
		 TI1[j++][1]=Double.parseDouble(s);
		 s=scr.next();
		 }
		this.length=j;
		 }
			
		 this.TI=new double[this.length][2];
		 for(int k=0;k<this.length;k++)
		 this.TI[k]=TI1[k];
		 
		 scr.close();

		}
	
	
	public CurrentWaveForm(double[][] ti) {
		TI=new double[ti.length][2];
		this.length=ti.length;

		 for(int k=0;k<this.length;k++){
		 this.TI[k][0]=ti[k][0];
		 this.TI[k][1]=ti[k][1];
		 }
		 

		}
	
	public CurrentWaveForm(Vect arg, Vect value) {
		TI=new double[arg.length][2];
		this.length=arg.length;

		 for(int k=0;k<this.length;k++){
		 this.TI[k][0]=arg.el[k];
		 this.TI[k][1]=value.el[k];
		 }
		 

		}
	
	public void setPeriodic(boolean b) {
		this.periodic=b;
	}
	
	public static void main(String[] args) throws Exception{
		
/*		BHCurve BH=new BHCurve("H350");
		
		Curve cv=new Curve(BH,600,600);
		cv.show(true);*/

CurrentWaveForm Ia=new CurrentWaveForm("emf//IaMath.txt");

CurrentWaveForm Ib=new CurrentWaveForm("emf//Ibn.txt");

CurrentWaveForm Ic=new CurrentWaveForm("emf//Icn.txt");
/*Mat M=new Mat(3,3);
double f=100;

M.el[0][0]=1;
M.el[0][1]=-0.5;
M.el[0][2]=0.5;
M.el[1][0]=0;
M.el[1][1]=0.5*Math.sqrt(3);
M.el[1][2]=-0.5*Math.sqrt(3);
M.el[2][0]=.5;
M.el[2][1]=.5;
M.el[2][2]=.5;

M=M.times(2.0/3);

Mat[] Io=new Mat[3];
for(int i=0;i<3;i++){
	Io[i]=new Mat(Ia.length,2);
	for(int k=0;k<Ia.length;k++)
		Io[i].el[k][0]=Ia.TI[k][0];
}

for(int k=0;k<Ia.length;k++){
	double t=Ia.TI[k][0];
	

	double tm=1*2*Math.PI*f*t;

	
	M.el[0][0]=Math.cos(tm);
	M.el[0][1]=Math.cos(tm-2*Math.PI/3);
	M.el[0][2]=Math.cos(tm+2*Math.PI/3);
	M.el[1][0]=-Math.sin(tm);
	M.el[1][1]=-Math.sin(tm-2*Math.PI/3);
	M.el[1][2]=-Math.sin(tm+2*Math.PI/3);
	M.el[2][0]=0.5*Math.sqrt(2);
	M.el[2][1]=M.el[2][0];
	M.el[2][2]=M.el[2][0];
	
	M=M.times(Math.sqrt(2.0/3));

	
	double wt=2*Math.PI*f*t;
	
	double iat=Ia.TI[k][1]*0+Math.cos(wt);
	double ibt=Ib.TI[k][1]*0+Math.cos(wt-2*Math.PI/3);
	double ict=Ic.TI[k][1]*0+Math.cos(wt+2*Math.PI/3);
	
//t+=2.5/360*.01;
	double iat=Ia.getI(t);
	double ibt=Ib.getI(t);
	double ict=Ic.getI(t);
	
	Vect ii=new Vect(iat,ibt,ict);
	
	Vect iout=M.mul(ii);
	
	iout.el[0]=iat*Math.cos(Math.PI/6)+ibt*Math.cos(Math.PI/6+40*Math.PI/180)+ict*Math.cos(Math.PI/6-40*Math.PI/180);
	iout.el[1]=iat*Math.sin(Math.PI/6)+ibt*Math.sin(Math.PI/6+40*Math.PI/180)+ict*Math.sin(Math.PI/6-40*Math.PI/180);
	for(int i=0;i<3;i++)
	Io[i].el[k][1]=iout.el[i];
	

	Io[0].el[k][1]=iat;
	Io[1].el[k][1]=ibt;
	Io[2].el[k][1]=ict;
	

}
util.plotBunch(Io);
Vect v=new Vect(Io[0].nRow);

Mat arr=new Mat(Io[0].nRow,4);
for(int k=0;k<Ia.length;k++){
	arr.el[k][0]=Io[0].el[k][0];
	
	arr.el[k][1]=Io[0].el[k][1];
	arr.el[k][2]=Io[1].el[k][1];
	arr.el[k][3]=Io[2].el[k][1];
	
	v.el[k]=util.getAng(new Vect(arr.el[k][1],arr.el[k][2]));
}


util.plot(v);
String file2 = "C:\\Works\\MeshView\\Iclark.txt";

//util.plot(arr.el);
	
new Writer().writeArray(arr.el, file2);
*/

/*CurrentWaveForm Iy=new CurrentWaveForm("vb0.txt");
CurrentWaveForm Iz=new CurrentWaveForm("vc0.txt");
		*/
/*		Curve cv2=new Curve(Ix,600,600);
		cv2.show(true);
		//util.show(Ix.TI);
		for(int k=0;k<Ix.TI.length;k++)
			util.pr(Ix.TI[k][1]);

		double[][] iy=new double[Ix.length][2];
		
		Vect  dv=new Vect(1800);

		for(int k=0;k<0*dv.length;k++){
		
			double t=.0025*0+k*5.55555e-6;
			//iy[k][1]=Ix.getI(iy[k][0])+Iy.getI(iy[k][0])+Iz.getI(iy[k][0]);
			dv.el[k]=Ix.getI(t);
			util.pr(k*5.55555e-6+"\t"+dv.el[k]);
			//util.pr(Ix.TI[k+1800][1]);
		}*/
		
		//dv.show();
	//util.plot(dv);
		
		//util.plot(dv);


		/*String bh = System.getProperty("user.dir") + "\\bh.xls";
		BH2.writexls(bh);
	*/
	}
	
	public  double getI(double t){


		if(t<this.TI[0][0]|| t>this.TI[length-1][0]){
			if(!this.periodic) throw new IllegalArgumentException("out of time span.");
			else{
				 t=t-Math.floor(t/this.TI[length-1][0])*this.TI[length-1][0];
			}
		}

		int j=getj(t);
		 double i=(this.TI[j][1]*(this.TI[j+1][0]-t)+this.TI[j+1][1]*(t-this.TI[j][0]))/(this.TI[j+1][0]-this.TI[j][0]);
	
		 return i;
	}
	
	public int getj(double t){
		
		int j=0;
		if(this.TI.length==1) return j;
		while(this.TI[j+1][0]<t){j++;}
		return j;
	}
	

}
