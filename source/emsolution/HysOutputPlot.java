package emsolution;

import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

import math.Mat;
import math.Vect;
import math.util;
import fem.Model;



public class HysOutputPlot {
	
Mat BB,HH;

String regex="[:; ,\\t]+";


public static void main(String[] args){

	HysOutputPlot pg=new HysOutputPlot();
	
	//pg.loadOutputXY();
	//pg.loadOutputRT();
	
	pg.loadData();
	
}
	
public  HysOutputPlot(){


		}

public void loadData(){


	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
//	String file="C:\\Works\\HVID\\folder1\\data\\A_B窶愿ｼ窶氾坂�佚趣ｿｽﾃ姑停�ｹ�ｿｽ[ﾆ致hts_data\\hys_data";

	//String file="C:\\Users\\Hassan Ebrahimi\\JavaWorks\\MagFem\\hys_data";
	//String file="C:\\Works\\HVID\\output";
	String file="C:\\Works\\PlayModel\\output";
		   file="C:\\Works\\HVIDConv\\output";
	

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			boolean generateData=false;
			Mat[][] BHforData=null;
			if(generateData){
				BHforData=new Mat[18][17];
			}
int Nmax=1000;

int ia=-1;
Mat[] BB=new Mat[Nmax];
Mat[] HH=new Mat[Nmax];

Mat[] BHs1=new Mat[Nmax];


Mat[] XX=new Mat[2*Nmax];

int numbCurves=0;
for(int i=0;i<Nmax;i++){
	
		line=br.readLine();
		
		if(line==null) break;
		
			line=br.readLine();
			line=br.readLine();

int L=Integer.parseInt(line);
numbCurves++;



Mat bbhh=new Mat(L,4);
	
	for( int p=0;p<L;p++)
	{
		line=br.readLine();
		double[] ar=this.getCSV(line);
		bbhh.el[p]=ar;
	}
		
			
			BB[i]=new Mat(L,2);
			HH[i]=new Mat(L,2);
			
		
		BB[i].setCol(bbhh.getColVect(0).times(1), 0);
		BB[i].setCol(bbhh.getColVect(1).times(1), 1);
			
			HH[i].setCol(bbhh.getColVect(2), 0);
			HH[i].setCol(bbhh.getColVect(3), 1);
	

	
			XX[2*i]=BB[i].times(1);
			XX[2*i+1]=HH[i];
			
			BHs1[i]=new Mat(L,2);
			
		
	//	HH[0].show();
			
				Vect Hr=new Vect(HH[i].nRow);
			Vect Br=new Vect(BB[i].nRow);
			
						
		double ang=Math.atan(BB[i].el[1][1]/BB[i].el[1][0]);

		//util.pr(ang/Math.PI*180);
		//util.pr(ang/Math.PI*180);
		
		Vect er=new Vect(Math.cos(ang),Math.sin(ang));
		
		//er.hshow();
		for(int j=0;j<Hr.length;j++){
			Hr.el[j]=new Vect(HH[i].el[j][0],HH[i].el[j][1]).dot(er);
			Br.el[j]=new Vect(BB[i].el[j][0],BB[i].el[j][1]).dot(er);
			//util.pr(Hr.el[i]+"\t"+Br.el[i]);
		}
		

		BHs1[i].setCol(Hr,0);
		BHs1[i].setCol(Br,1);
		
		//BHs1[i].show();
		
		line=br.readLine();
		if(line==null) break;
		
		if(generateData){
			int n1=0;
			while( Hr.el[n1+1]<Hr.el[n1]){n1++;}

			int n2=n1;
			while( Hr.el[n2+1]>Hr.el[n2]){n2++;}
			int n3=n2;
			while(n3<Hr.length-1 && Hr.el[n3+1]<Hr.el[n3]){n3++;}
			
			int Lds=n3-n2+1;
		
	
			int j=(numbCurves-1)%BHforData[0].length;
			if( j==0) {
				ia++;
			if(ia>=BHforData.length) break;
			}
			ia=(numbCurves-1)/17;
			j=(numbCurves-1)%17;
			 BHforData[ia][j]=new Mat(Lds,2);

			 for(int k=0;k<Lds;k++){
				 BHforData[ia][j].el[k][0]=Hr.el[k+n2];
				 BHforData[ia][j].el[k][1]=Br.el[k+n2];
			 }
			
		}
		

}

if(generateData){
	file="C:\\Works\\PlayModel\\hys_data_generated";
	util.plotBunch(BHforData[14]);
//util.pr(BHforData.length);
//util.pr(BHforData[0].length);
	new PlayModel2D().writeHystData(BHforData, file,0);	
}

//util.pr(numbCurves);

//BHforData[0][0].show();
			//BH[1].show();
util.plotBunch(BHs1,numbCurves);

//util.plotBunch(BHforData[0],17);

//file="C:\\Works\\PlayModel\\hys_data_generated";

//new PlayModel2D().writeHystData(BHforData, file, 0);	

Mat BH=BHs1[0];
Mat HB=new Mat(BH.size());
HB.setCol(BH.getColVect(1), 0);
HB.setCol(BH.getColVect(0), 1);
//HB.show();
//BHs1[0].show();
			br.close();
			fr.close();
			
		util.plotBunch(HH,1);
		
	//	HH[0].show();

			int L=HH[0].nRow/2;
			Mat M=new Mat(L,2);
			for(int i=0;i<L;i++){
				M.el[i]=HH[0].el[i+L];
			}
			
			//M.show();	
			int nangs=9;
			Mat BHall=new Mat(BHs1[0].nRow,nangs+1);
			for(int i=0;i<-BHs1[0].nRow;i++){
				BHall.el[i][0]=BHs1[0].el[i][1];
				for(int j=1;j<nangs+1;j++)
					BHall.el[i][j]=BHs1[j-1].el[i][0];
			
			}
			//BHall.show();
/*			for(int i=0;i<numbCurves;i++)
				BHs[i]=BHs1[i];*/
			//M.show();
		//	util.plot(M);
			//util.plotBunch(BHs);
		//	BHs[0].show();
			String file1="C:\\Works\\HVID\\b_times";

				
			}
		
			catch(IOException e){System.err.println("Error in loading BH data file.");

			}

		}

public void loadOutputRT(){

	

	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
//	String file="C:\\Works\\HVID\\folder1\\data\\A_B窶愿ｼ窶氾坂�佚趣ｿｽﾃ姑停�ｹ�ｿｽ[ﾆ致hts_data\\hys_data";

	//String file="C:\\Users\\Hassan Ebrahimi\\JavaWorks\\MagFem\\hys_data";
	//String file="C:\\Works\\HVID\\output";
	String file="C:\\Works\\PlayModel\\output";
	//String file="C:\\Works\\HVIDConv\\output";

	int Nmax=20;

	Mat[] BB=new Mat[Nmax];
	Mat[] HH=new Mat[Nmax];


		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;


int numbCurves=0;
for(int i=0;i<Nmax;i++){
	
		line=br.readLine();
		
		if(line==null) break;
		
			line=br.readLine();
			line=br.readLine();

int L=Integer.parseInt(line);

numbCurves++;

Mat bbhh=new Mat(L,4);
	
	for( int p=0;p<L;p++)
	{
		line=br.readLine();
		double[] ar=this.getCSV(line);
		bbhh.el[p]=ar;
	}
		
			
			BB[i]=new Mat(L,2);
			HH[i]=new Mat(L,2);
			
		
		BB[i].setCol(bbhh.getColVect(0).times(1), 0);
		BB[i].setCol(bbhh.getColVect(1).times(1), 1);
			
			HH[i].setCol(bbhh.getColVect(2), 0);
			HH[i].setCol(bbhh.getColVect(3), 1);
	
		line=br.readLine();
		

}
		}
		catch(IOException e){System.err.println("Error in loading BH data file.");

		}


 file="C:\\Works\\PlayModel\\outputn";
//String file="C:\\Works\\HVIDConv\\output";



Mat[] BBn=new Mat[Nmax];
Mat[] HHn=new Mat[Nmax];


	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;


int numbCurves=0;
for(int i=0;i<Nmax;i++){

	line=br.readLine();
	
	if(line==null) break;
	
		line=br.readLine();
		line=br.readLine();

int L=Integer.parseInt(line);

numbCurves++;

Mat bbhh=new Mat(L,4);

for( int p=0;p<L;p++)
{
	line=br.readLine();
	double[] ar=this.getCSV(line);
	bbhh.el[p]=ar;
}
	
		
		BBn[i]=new Mat(L,2);
		HHn[i]=new Mat(L,2);
		
	
	BBn[i].setCol(bbhh.getColVect(0).times(1), 0);
	BBn[i].setCol(bbhh.getColVect(1).times(1), 1);
		
		HHn[i].setCol(bbhh.getColVect(2), 0);
		HHn[i].setCol(bbhh.getColVect(3), 1);

	line=br.readLine();
	

}
	}	catch(IOException e){System.err.println("Error in loading BH data file.");

	}

	
	double ang=util.getAng(new Vect(BB[0].el[1]));
	
	for(int i=0;i<HH[0].nRow;i++){
		HH[0].el[i][0]+=-HHn[0].el[i][1];
		HH[0].el[i][1]+=HHn[0].el[i][0];
	}
	
	util.plot(HH[0]);
	
	HH[0].show();

		}

public void loadOutputXY(){

	

	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
//	String file="C:\\Works\\HVID\\folder1\\data\\A_B窶愿ｼ窶氾坂�佚趣ｿｽﾃ姑停�ｹ�ｿｽ[ﾆ致hts_data\\hys_data";

	//String file="C:\\Users\\Hassan Ebrahimi\\JavaWorks\\MagFem\\hys_data";
	//String file="C:\\Works\\HVID\\output";
	String file="C:\\Works\\PlayModel\\outputx";
	//String file="C:\\Works\\HVIDConv\\output";

	int Nmax=20;

	Mat[] BB=new Mat[Nmax];
	Mat[] HH=new Mat[Nmax];


		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;


int numbCurves=0;
for(int i=0;i<Nmax;i++){
	
		line=br.readLine();
		
		if(line==null) break;
		
			line=br.readLine();
			line=br.readLine();

int L=Integer.parseInt(line);

numbCurves++;

Mat bbhh=new Mat(L,4);
	
	for( int p=0;p<L;p++)
	{
		line=br.readLine();
		double[] ar=this.getCSV(line);
		bbhh.el[p]=ar;
	}
		
			
			BB[i]=new Mat(L,2);
			HH[i]=new Mat(L,2);
			
		
		BB[i].setCol(bbhh.getColVect(0).times(1), 0);
		BB[i].setCol(bbhh.getColVect(1).times(1), 1);
			
			HH[i].setCol(bbhh.getColVect(2), 0);
			HH[i].setCol(bbhh.getColVect(3), 1);
	
		line=br.readLine();
		

}
		}
		catch(IOException e){System.err.println("Error in loading BH data file.");

		}


 file="C:\\Works\\PlayModel\\outputy";
//String file="C:\\Works\\HVIDConv\\output";



Mat[] BBn=new Mat[Nmax];
Mat[] HHn=new Mat[Nmax];


	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;


int numbCurves=0;
for(int i=0;i<Nmax;i++){

	line=br.readLine();
	
	if(line==null) break;
	
		line=br.readLine();
		line=br.readLine();

int L=Integer.parseInt(line);

numbCurves++;

Mat bbhh=new Mat(L,4);

for( int p=0;p<L;p++)
{
	line=br.readLine();
	double[] ar=this.getCSV(line);
	bbhh.el[p]=ar;
}
	
		
		BBn[i]=new Mat(L,2);
		HHn[i]=new Mat(L,2);
		
	
	BBn[i].setCol(bbhh.getColVect(0).times(1), 0);
	BBn[i].setCol(bbhh.getColVect(1).times(1), 1);
		
		HHn[i].setCol(bbhh.getColVect(2), 0);
		HHn[i].setCol(bbhh.getColVect(3), 1);

	line=br.readLine();
	

}
	}	catch(IOException e){System.err.println("Error in loading BH data file.");

	}

	
	double ang=util.getAng(new Vect(BB[0].el[1]));
	
	for(int i=0;i<HH[0].nRow;i++){
		double Hx1=HH[0].el[i][0]*Math.cos(ang)+HH[0].el[i][1]*Math.sin(ang);
		double Hy1=-HH[0].el[i][0]*Math.sin(ang)+HH[0].el[i][1]*Math.cos(ang);
		
		double Hx2=HHn[0].el[i][0]*Math.cos(ang)-HHn[0].el[i][1]*Math.sin(ang);
		double Hy2=HHn[0].el[i][0]*Math.sin(ang)+HHn[0].el[i][1]*Math.cos(ang);
		HH[0].el[i][0]=Hx1+Hx2;
		HH[0].el[i][1]=Hy1+Hy2;
	
	}
	
	util.plot(HH[0]);
	
	HH[0].show();

		}



		private double[] getCSV(String line){
			
			String[] sp=line.split(regex);	
		
			int p0=0;
			if(sp[0].equals(""))
			{
				p0=1;
			}
			int L=sp.length-p0;

			double[] v=new double[L];

			for( int p=0;p<L;p++){

				v[p]=Double.parseDouble(sp[p+p0]);
			}

			return v;
		}
	
	
}
