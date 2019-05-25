package emsolution;

import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Scanner;

import math.Mat;
import math.Vect;
import math.util;
import fem.Model;



public class HysEnergyPlot {
	
Mat BE[];

String regex="[:; ,\\t]+";


public static void main(String[] args){

	HysEnergyPlot pg=new HysEnergyPlot();
	
	pg.loadData();
	
}
	
public  HysEnergyPlot(){}

public void loadData(){


	String file="C:\\Works\\PlayModel\\energy";
	
	

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

int Nmax=20;

BE=new Mat[Nmax];


int numbCurves=0;
for(int i=0;i<Nmax;i++){
	
		line=br.readLine();
		
		if(line==null) break;
		
		sp=line.split(regex);

int nRow=Integer.parseInt(sp[0]);
int nCol=Integer.parseInt(sp[1]);

numbCurves++;

Mat be=new Mat(nRow,nCol);
	
	for( int p=0;p<nRow;p++)
	{
		line=br.readLine();
		double[] ar=this.getCSV(line);
		be.el[p]=ar;
	}
		

			
			BE[i]=new Mat(nRow,2);
			
		
		BE[i].setCol(be.getColVect(0).times(1), 0);
		BE[i].setCol(be.getColVect(nCol-1).times(1), 1);


		
		line=br.readLine();
}

			br.close();
			fr.close();
	
			util.plotBunch(BE,1);
			
			util.plot(BE[0].getColVect(1));

		}
		
			catch(IOException e){System.err.println("Error in loading BH data file.");

			}
		
		

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

			pwBun.println("* 初磁化曲線数 * メジャーループ数 * 対称ループ数 * 下降曲線数 * 上昇曲線数 *");
			pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);

			for(int i=0;i<nTot;i++){
				pwBun.println("*xxx");
				pwBun.println(BHs[i].nRow);
				for(int j=0;j<BHs[i].nRow;j++)
					pwBun.println(BHs[i].el[j][0]+"\t"+BHs[i].el[j][1]);
			}

			pwBun.println("* ----- 回転ヒステリシス損");
			pwBun.println("* B数 *");
			pwBun.println("0");
			pwBun.println("* B * 損失");
			pwBun.println("* ----- 異方性");
			pwBun.println("* B数 * 角度数 *");
			pwBun.println(Lani+"\t"+nAni); 		//	pwBun.println(Lani+"\t"+nAni);
			pwBun.println("* B * H ･････ *　磁化容易軸");
			
			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][0])+"\t");
				}
				pwBun.println();
			}
			pwBun.println("* B * H ･････ *　磁化困難軸");
			for(int i=0;i<Lani;i++){
				pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
				}
				pwBun.println();
			}
				

			util.pr("Simulated angle-dependent hysteresis data was written to "+file+".");

			pwBun.close();
		}
		catch(IOException e){}
		

}

public void writeHystData(Mat[][] BHs, String file){
	
	int nSet=BHs.length;
	
	double Bseff= BHs[0][0].el[BHs[0][0].nRow-1][1];
	
	Mat[] BHani=new Mat[1];
	int nAni=0;
	int Lani=0;
	int nInit=1;
	int nMajor=1;
	int nTot=BHs[0].length;
	int nSymLoops=nTot-2;
	int nDescending=0;
	int nAscending=0;
	
	double Hseff=BHs[0][0].el[BHs[0][0].nRow-1][0];

	DecimalFormat dfB=new DecimalFormat("#.00");
	DecimalFormat dfH=new DecimalFormat("#.0");


		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		

for(int ia=0;ia<nSet;ia++){

			pwBun.println(1+"\t"+1+"\t"+nSet+"\t"+ia);
			pwBun.println("*Bs*Hs*");
			pwBun.println(Bseff+"\t"+Hseff);

			pwBun.println("* 初磁化曲線数 * メジャーループ数 * 対称ループ数 * 下降曲線数 * 上昇曲線数 *");
			pwBun.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);

			for(int i=0;i<nTot;i++){
				pwBun.println("*xxx");
				pwBun.println(BHs[ia][i].nRow);
				for(int j=0;j<BHs[ia][i].nRow;j++)
					pwBun.println(BHs[ia][i].el[j][0]+"\t"+BHs[ia][i].el[j][1]);
			}

			pwBun.println("* ----- 回転ヒステリシス損");
			pwBun.println("* B数 *");
			pwBun.println("0");
			pwBun.println("* B * 損失");
			pwBun.println("* ----- 異方性");
			pwBun.println("* B数 * 角度数 *");
			pwBun.println(0+"\t"+0); 		//	pwBun.println(Lani+"\t"+nAni);
			pwBun.println("* B * H ･････ *　磁化容易軸");
			
			for(int i=0;i<0*Lani;i++){
				pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][0])+"\t");
				}
				pwBun.println();
			}
			pwBun.println("* B * H ･････ *　磁化困難軸");
			for(int i=0;i<0*Lani;i++){
				pwBun.print(dfB.format(BHani[0].el[i][0])+"\t");
				for(int j=0;j<nAni;j++){
					pwBun.print(dfH.format(BHani[j].el[i][1])+"\t");
				}
				pwBun.println();
			}
/*			pwBun.println();
			pwBun.println("End of hysteresis data set "+ia);
			pwBun.println();*/
}
				

			util.pr("Simulated angle-dependent hysteresis data was written to "+file+".");

			pwBun.close();
		}
		catch(IOException e){}
		

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
