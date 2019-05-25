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
import math.util;
import fem.Model;



public class PlayModelShapeGraph {
	
Mat[] BH;
int numb;
double Bs, Hs;
int numb2,numb3,numb4,numb5;
String regex="[:; ,\\t]+";


public static void main(String[] args){

	PlayModelShapeGraph pg=new PlayModelShapeGraph();
	
	pg.loadShapeFunc();
}
	
public boolean loadShapeFunc(){

	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/

//	String file="C:\\Works\\PlayModel\\shape";
	String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Small model\\shapeAve\\shape";

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			
			line=br.readLine();
			
			line=br.readLine();
			sp=line.split(regex);	
			
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
			
			this.numb=Integer.parseInt(sp[0]);
			this.Bs=Double.parseDouble(sp[1]);
			this.Hs=Double.parseDouble(sp[2]);
			this.numb2=Integer.parseInt(sp[3]);
			if(sp.length>4){
			this.numb3=Integer.parseInt(sp[4]);
			this.numb4=Integer.parseInt(sp[5]);
			}
			//this.numb5=Integer.parseInt(sp[6]);
			
		
			//util.show(sp);
			
			
			Mat[] BH1=new Mat[250];
			int[] LBH1=new int[250];
			
			int Lx=10000;
			for( int p=0;p<BH1.length;p++)
			BH1[p]=new Mat(Lx,2);
		
			int iloop=0;
			int jx=0;
			while(br.ready()){
				line=br.readLine();
				sp=line.split(regex);	
				if(sp.length==1) break;
				double[] bh=getCSV(line);
				if(jx>0 && bh[0]<BH1[iloop].el[jx-1][0])
				{
					iloop++;
					jx=0;
				}

				BH1[iloop].el[jx][0]=bh[0];
				BH1[iloop].el[jx][1]=bh[1];
				
				LBH1[iloop]++;
				jx++;
			}
			
			this.BH=new Mat[iloop+1];
			for( int p=0;p<BH.length;p++){
				BH[p]=new Mat(LBH1[p],2);
				for( int i=0;i<BH[p].nRow;i++)
				{
					BH[p].el[i][0]=BH1[p].el[i][0];
					BH[p].el[i][1]=BH1[p].el[i][1];
				}

			}
			
			
			util.plotBunch(BH);
			
			Mat[][] VV=new Mat[7][BH.length];
			for( int t=0;t<7;t++){
			
				for( int p=0;p<BH.length;p++){
					VV[t][p]=new Mat(BH[p].nRow,2);
					for( int i=0;i<BH[p].nRow;i++)
					{
						VV[t][p].el[i][0]=BH[p].el[i][0];
						VV[t][p].el[i][1]=BH[p].el[i][1];
					}
				}
			}
			
			for( int p=0;p<BH.length;p++){
				util.pr(BH[p].el[0][0]);
			}
			
			String fileshape="C:\\Works\\EMSolBuild_C\\EMSolBatch\\Small model\\shapeAve\\shape7iso";
		
			this.writeShape(VV, Bs, Hs, numb2, 7, 1, fileshape);
			
		//	BH1[0].show();
			return true;
		
			}
			catch(IOException e){System.err.println("Error in loading BH data file.");
			return false;
			}

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
		
		
		public void writeShape(Mat[][] shapes, double Bs,double Hs,int nHysteron, int nAngs, int symMode,String file){


				try{
					PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(file)));		

					pwBun.println(1);
					pwBun.println(1+"\t"+Bs+"\t"+Hs+"\t"+nHysteron+"\t"+0+"\t"+0+"\t"+nAngs+"\t"+symMode);
		
					for(int i=0;i<shapes[0].length;i++){
						for(int j=0;j<shapes[0][i].nRow;j++){
							pwBun.print(shapes[0][i].el[j][0]+"\t");
							for(int ia=0;ia<nAngs;ia++)
								pwBun.print(shapes[ia][i].el[j][1]+"\t");
							pwBun.println();
						}
							
					}
					pwBun.println(0);
					pwBun.println(0);
						

					util.pr("Simulated angle-dependent shape data was written to "+file+".");

					pwBun.close();
				}
				catch(IOException e){}
				

		}
	
	
}
