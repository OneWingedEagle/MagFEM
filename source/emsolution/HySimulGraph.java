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
import math.util;
import fem.Model;



public class HySimulGraph {
	
Mat[] BH;
int numb1;
double Bs, Hs;
int numb2,numb3,numb4,numb5;
String regex="[:; ,\\t]+";


public static void main(String[] args){

	HySimulGraph pg=new HySimulGraph();
	
	pg.loadHysData();
}
	
public boolean loadHysData(){

	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
//	String file="C:\\Works\\HVID\\folder1\\data\\A_B“ü—Í‘ÎÌƒ‹[ƒvhts_data\\hys_data";

	String file="C:\\Works\\HVID\\hys_data_p";

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			this.numb4=10;
			
			BH=new Mat[numb4];
			

			for( int ip=0;ip<numb4;ip++){
				line=br.readLine();
				
					int L1=Integer.parseInt(line);
						
		
				BH[ip]=new Mat(L1,2);

				for( int i=0;i<L1;i++){
					line=br.readLine();
		
					double[] bh=getCSV(line);

					BH[ip].el[i][0]=bh[1];
					BH[ip].el[i][1]=bh[0];
				}
			
						
			}
		
			util.plotBunch(BH,2);
			//BH[1].show();
			
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
	
	
}
