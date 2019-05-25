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



public class ShapeGen {
	
Mat[] BH;
int numb;
double Bs, Hs;
int numb2,numb3,numb4,numb5;
String regex="[:; ,\\t]+";


public static void main(String[] args){

	ShapeGen pg=new ShapeGen();
	
	pg.genetate();
}
	
public void genetate(){

	/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/

	double dang=15;
	double angMax=90;
	int nAng;
	int symMode=0;
	int nHys=36;
	int sym=1;
	double Bs=1.8;
	double dB=Bs/nHys;
	
	nAng=(int)Math.floor(angMax/dang);
	if(sym==1){
		nAng++;
	}
	

Vect ang=new Vect(nAng);
for(int ia=0;ia<nAng;ia++)
	ang.el[ia]=ia*dang;


	Mat[][] shape=new Mat[nAng][nHys];
	
	for(int ia=0;ia<nAng;ia++)
		for(int i=0;i<nHys;i++){
			double Bp=Bs-i*dB;
			int ndiv=nHys+1-i;
		shape[ia][i]=new Mat(ndiv,2);	
		Vect bb=new Vect().linspace(-Bp, Bp, ndiv);
		for(int j=0;j<ndiv;j++){
			shape[ia][i].el[j][0]=bb.el[j];
			shape[ia][i].el[j][1]=Math.pow(10*bb.el[j],3)/(1+i*.2)*(1+2*Math.sin(ia*Math.PI/nAng));
		}
	}

	
	double Hs=shape[0][0].el[shape[0][0].nRow-1][1];
	
	String file="C:\\Works\\EMSolBuild_C\\EMSolBatch\\ThinDisk\\Small model\\shape";

		try{
			
			DecimalFormat formatter= new DecimalFormat("0.00");
			
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));		
			
			pw.println(1);
			pw.print(1+"\t");
			pw.print(formatter.format(Bs)+"\t");
			pw.print(formatter.format(Hs)+"\t");
			pw.print(nHys+"\t");
			pw.print(0+"\t");
			pw.print(0+"\t");
			pw.print(nAng+"\t");
			pw.print(sym);
			pw.println();
			
			for(int ia=0;ia<nAng;ia++){
				pw.println(formatter.format(ang.el[ia]));
				for(int i=0;i<nHys;i++){
				for(int j=0;j<shape[ia][i].nRow;j++){
					pw.println(formatter.format(shape[ia][i].el[j][0])+"\t"+formatter.format(shape[ia][i].el[j][1]));
				}
				}
				
			}
			
			pw.println(0);
			pw.println(0);
			
			pw.close();
			
			util.pr(" sahe functions was written to "+file);
		
			}
			catch(IOException e){System.err.println(" failed.");
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
