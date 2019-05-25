package emsolution;

import java.io.BufferedReader;

import java.io.FileReader;
import java.io.IOException;

import java.util.Arrays;

import math.Mat;
import math.Vect;
import math.util;
import fem.*;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 15, 2012.
 */
public class HystDataLoader {

	private String regex="[:; ,=\\t]+";
	private String regex2="[\\[\\]\\s: ,=\\t]+";

	
public static void main2(String[] args){
	
	
}
	

	
	
	public Mat[][] loadDataSym( String file){

		Mat[][] syms=new Mat[2][1];

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;


			line=br.readLine();

			sp=line.split(regex);	

			int nc=Integer.parseInt(sp[0]);
			
			syms=new Mat[2][nc];
			
			line=br.readLine();
			
			for(int i=0;i<nc;i++){
				int jx=0;
				double[][] data=new double[1000][3];
				
				line=br.readLine();

			while(line!=null && line.length()>0){
				data[jx]=this.getCSV(line);
				jx++;
				line=br.readLine();

			}
			
			syms[0][i]=new Mat(jx,2);
			syms[1][i]=new Mat(jx,2);
			for(int j=0;j<jx;j++){
				syms[0][i].el[j][0]=data[jx-1-j][1];
				syms[0][i].el[j][1]=data[jx-1-j][0];
				
				syms[1][i].el[j][0]=data[jx-1-j][2];
				syms[1][i].el[j][1]=data[jx-1-j][0];
	
				
			}
			
		}

		
			}
		
			catch(IOException e){System.err.println("Error in loading BH data file.");
		
			}
		
		return syms;

}	



	private Vect getVectData(String line, int dim){

		String[] sp=line.split(regex2);	
	
		Vect v=new Vect(dim);
	
		int k=0;
		while(sp[k].equals("")){k++;}	

		k++;
			
		for( int p=k;p<sp.length;p++){
		v.el[p-k]=Double.parseDouble(sp[p]);
		}

		return v;
	}



	
	private double[] getTabedData(String line){
		String[] sp=line.split(regex);	
		int L=sp.length;
		double[] v=new double[L];
		for( int p=0;p<L;p++)
			v[p]=Double.parseDouble(sp[p]);

		return v;
	}
	
	private double[] getPair(String line){
		String[] sp=line.split(regex);	
		double[] v=new double[2];
		int k=0;
		if(sp[k].equals("")) k++;

		for( int p=0;p<2;p++)
			v[p]=Double.parseDouble(sp[k+p]);

		return v;
	}
	
	private double getScalarData(String line){
		String[] sp=line.split(regex);	
		return Double.parseDouble(sp[sp.length-1]);
	}

	private int getIntData(String line){
		String[] sp=line.split(regex);	
		return Integer.parseInt(sp[sp.length-1]);
	}

	private boolean getBooleanData(String line){
		boolean b=false;
		String[] sp=line.split(regex);	
		
		if(sp[sp.length-1].startsWith("t"))	
			b=true;
		
		return b;

	}

	private String getStringData(String line){
		String[] sp=line.split(regex);	
		String[] sp2=sp[sp.length-1].split(" ");
		return sp2[sp2.length-1];
	}



public double[] loadArray(){
	String file=util.getFile();
	if(file==null || file.equals("") )  throw new NullPointerException("file not found.");
	return loadArray(file);
}

public double[] loadArray(String arrayPath){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		int N=100000;
		
		double[] x1=new double[N];
		
		int i=0;
		line=br.readLine();
		while(line!=null){
			if(i>N) break;
			x1[i++]=Double.parseDouble(line);
			line=br.readLine();
			
		}

		double[] x=Arrays.copyOf(x1, i);
		
			return x;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}	

public double[][] loadArrays(int n, int m,String arrayPath){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

	
		
		double[][] A=new double[n][m];
		
		for(int i=0;i<n;i++){
			line=br.readLine();
			if(line==null) continue;
			double[] x=getCSV(line);
			for(int j=0;j<m;j++)
				A[i][j]=x[j];
	
			
			
		}

		
			return A;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
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

private int[] getCSInt(String line){
	String[] sp=line.split(regex);	
	int L=sp.length;
	int[] v=new int[L];
	for( int p=0;p<L;p++)
				v[p]=Integer.parseInt(sp[p]);

	return v;
}


}
