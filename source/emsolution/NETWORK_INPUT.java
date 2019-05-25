package emsolution;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import math.Vect;
import meshFactory.MeshManipulator;


public class NETWORK_INPUT {
	String regex="[:; ,\\t]+";
	
	
	public static void main(String[] args){

		NETWORK_INPUT network=new NETWORK_INPUT();
		
		network.parallelFEM();
	}


	
	
public void parallelFEM(){

	
	int nFEMs=250;
		
	//DecimalFormat formatter= new DecimalFormat("0.000E00");


	try{


		String fout=System.getProperty("user.dir")+"\\EMSol\\netwrok.txt";
		PrintWriter networkWriter = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		


		int ng=1000;
		int nNodes1=1001;

		int nNodes=nNodes1;
		int nElem=1000;

		double R=1e-4;
		int nFEM=1;
		for(int i=0;i<nFEMs;i++){
			networkWriter.print("R"+"\t"+nElem+"\t"+nNodes+"\t"+(nNodes+1+i)+"\t"+R);
			networkWriter.println();
			nElem++;
			//nNodes++;
			networkWriter.print("FEM"+"\t"+nFEM+"\t"+(nNodes+1+i)+"\t"+ng+"\t"+nFEM);
			networkWriter.println();
			
	

			nFEM++;
		}
		
		
		networkWriter.close();
		
		System.out.println();
		System.out.println(" Network data was written to:");
		System.out.println("    "+fout);
	}
	catch(IOException e){ System.err.println("error");e.printStackTrace(); }



}


}
