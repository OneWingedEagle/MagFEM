package io;

import static java.lang.Math.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Scanner;

import fem.Edge;
import fem.EdgeSet;
import fem.Model;
import fem.Node;
import math.Complex;
import math.Mat;
import math.Vect;
import math.util;

public class Writer {

	public static void main2(String[] args) throws IOException{
		Writer wr=new Writer();
		String bunQ=System.getProperty("user.dir")+"\\bunQ.txt";
		String nodeFile ="C:\\Users\\Hassan\\Desktop\\Triangle\\model.node";
		String elFile = "C:\\Users\\Hassan\\Desktop\\Triangle\\model.1.ele";
		//wr.writeNodeX(nodeFile);

		//	wr.writeMesh(nodeFile,elFile);

	}

	public void writeMesh(Model model,String bunFilePath){
		writeMesh(model,bunFilePath,false);

	}

	public void writeMesh(Model model,String bunFilePath, boolean modeShape){

		DecimalFormat formatter;
		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.00000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		
	
			pwBun.println(model.elType);
			pwBun.println("//Number_of_Node");
			pwBun.println(model.numberOfNodes);

			pwBun.println("//Number_of_Element");
			pwBun.println(model.numberOfElements);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=model.numberOfElements;i++){
				int[] vertNumb=model.element[i].getVertNumb();
				for(int j=0;j<model.nElVert;j++)
					pwBun.print(vertNumb[j]+",");
				pwBun.println();
			}
			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect xyz;

				xyz=model.node[i].getCoord();

				if(modeShape && model.node[i].isDeformable()){ xyz=xyz.add(model.node[i].getU()); }
				for(int j=0;j<model.dim;j++){
					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(model.region[i].getFirstEl()+","+model.region[i].getLastEl()+","+model.region[i].getName());

			}

			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();
			File file=new File(bunFilePath);
			long filesize = file.length();
			double filesizeInMB=(double)filesize /(1024*1024);

			pwBun.close();
		}
		catch(IOException e){}

	}
	
	public void writeMeshOneNode(Model model,String bunFilePath){

		DecimalFormat formatter;
		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.00000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		
			pwBun.println(model.elType);
			pwBun.println("//Number_of_Node");
			pwBun.println(model.numberOfNodes);

			pwBun.println("//Number_of_Element");
			pwBun.println(model.numberOfElements);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=model.numberOfElements;i++){
				int[] vertNumb=model.element[i].getVertNumb();
				for(int j=0;j<model.nElVert;j++){
					if(model.node[vertNumb[j]].common && model.node[vertNumb[j]].getMap()>0)
						pwBun.print(model.node[vertNumb[j]].getMap()+",");				
					else	pwBun.print(vertNumb[j]+",");
				}
				pwBun.println();
			}
			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect xyz;;

				xyz=model.node[i].getCoord();
				for(int j=0;j<model.dim;j++){
					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(model.region[i].getFirstEl()+","+model.region[i].getLastEl()+","+model.region[i].getName());

			}

			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();
			File file=new File(bunFilePath);
			long filesize = file.length();
			double filesizeInMB=(double)filesize /(1024*1024);

			pwBun.close();

		}
		catch(IOException e){}

	}

	public void writeModeShape(Model model,String bunFilePath, double a){

		int dim=model.dim;
		DecimalFormat formatter;
		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.00000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		
			pwBun.println(model.elType);
			pwBun.println("//Number_of_Node");
			pwBun.println(model.numberOfNodes);

			pwBun.println("//Number_of_Element");
			pwBun.println(model.numberOfElements);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			double r0=model.node[1].getCoord().norm();

			for(int i=1;i<=model.numberOfElements;i++){
				int[] vertNumb=model.element[i].getVertNumb();
				for(int j=0;j<model.nElVert;j++)
					pwBun.print(vertNumb[j]+",");
				pwBun.println();
			}
			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect xyz;

				xyz=model.node[i].getCoord();
				
				if(model.node[i].isDeformable()){
				Vect v=model.node[i].getU().times(a);
					if(model.coordCode==1){

						Mat R=new Mat();

						double ang=util.getAng(model.node[i].getCoord().v2());
						
						R=util.rotMat2D(ang);
						
						if(model.dim==2){
							
							v=R.mul(v);		

						}
						else {

							Mat R3=new Mat(dim,dim);
							for(int p=0;p<2;p++)
								for(int q=0;q<2;q++)
									R3.el[p][q]=R.el[p][q];

							R3.el[2][2]=1;		
							v=R3.mul(v);	
						

						
						}

					
					}
					 xyz=xyz.add(v);
				}
			
				for(int j=0;j<model.dim;j++){
					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(model.region[i].getFirstEl()+","+model.region[i].getLastEl()+","+model.region[i].getName());

			}

			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();
			File file=new File(bunFilePath);
			long filesize = file.length();
			double filesizeInMB=(double)filesize /(1024*1024);

			pwBun.close();
		}
		catch(IOException e){}

	}


	public void writeMeshTriToQuad2(Model model,String bunFilePath , boolean deformed){
		DecimalFormat formatter;
		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000");
		else 
			formatter= new DecimalFormat("0.00");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("quadrangle");
			pwBun.println("//Number_of_Node");
			pwBun.println(model.numberOfNodes/2);

			pwBun.println("//Number_of_Element");
			pwBun.println(model.numberOfElements);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=model.numberOfElements;i++){ 
				int[] vertNumb=model.element[i].getVertNumb();
				for(int j=4;j<model.nElVert;j++)
					pwBun.print(vertNumb[j]+",");
				pwBun.println();
			}

			int m=0;

			for(int i=1;i<=model.numberOfNodes/2;i++){ 
				Vect xyz;;

				Vect v=new Vect(model.dim);
				v.rand();
				v=v.add(-.5).times(.1);


				xyz=model.node[i].getCoord();
				if(m==1)
					if(xyz.norm()<.99) xyz=xyz.add(v);

				for(int j=0;j<2;j++){

					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(model.region[i].getFirstEl()+","+model.region[i].getLastEl()+","+model.region[i].getName());

			}

			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();
			File file=new File(bunFilePath);
			long filesize = file.length();
			double filesizeInMB=(double)filesize /(1024*1024);
			//System.out.println(" Size of the generated mesh file :"+trunc(filesizeInMB,3)+" MB");

			pwBun.close();
		}
		catch(IOException e){}
	}

	public void writeMeshq2h(Model model,String bunFilePath , boolean deformed){

		double h=.02;

		DecimalFormat formatter;
		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		
			pwBun.println("hexahedron");
			pwBun.println("//Number_of_Node");
			pwBun.println(2*model.numberOfNodes);

			pwBun.println("//Number_of_Element");
			pwBun.println(model.numberOfElements);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=model.numberOfElements;i++){
				int[] vertNumb=model.element[i].getVertNumb();

				for(int j=4;j<8;j++)
					pwBun.print((model.numberOfNodes+vertNumb[j-4])+",");
				for(int j=0;j<model.nElVert;j++)
					pwBun.print(vertNumb[j]+",");

				pwBun.println();
			}
			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect xyz;;

				xyz=model.node[i].getCoord();
				for(int j=0;j<model.dim;j++)
					pwBun.print(formatter.format(xyz.el[j]*model.scaleFactor)+" ,");
				pwBun.print(formatter.format(0*model.scaleFactor)+" ,");

				pwBun.println();	

			}

			for(int i=model.numberOfNodes+1;i<=2*model.numberOfNodes;i++){ 
				Vect xyz;;

				xyz=model.node[i-model.numberOfNodes].getCoord();
				for(int j=0;j<model.dim;j++)
					pwBun.print(formatter.format(xyz.el[j]*model.scaleFactor)+" ,");
				pwBun.print(formatter.format(h*model.scaleFactor)+" ,");

				pwBun.println();	

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(model.region[i].getFirstEl()+","+model.region[i].getLastEl()+","+model.region[i].getName());

			}

			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();
			File file=new File(bunFilePath);
			long filesize = file.length();
			double filesizeInMB=(double)filesize /(1024*1024);
			//	System.out.println(" Size of the generated mesh file :"+trunc(filesizeInMB,3)+" MB");

			pwBun.close();
		}
		catch(IOException e){}
		if(model.dim==3){
			String bun2D=System.getProperty("user.dir")+"\\bun2D.txt";
			model.writeMeshTriToQuad2(bun2D,false);
		}
	}


	public void writeMeshTriToTri(Model model,String bunFilePath,double r1, double r2){

		Vect[] nodep=new Vect[model.numberOfNodes+model.numberOfElements+1];
		for(int i=1;i<=model.numberOfNodes;i++)
			nodep[i]=model.node[i].getCoord().deepCopy();




		int[][] elv=new int[3*model.numberOfElements+1][3];

		int[] ixr=new int[model.numberOfRegions+1];
		int in=0;
		int ix=0;

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				ixr[ir]++;
				ix++;
				int[] vertNumb=model.element[i].getVertNumb();
				Vect v=model.getElementCenter(i);
				if(v.norm()<r1 || v.norm()>r2) {
					elv[ix]=vertNumb;

					continue;
				}

				in++;
				nodep[in+model.numberOfNodes]=v.deepCopy();
				elv[ix][0]=vertNumb[0];
				elv[ix][1]=model.numberOfNodes+in;
				elv[ix][2]=vertNumb[2];

				ixr[ir]++;
				ix++;	

				elv[ix][0]=vertNumb[0];
				elv[ix][1]=vertNumb[1];
				elv[ix][2]=model.numberOfNodes+in;

				ixr[ir]++;
				ix++;
				elv[ix][0]=model.numberOfNodes+in;
				elv[ix][1]=vertNumb[1];
				elv[ix][2]=vertNumb[2];

			}


		int nNodeNew=model.numberOfNodes+in;;
		int nElNew=ix;
		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("triangle");
			pwBun.println("//Number_of_Node");
			pwBun.println(nNodeNew);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=nElNew;i++){ 
				for(int j=0;j<3;j++)
					pwBun.print(elv[i][j]+",");
				pwBun.println();
			}
			for(int i=1;i<=nNodeNew;i++){ 
				Vect xyz;;

				xyz=nodep[i].deepCopy();
				for(int j=0;j<2;j++){

					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			int[][] regEnds=new int[model.numberOfRegions+1][2];
			regEnds[1][0]=1;
			regEnds[1][1]=ixr[1];
			for(int i=2;i<=model.numberOfRegions;i++){
				regEnds[i][0]=regEnds[i-1][1]+1;
				regEnds[i][1]=regEnds[i][0]-1+ixr[i];

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}



			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();

			pwBun.close();
		}
		catch(IOException e){}
	}


	public void writeMeshTriToTriW(Model model,String bunFilePath){


		int[][] elEdge=new int[model.numberOfElements+1][3];

		EdgeSet eds=new EdgeSet();
		int[][] edgeNodes=eds.getEdgeNodesXY(model,elEdge);

		/*		int[] nodeEdge1=new int[model.numberOfEdges+1];
		int[] nodeEdge2=new int[model.numberOfEdges+1];
		for(int i=1;i<=model.numberOfNodes;i++){
			nodeEdge1[edgeNodes[i][0]]=i;
			nodeEdge2[edgeNodes[i][1]]=i;
		}
		 */                   
		int nEdge=edgeNodes.length-1;

		Vect[] nodep=new Vect[nEdge+model.numberOfNodes+1];
		for(int i=1;i<=model.numberOfNodes;i++)
			nodep[i]=model.node[i].getCoord().deepCopy();

		for(int i=1;i<=nEdge;i++){
			nodep[model.numberOfNodes+i]=model.node[edgeNodes[i][0]].getCoord().add(model.node[edgeNodes[i][1]].getCoord()).times(.5);
		}



		int[][] elv=new int[4*model.numberOfElements+1][3];

		int ne=0;

		for(int i=1;i<=model.numberOfElements;i++){

			int[] vertNumb=model.element[i].getVertNumb();

			ne++;
			elv[ne][0]=vertNumb[0];
			elv[ne][1]=model.numberOfNodes+elEdge[i][0];
			elv[ne][2]=model.numberOfNodes+elEdge[i][2];

			ne++;
			elv[ne][0]=vertNumb[1];
			elv[ne][1]=model.numberOfNodes+elEdge[i][1];
			elv[ne][2]=model.numberOfNodes+elEdge[i][0];

			ne++;
			elv[ne][0]=vertNumb[2];
			elv[ne][1]=model.numberOfNodes+elEdge[i][2];
			elv[ne][2]=model.numberOfNodes+elEdge[i][1];

			ne++;
			elv[ne][0]=model.numberOfNodes+elEdge[i][0];
			elv[ne][1]=model.numberOfNodes+elEdge[i][1];
			elv[ne][2]=model.numberOfNodes+elEdge[i][2];


		}

		int nNodeNew=nEdge+model.numberOfNodes;
		int nElNew=4*model.numberOfElements;
		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("triangle");
			pwBun.println("//Number_of_Node");
			pwBun.println(nNodeNew);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=nElNew;i++){ 
				for(int j=0;j<3;j++)
					pwBun.print(elv[i][j]+",");
				pwBun.println();
			}


			for(int i=1;i<=nNodeNew;i++){ 
				Vect xyz;;

				xyz=nodep[i].deepCopy();
				for(int j=0;j<2;j++){

					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			int[][] regEnds=new int[model.numberOfRegions+1][2];
			regEnds[1][0]=1;
			regEnds[1][1]=4*model.region[1].getLastEl();
			for(int i=2;i<=model.numberOfRegions;i++){
				regEnds[i][0]=regEnds[i-1][1]+1;
				regEnds[i][1]=regEnds[i][0]-1+4*model.region[i].getNumbElements();

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}



			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();


			pwBun.close();
		}
		catch(IOException e){}
	}

	public void writeMeshTriToQuad(Model model,String bunFilePath){


		int[][] elEdge=new int[model.numberOfElements+1][3];

		EdgeSet eds=new EdgeSet();
		int[][] edgeNodes=eds.getEdgeNodesXY(model,elEdge);

		int nEdge=edgeNodes.length-1;

		Vect[] nodep=new Vect[model.numberOfNodes+nEdge+model.numberOfElements+1];
		for(int i=1;i<=model.numberOfNodes;i++)
			nodep[i]=model.node[i].getCoord().deepCopy();

		for(int i=1;i<=nEdge;i++){
			nodep[model.numberOfNodes+i]=model.node[edgeNodes[i][0]].getCoord().add(model.node[edgeNodes[i][1]].getCoord()).times(.5);
		}


		int nx=model.numberOfNodes+nEdge;

		int[][] elv=new int[3*model.numberOfElements+1][4];

		for(int i=1;i<=model.numberOfElements;i++){

			int[] vertNumb=model.element[i].getVertNumb();

			nodep[++nx]=model.getElementCenter(i);

			for(int j=0;j<model.nElVert;j++){
				int ne=3*(i-1)+j+1;
				elv[ne][0]=vertNumb[j];
				elv[ne][1]=model.numberOfNodes+elEdge[i][j];
				elv[ne][2]=nx;
				elv[ne][3]=model.numberOfNodes+elEdge[i][(j+2)%3];
			}
		}

		int nNodeNew=nx;
		int nElNew=3*model.numberOfElements;
		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("quadrangle");
			pwBun.println("//Number_of_Node");
			pwBun.println(nNodeNew);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=nElNew;i++){ 
				for(int j=0;j<4;j++)
					pwBun.print(elv[i][j]+",");
				pwBun.println();
			}
			for(int i=1;i<=nNodeNew;i++){ 
				Vect xyz;;

				xyz=nodep[i].deepCopy();
				for(int j=0;j<2;j++){

					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			int[][] regEnds=new int[model.numberOfRegions+1][2];
			regEnds[1][0]=1;
			regEnds[1][1]=3*model.region[1].getLastEl();
			for(int i=2;i<=model.numberOfRegions;i++){
				regEnds[i][0]=regEnds[i-1][1]+1;
				regEnds[i][1]=regEnds[i][0]-1+3*model.region[i].getNumbElements();

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}



			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();


			pwBun.close();
		}
		catch(IOException e){}
	}


	public void writeMeshTriToQuadCoarse(Model model,String bunFilePath){


		int[][] elEdge=new int[model.numberOfElements+1][3];

		EdgeSet eds=new EdgeSet();
		int[][] edgeNodes=eds.getEdgeNodesXY(model,elEdge);


		int nEdges=edgeNodes.length-1;
		int[] ed=new int[3];
		int[][] edgeEls=new int[nEdges+1][2];
		for(int i=1;i<=model.numberOfElements;i++){
			for(int j=0;j<3;j++)
				ed[j]=elEdge[i][j];

			Vect c=model.getElementCenter(i);
			Vect v1=model.node[edgeNodes[ed[0]][1]].getCoord().sub(model.node[edgeNodes[ed[0]][0]].getCoord());
			Vect v2=c.sub(model.node[edgeNodes[ed[0]][0]].getCoord());
			int k=0;

			if(v1.cross(v2).el[2]<0) k=1;

			for(int j=0;j<3;j++)	
				edgeEls[elEdge[i][j]][k]=i;

		}

		//util.show(edgeEls);

		boolean[] edgeCommon=new boolean[nEdges+1];

		for(int i=1;i<=nEdges;i++)
			if(edgeEls[i][0]>0 && edgeEls[i][1]>0 )edgeCommon[i]=true;

		boolean[] ner=new boolean[model.numberOfNodes+1];
		boolean[] removeEdge=new boolean[nEdges+1];
		int ix=0;
		int iy=0;
		int[][] regEnds=new int[model.numberOfRegions+1][2];
		regEnds[1][0]=1;
		int[] dd=new int[3];
		for(int ir=1;ir<=model.numberOfRegions;ir++){
			iy=0;
			if(ir>1)
				regEnds[ir][0]=regEnds[ir-1][1]+1;

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				for( int j=0;j<3;j++)
					dd[j]=elEdge[i][j];

				for( int j=0;j<3;j++){
					if(!ner[edgeNodes[dd[j]][0]] && !ner[edgeNodes[dd[j]][1]] && edgeCommon[dd[j]]) {

						removeEdge[dd[j]]=true;
						ner[edgeNodes[dd[j]][0]]=true;
						ner[edgeNodes[dd[j]][1]]=true;
						ix++;
						iy++;
					}

				}


			}
			regEnds[ir][1]=regEnds[ir][0]+iy-1;
		}
		int nElNew=ix;

		int[][] elv=new int[nElNew+1][4];

		ix=0;
		for(int i=1;i<=nEdges;i++){


			if(!removeEdge[i]  ) continue;

			int[] vn1=model.element[edgeEls[i][0]].getVertNumb();
			int[] vn2=model.element[edgeEls[i][1]].getVertNumb();

			int n1=edgeNodes[i][0];
			int n2=edgeNodes[i][1];

			for(int j=0;j<3;j++)
				if(vn2[j]!=n1){}

			int[] vnq=new int[4];

			vnq[0]=n1;
			vnq[2]=n2;

			for(int j=0;j<3;j++)
				if(vn1[j]!=n1 && vn1[j]!=n2){vnq[1]=vn1[j]; break;}

			for(int j=0;j<3;j++)
				if(vn2[j]!=n1 && vn2[j]!=n2){vnq[3]=vn2[j]; break;}


			++ix;
			for(int j=0;j<4;j++){
				elv[ix][j]=vnq[j];
			}

		}



		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("quadrangle");
			pwBun.println("//Number_of_Node");
			pwBun.println(model.numberOfNodes);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=nElNew;i++){ 
				for(int j=0;j<4;j++)
					pwBun.print(elv[i][j]+",");
				pwBun.println();
			}
			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect xyz;;

				xyz=model.node[i].getCoord();
				for(int j=0;j<2;j++){
					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}


			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}



			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();


			pwBun.close();
		}
		catch(IOException e){}
	}

	public void writeMeshQuadToTri(Model model,String bunFilePath){


		int nElNew=2*model.numberOfElements;

		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("triangle");
			pwBun.println("//Number_of_Node");
			pwBun.println(model.numberOfNodes);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=model.numberOfElements;i++){ 
				int[] vertNumb=model.element[i].getVertNumb();
				
				boolean is3ang=false;
				int jx=-1,kx=-1;
				for(int j=0;j<vertNumb.length;j++)
					for(int k=0;k<j;k++)
						if(vertNumb[k]==vertNumb[j]) {
							jx=j;
							kx=k;
							is3ang=true;
						}
				if(is3ang){
					int cc=0;
					int[] v3ang=new int[3];
					int px=0;
					int tx=-1;
					for(int j=0;j<vertNumb.length;j++){
						util.pr(j+"  "+jx+"  "+kx+"  "+cc);
						if(j==jx || j==kx){
							cc++;
						}
						if(cc<2){
							v3ang[px++]=vertNumb[j];
						}
						else{
							tx=j;
						}
					}
					pwBun.println(v3ang[0]+","+v3ang[1]+","+v3ang[2]);
					pwBun.println(vertNumb[tx]+","+vertNumb[tx]+","+vertNumb[tx]);
				}
					
				
				else{
				double d1=model.node[vertNumb[0]].getCoord().sub(model.node[vertNumb[2]].getCoord()).norm();
				double d2=model.node[vertNumb[1]].getCoord().sub(model.node[vertNumb[3]].getCoord()).norm();
				if(d1>d2){

					pwBun.println(vertNumb[0]+","+vertNumb[1]+","+vertNumb[3]);
					pwBun.println(vertNumb[1]+","+vertNumb[2]+","+vertNumb[3]);
				}
				else{

					pwBun.println(vertNumb[0]+","+vertNumb[1]+","+vertNumb[2]);
					pwBun.println(vertNumb[2]+","+vertNumb[3]+","+vertNumb[0]);
				}
				}
			}
		
			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect xyz;;

				xyz=model.node[i].getCoord();
				for(int j=0;j<2;j++){

					pwBun.print(formatter.format(xyz.el[j]*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			int[][] regEnds=new int[model.numberOfRegions+1][2];
			regEnds[1][0]=1;
			regEnds[1][1]=2*model.region[1].getLastEl();
			for(int i=2;i<=model.numberOfRegions;i++){
				regEnds[i][0]=regEnds[i-1][1]+1;
				regEnds[i][1]=regEnds[i][0]-1+2*model.region[i].getNumbElements();

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}


			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();
			File file=new File(bunFilePath);
			long filesize = file.length();
			double filesizeInMB=(double)filesize /(1024*1024);
			System.out.println(" Size of the generated mesh file :"+trunc(filesizeInMB,3)+" MB");

			pwBun.close();
		}
		catch(IOException e){}
	}
	
	public void writeMeshHexaToTetra(Model model,String bunFilePath){

		int nSubEl=5;
	
		
		
		int[][] tetNodes={{1,3,4,0},{1,4,6,5},{1,6,3,2},{4,3,6,7},{1,4,3,6}};
		int[][] tetNodes2={{0,2,7,3},{5,0,7,4},{5,2,0,1},{5,7,2,6},{0,7,2,5}};


		int[][] elv=new int[nSubEl*model.numberOfElements+1][4];

		int ix=0;
		for(int i=1;i<=model.numberOfElements;i++){

			int[] vertNumb=model.element[i].getVertNumb();

		

			
			for(int ne=0;ne<nSubEl;ne++){

				ix++;
				
				for(int j=0;j<4;j++)
					if(i%2==1)
					elv[ix][j]=vertNumb[tetNodes[ne][j]];
					else
						elv[ix][j]=vertNumb[tetNodes2[ne][j]];

			}

		}

		int nElNew=nSubEl*model.numberOfElements;
		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("tetrahedron");
			pwBun.println("//Number_of_Node");
			pwBun.println(model.numberOfNodes);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=nElNew;i++){ 
				for(int j=0;j<4;j++)
					pwBun.print(elv[i][j]+",");
				pwBun.println();
			}
			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect v=model.node[i].getCoord();
				for(int j=0;j<3;j++){

					pwBun.print(formatter.format((v.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			int[][] regEnds=new int[model.numberOfRegions+1][2];
			regEnds[1][0]=1;
			regEnds[1][1]=nSubEl*model.region[1].getLastEl();
			for(int i=2;i<=model.numberOfRegions;i++){
				regEnds[i][0]=regEnds[i-1][1]+1;
				regEnds[i][1]=regEnds[i][0]-1+nSubEl*model.region[i].getNumbElements();

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}



			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();


			pwBun.close();
		}
		catch(IOException e){}
	}

	public void writeMeshHexaToPyramid(Model model,String bunFilePath){

		int nSubEl=6;
		int nNodes1=model.numberOfNodes;
		Vect[] nodep=new Vect[nNodes1+model.numberOfElements+1];
		for(int i=1;i<=model.numberOfNodes;i++)
			nodep[i]=model.node[i].getCoord().deepCopy();
		int[][] baseNodes={{1,2,6,5},{3,0,4,7},{3,7,6,2},{0,1,5,4},{6,7,4,5},{3,2,1,0}};

		int[][] elv=new int[nSubEl*model.numberOfElements+1][5];

		int ix=0;
		for(int i=1;i<=model.numberOfElements;i++){

			int[] vertNumb=model.element[i].getVertNumb();

			nodep[nNodes1+i]=model.getElementCenter(i);
			int inx=i+nNodes1;

			for(int ne=0;ne<nSubEl;ne++){

				ix++;

				for(int j=0;j<4;j++)
					elv[ix][j]=vertNumb[baseNodes[ne][j]];
				elv[ix][4]=inx;
			}

		}

		int nNodeNew=nNodes1+model.numberOfElements;
		int nElNew=nSubEl*model.numberOfElements;
		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("pyramid");
			pwBun.println("//Number_of_Node");
			pwBun.println(nNodeNew);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=nElNew;i++){ 
				for(int j=0;j<5;j++)
					pwBun.print(elv[i][j]+",");
				pwBun.println();
			}
			for(int i=1;i<=nNodeNew;i++){ 
				for(int j=0;j<3;j++){

					pwBun.print(formatter.format((nodep[i].el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			int[][] regEnds=new int[model.numberOfRegions+1][2];
			regEnds[1][0]=1;
			regEnds[1][1]=nSubEl*model.region[1].getLastEl();
			for(int i=2;i<=model.numberOfRegions;i++){
				regEnds[i][0]=regEnds[i-1][1]+1;
				regEnds[i][1]=regEnds[i][0]-1+nSubEl*model.region[i].getNumbElements();

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}



			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();


			pwBun.close();
		}
		catch(IOException e){}
	}
	
	public void writeMeshHexaToPrism(Model model,String bunFilePath){
		writeMeshHexaToPrism(model,bunFilePath,2);
	}
	
	public void writeMeshHexaToPrism(Model model,String bunFilePath, int dir){

		int nSubEl=2;
		
		int[][][] prismNodesAll={{{0,3,7,1,2,6},{0,7,4,1,6,5}},
				{{0,1,2,4,5,6},{2,3,0,6,7,4}},
				{{0,1,2,4,5,6},{2,3,0,6,7,4}}};
				
				
				int[][] prismNodes=prismNodesAll[dir];
		
		int[][] elv=new int[nSubEl*model.numberOfElements+1][6];

		int ix=0;
		for(int i=1;i<=model.numberOfElements;i++){

			int[] vertNumb=model.element[i].getVertNumb();

			for(int ne=0;ne<nSubEl;ne++){

				ix++;

				for(int j=0;j<6;j++)
					elv[ix][j]=vertNumb[prismNodes[ne][j]];
			}

		}

		int nElNew=nSubEl*model.numberOfElements;
		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("prism");
			pwBun.println("//Number_of_Node");
			pwBun.println(model.numberOfNodes);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=nElNew;i++){ 
				for(int j=0;j<6;j++)
					pwBun.print(elv[i][j]+",");
				pwBun.println();
			}
			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect v=model.node[i].getCoord();
				for(int j=0;j<3;j++){
					pwBun.print(formatter.format(v.el[j]*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			int[][] regEnds=new int[model.numberOfRegions+1][2];
			regEnds[1][0]=1;
			regEnds[1][1]=nSubEl*model.region[1].getLastEl();
			for(int i=2;i<=model.numberOfRegions;i++){
				regEnds[i][0]=regEnds[i-1][1]+1;
				regEnds[i][1]=regEnds[i][0]-1+nSubEl*model.region[i].getNumbElements();

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}



			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();


			pwBun.close();
		}
		catch(IOException e){}
	}

	
	public void writeMeshPrismToHexa(Model model,String bunFilePath){

		
		// has problems **********************
		
		model.setEdge();
		
		int ix=0;

		for(int i=1;i<=model.numberOfEdges;i++)
			if(model.edge[i].direction!=2)
			{
				ix++;
			}
		
		
		
		int nNodex=ix;
		
		Vect[] nodep=new Vect[model.numberOfNodes+nNodex+2*model.numberOfElements+1];
		for(int i=1;i<=model.numberOfNodes;i++)
			nodep[i]=model.node[i].getCoord().deepCopy();

		ix=model.numberOfNodes;
		int[] mapEdNd=new int[1+model.numberOfEdges];
		for(int i=1;i<=model.numberOfElements;i++){

			int[] edgeNumb=model.element[i].getEdgeNumb();
			for(int j=0;j<model.nElEdge;j++){
				
				int ned=edgeNumb[j];
				
				if(mapEdNd[ned]>0) continue;
			
			//int[] end=model.edge[ned].endNodeNumber;
			
			if(model.edge[ned].direction!=2){
				ix++;
				nodep[ix]=model.edge[ned].node[0].getCoord().add(model.edge[ned].node[1].getCoord()).times(.5);
				mapEdNd[ned]=ix;
			}
			}
		}
		
		

		int[][] mapElNd=new int[1+model.numberOfElements][2];
		int[][] fc=new int[1+model.numberOfNodes][10];
		int[] ic=new int[1+model.numberOfNodes];
		
		


		int nx=model.numberOfNodes+nNodex;
		for(int i=1;i<=model.numberOfElements;i++){

			int[] vertNumb=model.element[i].getVertNumb();
			
			int[] vn1=new int[3];
			for(int k=0;k<3;k++)
				vn1[k]=vertNumb[k];

		
			int px=hasCommon(fc[vn1[0]],fc[vn1[1]],fc[vn1[2]]);
		
			if(px==-1){
			nodep[++nx]=model.node[vertNumb[0]].getCoord().add(model.node[vertNumb[1]].getCoord()).
			add(model.node[vertNumb[2]].getCoord()).times(1.0/3);
			
			mapElNd[i][0]=nx;
			
			for(int k=0;k<3;k++){
			fc[vn1[k]][ic[vn1[k]]++]=nx;
		
			}
			}
			else
				mapElNd[i][0]=px;
			
			int[] vn2=new int[3];
			for(int k=0;k<3;k++)
				vn2[k]=vertNumb[k+3];
			
			px=hasCommon(fc[vn2[0]],fc[vn2[1]],fc[vn2[2]]);

			if(px==-1){
			nodep[++nx]=model.node[vertNumb[3]].getCoord().add(model.node[vertNumb[4]].getCoord()).
			add(model.node[vertNumb[5]].getCoord()).times(1.0/3);
			
			mapElNd[i][1]=nx;
			

			for(int k=0;k<3;k++){
				fc[vn2[k]][ic[vn2[k]]++]=nx;
			}
			}
			else
				mapElNd[i][1]=px;
		}


			
		


		
		int[][] elv=new int[3*model.numberOfElements+1][8];

		for(int i=1;i<=model.numberOfElements;i++){

			int[] vertNumb=model.element[i].getVertNumb();
			int[] edgeNumb=model.element[i].getEdgeNumb();
			
		
			for(int j=0;j<3;j++){
				int ne=3*(i-1)+j+1;
				elv[ne][0]=vertNumb[j];
				elv[ne][1]=mapEdNd[edgeNumb[j]];
				elv[ne][2]=mapElNd[i][0];
				elv[ne][3]=mapEdNd[edgeNumb[(j+2)%3]];
				
				elv[ne][4]=vertNumb[j+3];
				elv[ne][5]=mapEdNd[edgeNumb[j+3]];
				elv[ne][6]=mapElNd[i][1];
				elv[ne][7]=mapEdNd[edgeNumb[(j+2)%3+3]];
			
			}
			
		
			
		}
		

		int nNodeNew=nx;
		
		util.pr(nNodeNew+"  "+nodep.length+"  "+model.numberOfElements);
		int nElNew=3*model.numberOfElements;
		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("hexahedron");
			pwBun.println("//Number_of_Node");
			pwBun.println(nNodeNew);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=nElNew;i++){ 
				for(int j=0;j<8;j++)
					pwBun.print(elv[i][j]+",");
				pwBun.println();
			}
			for(int i=1;i<=nNodeNew;i++){ 
				Vect xyz;;
				xyz=nodep[i].deepCopy();
				for(int j=0;j<3;j++){

					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			int[][] regEnds=new int[model.numberOfRegions+1][2];
			regEnds[1][0]=1;
			regEnds[1][1]=3*model.region[1].getLastEl();
			for(int i=2;i<=model.numberOfRegions;i++){
				regEnds[i][0]=regEnds[i-1][1]+1;
				regEnds[i][1]=regEnds[i][0]-1+3*model.region[i].getNumbElements();

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}


			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();


			pwBun.close();
		}
		catch(IOException e){}
	}
	
	private int hasCommon(int[] x1, int[] x2, int[]x3){
		int px=-1;
		
		for(int i=0;i<x1.length;i++)
			for(int j=0;j<x2.length;j++)
				if(x1[i]!=0 && x1[i]==x2[j])
					for(int k=0;k<x3.length;k++)
						if(x1[i]==x3[k]) return x1[i];
		
		return px;
	}

	public void writeMeshPyramidToTetra(Model model,String bunFilePath){

		int nSubEl=2;


		int[][] tetNodes={{0,1,2,4},{2,3,0,4}};

		int[][] elv=new int[nSubEl*model.numberOfElements+1][4];

		int ix=0;
		for(int i=1;i<=model.numberOfElements;i++){

			int[] vertNumb=model.element[i].getVertNumb();

			for(int ne=0;ne<nSubEl;ne++){
				ix++;

				for(int j=0;j<4;j++)
					elv[ix][j]=vertNumb[tetNodes[ne][j]];
			}

		}

		int nNodeNew=model.numberOfNodes;
		int nElNew=nSubEl*model.numberOfElements;
		DecimalFormat formatter;

		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.0000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		

			pwBun.println("tetrahedron");
			pwBun.println("//Number_of_Node");
			pwBun.println(nNodeNew);

			pwBun.println("//Number_of_Element");
			pwBun.println(nElNew);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			for(int i=1;i<=nElNew;i++){ 
				for(int j=0;j<4;j++)
					pwBun.print(elv[i][j]+",");
				pwBun.println();
			}
			for(int i=1;i<=nNodeNew;i++){ 
				for(int j=0;j<3;j++){

					pwBun.print(formatter.format(model.node[i].getCoord(j)*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			int[][] regEnds=new int[model.numberOfRegions+1][2];
			regEnds[1][0]=1;
			regEnds[1][1]=nSubEl*model.region[1].getLastEl();
			for(int i=2;i<=model.numberOfRegions;i++){
				regEnds[i][0]=regEnds[i-1][1]+1;
				regEnds[i][1]=regEnds[i][0]-1+nSubEl*model.region[i].getNumbElements();

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(regEnds[i][0]+","+regEnds[i][1]+","+model.region[i].getName());

			}



			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();


			pwBun.close();
		}
		catch(IOException e){}
	}



	public void writeB(Model model,String fluxFile){
		int dim=model.dim;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(fluxFile)));

			
				pw.println("flux");
				pw.println(dim);
				pw.println(model.numberOfElements);
				for(int i=1;i<=model.numberOfElements;i++){
					Vect B=model.element[i].getB();
					/*if(model.element[i].hasJ())*/
					//	B=new Vect(model.getElement2DA(i),0);
							
					
				//B=new Vect(model.edge[model.element[i].getEdgeNumb(0)].A,0);
					
					for(int k=0;k<dim;k++)					
						pw.format("%E\t",B.el[k]);

					pw.println();
				}
			

			pw.close();
		} catch(IOException e){System.out.println("writing flux file failed.");}

		System.out.println(" Magnetic flux density was written to "+fluxFile);

	}
	
	public void writeB1Layer(Model model,int nr,String fluxFile){
		int dim=2;
		int nEls=model.region[nr].getNumbElements()/6;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(fluxFile)));

			
				pw.println("flux");
				pw.println(dim);
				pw.println(nEls);
				for(int ir=1;ir<=model.numberOfRegions;ir++){
					if(ir!=nr) continue;
					for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
						if(i>nEls) break;
					Vect B=model.element[i].getB();

					for(int k=0;k<dim;k++)					
						pw.format("%E\t",B.el[k]);

					pw.println();
				}
			}

			pw.close();
		} catch(IOException e){System.out.println("writing flux file failed.");}

		System.out.println(" Magnetic flux density was written to "+fluxFile);

	}


	

	public void writeA(Model model,String vPotFile){
		int dim=model.dim;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(vPotFile)));
			pw.println("vPot");
			pw.println(dim);
			pw.println(model.numberOfEdges);
			for(int i=1;i<=model.numberOfEdges;i++){
				double A=model.edge[i].A;
				pw.println(A);

			}
			pw.close();
		} catch(IOException e){System.out.println("writing vector potential file failed.");}

		System.out.println(" Magnetic vector potential was written to "+vPotFile);

	}	
	
	public void writeA_as_flux(Model model,String a_file){
		int dim=model.dim;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(a_file)));

			
				pw.println("flux");
				pw.println(dim);
				pw.println(model.numberOfElements);
				for(int i=1;i<=model.numberOfElements;i++){

					Vect A=model.getElementA(i);

					for(int k=0;k<dim;k++)					
						pw.format("%E\t",A.el[k]);

					pw.println();
				}
			

			pw.close();
		} catch(IOException e){System.out.println("writing flux file failed.");}

		System.out.println(" A field density was written to "+a_file);
	}	

	public void writeJe(Model model,String eddyFile){
		//int dim=model.dim;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(eddyFile)));
			pw.println("fluxJ");
			pw.println(3);
			pw.println(model.numberOfElements);

			for(int ir=1;ir<=model.numberOfRegions;ir++){

				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				if(model.element[i].isConductor()){

					Vect J=new Vect(3);
					if(model.region[ir].isConductor) 	{				
							J=model.element[i].getJe();
							pw.format("%d\t",i);
							for(int k=0;k<3;k++)
								pw.format("%12.6e\t",J.el[k]);
							pw.println();
					}else{
						//for(int k=0;k<3;k++)
							//pw.format("%20.4f",0.0);
						//pw.println();
					}
				}
				}
			}
			pw.close();
		} catch(IOException e){System.out.println("writing flux file failed.");}

		System.out.println(" Eddy current was written to "+eddyFile);
	}


	public void writeJ0(Model model,String J0File){
		int dim=3;
		try{

			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(J0File)));
			pw.println("fluxJ");
			pw.println(dim);
			pw.println(model.numberOfElements);
			for(int ir=1;ir<=model.numberOfRegions;ir++){

				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
					Vect J=new Vect(3);
		
				if(model.region[ir].hasJ) 	{				
						J=model.element[i].getJ();
						pw.format("%d\t",i);
						for(int k=0;k<3;k++)
							pw.format("%12.6e\t",J.el[k]);
						pw.println();
				}else{
					//for(int k=0;k<3;k++)
						//pw.format("%20.4f",0.0);
					//pw.println();
				}
			
				
				}
				
			}
			pw.close();
		} catch(IOException e){System.out.println("writing flux file failed.");}

		System.out.println(" J0 density was written to "+J0File);
	}




	public void writeNodalField(Model model,String nodalForceFile,int mode){
		int dim=model.dim;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(nodalForceFile)));
			if(mode==1)
				pw.println("nodal");
			else if(mode==2)
				pw.println("nodal");
			else if(mode==3)
				pw.println("nodal");
			else if(mode==-1)
				pw.println("displacement");


			pw.println(dim);

			pw.println(model.numberOfNodes);

			for(int n=1;n<=model.numberOfNodes;n++)
			{

	Vect r=model.node[n].getCoord();
	if(r.norm()<model.rm) continue;

				Vect v=model.node[n].getNodalVect(mode);

				if(v==null) continue;
			
				if(v.norm()==0) continue;

		
				if(model.coordCode==1){

					Mat R=new Mat();

					double ang=util.getAng(model.node[n].getCoord().v2());
					
					R=util.rotMat2D(ang);
					
					if(model.dim==2){
						
						v=R.mul(v);		

					}
					else {

						Mat R3=new Mat(dim,dim);
						for(int p=0;p<2;p++)
							for(int q=0;q<2;q++)
								R3.el[p][q]=R.el[p][q];

						R3.el[2][2]=1;		
						v=R3.mul(v);	
					

					
					}

				
				}
				
				if(mode==-1) v=v.times(1e9);
		
				pw.format("%d\t",n);

				for(int k=0;k<dim;k++)
					pw.format("%E\t",v.el[k]);

				pw.println();

			}
			
			pw.close();
		} catch(IOException e){System.out.println("writing nodal file failed.");}
		
		String str="force";
		if(mode==-1) str="displacement";

		System.out.println(" Magnetic nodal "+str+" was written to "+nodalForceFile);
	}
	
	
	
	public void writeNodalScalar(Model model,String file){
		int dim=model.dim;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(file)));
			
				pw.println("temperature");


			pw.println(dim);

			pw.println(model.numberOfNodes);

			for(int n=1;n<=model.numberOfNodes;n++)
			{


				
				if(model.node[n].T==0) continue;

				pw.format("%d\t%E",n,model.node[n].T);


				pw.println();
			}
			
			pw.close();
		} catch(IOException e){}

		System.out.println(" Temperature was written to "+file);
	}
	
	
	public void writePhi(Model model,String file){
		int dim=model.dim;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(file)));
			
				pw.println("temperature");


			pw.println(dim);

			pw.println(model.numberOfNodes);

			for(int n=1;n<=model.numberOfNodes;n++)
			{
				if(!model.node[n].isPhiVar()) continue;
				if(model.node[n].getPhi()==0) continue;

				pw.format("%d\t%E",n,model.node[n].getPhi());


				pw.println();
			}
			
			pw.close();
		} catch(IOException e){}

		System.out.println(" Temperature was written to "+file);
	}
	
	


	public void writeVects(Vect[] v,String vectFile){
		int N=v.length;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(vectFile)));

				for(int i=0;i<N;i++){
						for(int k=0;k<v[i].length;k++)					
						pw.format("%E\t",v[i].el[k]);

					pw.println();
				}
			

			pw.close();
		} catch(IOException e){e.printStackTrace();}

		System.out.println(" vectors were written to "+vectFile);

	}
	
	public void writeArray(double[] x,String file){
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(file)));

		
		
						for(int k=0;k<x.length;k++)					
						pw.println(x[k]);

			

			pw.close();
		} catch(IOException e){e.printStackTrace();}

		System.out.println(" array was vertically written to "+file);

	}
	
	public void writeFFT(Complex[] x,String file){
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(file)));

		
						pw.println(x.length);
				
						for(int k=0;k<x.length;k++)					
						pw.println(x[k].re+","+x[k].im);

			

			pw.close();
		} catch(IOException e){e.printStackTrace();}

		System.out.println(" FFT was  written to "+file);

	}
	
	
	public void writeMat(Mat A,String file){

	
		writeArray(A.el,file);
				
			}
	
	
	public void writeArray(double[][] A,String file){

Vect[] v=new Vect[A.length];

		int N=v.length;
		for(int i=0;i<N;i++){
			v[i]=new Vect(A[i]);
		}
		writeVects(v,file);
		
	}

	/*public void writeTemper(Model model,String temperFile){
		model.setTemp(1);
		int dim=model.dim;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(temperFile)));
			pw.println("temperature");
			pw.println(dim);		
			pw.println(model.numberOfNodes);

				for(int n=1;n<=model.numberOfNodes;n++)
					{

						double t=model.node[n].getTemp();
						pw.println(t);
					}

			pw.close();
		} catch(IOException e){System.out.println("writing temperature file failed.");}

		System.out.println(" temperature was written to "+temperFile);
	}*/

	private double trunc(double a,int n){
		return floor(a*pow(10,n))/pow(10,n);
	}
	
	@SuppressWarnings("javadoc")
	public void copyFile(String ss, String sd)
		        throws IOException {
		
		File source = new File(ss);
		File dest = new File(sd);	
		    InputStream input = null;
		    OutputStream output = null;
		    try {
		        input = new FileInputStream(source);
		        output = new FileOutputStream(dest);
		        byte[] buf = new byte[1024];
		        int bytesRead;
		        while ((bytesRead = input.read(buf)) > 0) {
		            output.write(buf, 0, bytesRead);
		        }
		    } finally {
		        input.close();
		        output.close();
		    }
		}
	
	
	
	public void writeMeshNeu(Model model,String bunFilePath){

		String line=null;

		int[] numVerts=new int[8];
		String lineZeros="0,0,0,0,0,0,0,0,0,0,";
		String lineZerosLong="0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,";
		
		DecimalFormat formatter;
		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.00000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		
			pwBun.println("   -1");
			pwBun.println("   100");
			pwBun.println("<NULL>");
			pwBun.println("8.,");
			pwBun.println("   -1");
			pwBun.println("   -1");
			pwBun.println("   403");
			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect xyz;

				xyz=model.node[i].getCoord().v3();
				pwBun.print(i+",0,0,1,46,0,0,0,0,0,0,");
				
				
				for(int j=0;j<3;j++){
					pwBun.print(formatter.format((xyz.el[j]))+" ,");
				}
				pwBun.println("0,");	

			}

			pwBun.println("   -1");
			pwBun.println("   -1");
			pwBun.println("   404");
			for(int i=1;i<=model.numberOfElements;i++){
				
				//if(i>10) break;
				int[] vertNumb=model.element[i].getVertNumb();
				int nz=0;
				for(int j=0;j<vertNumb.length;j++){
					numVerts[j]=vertNumb[j];
					if(j>0 && (numVerts[j]==numVerts[j-1] || numVerts[j]==0)) nz++;
				}
				boolean tet=false;
				if(nz==4 && numVerts[3]==numVerts[4]) tet=true;
				
				boolean wedge=false;
				
		
				if(model.elCode==2) tet=true;
				if(model.elCode==3) wedge=true;
				
				boolean pyramid=(nz==3);
				
				int xx=124;
				int type=4;
				
				if(numVerts[3]==numVerts[2]) numVerts[3]=0;
				
				if(tet) type=6;
				else if(wedge) type=14;
				
				else if(pyramid) type=7;
				
				else 	if(model.elCode==0) type=2;
				
				else 	if(model.elCode==1) type=4;
				
				else if(nz==0) type=8;

				
				else if(numVerts[3]==0) type=2;
				

				if(tet || wedge)
				 line=i+","+xx+","+(model.element[i].getRegion())+",25,"+type+",1,0,0,9,0,0,0,0,";
				else if(model.elCode==4)
					 line=i+","+xx+","+(model.element[i].getRegion())+",25,"+type+",1,0,0,9,0,0,0,0,";
				else
					 line=i+","+xx+","+(model.element[i].getRegion())+",17,"+type+",1,0,0,9,0,0,0,0,";
				 pwBun.println(line);

				 if(wedge)
					 line=numVerts[0]+","+numVerts[1]+","+numVerts[2]+",0,"+numVerts[3]+","+numVerts[4]+","+ numVerts[5]+",0,0,0,000000000";
				 else if(tet)
					 line=numVerts[0]+","+numVerts[1]+","+numVerts[2]+",0,"+numVerts[3]+",0,0,0,0,0,";
				 else if(pyramid)
					 line=numVerts[0]+","+numVerts[1]+","+numVerts[2]+","+numVerts[3]+","+numVerts[4]+",0,0,0,0,0,";
				// else	 if(nz==0)						 line=numVerts[0]+","+numVerts[1]+","+numVerts[2]+numVerts[3]+","+numVerts[4]+","+numVerts[5]+","+numVerts[6]+numVerts[7]+",0,0,";
				 else
				 line=numVerts[0]+","+numVerts[1]+","+numVerts[2]+","+numVerts[3]+","+numVerts[4]+","+ numVerts[5]+","+numVerts[6]+","+numVerts[7]+",0,0,";
				 pwBun.println(line);
				 pwBun.println(lineZeros);
				 pwBun.println("0.,0.,0.,");
				 pwBun.println("0.,0.,0.,");
				 pwBun.println("0.,0.,0.,");
				 pwBun.println(lineZerosLong);

			}
			pwBun.println("   -1");
			
			pwBun.flush();

			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);
			System.out.println();


			pwBun.close();
		}
		catch(IOException e){}

	}
	

	public void reportData(Model model){
		System.out.println();

		DecimalFormat formatter = new DecimalFormat("0.00E0");

		System.out.println(" Element type: "+model.elType);
		System.out.println(" Number of regions: "+model.numberOfRegions);
		System.out.println(" Number of elements: "+model.numberOfElements);
		System.out.println(" Number of nodes   : "+model.numberOfNodes+"    known: "+model.numberOfKnownPhis+" , unknown: "+model.numberOfVarNodes);
		System.out.println(" Number of edges   : "+model.numberOfEdges+"    known: "+model.numberOfKnownEdges+" , unknown: "+model.numberOfUnknownEdges);	
		System.out.println(" Number of unknown currents   : "+model.network.no_unknown_currents);	
		System.out.println(" Total number of unknows   : "+model.numberOfUnknowns);	
		System.out.println();



	}
	

	public double outputLoss(Model model,String file,int step,double phase_time,boolean append){
		

		double totalLoss=0;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(file,append)));
			pw.println("Step: "+step+"\t phase/time: "+phase_time);
			pw.println();
			util.pr("step "+step);
			pw.println(" ****** ============= ******** ==================== ******");
				pw.println("Joule Losses [W]");
				pw.println();
				util.pr("Joule Losses [W]");
				for(int ir=1;ir<=model.numberOfRegions;ir++){
					if((model.region[ir].isConductor && model.analysisMode>0) || model.coilInicesByRegion[ir]>=0){
				double  loss=0;
				
				loss=model.obtainLoss(ir);
				
				totalLoss+=loss;
				pw.format("%5d %12.5e",ir,loss);
				pw.println();
				System.out.format("%5d %25.12e\n",ir,loss);

					}
				}
				pw.format("total %12.5e",totalLoss);
				pw.println();
		
			pw.close();
		} catch(IOException e){System.err.println("IOException: " + e.getMessage());}

		
		//System.out.println(" Temperature was written to "+file);

		
		return totalLoss;

	}


	public double outputEnergies(Model model,String file,int step,double phase_time,boolean append){
	


		double totalEnergy=0;
		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(file,append)));
				
				pw.println("Energies [VAR]");
				pw.println(" ****** -------------------------------------- ******");
				pw.println();
				util.pr("Energies ");
				for(int ir=1;ir<=model.numberOfRegions;ir++){

				double  energy=model.obtainEnergies(ir);
				totalEnergy+=energy;
				pw.format("%5d %12.5e",ir,energy);
				pw.println();
				System.out.format("%5d %12.5e\n",ir,energy);
					
				}
				pw.format("total %12.5e",totalEnergy);
				pw.println();
				pw.println(" ****** -------------------------------------- ******");
				pw.println();
			pw.close();
		} catch(IOException e){System.err.println("IOException: " + e.getMessage());}


		
		return totalEnergy;

	}



	
public void writeR_L(Model model,String file,Vect freqs, Vect[] R_L){
		
		boolean b=true;

		try{
			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(file,b)));
			pw.println("R and L versus frequency:");
			pw.println(" ****** -------------------------------------- ******");
			pw.println();
			for(int i=0;i<freqs.length;i++){
				pw.format("%12.5e ",freqs.el[i]);
				pw.format("%12.5e ",R_L[i].el[0]);
				pw.format("%12.5e ",R_L[i].el[1]);
				pw.println();
			}	
			pw.println(" ****** -------------------------------------- ******");
				pw.println();
			pw.close();
		} catch(IOException e){System.err.println("IOException: " + e.getMessage());}

	

	}
}
