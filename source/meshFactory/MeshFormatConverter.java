package meshFactory;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import fem.BoundarySet;
import fem.Element;
import fem.Model;
import fem.Node;
import fem.Region;
import io.Loader;
import io.Writer;
import math.Complex;
import math.Mat;
import math.SpBlockMat;
import math.SpMat;
import math.Vect;
import math.util;



public class MeshFormatConverter {

	String regex="[ : ,=\\t]+";
	BoundarySet bs=new BoundarySet();
	Writer writer=new Writer();

	public static void main(String[] args){
		MeshFormatConverter mfc=new MeshFormatConverter();
		//mfc.getPreHexAtlasOrig();
		
	//mfc.convertToNeu();

		//String file=System.getProperty("user.dir")+"\\ladder_current";
	//	String file=System.getProperty("user.dir")+"\\magnetic";
		//mfc.getEMSolFlux(file,3,102700,"BMAG");
		//mfc.getPostHexAtlasOrig();
		//mfc.getPostHexaNeu();
	//	mfc.getPostHexAtlas();
		//	mfc.getNeuMeshQ();
	//	mfc.getNeuMeshHexa();
	//	mfc.getPostHexaNeu(8);
	//	mfc.get2DFormNeu(4);
		
		mfc.getFromIr3(8);
		//mfc.convertTetraNeu();
/*int K=1000000;
double dr=.1/K;
		for(int k=0;k<K;k++){
		double R=3.29+k*dr;
		double fR=R*R*PI-2*R*R*(Math.acos(1-4/R))+2*(R-4)*Math.sqrt(4*(2*R-4))+3*PI-(R-3)*(R-3)*PI;
		if(abs(fR)<.0000001){
		util.pr(R+"  "+fR+"  "+2*(R-3));
		}
		}*/
		
			
		boolean anim=false;
		if(anim){
		//String bun="C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\viewer\\board\\bun.txt";
		String bun="C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\viewer\\2dFourier\\bun.txt";
		Model model=new Model(bun);
		
		Vect unitB=new Vect(0,0,1);
		Vect zeroB=new Vect(0,0,.5);
		
		for(int k=1;k<=model.numberOfNodes;k++){
			model.node[k].setDeformable(true);
		//	if(model.node[k].getCoord(2)<-.0001)
			//	model.node[k].setCoord(new Vect(2,2,model.node[k].getCoord(2)-1));
		}
		//model.writeMesh("C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\viewer\\2dFourier\\bun1.txt");
		Vect ff=new Vect(300);
		double d=2;
		double a=.8;
		double fill=PI*a*a/d/d;
		boolean cmplx=false;
		double dc=.0;
		double peak=2;
		double deltaz=peak-dc;
		boolean disp=true;
		
		int nGmax=3;
	
		for(int nG=nGmax;nG<=nGmax;nG++){
		
			//if(nG<10) continue;
		for(int i=0;i<1;i++)
			for(int j=0;j<1;j++){
				
				String file1="C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\viewer\\2dFourier\\disp"+nG+".txt";
			
				
				if(nG==-10){
					double hh=1;
				int ir=1;
				for(ir=1;ir<=2;ir++)
					for(int k=model.region[ir].getFirstEl();k<=model.region[ir].getLastEl();k++){
						int vn[] =model.element[k].getVertNumb();
						for(int m=0;m<4;m++)
							model.node[vn[m]].setCoord(2,1.4);
						//	model.node[vn[m]].setU(2, hh*1e-7);
					}

				}
				else{
				
					
				for(int k=1;k<=model.numberOfNodes;k++){
					//model.node[k].setU(2, 0);
					Vect coord=	model.node[k].getCoord();
					for(int px=0;px<3;px++)
						for(int py=0;py<3;py++){
							Vect coord_shifted=	coord.sub(new Vect(2*px,2*py,0));
							double r=coord_shifted.v2().norm();
							//util.pr(r);
						
					if(r<100.72999 &&coord_shifted.el[2]>1e-9){
					
						double x=coord_shifted.el[0];
						double y=coord_shifted.el[1];
						//util.pr(nG);
						double cc=0;
						double hh=0;
						double sumR=0;
						
						if(cmplx){
							Complex sum=new Complex(0,0);
						for(int u=-nG;u<=nG;u++){
							for(int v=-nG;v<=nG;v++){
								double rho=new Vect(2*u*PI/d,2*v*PI/d).norm();
								if(rho==0)
									cc=fill*deltaz;
								else
									cc=2*deltaz*fill*util.J1(rho*a)/(rho*a);
								//util.pr(cc);
								Complex tt=new Complex(Math.cos(2*PI*u/d*x+2*PI*v/d*y),Math.sin(2*PI*u/d*x+2*PI*v/d*y));
								
								sum =sum.add(tt.times(cc));

							}
							}
						hh=sum.re;
						}else{

							double sum=0;

							for(int u=0;u<=nG;u++){
								double partialsumU=Math.cos(2*PI*u/d*x);
								
								for(int v=0;v<=nG;v++){
									double partialsumV=Math.cos(2*PI*v/d*y);
								
									double rho=new Vect(2*u*PI/d,2*v*PI/d).norm();
									if(rho==0){
										cc=fill*deltaz;
										sum+=cc;
									}
									else{
										cc=2*deltaz*fill*util.J1(rho*a)/(rho*a);
										if(u==0)
											sum+=2*cc*partialsumV;
										else if(v==0)
											sum+=2*cc*partialsumU;
										else
											sum+=4*cc*partialsumU*partialsumV;//  coeff must be 4 times to represent cylinder
									}
								}
								}
							hh=sum;
						}
					
						if(x>=0 &&y==0) ff.el[(int)((x+1e-4)/.025)]=hh;
						//ff.el[nG]=cc;
						if(!disp)
						model.node[k].setCoord(2,dc+hh);
						else
							model.node[k].setU(2, hh*1e-7);
					}
						}
				}
					
				}
				
				for(int k=1;k<=0*model.numberOfElements;k++){
					Vect cc=model.getElementCenter(k);	
					int[] vn=model.element[k].getVertNumb();
					//cc.hshow();
					if(cc.el[0]>i&& cc.el[0]<i+1 &&cc.el[1]>j&& cc.el[1]<j+1){
						//model.element[k].setB(unitB);
					//}
				
					for(int p=0;p<4;p++){
						model.node[vn[p+4]].setU(2, 1e-7);
					}
					}else{
						//model.element[k].setB(unitB);
					//}
				
					for(int p=0;p<4;p++){
						model.node[vn[p+4]].setU(2, 1e-10);
					}
					}
				}
					/*else {if(ix==95)
						model.element[k].setB(zeroB.times(0));
					else 
						model.element[k].setB(zeroB);
					}*/
				if(disp){
				String filex="C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\viewer\\2dFourier\\disp"+nG+".txt";
				model.writeNodalField(filex, -1);
				}else{
				String file2="C:\\Users\\Me\\Desktop\\Farsi\\2018\\darbareh\\viewer\\2dFourier\\bun"+nG+".txt";
				model.writeMesh(file2);
				}
			}
				
		}
				//model.writeB(file1);
		util.plot(ff);

		}

	}




	public void getNeuMeshQ(){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];



			for(int i=1;i<100000;i++)
			{

				line=br.readLine();
				sp=line.split(regex);
				if(sp.length==15 && !sp[0].equals("0")) break;


			}



			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				sp=line.split(regex);

				line=br.readLine();

				if(sp.length!=15) break;


				nn=Integer.parseInt(sp[0]);
				nx++;

				map[nn]=nx;


				coord1[nx]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));


			}

			nnMax=nx;


			int[][] vernumb=new int[10*nnMax+1][4];
			int[] nReg=new int[1000000+1];

			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;

			int ix=0;

			line=br.readLine();
			//	line=br.readLine();


			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				line=br.readLine();

				sp=line.split(regex);
				if(sp.length<5) break;
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
				line=br.readLine();

				sp=line.split(regex);

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				vernumb[ix][0]=map[Integer.parseInt(sp[0])];
				vernumb[ix][1]=map[Integer.parseInt(sp[1])];
				vernumb[ix][2]=map[Integer.parseInt(sp[2])];
				vernumb[ix][3]=map[Integer.parseInt(sp[3])];

				for(int j=0;j<4;j++){
					if(vernumb[ix][j]==0) {
						if(j>0) vernumb[ix][j]=vernumb[ix][j-1];
						//ix--;
						break;
					}
				}

			}

			int nEl=ix;



			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<nReg.length;i++){
				if(nReg[i]!=-1)
					list1.add(nReg[i]);
			}

			Set<Integer> set = new HashSet<Integer>(list1);

			ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);

			int nRegions=regNumbs.size();


			int[] regNumber=new int[nRegions+1];
			for(int ir=1;ir<=nRegions;ir++)
				regNumber[ir]=regNumbs.get(ir-1);



			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];

			nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir])
						elOrd[++nx]=i;
				regEnd[ir][1]=nx;

			}


			int nNodes=nnMax;
			double scaleFactor=1;
			double scx=1;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");


			try{


				String fout=System.getProperty("user.dir")+"\\EMSol\\quad.txt";
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("quadrangle");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scx);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<4;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}

						pwBun.println();
					}

				Vect v;
				for(int i=1;i<=nnMax;i++){ 
					if(coord1[i]==null)
						v=new Vect(0,0);
					else
						v=coord1[i].deepCopy();
					for(int j=0;j<2;j++){
						pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
					}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);

				pwBun.close();
				br.close();
				fr.close();

				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");e.printStackTrace(); }



		}



		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}

	public void getNeuMeshHexa(){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];



			for(int i=1;i<100000;i++)
			{

				line=br.readLine();
				util.pr(line);
				sp=line.split(regex);
				if(sp.length==15 && !sp[0].equals("0")) break;


			}



			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				sp=line.split(regex);

				line=br.readLine();

				if(sp.length!=15) break;


				nn=Integer.parseInt(sp[0]);
				nx++;

				map[nn]=nx;


				coord1[nx]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));


			}

			nnMax=nx;


			int[][] vernumb=new int[10*nnMax+1][4];
			int[] nReg=new int[1000000+1];

			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;

			int ix=0;

			line=br.readLine();
			//	line=br.readLine();


			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				line=br.readLine();

				sp=line.split(regex);
				if(sp.length<5) break;
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
				line=br.readLine();

				sp=line.split(regex);
				int ib=0;
				if(sp[0].equals("")) ib++;
				

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				vernumb[ix][0]=map[Integer.parseInt(sp[ib+4])];
				vernumb[ix][1]=map[Integer.parseInt(sp[ib+5])];
				vernumb[ix][2]=map[Integer.parseInt(sp[ib+6])];
				vernumb[ix][3]=map[Integer.parseInt(sp[ib+7])];
				
				vernumb[ix][4]=map[Integer.parseInt(sp[ib+0])];
				vernumb[ix][5]=map[Integer.parseInt(sp[ib+1])];
				vernumb[ix][6]=map[Integer.parseInt(sp[ib+2])];
				vernumb[ix][7]=map[Integer.parseInt(sp[ib+3])];

				for(int j=0;j<8;j++){
					if(vernumb[ix][j]==0) {
						if(j>0) vernumb[ix][j]=vernumb[ix][j-1];
						//ix--;
						break;
					}
				}

			}

			int nEl=ix;



			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<nReg.length;i++){
				if(nReg[i]!=-1)
					list1.add(nReg[i]);
			}

			Set<Integer> set = new HashSet<Integer>(list1);

			ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);

			int nRegions=regNumbs.size();

			util.pr(nRegions);


			int[] regNumber=new int[nRegions+1];
			for(int ir=1;ir<=nRegions;ir++)
				regNumber[ir]=regNumbs.get(ir-1);



			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];

			nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir])
						elOrd[++nx]=i;
				regEnd[ir][1]=nx;

			}


			int nNodes=nnMax;
			double scaleFactor=1;
			double scx=1;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");


			try{


				String fout=System.getProperty("user.dir")+"\\EMSol\\quad.txt";
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("quadrangle");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scx);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<4;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}

						pwBun.println();
					}

				Vect v;
				for(int i=1;i<=nnMax;i++){ 
					if(coord1[i]==null)
						v=new Vect(0,0);
					else
						v=coord1[i].deepCopy();
					for(int j=0;j<2;j++){
						pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
					}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);

				pwBun.close();
				br.close();
				fr.close();

				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");e.printStackTrace(); }



		}



		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
	}
	
	
	public void getNeuMeshHexabbbbb(int mode){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String linep="";
			String[] sp=new String[15];
			
			

			if(mode==0)
				while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){		}

			if(mode==2){
				for(int k=0;k<8;k++)
					line=br.readLine();
			}else
	for(int i=1;i<100000;i++)
			{

				line=br.readLine();
				sp=line.split(regex);
				if(sp.length==14 && !sp[0].equals("0")) break;
			
				
			}
			
	

			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				sp=line.split(regex);
				
				line=br.readLine();
				if(sp.length!=14) break;

			
				nn=Integer.parseInt(sp[0]);
				nx++;
			
				map[nn]=nx;

				coord1[nx]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
		

			}
			nnMax=nx;
			int[][] vernumb=new int[10*nnMax+1][8];
			int[] nReg=new int[1000000+1];
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		
				line=br.readLine();
			//	line=br.readLine();
/*	for(int i=0;i<7;i++)
		line=br.readLine();*/

				
			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				linep=br.readLine();
				line=br.readLine();

				if(line==null) break;
				
				sp=line.split(regex);
	
				if(sp.length<5) break;
				
				if(sp[6].equals("0")) {
					for(int j=0;j<5;j++)
						line=br.readLine();
					continue;
				}
				

				if(i>-1500) {
				
				sp=linep.split(regex);
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
			
			
				sp=line.split(regex);
				
				for(int j=0;j<8;j++){
					if(vernumb[ix][j]==0) {
						
						nnMax++;

						vernumb[ix][j]=nnMax;
		
						if(j>0 && j!=4) coord1[nnMax]=coord1[vernumb[ix][j-1]];
						else if(j==1) coord1[nnMax]=coord1[vernumb[ix][1]];
						else if(j==4) coord1[nnMax]=coord1[vernumb[ix][5]];
						//ix--;
						//break;
					}
				}

				for(int j=0;j<8;j++){
					if(j<4)
					vernumb[ix][j]=map[Integer.parseInt(sp[j+4])];
					else
						vernumb[ix][j]=map[Integer.parseInt(sp[j-4])];
					
					//if(vernumb[ix][j]==1) util.pr(line);;//{ util.pr(ix);util.hshow(vernumb[ix]);}
					
				}
				
			

				for(int j=0;j<8;j++){
					if(vernumb[ix][j]==0) {
						
						nnMax++;

						vernumb[ix][j]=nnMax;
		
						if(j>0 && j!=4) coord1[nnMax]=coord1[vernumb[ix][j-1]];
						else if(j==1) coord1[nnMax]=coord1[vernumb[ix][1]];
						else if(j==4) coord1[nnMax]=coord1[vernumb[ix][5]];
						ix--;
						break;
					}
				}
				
				}
	
				
				for(int j=0;j<5;j++)
					line=br.readLine();

			}
			
			int nEl=ix;

		
			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<nReg.length;i++){
				if(nReg[i]!=-1)
					list1.add(nReg[i]);
			}
		
				Set<Integer> set = new HashSet<Integer>(list1);
				
				ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);
				
				int nRegions=regNumbs.size();
				
				util.pr(nRegions);

				
				int[] regNumber=new int[nRegions+1];
				for(int ir=1;ir<=nRegions;ir++)
					regNumber[ir]=regNumbs.get(ir-1);
				

					
			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];
			
			 nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir])
						elOrd[++nx]=i;
				regEnd[ir][1]=nx;
				
			}
			
			boolean[] nodeUsed=new boolean[vernumb.length];
			for(int i=1;i<=nEl;i++){
				for(int j=0;j<8;j++){
					nodeUsed[vernumb[elOrd[i]][j]]=true;
				}
			}

			int nNodes=0;
			for(int i=1;i<vernumb.length;i++)
				{
					if(nodeUsed[i]) nNodes++;
				}
			
			int[] nodeMapFinal=new int[vernumb.length];
			ix=1;
			for(int i=1;i<vernumb.length;i++)
			{
				if(nodeUsed[i]) nodeMapFinal[i]=ix++;
			}

			double scaleFactor=1000;
			double scx=1;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");
			

			try{


				String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";

				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("hexahedron");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scx);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<8;j++){
							int mpn=nodeMapFinal[vernumb[elOrd[i]][j]];
							pwBun.print(mpn+",");
						}
					
						pwBun.println();
					}

			Vect v;
			for(int i=1;i<=nnMax;i++){ 
				if(!nodeUsed[i]) continue;
				if(coord1[i]==null)
					v=new Vect(0,0,0);
				else
					v=coord1[i];
						for(int j=0;j<3;j++){
							pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
						}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);
				
				pwBun.close();
				br.close();
				fr.close();
				
				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");e.printStackTrace(); }

		
		
		}
		
		

		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	
	public void getPostMeshQ(){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
			
/*	while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){
	
			}*/

	for(int i=1;i<1000;i++)
			{

				line=br.readLine();
				sp=line.split(regex);
				//if(sp.length==15 && !sp[0].equals("0")) break;
				if(sp.length>5) break;
			
				
			}
			
	

			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			//int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				sp=line.split(regex);
				
				line=br.readLine();

				if(sp.length!=14) break;

			
				nn=Integer.parseInt(sp[0]);
				nx++;
			
			//	map[nn]=nx;
		

				coord1[nn]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
				

			}
			
			nnMax=nx;


			int[][] vernumb=new int[10*nnMax+1][4];
			int[] nReg=new int[1000000+1];
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		
				line=br.readLine();
			//	line=br.readLine();
			

			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				line=br.readLine();

				sp=line.split(regex);
				if(sp.length<5) break;
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
				line=br.readLine();
		
				sp=line.split(regex);
				
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

			/*	vernumb[ix][0]=map[Integer.parseInt(sp[0])];
				vernumb[ix][1]=map[Integer.parseInt(sp[1])];
				vernumb[ix][2]=map[Integer.parseInt(sp[2])];
				vernumb[ix][3]=map[Integer.parseInt(sp[3])];*/
				
				vernumb[ix][0]=Integer.parseInt(sp[0]);
				vernumb[ix][1]=Integer.parseInt(sp[1]);
				vernumb[ix][2]=Integer.parseInt(sp[2]);
				vernumb[ix][3]=Integer.parseInt(sp[3]);

				for(int j=0;j<4;j++){
					if(vernumb[ix][j]==0) {
						if(j>0) vernumb[ix][j]=vernumb[ix][j-1];
						//ix--;
						break;
					}
				}

			}
			
			int nEl=ix;

		
			
			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<nReg.length;i++){
				if(nReg[i]!=-1)
					list1.add(nReg[i]);
			}
		
				Set<Integer> set = new HashSet<Integer>(list1);
				
				ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);
				
				int nRegions=regNumbs.size();
				
				util.pr(nRegions);

				
				int[] regNumber=new int[nRegions+1];
				for(int ir=1;ir<=nRegions;ir++)
					regNumber[ir]=regNumbs.get(ir-1);
				

					
			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];
			
			 nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir])
						elOrd[++nx]=i;
				regEnd[ir][1]=nx;
				
			}


			int nNodes=nnMax;
			double scaleFactor=1;
			double scx=1;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");
			

			try{


				String fout=System.getProperty("user.dir")+"\\quadExtracted.txt";
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("quadrangle");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scx);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<4;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}
					
						pwBun.println();
					}

			Vect v;
			for(int i=1;i<=nnMax;i++){ 
				if(coord1[i]==null)
					v=new Vect(0,0);
				else
					v=coord1[i].deepCopy();
						for(int j=0;j<2;j++){
							pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
						}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);
				
				pwBun.close();
				br.close();
				fr.close();
				
				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");e.printStackTrace(); }

		
		
		}
		
		

		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	
	public void getNeuMeshHexOld(){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String linep="";
			String[] sp=new String[15];

			int numbAddedNodes=0;
	//if(mode==0)
	//while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){		}

			for(int i=0;i<7;i++)
				line=br.readLine();
			
			line=br.readLine();

			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				
				
				
				line=br.readLine();

				sp=line.split(regex);
			
				if(sp.length<14) break;

			
				nn=Integer.parseInt(sp[0]);
			//	nx++;
			nx=nn;
				
			map[nn]=nx;
			
		

				coord1[nx]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
		

			}
			
			
			nnMax=0;
			for(int k=0;k<map.length;k++)
				if(map[k]>nnMax) nnMax=map[k];


			int[][] vernumb=new int[10*nnMax+1][8];
			int[] nReg=new int[1000000+1];
			int firstElement=0;
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		int nEl=0;
	for(int i=0;i<2;i++)
		line=br.readLine();

			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				
				linep=br.readLine();

				line=br.readLine();
				if(line==null) break;
					
				sp=line.split(regex);
	
				if(sp.length<5) break;
				
				if(sp[6].equals("0")) {
					for(int j=0;j<5;j++)
						line=br.readLine();
					continue;
				}
				

				
				sp=linep.split(regex);
		
				ix=Integer.parseInt(sp[0]);
				if(firstElement==0) firstElement=ix;
				
				if(ix>nEl) nEl=ix;

				nReg[ix]=Integer.parseInt(sp[2]);
			
				sp=line.split(regex);


				for(int j=0;j<8;j++){

					vernumb[ix][j]=map[Integer.parseInt(sp[j])];
				}
				
			

				for(int j=0;j<8;j++){
					if(vernumb[ix][j]==0) {/*

						numbAddedNodes++;
				
						 vernumb[ix][j]=nnMax+numbAddedNodes;//vernumb[ix][j-1];
							if(j>=0){
						 coord1[numbAddedNodes]=coord1[vernumb[ix][j-1]].deepCopy().times(0);
							}
						//ix--;
						//break;
					*/}
				}
				
				// change order
				
				int[] tmp=new int[8];
				for(int j=0;j<8;j++){
					tmp[j]=vernumb[ix][j];
				
				}
				vernumb[ix][0]=tmp[0];
				vernumb[ix][1]=tmp[3];
				vernumb[ix][2]=tmp[2];
				vernumb[ix][3]=tmp[1];
				
				vernumb[ix][4]=tmp[4];
				vernumb[ix][5]=tmp[7];
				vernumb[ix][6]=tmp[6];
				vernumb[ix][7]=tmp[5];

				//======
				
				for(int j=0;j<5;j++)
					line=br.readLine();

			}
			
		

		
			
			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<=nEl;i++){
				if(nReg[i]==-1){
					nReg[i]=1;
					for(int j=0;j<8;j++)
						vernumb[i][j]=vernumb[firstElement][j];
				}
				else
					nReg[i]+=1;
				
					list1.add(nReg[i]);

			}

				Set<Integer> set = new HashSet<Integer>(list1);
				
				ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);
				
				int nRegions=regNumbs.size();
				
				util.pr(nRegions);

				
				int[] regNumber=new int[nRegions+1];
				for(int ir=1;ir<=nRegions;ir++)
					regNumber[ir]=regNumbs.get(ir-1);

					
			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];
			
			 nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
	
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir]){
						elOrd[i]=i;
						//elOrd[++nx]=i;
			
				nx++;
					}
				regEnd[ir][1]=nx;
				
				util.hshow(regEnd[ir]);
			}

			int nNodes=nnMax+numbAddedNodes;
			double scaleFactor=1;
			double scx=1;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");
			

			try{


				String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("hexahedron");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scx);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<8;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}
					
						pwBun.println();
					}

			Vect v;
			for(int i=1;i<=nnMax+numbAddedNodes;i++){ 
				if(coord1[i]==null)
					v=new Vect(0,0,0);
				else
					v=coord1[i].deepCopy();
						for(int j=0;j<3;j++){
							pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
						}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);
				
				pwBun.close();
				br.close();
				fr.close();
				
				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");e.printStackTrace(); }

		
		
		}
		
		

		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	public void getPreHexAtlas(){

		String s=util.getFile();
		
		String line;
		int max=2000000;
		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeMap=new int[max];

		int[] elMap=new int[max];
		int[] regNumb=new int[max];
		
		int[] nRegEls=new int[15];



		int nnx=1;
		int nex=1;
		
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
			while(true){
				line=br.readLine();
				
				if(line==null) break;
				if(line.startsWith("GRID")){
					nnx=readAtlasNodes(br,nodeMap,coord1,nnx);


				}
				
				
				if(line.startsWith("CONC")){	
		
					nex=readAtlasElements(br,elMap,regNumb,nRegEls,vertNumb,nex);
				}

			}
		
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }

			int nNodes=0;
			for(int i=0;i<nodeMap.length;i++)
				if(nodeMap[i]>0) nNodes++;


			int nEls=0;
		
			int nRegions=0;
			
			for(int i=0;i<nRegEls.length;i++)
				if(nRegEls[i]>0) {
					nRegions++;
					nEls+=nRegEls[i];
				}
			
		
			
			int[] regMap=new int[nRegEls.length];
			int nr=0;
			int[][]regEnd=new int[nRegions+1][2];
			regEnd[0][0]=1;
					
			for(int i=0;i<nRegEls.length;i++){
				
				if(nRegEls[i]>0) {
					nr++;
					regMap[i]=nr;
					regEnd[nr][0]=regEnd[nr-1][1]+1;
				
					regEnd[nr][1]=regEnd[nr][0]+nRegEls[i]-1;
				}
			}
			
			
			Model model=new Model(nRegions,nEls,nNodes,"hexahedron");
				
			for(int ir=1;ir<=nRegions;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+regMap[ir]);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				for(int j=0;j<8;j++){
					int mpn=nodeMap[vertNumb[i][j]];
					model.element[i].setVertNumb(j,mpn);
					
				}

			}
			
		
			}
			
			for(int i=1;i<=nNodes;i++){
				model.node[i].setCoord(coord1[i]);
				}
		
			model.scaleFactor=1000;
			
			for(int i=1;i<=nEls;i++)
			model.element[i].setRegion(regMap[regNumb[i]]);
			
			reRegionGroupEls(model);
			
			String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";
			
			model.writeMesh(fout);



	}
	
	
	public void getPostHexAtlas(){

		String s=util.getFile();
		
		String line;
		int max=1000000;
		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeMap=new int[max];

		int[] elMap=new int[max];
		int[] regNumb=new int[max];
		
		int[] nRegEls=new int[1000];


		int nElMax=0;

		int nnx=1;
		int nex=1;
		
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
			while(true){
				line=br.readLine();
				
				if(line==null) break;
				if(line.startsWith("GRID")){
					nnx=readAtlasNodes(br,nodeMap,coord1,nnx);


				}
				
				
				if(line.startsWith("CONC")){	
		
					nex=readAtlasElements(br,elMap,regNumb,nRegEls,vertNumb,nex);
				}

			}
		
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }

		
			int nNodes=0;
			for(int i=0;i<nodeMap.length;i++)
				if(nodeMap[i]>0) nNodes++;
			
			for(int i=1;i<elMap.length;i++)
				if(elMap[i]>nElMax) nElMax=elMap[i];


			int nEls=0;
		
			int nRegions=0;
			
			for(int i=0;i<nRegEls.length;i++)
				if(nRegEls[i]>0) {
					nRegions++;
					nEls+=nRegEls[i];
				}
			
		
			
			int[] regMap=new int[nRegEls.length];
			int nr=0;
			int[][]regEnd=new int[nRegions+1][2];
			regEnd[0][0]=1;
					
			for(int i=0;i<nRegEls.length;i++){
				
				if(nRegEls[i]>0) {
					nr++;
					regMap[i]=nr;
					regEnd[nr][0]=regEnd[nr-1][1]+1;
				
					regEnd[nr][1]=regEnd[nr][0]+nRegEls[i]-1;
				}
			}
			
			
			Model model=new Model(nRegions,nEls,nNodes,"hexahedron");
				
			for(int ir=1;ir<=nRegions;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+regMap[ir]);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				for(int j=0;j<8;j++){
					int mpn=nodeMap[vertNumb[i][j]];
					model.element[i].setVertNumb(j,mpn);
					
				}

			}
			
		
			}
			
			for(int i=1;i<=nNodes;i++){
				model.node[i].setCoord(coord1[i]);
				}
		
			model.scaleFactor=1;
			
			for(int i=1;i<=nEls;i++){
			model.element[i].setRegion(regMap[regNumb[i]]);
			}
			
			int[] elMapReorderd=reRegionGroupEls(model);
			

			
			int[] elMapReorderdRev=new int[elMapReorderd.length];
			for(int i=1;i<elMapReorderd.length;i++){
				elMapReorderdRev[elMapReorderd[i]]=i;
	
			}
	
			
			String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";

			model.writeMesh(fout);
			
			String map=System.getProperty("user.dir")+"\\EMSol\\map.txt";
			
			try {
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(map)));
				
				for(int i=1;i<elMap.length;i++)
					if(elMap[i]>0)
						pwBun.println(i+"\t"+elMapReorderdRev[elMap[i]]);
				
				pwBun.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	
			




	}
	

	
	public void getPostHexaNeu(int numElemNodes){


		String file=util.getFile();
		
		String line;
		int max=5000*1000;

		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeMap=new int[max];

		int[] elMap=new int[max];
		int[] regNumb=new int[max];
		
		int[] nRegEls=new int[10000];



		int nElMax=0;

		int nnx=1;
		int nex=1;
		
		try{
			File f=new File(file);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
	
			while(true){
	
				line=br.readLine();
			
				if(line==null) break;
		
				if(util.first(line).equals("403")){
		
					nnx=readNeuNodes(br,nodeMap,coord1,nnx);

				}

				
				if(util.first(line).equals("404")){

					nex=readNeuElements(br,elMap,regNumb,nRegEls,vertNumb,nex,numElemNodes);

			}

			}

				
	
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }


			int nNodes=nnx-1;

	
			
			for(int i=1;i<elMap.length;i++)
				if(elMap[i]>nElMax) nElMax=elMap[i];


			int nEls=0;
		
			int nRegions=0;
			
			for(int i=0;i<nRegEls.length;i++)
				if(nRegEls[i]>0) {
					nRegions++;
					nEls+=nRegEls[i];
				}
			
		
			
			int[] regMap=new int[nRegEls.length];
			int nr=0;
			int[][]regEnd=new int[nRegions+1][2];
			regEnd[0][0]=1;
					
			for(int i=0;i<nRegEls.length;i++){
		
				if(nRegEls[i]>0) {
				
					nr++;
					regMap[i]=nr;

					regEnd[nr][0]=regEnd[nr-1][1]+1;
				
					regEnd[nr][1]=regEnd[nr][0]+nRegEls[i]-1;
				}
			}

			int[] regMapRev=new int[regMap.length];
			for(int ir=1;ir<regMap.length;ir++)
				regMapRev[regMap[ir]]=ir;

			String type="hexahedron";
			if(numElemNodes==6) type="prism";
			else if(numElemNodes==4) type="tetrahedron";
			
			Model model=new Model();
			model.alloc(nRegions,nEls,nNodes,type);
			
			for(int ir=1;ir<=nRegions;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+regMapRev[ir]);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				for(int j=0;j<numElemNodes;j++){
					int mpn=nodeMap[vertNumb[i][j]];
					model.element[i].setVertNumb(j,mpn);
					
				}

			}
			
		
			}
			

			double unit=1;
			
			for(int i=1;i<=model.numberOfNodes;i++){

				model.node[i].setCoord(coord1[i].times(unit));
				}
			
			
		
			
			model.scaleFactor=1;
						
			for(int i=1;i<=nEls;i++){
			model.element[i].setRegion(regMap[regNumb[i]]);
			}
			
			int[] elMapReorderd=reRegionGroupEls(model);
			
			int[] elMapReorderdRev=new int[elMapReorderd.length];
			for(int i=1;i<elMapReorderd.length;i++){
				elMapReorderdRev[elMapReorderd[i]]=i;

			}
			
			
			int[] nodeMapRev=new int[nodeMap.length];
			
			for(int i=1;i<=model.numberOfElements;i++){
				int[] vn=model.element[i].getVertNumb();
				for(int j=0;j<numElemNodes;j++)
					nodeMapRev[nodeMap[i]]=vn[j];
		
			}

	
	
			
			String folder=new File(file).getParentFile().getPath();
			String fout=folder+"\\hexa.txt";

			model.writeMesh(fout);

			
			String map=folder+"\\elMap.txt";
			
			try {
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(map)));
				
				for(int i=1;i<elMap.length;i++)
					if(elMap[i]>0)
						pwBun.println(i+"\t"+elMapReorderdRev[elMap[i]]);
				
						
				pwBun.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	
			
			map=folder+"\\nodeMap.txt";
			

			
			try {
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(map)));
				
				for(int i=1;i<nodeMap.length;i++)
					if(nodeMap[i]>0)
						pwBun.println(i+"\t"+nodeMap[i]);
				
				pwBun.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	



	
	}
	

	public void get2DFormNeu(int numElemNodes){


		String file=util.getFile();
		
		String line;
		int max=1000000;

		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeMap=new int[max];

		int[] elMap=new int[max];
		int[] regNumb=new int[max];
		
		int[] nRegEls=new int[10000];



		int nElMax=0;

		int nnx=1;
		int nex=1;
		
		try{
			File f=new File(file);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
	
			while(true){
	
				line=br.readLine();
			
				if(line==null) break;
		
				if(util.first(line).equals("403")){
		
					nnx=readNeuNodes(br,nodeMap,coord1,nnx);

				}

				
				if(util.first(line).equals("404")){

					nex=readNeuElements(br,elMap,regNumb,nRegEls,vertNumb,nex,numElemNodes);
				
					

			}

			}

				
	
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }


			int nNodes=nnx-1;

	
			
			for(int i=1;i<elMap.length;i++)
				if(elMap[i]>nElMax) nElMax=elMap[i];


			int nEls=0;
		
			int nRegions=0;
			
			for(int i=0;i<nRegEls.length;i++)
				if(nRegEls[i]>0) {
					nRegions++;
					nEls+=nRegEls[i];
				}
			
		
			
			int[] regMap=new int[nRegEls.length];
			int nr=0;
			int[][]regEnd=new int[nRegions+1][2];
			regEnd[0][0]=1;
					
			for(int i=0;i<nRegEls.length;i++){
		
				if(nRegEls[i]>0) {
				
					nr++;
					regMap[i]=nr;

					regEnd[nr][0]=regEnd[nr-1][1]+1;
				
					regEnd[nr][1]=regEnd[nr][0]+nRegEls[i]-1;
				}
			}

			int[] regMapRev=new int[regMap.length];
			for(int ir=1;ir<regMap.length;ir++)
				regMapRev[regMap[ir]]=ir;

			String type="quadrangle";
			if(numElemNodes==3) type="triangle";
			
			Model model=new Model(nRegions,nEls,nNodes,type);
			
			for(int ir=1;ir<=nRegions;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+regMapRev[ir]);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				for(int j=0;j<numElemNodes;j++){
					int mpn=nodeMap[vertNumb[i][j]];
					model.element[i].setVertNumb(j,mpn);
					
				}

			}
			
		
			}
			

			double unit=1;
			
			for(int i=1;i<=model.numberOfNodes;i++){

				model.node[i].setCoord(coord1[i].times(unit));
				}
			
			
		
			
			model.scaleFactor=1;
						
			for(int i=1;i<=nEls;i++){
			model.element[i].setRegion(regMap[regNumb[i]]);
			}
			
			int[] elMapReorderd=reRegionGroupEls(model);
			
			int[] elMapReorderdRev=new int[elMapReorderd.length];
			for(int i=1;i<elMapReorderd.length;i++){
				elMapReorderdRev[elMapReorderd[i]]=i;

			}
			
			
			int[] nodeMapRev=new int[nodeMap.length];
			
			for(int i=1;i<=model.numberOfElements;i++){
				int[] vn=model.element[i].getVertNumb();
				for(int j=0;j<numElemNodes;j++)
					nodeMapRev[nodeMap[i]]=vn[j];
		
			}

	
	
			String folder=new File(file).getParentFile().getPath();
			String fout=folder+"\\quad.txt";

			model.writeMesh(fout);

			
			String map=folder+"\\elMap.txt";
			
			try {
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(map)));
				
				for(int i=1;i<elMap.length;i++)
					if(elMap[i]>0)
						pwBun.println(i+"\t"+elMapReorderdRev[elMap[i]]);
				
						
				pwBun.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	
			
			map=folder+"\\nodeMap.txt";
			

			
			try {
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(map)));
				
				for(int i=1;i<nodeMap.length;i++)
					if(nodeMap[i]>0)
						pwBun.println(i+"\t"+nodeMap[i]);
				
				pwBun.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	



	
	}
	

	
	
	public int readAtlasNodes(BufferedReader br,int[] nodeMap,Vect[] coord1,int nNode){
		String line;
		String[] sp;
		try{
		line=br.readLine();
		int nn;
		int nIdent=1;
	
		while(true){
			if(line==null) break;

			sp=line.split(regex);
		
			if(sp.length<3) return nNode;

			nn=Integer.parseInt(sp[nIdent]);

			if(nn==-1) break;
			
			nodeMap[nn]=nNode;
		

			coord1[nNode]=new Vect(Double.parseDouble(sp[nIdent+1]),Double.parseDouble(sp[nIdent+2]),Double.parseDouble(sp[nIdent+3]));
	
			line=br.readLine();
			
			nNode++;

		}
		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nNode;
		
	}
	
	public int readNeuNodes(BufferedReader br,int[] nodeMap,Vect[] coord1,int nNode){
		String line;
		String[] sp;
		try{
	
			int nmax1=4*coord1.length/5;
			int nmax2=coord1.length-nmax1;

		line=br.readLine();
		
		if(line==null) return nNode;

		sp=line.split(regex);
		
		int nn;
		int nIdent=0;
		
		while(!util.first(line).equals("-1")){
		
			nn=Integer.parseInt(sp[nIdent]);
			
			if(nn>nmax1) {
	
				nn=nmax1+nn%nmax2;

			}
			
			nodeMap[nn]=nNode;
		
			coord1[nNode]=new Vect(Double.parseDouble(sp[nIdent+11]),Double.parseDouble(sp[nIdent+12]),Double.parseDouble(sp[nIdent+13]));
	
			line=br.readLine();
			sp=line.split(regex);
			
			nNode++;


		}
		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nNode;
		
	}
	
	
	
public int readNeuElements(BufferedReader br,int[] elMap,int[]regNumb,int[] nRegEls,int[][] vertNumb,int nEl,int numElemNodes){
		
	
	int nmax1=4*elMap.length/5;
	int nmax2=elMap.length-nmax1;
	
	
		String line1="",line2;
		String[] sp1,sp2;
		try{

		int nIdent=0;

		
		int ne;
		line1=br.readLine();
		line2=br.readLine();
		
		while(!util.first(line1).equals("-1")){
		
			sp1=line1.split(regex);
			sp2=line2.split(regex);

			int[] vertNumb1=new int[8];
			 boolean def=false;
			boolean tetra=true;
			 int nz=0;

			for(int j=0;j<8;j++){

				vertNumb1[j]=Integer.parseInt(sp2[nIdent+j]);
			
				
				if(vertNumb1[j]==0){
					nz++;
					//vertNumb1[j]=vertNumb1[j-1];
	
			
				}
				

				
				if(vertNumb1[j]>nmax1){

					vertNumb1[j]=nmax1+vertNumb1[j]%nmax2;
				}
			
			}
			
			if(numElemNodes==8 && nz>0){
				def=true;
			}
		else 	if(numElemNodes==6 && nz!=2){
			
			def=true;
		}
else 	if(numElemNodes==4 && nz!=4){
			
			def=true;
		}
else 	if(numElemNodes==3 && nz!=5){
	
	def=true;
}
		
			
		if(!def){

			for(int j=0;j<8;j++)
				vertNumb[nEl][j]=vertNumb1[j];
			
			
				ne=Integer.parseInt(sp1[nIdent]);
				if(ne>nmax1) ne=nmax1+ne%nmax2;
				
				
				elMap[ne]=nEl;
							
				regNumb[nEl]=Integer.parseInt(sp1[nIdent+2]);
			

				nRegEls[regNumb[nEl]]++;

			
			// change order
			
			int[] tmp=new int[8];
			for(int j=0;j<8;j++){
				tmp[j]=vertNumb[nEl][j];
		
			}
			if(nz==0){
			vertNumb[nEl][0]=tmp[4];
			vertNumb[nEl][1]=tmp[7];
			vertNumb[nEl][2]=tmp[6];
			vertNumb[nEl][3]=tmp[5];
			
/*			vertNumb[nEl][4]=tmp[4];
			vertNumb[nEl][5]=tmp[7];
			vertNumb[nEl][6]=tmp[6];
			vertNumb[nEl][7]=tmp[5];*/
			vertNumb[nEl][4]=tmp[0];
			vertNumb[nEl][5]=tmp[3];
			vertNumb[nEl][6]=tmp[2];
			vertNumb[nEl][7]=tmp[1];
			}
			else if(nz==2){
				vertNumb[nEl][0]=tmp[4];
				vertNumb[nEl][1]=tmp[5];
				vertNumb[nEl][2]=tmp[6];	
				vertNumb[nEl][3]=tmp[0];
				vertNumb[nEl][4]=tmp[1];
				vertNumb[nEl][5]=tmp[2];
		
			
			}else if(nz==4 ){
				if(tmp[4]>0){
				vertNumb[nEl][0]=tmp[0];
				vertNumb[nEl][1]=tmp[1];
				vertNumb[nEl][2]=tmp[2];
				vertNumb[nEl][3]=tmp[4];	
				}else{

					vertNumb[nEl][0]=tmp[0];
					vertNumb[nEl][1]=tmp[1];
					vertNumb[nEl][2]=tmp[2];
					vertNumb[nEl][3]=tmp[3];	
				}
			}
			else if(nz==5 ){

					vertNumb[nEl][0]=tmp[0];
					vertNumb[nEl][1]=tmp[1];
					vertNumb[nEl][2]=tmp[2];
				
			}
			
			for(int j=0;j<8;j++){
				if(vertNumb[nEl][j]==0) {
				//	vertNumb[nEl][j]=vertNumb[nEl][j-1];
				//	if(j>0 && j!=4) vertNumb[nEl][j]=vertNumb[nEl][j-1];
				//	else if(j==0) vertNumb[nEl][j]=vertNumb[nEl][1];
					//else if(j==4) vertNumb[nEl][j]=vertNumb[nEl][5];
				//	ix--;
				//	break;
				}
			}
			
			nEl++;

		
	}

		
			br.readLine();
			br.readLine();
			br.readLine();
			br.readLine();
			br.readLine();
			
			line1=br.readLine();
			line2=br.readLine();

	
		}


		}
		
	
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nEl;
	}
	
	public int readAtlasElements(BufferedReader br,int[] elMap,int[]regNumb,int[] nRegEls,int[][] vertNumb,int nEl){
		
		String line1,line2;
		String[] sp1,sp2;
		try{
	//	line=br.readLine();
		int nIdent=1;

		int ne;
		while(true){
			
			line1=br.readLine();
		
			sp1=line1.split(regex);
			if(sp1.length<7) {return nEl;}
			

			line2=br.readLine();
			if(line2==null) return nEl;
			
			sp2=line2.split(regex);


			ne=Integer.parseInt(sp1[nIdent]);
			elMap[ne]=nEl;
			
			regNumb[nEl]=Integer.parseInt(sp1[nIdent+2]);


			nRegEls[regNumb[nEl]]++;

			for(int j=0;j<5;j++){

				vertNumb[nEl][j]=Integer.parseInt(sp1[nIdent+j+3]);
			}
			
	
		
			for(int j=0;j<3;j++){

				vertNumb[nEl][j+5]=Integer.parseInt(sp2[nIdent+j]);
			}
			
/*
			for(int j=0;j<8;j++)
			{



					numbAddedNodes++;
			
					 vernumb[ix][j]=nnMax+numbAddedNodes;//vernumb[ix][j-1];
						if(j>=0){
					 coord1[numbAddedNodes]=coord1[vernumb[ix][j-1]].deepCopy().times(0);
						}
					//ix--;
					//break;
				}
			}*/
	
			
			// change order
			
			int[] tmp=new int[8];
			for(int j=0;j<8;j++){
				tmp[j]=vertNumb[nEl][j];
			
			}
			vertNumb[nEl][0]=tmp[0];
			vertNumb[nEl][1]=tmp[3];
			vertNumb[nEl][2]=tmp[2];
			vertNumb[nEl][3]=tmp[1];
			
			vertNumb[nEl][4]=tmp[4];
			vertNumb[nEl][5]=tmp[7];
			vertNumb[nEl][6]=tmp[6];
			vertNumb[nEl][7]=tmp[5];
			
			nEl++;

		}

		}
		
	
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nEl;
	}
	
	public int readAtlasNodesOrig(BufferedReader br,int[] nodeNumb,Vect[] coord1, int nNode){
		String line;
		String[] sp;
		try{
		line=br.readLine();
		int nn;
		int nIdent=1;
	
		while(true){
			if(line==null) return nNode;

			sp=line.split(regex);

			

			nn=Integer.parseInt(sp[nIdent]);

			if(nn==-1) return nNode;
			
			nNode++;
			
			nodeNumb[nNode]=nn;

			coord1[nn]=new Vect(Double.parseDouble(sp[nIdent+1]),Double.parseDouble(sp[nIdent+2]),Double.parseDouble(sp[nIdent+3]));
	
			line=br.readLine();
			
		
		}
		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		return nNode;
		
	}
	
	public int readAtlasElementsOrig(BufferedReader br,int[] elNumb,int[]regNumb,int[][] vertNumb, int nEl){
		
		String line1,line2;
		String[] sp1,sp2;
		try{
	//	line=br.readLine();
		int nIdent=1;

		int ne;
		while(true){
			
			line1=br.readLine();
		
			sp1=line1.split(regex);
			if(sp1.length<7) { return nEl;}

			line2=br.readLine();
			if(line2==null) return nEl;
			
			sp2=line2.split(regex);


			ne=Integer.parseInt(sp1[nIdent]);
			
			nEl++;
			
			elNumb[nEl]=ne;
			

			
			regNumb[ne]=Integer.parseInt(sp1[nIdent+2]);


			for(int j=0;j<5;j++){

				vertNumb[ne][j]=Integer.parseInt(sp1[nIdent+j+3]);
			}
			
	
		
			for(int j=0;j<3;j++){

				vertNumb[ne][j+5]=Integer.parseInt(sp2[nIdent+j]);
			}
			
/*
			for(int j=0;j<8;j++)
			{



					numbAddedNodes++;
			
					 vernumb[ix][j]=nnMax+numbAddedNodes;//vernumb[ix][j-1];
						if(j>=0){
					 coord1[numbAddedNodes]=coord1[vernumb[ix][j-1]].deepCopy().times(0);
						}
					//ix--;
					//break;
				}
			}*/
	
			
			// change order
			
			int[] tmp=new int[8];
			for(int j=0;j<8;j++){
				tmp[j]=vertNumb[ne][j];
			
			}
			vertNumb[ne][0]=tmp[0];
			vertNumb[ne][1]=tmp[3];
			vertNumb[ne][2]=tmp[2];
			vertNumb[ne][3]=tmp[1];
			
			vertNumb[ne][4]=tmp[4];
			vertNumb[ne][5]=tmp[7];
			vertNumb[ne][6]=tmp[6];
			vertNumb[ne][7]=tmp[5];
				
		}

		}
		
	
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
		return nEl;
	}
	
	
	public void getEMSolFluxOld(String bbf,int dim, int numb,String type){

		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
		
			Mat BB=new Mat(1000*1000,dim);
			
			for(int tt=0;tt<numb;tt++){

				line="";
				while(!line.startsWith("STEP")){
					line=br.readLine();
					
					}
				int[] nx=new int[dim];;

				String fout=System.getProperty("user.dir")+"\\EMSol\\flux"+tt+".txt";
				
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));	
				
	while(!line.startsWith(type)){
	line=br.readLine();
			}
	
	for(int k=0;k<6;k++)
		line=br.readLine();
	
	int k=0;

	
			for(int i=1;i<1800000;i++){
				
	
		if(line.startsWith("-1")){
			k++;
			for(int j=0;j<8;j++)
				line=br.readLine();

		}
		else

		sp=line.split(regex);
			

		double Bu=Double.parseDouble(sp[1]);
		BB.el[nx[k]][k]=Bu;

				line=br.readLine();
	
				nx[k]++;
				
			

				if(k==dim-1 && nx[k]==nx[k-1]){
					break;
				}
				
				

	}



	int Ne=nx[0];
	pwBun.println("flux");
	pwBun.println(dim);
	pwBun.println(Ne);
	
	for(int j=0;j<Ne;j++){
		for(int p=0;p<dim;p++)
		pwBun.print(BB.el[j][p]+"\t");
		
		pwBun.println();
	}
	
	System.out.println("Flux was written to "+fout);
	pwBun.close();
	
			}
		
	br.close();
	fr.close();

	
	
		}


		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		
		

	
	
	
	}
	
	
	
	
	
	public void getNeuMeshTri(){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
			
			while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){
				
			}
		
			for(int i=1;i<10000;i++)
			{

				line=br.readLine();

				sp=line.split(regex);

				if(sp.length==15 && !sp[0].equals("0")) break;
			
				
			}
			


			int nnMax=0,nn=1;
			Vect[] coord1=new Vect[1000000];
			int[] map=new int[1000000];
			int nx=0;
			for(int i=1;i<1000000;i++)
			{
				sp=line.split(regex);
		
				line=br.readLine();

				if(sp.length!=15) break;

			
				nn=Integer.parseInt(sp[0]);
				nx++;
			
				map[nn]=nx;
		

				coord1[nx]=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
				

			}
			
			nnMax=nx;


			int[][] vernumb=new int[10*nnMax+1][3];
			int[] nReg=new int[1000000+1];
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		
				line=br.readLine();
			//	line=br.readLine();
			

			sp=new String[12];
			for(int i=1;i<=vernumb.length;i++)
			{
				line=br.readLine();

				sp=line.split(regex);
				if(sp.length<5) break;
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
				line=br.readLine();
		
				sp=line.split(regex);
				
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				vernumb[ix][0]=map[Integer.parseInt(sp[0])];
				vernumb[ix][1]=map[Integer.parseInt(sp[1])];
				vernumb[ix][2]=map[Integer.parseInt(sp[2])];

				

				for(int j=0;j<3;j++){
					if(vernumb[ix][j]==0) {
						ix--;
						break;
					}
				}

			}
			
			int nEl=ix;


			
			List<Integer> list1=new ArrayList<Integer>();
			for(int i=1;i<nReg.length;i++){
				if(nReg[i]!=-1)
					list1.add(nReg[i]);
			}
		
				Set<Integer> set = new HashSet<Integer>(list1);
				
				ArrayList<Integer> regNumbs = new ArrayList<Integer>(set);
				
				int nRegions=regNumbs.size();
				
				
				int[] regNumber=new int[nRegions+1];
				for(int ir=1;ir<=nRegions;ir++)
					regNumber[ir]=regNumbs.get(ir-1);
				

					
			int[] elOrd=new int[nEl+1];
			int[][]regEnd=new int[nRegions+1][2];
			
			 nx=0;
			for(int ir=1;ir<=nRegions;ir++)
			{
				regEnd[ir][0]=nx+1;
				for(int i=1;i<=nEl;i++)
					if(nReg[i]==regNumber[ir])
						elOrd[++nx]=i;
				regEnd[ir][1]=nx;
				
			}


			int nNodes=nnMax;
			double scaleFactor=1000.0;
			DecimalFormat formatter;
			if(scaleFactor==1)
				formatter= new DecimalFormat("0.000000000");
			else 
				formatter= new DecimalFormat("0.000000");
			

			try{


				String fout=System.getProperty("user.dir")+"\\\\EMSol\\tri.txt";
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

				pwBun.println("triangle");
				pwBun.println("//Number_of_Node");
				pwBun.println(nNodes);

				pwBun.println("//Number_of_Element");
				pwBun.println(nEl);

				pwBun.println("//Number_of_Region");
				pwBun.println(nRegions);
				pwBun.println("//Factor");
				pwBun.println(scaleFactor);

				for(int ir=1;ir<=nRegions;ir++)
					for(int i=regEnd[ir][0];i<=regEnd[ir][1];i++){
						for(int j=0;j<3;j++){
							int mpn=vernumb[elOrd[i]][j];
							pwBun.print(mpn+",");
						}
					
						pwBun.println();
					}

			Vect v;
			for(int i=1;i<=nnMax;i++){ 
				if(coord1[i]==null)
					v=new Vect(0,0);
				else
					v=coord1[i].deepCopy();
						for(int j=0;j<2;j++){
							pwBun.print(formatter.format(v.el[j]*scaleFactor)+" ,");
						}
					pwBun.println();	

				}

				for(int ir=1;ir<=nRegions;ir++)
					pwBun.println(regEnd[ir][0]+","+regEnd[ir][1]+","+"region"+ir);
				
				pwBun.close();
				br.close();
				fr.close();
				
				System.out.println();
				System.out.println(" Bun data was written to:");
				System.out.println("    "+fout);
			}
			catch(IOException e){ System.err.println("error");}

		
		
		}
		
		

		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	


public void getFluxAtlas(int dim, int nElMax){

	String regex="[ ,\\t]+";
	String s=util.getFile();
	try{
		
		String fout=System.getProperty("user.dir")+"\\EMSolFlux.txt";
		PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));		

		File f=new File(s);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		String line="";
		String[] sp=new String[15];
		
		
		Mat BB=new Mat(nElMax,dim);
		
		int[] nx=new int[dim];;
		
		
while(!line.startsWith("STEP")){
	line=br.readLine();
	pwBun.println(line);
		}
for(int k=0;k<6;k++){
	line=br.readLine();
	pwBun.println(line);
}

int k=0;


		for(int i=1;i<100000;i++){
			

/*		if(line.startsWith("-1")){
		k++;
		for(int j=0;j<8;j++){
			line=br.readLine();
			pwBun.println(line);
		}

	}
	else*/
	
	sp=line.split(regex);
		

	double Bu=Double.parseDouble(sp[1]);
	String line2=sp[0]+",\t"+sp[1]+",";
//	BB.el[nx[k]][k]=Bu;*/
	pwBun.println(line);
			line=br.readLine();
		
		//	util.pr(nx[0]+" - "+Bu);

		//	if(nx[k]<5) BB.el[nx[k]][k]=nx[k];
			
		//	nx[k]++;
			
		
/*				if(line.startsWith("-1")){
				pwBun.println(line);
				
				pwBun.println("  -1");
				break;
			}
			*/
			if(line==null || line.length()==0) break;
			

}


	
br.close();
fr.close();
pwBun.close();

System.out.println("Flux was written to "+fout);

	}


	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	

}

public int[] reRegionGroupEls(Model model){
	
	MeshFactory mf=new MeshFactory();
	
	return mf.reRegionGroupEls(model);
}


public void convertTetraNeu(){
	
	int nElems=0;

	int nNodes=0;
	int nnMax=0;
	String file=util.getFile();
	String line="";


	int[][] hexVert=new int[1000000][8];
	int[] regno=new int[1000000];
	
	int[] NodeIndex=new int[1000000];
	//int[] reverseIndex=new int[1000000];
	boolean[] nodeInuse=new boolean[1000000];
	Vect[] coords=new Vect[1000000];

	String[] sp;
	try{
		File f=new File(file);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		while(!line.startsWith("   403")){

			line=br.readLine();
			if(line==null) break;

		}
		

		
		while(true){
			line=br.readLine();

			sp=line.split(regex);


			if(!util.first(line).equals("-1")) {
				int ib=0;
				if(sp[0].equals("")) ib=1;
				int nn=Integer.parseInt(sp[ib]);
				coords[nn]=new Vect(Double.parseDouble(sp[ib+11]),Double.parseDouble(sp[ib+12]),Double.parseDouble(sp[ib+13]));
				if(nn>nnMax) nnMax=nn;
				
			} 
			else break;
		}
		while(!line.startsWith("   404")){
			line=br.readLine();
			if(line==null) break;

		}

		String[] lines=new String[7];
		String[] sp1;
		while(true){

			for(int i=0;i<7;i++){
				lines[i]=br.readLine();

			}

			if(util.first(lines[0]).equals("-1")) {
				break;
			}

			sp1=lines[0].split(regex);

			int ib=0;
			if(sp1[0].equals("")) ib=1;
			int elNum=Integer.parseInt(sp1[ib]);
		
			int regNo=Integer.parseInt(sp1[ib+2]);

			sp=lines[1].split(regex);

			ib=0;
			if(sp[0].equals("")) ib=1;
			for(int i=0;i<8;i++){
				hexVert[nElems][i]=Integer.parseInt(sp[i+ib]);
			if(hexVert[nElems][i]>0)
					nodeInuse[hexVert[nElems][i]]=true;
			}
			regno[nElems]=regNo;
			nElems++;					

			//if(nElems>0) break;

		}

		br.close();
		fr.close();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }

/*	int[] regno2=Arrays.copyOf(regno, regno.length);
	Arrays.sort(regno2);
	*/
	//int[] regMap=new int[regno.length];

	int nRegions=1;
/*	for(int i=regno2.length-2;i>0 && regno2[i]>0;i--){
		if(regno2[i]!=regno2[i+1]) nRegions++;
	}*/
	
	int[] map=new int[nodeInuse.length];
	int[] revMap=new int[nodeInuse.length];
	 nNodes=0;
	for(int i=0;i<nodeInuse.length;i++){
		if(nodeInuse[i]) {
			map[nNodes+1]=i;
			revMap[i]=nNodes+1;
			nNodes++;
			}
	}
	
	util.pr(nElems);
	util.pr(nNodes);
	
	Model model=new Model(nRegions,nElems,nNodes,"tetrahedron");

	model.region[1].setFirstEl(1);
	model.region[1].setLastEl(nElems);

	for(int i=1;i<=nNodes;i++){
		model.node[i].setCoord(coords[map[i]]);
		}
	
	
	for(int i=0;i<nElems;i++){
		int[] vn=new int[4];
		vn[0]=revMap[hexVert[i][0]];
		vn[1]=revMap[hexVert[i][2]];
		vn[2]=revMap[hexVert[i][1]];
		vn[3]=revMap[hexVert[i][4]];
		model.element[i+1].setVertNumb(vn);
		model.element[i+1].setRegion(1); //----
	}
	

	MeshFactory mf=new MeshFactory();
	
	mf.reRegionGroupEls(model);
	//model=mf.dropUnusedNodes(model);
	
	model.scaleFactor=1;
	String folder=new File(file).getParentFile().getPath();
	String tetFile=folder+"\\tet_from_neu.txt";
	model.writeMesh(tetFile);


}

public void convertToNeu(){
	String file=util.getFile();
	Model model=new Model(file);
	String folder=new File(file).getParentFile().getPath();
	String fileNeu=folder+"\\bun.neu";
	model.writer.writeMeshNeu(model, fileNeu);
}
	



public void getFromIr3(int numElemNodes){


	String file=util.getFile();
	
	String line;

	
	int nRegions=1000;
	int nNodes=0;
	int nVolElems=0;
	int nQuads=0;
	String type="hexahedron";
	if(numElemNodes==6) type="prism";
	else if(numElemNodes==4) type="tetrahedron";
	else if(numElemNodes==3) type="triangle";

	
	try{
		File f=new File(file);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);

		line=br.readLine();
		int nx=Integer.parseInt(line);
		
		line=br.readLine();
		util.pr(line);
		String[] sp=line.split(regex);
		int kx=0;
		nNodes=Integer.parseInt(sp[kx++]);
	
		nQuads=Integer.parseInt(sp[kx++]);

		if(numElemNodes>4)
		nVolElems=Integer.parseInt(sp[kx++]);
		if(numElemNodes==3 || numElemNodes==4 )
			nVolElems=nQuads;

			Model model1=new Model();
			model1.alloc(nRegions,nVolElems,nNodes,type);
			
			model1.region[1].setFirstEl(1);
			model1.region[1].setLastEl(nVolElems);
			for(int ir=2;ir<=model1.numberOfRegions;ir++){
				model1.region[ir].setFirstEl(1);
				model1.region[ir].setLastEl(0);
			}

			double unit=1;
			
			int dim=model1.dim;
			for(int i=1;i<=model1.numberOfNodes;i++){
	
					line=br.readLine();
					sp=line.split(regex);
					int nn=Integer.parseInt(sp[0]);
					Vect v=new Vect(dim);
					v.el[0]=Double.parseDouble(sp[1]);
					v.el[1]=Double.parseDouble(sp[2]);
					if(dim>2)
					v.el[2]=Double.parseDouble(sp[3]);
				model1.node[i].setCoord(v.times(unit));
				}
			
			if(numElemNodes>4 )
			for(int i=0;i<nQuads;i++){
				line=br.readLine();
			}
			
			for(int i=1;i<=model1.numberOfElements;i++){
				
				line=br.readLine();

				sp=line.split(regex);
				int ne=Integer.parseInt(sp[0]);
				int nReg=Integer.parseInt(sp[1]);
				int nVerts=Integer.parseInt(sp[2]);
				int[] nv1=new int[nVerts];
				for(int j=0;j<nVerts;j++){
					nv1[j]=Integer.parseInt(sp[3+j]);
				}
				int[] nv=new int[nVerts];	
				if(nVerts==8){
				nv[0]=nv1[0];
				nv[1]=nv1[3];
				nv[2]=nv1[2];
				nv[3]=nv1[1];
				nv[4]=nv1[4];
				nv[5]=nv1[7];
				nv[6]=nv1[6];
				nv[7]=nv1[5];
				}else{
					for(int j=0;j<nVerts;j++)
						nv[j]=nv1[j];	
				}

				model1.element[i].setVertNumb(nv);
				model1.element[i].setRegion(nReg);

			}
			
			reRegionGroupEls(model1);
			
			int nrx=0;
			
			int regMap[]=new int[model1.numberOfRegions+1];
	
			for(int ir=1;ir<=model1.numberOfRegions;ir++){
		
				if(model1.region[ir].getNumbElements()>0){
					
					regMap[ir]=++nrx;
				}
			}
			
util.pr(nrx);
			
			int nRegions2=nrx;
			
			Model model=new Model();
			model.alloc(nRegions2,nVolElems,nNodes,type);
			
			for(int i=1;i<=model1.numberOfNodes;i++){

			model.node[i].setCoord(model1.node[i].getCoord());
			}
			
			for(int i=1;i<=model.numberOfElements;i++){

			model.element[i].setVertNumb(model1.element[i].getVertNumb());
			int ir1=model1.element[i].getRegion();
			int ir=regMap[ir1];
			model.element[i].setRegion(ir);
			model.region[ir].setFirstEl(model1.region[ir1].getFirstEl());
			model.region[ir].setLastEl(model1.region[ir1].getLastEl());
			}

		
			reRegionGroupEls(model);
			
			String folder=new File(file).getParentFile().getPath();
			String fout=folder+"\\hexa.txt";

			model.writeMesh(fout);
	

	}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }

		






}

}
