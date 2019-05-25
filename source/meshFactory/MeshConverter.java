package meshFactory;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
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
import math.Mat;
import math.SpBlockMat;
import math.SpMat;
import math.Vect;
import math.util;



public class MeshConverter {

	String regex="[ : ,=\\t]+";
	BoundarySet bs=new BoundarySet();
	Writer writer=new Writer();

	public static void main(String[] args){
		MeshConverter mf=new MeshConverter();
	//	mf.hexaToTetraNeu();
		
	//mf.quadToTriang();
	//mf.triToQuad();

	//	mf.prismToHexa();
	
	//mf.hexaToTetra();
		
		//mf.hexToPrism();
	//	mf.addMidNodesNeu();
//	mf.refine2ndTo1stTetra();
//	mf.refineTetraTo3Tetra();
	//	mf.removeMidNodesNeu();
		
	//	mf.cutCorners();
		
	//mf.scaleNeu();
		
	//	mf.reRegionNeu();
	}

	


	public void hexToPrism()
	{
		hexToPrism(2);
	}
	
	public void hexToPrism(int dir)
	
	{
		String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model();
	model.loadMesh(bun);

	String prismMesh = System.getProperty("user.dir") + "//prismEl.txt";

	if(model.elCode==4)
		model.writeMeshHexaToPrism(prismMesh,dir);

	}
	
	public void hexaToPyramid()
	
	{
		String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model();
	model.loadMesh(bun);

	String pyrMesh = System.getProperty("user.dir") + "//pyrEl.txt";

	if(model.elCode==4)
		writer.writeMeshHexaToPyramid(model,pyrMesh);
	}
	
	public void prismToHexa(){
		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		String folder=new File(bun).getParentFile().getPath();
		String pyrMesh = folder + "//hexa.txt";

		if(model.elCode==3)
			writer.writeMeshPrismToHexa(model,pyrMesh);
		}
		
public void hexaToTetra()
	
	{
		String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model();
	model.loadMesh(bun);

	String folder=new File(bun).getParentFile().getPath();
	
	String pyrMesh = folder + "//tetEl.txt";

	if(model.elCode==4)
		writer.writeMeshHexaToTetra(model,pyrMesh);
	}
	
	public void reRegionGroupEls(Model model){
	
	MeshFactory mf=new MeshFactory();
	
	mf.reRegionGroupEls(model);
}

	public void hexaToTetraNeu(){


		//int[][] tetNodes={{1,3,4,0},{1,4,6,5},{1,6,3,2},{4,3,6,7},{1,4,3,6}};
		//int[][] tetNodes2={{0,2,7,3},{5,0,7,4},{5,2,0,1},{5,7,2,6},{0,7,2,5}};
		int[][] tetNodes={{5,7,0,4},{5,0,2,1},{5,2,7,6},{0,7,2,3},{5,0,7,2}};
		int[][] tetNodes2={{4,6,3,7},{1,4,3,0},{1,6,4,5},{1,3,6,2},{4,3,6,1}};
		int[][] tetNodes3={{4,6,3,7},{1,4,3,0},{1,6,4,5},{1,3,6,2},{4,3,6,1}};

		int nSubEl=5;

		int ix=0;

		String file=util.getFile();
		String folder=new File(file).getParentFile().getPath();
		String tetFile=folder+"\\tet_from_hex.neu";

		String line="";


		String[] sp;

		boolean[] isApex=new boolean[1000000];

		//	for(int i=0;i<canBeApex.length;i++)
		//	canBeApex[i]=true;

		try{
			File f=new File(file);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);


			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(tetFile)));

			while(!line.startsWith("   404")){

				line=br.readLine();
				pw.println(line);
				if(line==null) break;

			}



			int[] hexVert=new int[8];
			int[][] tetVert=new int[5][8];
			String[] lines=new String[7];
			String[] sp1;
			while(true){

				for(int i=0;i<7;i++){
					lines[i]=br.readLine();

				}

				if(util.first(lines[0]).equals("-1")) {
					pw.println(lines[0]);
					break;
				}

				sp1=lines[0].split(regex);

				int ib=0;
				if(sp1[0].equals("")) ib=1;
				int elNum=Integer.parseInt(sp1[ib]);
				if(!sp1[3].equals("25")) continue;
				//	sp1[1]=Integer.toString(elNum);
				sp1[3]="26";
				sp1[4]="6";
				String line1=sp1[ib+1];
				for(int i=2+ib;i<sp1.length;i++)
					line1=line1+","+sp1[ib+i];

				//	pw.println(lines[0]);

				sp=lines[1].split(regex);

				ib=0;
				if(sp[0].equals("")) ib=1;
				for(int i=0;i<8;i++)
					hexVert[i]=Integer.parseInt(sp[i+ib]);


				for(int ne=0;ne<nSubEl;ne++){
					for(int i=0;i<8;i++)
						tetVert[ne][i]=0;
					ix++;


					if(!isApex[hexVert[tetNodes[ne][0]]]  && !isApex[hexVert[tetNodes[ne][1]]] && !isApex[hexVert[tetNodes[ne][2]]]){
						for(int j=0;j<4;j++)
							tetVert[ne][j]=hexVert[tetNodes[ne][j]];
						if(ne<4)
							isApex[hexVert[tetNodes[ne][3]]]=true;
					}else if(!isApex[hexVert[tetNodes2[ne][0]]]  && !isApex[hexVert[tetNodes2[ne][1]]] && !isApex[hexVert[tetNodes2[ne][2]]]){
						for(int j=0;j<4;j++)
							tetVert[ne][j]=hexVert[tetNodes2[ne][j]];
						if(ne<4)
							isApex[hexVert[tetNodes2[ne][3]]]=true;
					}


					String line11=Integer.toString(ix)+","+line1;


					String line2=tetVert[ne][0]+","+tetVert[ne][1]+","+tetVert[ne][2]+",0,"+tetVert[ne][3]+",0,0,0,0,0";

					pw.println(line11);
					pw.println(line2);
					for(int i=2;i<7;i++)
						pw.println(lines[i]);


				}

			}

			pw.flush();

			pw.close();
			br.close();
			fr.close();
			
			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+tetFile);
			System.out.println();

		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}
	
	public void quadToTriang(){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		String folder=new File(bun).getParentFile().getPath();
		String triangMesh = folder + "//triangEl.txt";

		if(model.elCode==1)
			model.writeMeshQuadToTri(triangMesh);

	}

	
	public void triToTri(){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		String folder=new File(bun).getParentFile().getPath();
		String triangMesh = folder + "//triangEl.txt";

		if(model.elCode==0)
			model.writeMeshTriToTri(triangMesh,0,1e99);

	}

	public void triToQuad(){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		String folder=new File(bun).getParentFile().getPath();
		String triangMesh = folder + "//quad.txt";

		if(model.elCode==0)
			model.writeMeshTriToQuad(triangMesh);

	}


	
	public void refine2ndTo1stTetra(){



		//int[][] tetNodes={{8,12,10,0},{9,13,8,1},{10,14,9,2},{12,13,14,4}};
		int[][] tetNodes={{8,12,10,0},{9,13,8,1},{10,14,9,2},{12,13,14,4},{8,10,12,13},{10,9,14,13},{14,12,10,13},{8,9,10,13}};


		int nSubEl=8;

		int ix=0;

		String file=util.getFile();
		String folder=new File(file).getParentFile().getPath();
		String tetFile=folder+"\\2nd_refined_to_1st.neu";


		String line="";


		String[] sp;

		boolean[] isApex=new boolean[1000000];

		//	for(int i=0;i<canBeApex.length;i++)
		//	canBeApex[i]=true;

		try{
			File f=new File(file);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);


			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(tetFile)));

			while(!line.startsWith("   404")){

				line=br.readLine();
				pw.println(line);
				if(line==null) break;

			}


			int[] hexVert=new int[20];
			int[][] tetVert=new int[8][8];
			String[] lines=new String[7];
			String[] sp0,sp1;
			while(true){

				for(int i=0;i<7;i++){
					lines[i]=br.readLine();
				}

				if(util.first(lines[0]).equals("-1")) {
					pw.println(lines[0]);
					break;
				}

		
			//	int elNum=Integer.parseInt(sp0[ib]);
				//if(elNum>1) break;
				//if(!sp1[3].equals("25")) continue;
				//	sp1[1]=Integer.toString(elNum);


				//	pw.println(lines[0]);

				sp=lines[1].split(regex);

				int ib=0;
				if(sp[0].equals("")) ib=1;
				for(int i=0;i<10;i++)
					hexVert[i]=Integer.parseInt(sp[i+ib]);
				
				sp=lines[2].split(regex);

				ib=0;
				if(sp[0].equals("")) ib=1;
				for(int i=0;i<10;i++)
					hexVert[i+10]=Integer.parseInt(sp[i+ib]);


				for(int ne=0;ne<nSubEl;ne++){
					for(int i=0;i<8;i++)
						tetVert[ne][i]=0;
					ix++;

				for(int j=0;j<4;j++)
							tetVert[ne][j]=hexVert[tetNodes[ne][j]];
				
	
				sp0=lines[0].split(regex);

				 ib=0;
				if(sp0[0].equals("")) ib=1;
				//sp1[3]="26";
				//sp0[2]=Integer.toString(ix);
				sp0[4]="6";
				String line1=sp0[ib+1];
				for(int i=2+ib;i<sp0.length;i++)
					line1=line1+","+sp0[ib+i];

					String line11=Integer.toString(ix)+","+line1;

					String line2=tetVert[ne][0]+","+tetVert[ne][1]+","+tetVert[ne][2]+",0,"+tetVert[ne][3]+",0,0,0,0,0";
					String line3="0,0,0,0,0,0,0,0,0,0";

					pw.println(line11);
					pw.println(line2);
					pw.println(line3);
					
					for(int i=3;i<7;i++)
						pw.println(lines[i]);


				}

			}


			pw.close();
			br.close();
			fr.close();

		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	}


	

	private void addMidNodesNeu(){

		

		int[][] hexVert=new int[1000000][8];
		
		int[] NodeIndex=new int[1000000];
		Vect[] coords=new Vect[1000000];

		int nElems=0;

		int nNodes=0;
		int nnMax=0;
		String file=util.getFile();
		String folder=new File(file).getParentFile().getPath();
		String tetMidFile=folder+"\\midNodesAdded.neu";

		String line="";


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
					NodeIndex[nn]=nNodes;
					coords[nNodes]=new Vect(Double.parseDouble(sp[ib+11]),Double.parseDouble(sp[ib+12]),Double.parseDouble(sp[ib+13]));
					nNodes++;
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

				sp=lines[1].split(regex);

				ib=0;
				if(sp[0].equals("")) ib=1;
				for(int i=0;i<8;i++)
					hexVert[nElems][i]=Integer.parseInt(sp[i+ib]);

				nElems++;					


			}

			br.close();
			fr.close();

		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }




	//	byte[][] edgeLocalNodes={{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
		
		byte[][] edgeLocalNodes={{0,1},{1,2},{2,0},{0,4},{1,4},{2,4}}; // for tetra

		int[][][] nodeNode=new int[2][8*nNodes+1][5*6];
		byte[] nodeNodeIndex=new byte[8*nNodes+1];
		int nEdge=0;
		int n1,n2,m;
		int[][] edgeNodes=new int[20*nNodes+1][2];

		int nElEdges=edgeLocalNodes.length;
		
		int[][] elmEdges=new int[nElems][6];
		
		for(int i=0;i<nElems;i++){
			int[] vertNumb=hexVert[i];
			for(int j=0;j<nElEdges;j++){
				n1=vertNumb[edgeLocalNodes[j][0]];
				n2=vertNumb[edgeLocalNodes[j][1]];
				if(n1==0 || n2==0) continue;
				if(n2<n1) { int tmp=n1; n1=n2; n2=tmp;}
				m=util.search(nodeNode[0][n1],n2);

				if(m<0){

					nEdge++;
					edgeNodes[nEdge][0]=n1;
					edgeNodes[nEdge][1]=n2;
					nodeNode[0][n1][nodeNodeIndex[n1]]=n2;
					nodeNode[1][n1][nodeNodeIndex[n1]++]=nEdge;
					elmEdges[i][j]=nodeNode[1][n1][nodeNodeIndex[n1]-1];
				}
				else 
					for(int k=0;k<nodeNodeIndex[n1];k++)
					{
					if(nodeNode[0][n1][k]==n2){
						elmEdges[i][j]=nodeNode[1][n1][k];
					break;
					}
					}

			}

		}
		

		boolean[] edgeCounted=new boolean[nEdge+1];
		util.pr("numEdges: "+nEdge);
		
		util.pr(nnMax);
	
		//int[] map=new int[nnMax+1];
		
		int nnMax2=nnMax;
	
		int ix=0;
		try{
			util.pr(tetMidFile);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(tetMidFile)));
			File f=new File(file);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			while(!line.startsWith("   403")){

				line=br.readLine();
				pw.println(line);
				if(line==null) break;

			}
			
			String[] sp1=null;
			String line1;
			while(true){
				line=br.readLine();



				if(!util.first(line).equals("-1")) {
					ix++;
					if(ix==1)
					sp1=line.split(regex);
	
				} else{
					for(int i=0;i<nElems;i++){
						for(int j=0;j<6;j++)
						if(!edgeCounted[elmEdges[i][j]]){
							
						edgeCounted[elmEdges[i][j]]=true;
							
						nnMax2++;
						n1=edgeNodes[elmEdges[i][j]][0];
						n2=edgeNodes[elmEdges[i][j]][1];
						Vect mid=coords[NodeIndex[n1]].add(coords[NodeIndex[n2]]).times(.5);
						sp1[0]=Integer.toString(nnMax2);
						sp1[11]=Double.toString(mid.el[0]);
						sp1[12]=Double.toString(mid.el[1]);
						sp1[13]=Double.toString(mid.el[2]);
						
						int ib=0;
					//	if(sp1[0].equals("")) ib=1;	
						line1=sp1[ib];
						for(int k=ib+1;k<14;k++)
							line1=line1+","+sp1[k];
						
						pw.println(line1);
					}
			}
					
					pw.println(line);
					break;

				}
				pw.println(line);

			}
			
			while(!line.startsWith("   404")){

				line=br.readLine();
				pw.println(line);
				if(line==null) break;

			}

			String[] lines=new String[7];
			
			ix=0;
		
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
				//if(!sp1[3].equals("25")) continue;
				//	sp1[1]=Integer.toString(elNum);
				sp1[3]="26";
				sp1[4]="10";
				lines[0]=sp1[ib];
				for(int i=1+ib;i<sp1.length;i++)
					lines[0]=lines[0]+","+sp1[ib+i];


				sp=lines[1].split(regex);

				 ib=0;
				if(sp[0].equals("")) ib=1;
	//util.pr(ix);
				sp[ib+8]=Integer.toString(nnMax+elmEdges[ix][0]);
				sp[ib+9]=Integer.toString(nnMax+elmEdges[ix][1]);
				
				lines[1]=sp[0];
				for(int j=1;j<10;j++)
					lines[1]=	lines[1]+","+sp[j];
				
				sp=lines[2].split(regex);

				 ib=0;
				if(sp[0].equals("")) ib=1;
				
				//util.pr(elmEdges[ix][2]);
				sp[ib+0]=Integer.toString(nnMax+elmEdges[ix][2]);
				sp[ib+2]=Integer.toString(nnMax+elmEdges[ix][3]);
				sp[ib+3]=Integer.toString(nnMax+elmEdges[ix][4]);
				sp[ib+4]=Integer.toString(nnMax+elmEdges[ix][5]);
				
				ix++;
				lines[2]=sp[0];
				for(int j=1;j<10;j++)
					lines[2]=	lines[2]+","+sp[j];

				for(int i=0;i<7;i++)
					pw.println(lines[i]);		


			}
			pw.println("   -1");
			pw.close();
			br.close();
			fr.close();

		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }

	}
	
	private void removeMidNodesNeu(){


		String file=util.getFile();
	
		String folder=new File(file).getParentFile().getPath();
		String tetMidFile=folder+"\\2nd_to_1st.neu";

		String line="";


		String[] sp;



		int ix=0;
		try{
			
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(tetMidFile)));
			File f=new File(file);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			
			String[] sp1=null;
			String line1;
			
			while(!line.startsWith("   404")){

				line=br.readLine();
				pw.println(line);
				if(line==null) break;

			}

			String[] lines=new String[7];
			
			ix=0;
		
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
				int nReg=Integer.parseInt(sp1[ib+2]);
			//	if(nReg==1120){
				//if(!sp1[3].equals("25")) continue;
				//	sp1[1]=Integer.toString(elNum);
				sp1[3]="25";
				sp1[4]="6";
				lines[0]=sp1[ib];
				for(int i=1+ib;i<sp1.length;i++)
					lines[0]=lines[0]+","+sp1[ib+i];
				
				
				sp1=lines[1].split(regex);
				ib=0;
				if(sp1[0].equals("")) ib=1;
				for(int k=8;k<sp1.length;k++)
					sp1[ib+k]="0";
				
				lines[1]=sp1[ib];
				for(int i=1+ib;i<sp1.length;i++)
					lines[1]=lines[1]+","+sp1[ib+i];
				
				lines[2]="0,0,0,0,0,0,0,0,0,0,";
			//}

				for(int i=0;i<7;i++)
					pw.println(lines[i]);		


			}
			pw.println("   -1");
			pw.flush();
			pw.close();
			br.close();
			fr.close();
			
			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+tetMidFile);
			System.out.println();

		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }

	}
	
	

	
	


	
	public void scaleNeu(){

		String file=util.getFile();

		String folder=new File(file).getParentFile().getPath();
		String scaledFile=folder+"\\scaled.neu";

		String line="";


		String[] sp;


		try{
			
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(scaledFile)));
			File f=new File(file);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			while(!line.startsWith("   403")){

				line=br.readLine();

				pw.println(line);
				if(line==null) break;

			}

			
			while(!line.startsWith("   404")){

				line=br.readLine();
				sp=line.split(regex);
				
				boolean unused=false;
				if(sp.length>13){
					int ib=0;
					if(sp[0].equals("")) ib=1;
					double z=Double.parseDouble(sp[ib+13]);
					
					if(z<0) unused=true;
		
					if(z>.0 && z<.08){
						//util.pr(z);
						if(z>.0223 && z<.0224) {
							z=.03;
							util.pr("1----- "+z);
						}
						if(z>.025 && z<.027){
							z=.045;
							util.pr("2----- "+z);
						}
						if(z>.0285 && z<.03){
							z=.07;
							util.pr("3----- "+z);
						}
						if(z>.006 && z<.007){
							z=.019/4;
							util.pr("4----- "+z);
						}
						if(z>.012 && z<.013){
							z=.019/2;
							util.pr("4----- "+z);
						}
						if(z<0){
							//z*=.5;
							util.pr("5------------------------ "+z);
						}
						//z*=1.5*z/.019;
					sp[ib+13]=Double.toString(z);
					line=sp[ib];
					for(int i=1+ib;i<sp.length;i++)
						line=line+","+sp[ib+i];
					//util.pr(z);
					}
	}
				if(!unused)
				pw.println(line);
				if(line==null) break;

			}

			String[] lines=new String[7];
				
			while(true){

				for(int i=0;i<7;i++){
					lines[i]=br.readLine();

				}

				if(util.first(lines[0]).equals("-1")) {
					break;
				}


				for(int i=0;i<7;i++)
					pw.println(lines[i]);		


			}
			pw.println("   -1");
			pw.close();
			br.close();
			fr.close();

		}
		catch(Exception e){System.err.println("error");	e.printStackTrace(); }

	}
	

}
