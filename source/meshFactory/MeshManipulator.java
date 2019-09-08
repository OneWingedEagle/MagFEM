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
import math.Complex;

import math.Mat;

import math.SpBlockMat;
import math.SpMat;
import math.Vect;
import math.util;



public class MeshManipulator {

	String regex="[ : ,=\\t]+";
	BoundarySet bs=new BoundarySet();
	Writer writer=new Writer();

	public static void main(String[] args){

		MeshManipulator mf=new MeshManipulator();

		//Model model=new Model("D:\\JavaWorks\\FEM problems\\Hamed solver\\bun1elem.txt");
		//model.setEdge();
		//mf.reRegionb();
	//mf.connectivity(1e-6);	mf.dropUnusedNodes();
	//	mf.dropUnusedNodes();
	//mf.deform();
	//mf.rescale(1./500);
	//	mf.rotate(9*PI/180);
		
	//	mf.translate(new Vect(1,.0,1));
	//	Mat R=util.rotMat(new Vect(0,0,1), new Vect(1,1.4,0));
	//	Mat R=util.rotEuler(new Vect(0,0,1), 45*PI/180);
	//	mf.rotate(R);
	//	mf.pileRotate(10, 1*PI/180);

	//mf.reverseHexa();

	//Model model=mf.rotExtendNfold(29);
	//model.writeMesh(new File(model.meshFilePath).getParentFile().getPath()+"\\extenRot.txt");
		int Nr=10;
		int[] regs0=new  int[Nr];
		for(int i=1;i<=Nr;i++){
			regs0[i-1]=i;
		}
		//mf.revolveLine(new Vect().linspace(1, 10, Nr+1), regs0, 45, PI/2/45);
		int[] regs={1,2};
		mf.extractReg(regs);mf.dropUnusedNodes();
		String stat="D:\\JavaWorks\\FEM problems\\Solver Gaol\\noUnusedNodes.txt";
		 String rot="D:\\JavaWorks\\FEM problems\\Solver Gaol\\bunTranslated.txt";
		//mf.assemble(rot,stat);

		//mf.assemble();
	//mf.rescale(new Vect(1,1,0.54641016151377545870548926830117));
	//	mf.extractReg(-0.001,6.2,0,2*PI,0,.19e-3);
		//mf.copy(new Vect(0,0,2), 9);
		/*String bun=util.getFile();
		Model model1=new Model(bun);
		for(int j=1;j<10;j++){
			Model model2=new Model(bun);
		mf.translate(model2,new Vect(0,0,j),false);
		model1=mf.assemble(model2, model1, false);
		}
		String folder=new File(bun).getParentFile().getPath();
		model1.writeMesh(folder+"\\assembled.txt");*/

	//String outMesh = folder + "//quads2.txt";
	//	model.writeMeshQuadToTri(outMesh);
		//model.writeMeshTriToQuad(outMesh);
		/*		int[] map=new int[1+model.numberOfNodes];
		
		for(int i=1;i<=model.numberOfNodes;i++){
			Vect v=model.node[i].getCoord();
			double r=v.norm();
			double tt=util.getAng(v);
			int ir=(int )Math.round((r-.999)/.2);
			int it=(int )Math.round((tt-.5)/.175);
			
			util.pr(r+" "+tt);
		//	util.pr(ir+" "+it);
		//	map[i]=

		}*/
	//	model.writeMesh("D:\\JavaWorks\\FEM problems\\ipm_motor2D\\compacted.txt");*/
		//	mf.rotate(PI/2);
	//	mf.pileUpHexa(40, .05);
		
	//	mf.pileHelic(6*8, PI/4, .0125*18);
		
		//mf.pileRotate(10, PI/18);



	}
	


	public void setSliceBounds(Model model)
	{
		bs.setSliceBounds(model);
	}

	public void rotExtendNfoldW( int ns){
		String bun=util.getFile();
		if(bun==null || bun.equals("") )  throw new NullPointerException("file not found.");
		rotExtendNfold(bun,ns);
	}

	public void rotExtendNfold(String bun, int ns){

		Model model=rotExtendNfold(new Model(bun),ns);

		model.writeMesh(System.getProperty("user.dir") + "//extended.txt");

	}



	public void extendFlip(int hv){
		String bun=util.getFile();
		if(bun==null || bun.equals("") )  throw new NullPointerException("file not found.");
		Model model=new Model(bun);
		extendFlip(model,hv);

	}

	public void extendFlip(Model md0,int hv){

		boolean[] nodeIn=new boolean[md0.numberOfNodes+1];
		for(int i=1;i<=md0.numberOfElements;i++){
			int[] vertNum=md0.element[i].getVertNumb();
			for(int j=0;j<md0.nElVert;j++)
				nodeIn[vertNum[j]]=true;
		}



		int nRegions=md0.numberOfRegions;

		int nElements=2*md0.numberOfElements;

		int ix=0;
		boolean[] onAx=new boolean[md0.numberOfNodes+1];
		for(int i=1;i<=md0.numberOfNodes;i++)
		{
			if(!nodeIn[i]) continue;
			Vect z=md0.node[i].getCoord();
			if(Math.abs(z.el[hv])<1e-10) onAx[i]=true;

		}

		int nNodes=2*md0.numberOfNodes;


		Model md1=new Model(nRegions,nElements,nNodes,md0.elType);

		int n1=1,N;
		for(int i=1;i<=md0.numberOfRegions;i++){


			if(i>1)
				n1=md1.region[i-1].getLastEl()+1;
			N=md0.region[i].getNumbElements();
			md1.region[i].setFirstEl(n1);
			md1.region[i].setLastEl(n1+2*N-1);
			md1.region[i].setName(md0.region[i].getName());
			md1.region[i].setMaterial(md0.region[i].getMaterial());
			md1.region[i].setMur(md0.region[i].getMur());
			md1.region[i].setM(md0.region[i].getM());
			md1.region[i].setJ(md0.region[i].getJ());
			md1.region[i].setSigma(md0.region[i].getSigma());
			md1.region[i].setPois(md0.region[i].getPois());
			md1.region[i].setYng(md0.region[i].getYng());

		}

	/*	md1.BCtype=new int[4];
		md1.diricB=new Vect[4];
		for(int j=0;j<4;j++)
			md1.diricB[j]=new Vect(2);
*/


		for(int i=1;i<=md0.numberOfNodes;i++)	
			md1.node[i].setCoord(md0.node[i].getCoord());

		int[] map=new int[md0.numberOfNodes+1];
		for(int i=1;i<=md0.numberOfNodes;i++)	
		{

			if(onAx[i]) continue;

			Vect v=md0.node[i].getCoord();

			v.el[hv]*=-1;
			md1.node[i+md0.numberOfNodes].setCoord(v);
			map[i]=i+md0.numberOfNodes;
		}
		ix=0;

		for(int ir=1;ir<=md0.numberOfRegions;ir++)		{

			int ne=md0.region[ir].getNumbElements();
			ix=0;
			for(int i=md0.region[ir].getFirstEl();i<=md0.region[ir].getLastEl();i++)	{

				ix++;
				int p=1;	
				p=md1.region[ir].getFirstEl();
				int[] vertNumb=md0.element[i].getVertNumb();

				for(int j=0;j<2;j++){
					int ie=p-1+j*ne+ix;
					if(j==0)
						md1.element[ie].setVertNumb(vertNumb);	
					else{

						int[] vert1=new int[md0.nElVert];
						for(int k=0;k<md0.nElVert;k++)
							if(!onAx[vertNumb[k]] )
								vert1[k]=map[vertNumb[k]];
							else
								vert1[k]=vertNumb[k];

						int[] vert2=new int[md0.nElVert];
						
						if(md1.dim==2){
						for(int k=0;k<md0.nElVert;k++)				
							vert2[k]=vert1[md0.nElVert-1-k];	
						}
						else {

							for(int k=0;k<md0.nElVert;k++)	
								if(k<md0.nElVert/2)
								vert2[k]=vert1[md0.nElVert/2+k];	
								else
								vert2[k]=vert1[k-md0.nElVert/2];	
							}
						

						for(int k=0;k<md0.nElVert;k++){
							md1.element[ie].setVertNumb(k,vert2[k]);

						}


					}
				}
			}

		}


		md1.scaleFactor=md0.scaleFactor;
		String file = System.getProperty("user.dir") + "//reflected.txt";
		md1.writeMesh(file);

	}

	public void copy(Vect v, int N){
		String bun=util.getFile();
		Model model=new Model(bun);
		Model model1=model.deepCopy();

		for(int j=1;j<=N;j++){
			Model model2=model1.deepCopy();
		translate(model2,v.times(j),false);
		model=assemble(model2, model, false);
		}
		String folder=new File(bun).getParentFile().getPath();
		model.writeMesh(folder+"\\copies.txt");
	}


	public int[] getBorderNodeSorted(Model model,int nb){

		return bs.getBorderNodeSorted(model,nb);

	}


	public void setNodeOnMotorBound(Model model){

		bs.setNodeOnMotorBound(model);

	}

	public int[][] getBoundaryNodesMap(Model model,int nb1,int nb2){

		return  bs.mapBorderNodes(model,nb1,nb2);

	}



	public void triangToQuad(){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;

		Model model=new Model();
		model.loadMesh(bun);

		String quadMesh = System.getProperty("user.dir") + "//quadEl.txt";


		if(model.elCode==0)
			model.writeMeshTriToQuad(quadMesh);

	}

	public void triangToQuadCoarse(){



		//String bun=util.getFile(0);
		String bun = System.getProperty("user.dir") + "//trN.txt";

		if(bun==null || bun.equals("") )return;

		Model model=new Model();
		model.loadMesh(bun);

		String quadMesh = System.getProperty("user.dir") + "//quadEl.txt";


		if(model.elCode==0)
			writer.writeMeshTriToQuadCoarse(model,quadMesh);

	}
	
	public Model getModel(){
		String bun=util.getFile();
		if(bun==null || bun.equals("") )return null;

		Model model=new Model();
		model.loadMesh(bun);
		
		return model;
	}
	
	

	public Vect fetchNodeCoord(int n){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )  throw new NullPointerException("file not found.");
		Model model=new Model();
		model.loadMesh(bun);
		return fetchNodeCoord(model,n);

	}

	public Vect fetchNodeCoord(Model model, int n){

		Vect P=model.node[n].getCoord();
		util.pr("coodrinates of node "+n+":")	;	 
		P.show();
		return P;

	}
	
	public int getNodeOf(Vect v,double er){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )  throw new NullPointerException("file not found.");
		Model model=new Model();
		model.loadMesh(bun);
		return getNodeOf(model,v,er);

	}
	
	
	public int getNodeOf(Model model,Vect v,double er){

		int n=0;
		for(int i=1;i<=model.numberOfNodes;i++)
		{
			if(model.node[i].getCoord().sub(v).norm()<er){
				util.pr(i);
				n=i;
				model.node[i].getCoord().hshow();
						break;
			}
		}
		
		return n;
	}
	

	public void pileUpPrism(int K,double h){
		pileUpPrism(K,h,0);

	}

	public void pileUpPrism(int K,double h,double z0){
		double[] dh=new double[K];
		double dh0=h/K;
		for(int j=0;j<dh.length;j++)
			dh[j]=dh0;

		pileUpPrism(dh,z0);

	}


	public void pileUpPrism(double[] dh,double z0){

		int K=dh.length;


		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model(bun);

		if(model.elCode>1) return;
		String s="prism";
		int d=3,nElVert=6;
		if(model.elCode==1){
			s="hexahedron";
			d=4;
			nElVert=8;
		}

		int nReg=model.numberOfRegions;
		int nEls=model.numberOfElements*(K);
		int nNodes=model.numberOfNodes*(K+1);

		Model prismModel=new Model();
		prismModel.alloc(nReg, nEls, nNodes, s);

		double t=z0;
		for(int j=0;j<=K;j++){
			if(j>0)
				t+=dh[j-1];
			for(int i=1;i<=model.numberOfNodes;i++){
				Vect P=model.node[i].getCoord();
				{

					Vect Pp=new Vect(P.el[0],P.el[1],t);
					prismModel.node[j*model.numberOfNodes+i].setCoord(Pp);
				}
			}

		}

		int[] vertNumbP=new int[nElVert];
		int nn=model.numberOfNodes;
		util.pr(model.numberOfElements);
		util.pr(prismModel.numberOfElements);

		int net=0;
		for(int ir=1;ir<=model.numberOfRegions;ir++){
			for(int j=0;j<K;j++){
				for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

					net++;

					int[] vertNumb=model.element[i].getVertNumb();

					for(int k=0;k<d;k++)
						vertNumbP[k]=vertNumb[k]+(j+1)*nn;
					for(int k=d;k<nElVert;k++)
						vertNumbP[k]=vertNumb[k-d]+j*nn;
					

					prismModel.element[net].setVertNumb(vertNumbP);
				}
			}
		}

		prismModel.region[1].setFirstEl(1);
		prismModel.region[1].setLastEl(K*model.region[1].getLastEl());
		prismModel.region[1].setName(model.region[1].getName());

		for(int i=2;i<=model.numberOfRegions;i++){

			prismModel.region[i].setFirstEl(prismModel.region[i-1].getLastEl()+1);
			prismModel.region[i].setLastEl(prismModel.region[i].getFirstEl()+K*model.region[i].getNumbElements()-1);

			prismModel.region[i].setName(model.region[i].getName());
		}

		prismModel.scaleFactor=	model.scaleFactor;

		String folder=new File(bun).getParentFile().getPath();

		String prismMesh = folder + "//prismatic.txt";
		
		prismModel.writeMesh(prismMesh);

	}

	public void pileUpHexa(int K,double h){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model(bun);

		if(model.elCode!=1) return;
		int nReg=model.numberOfRegions;
		int nEls=model.numberOfElements*(K);
		int nNodes=model.numberOfNodes*(K+1);

		Model prismModel=new Model();
		prismModel.alloc(nReg, nEls, nNodes, "hexahedron");

		for(int j=0;j<=K;j++)
			for(int i=1;i<=model.numberOfNodes;i++){
				Vect P=model.node[i].getCoord();
				{
					Vect Pp=new Vect(P.el[0],P.el[1],j*h);
					prismModel.node[j*model.numberOfNodes+i].setCoord(Pp);
				}
			}

		int[] vertNumbP=new int[8];
		int nn=model.numberOfNodes;


		int net=0;
		for(int ir=1;ir<=model.numberOfRegions;ir++){
			int ner=model.region[ir].getNumbElements();
			for(int j=0;j<K;j++){
				for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

					int nf=model.region[ir].getFirstEl();
					net++;

					int[] vertNumb=model.element[i].getVertNumb();

					for(int k=0;k<4;k++)
						vertNumbP[k]=vertNumb[k]+(j+1)*nn;
					for(int k=4;k<8;k++)
						vertNumbP[k]=vertNumb[k-4]+j*nn;

					prismModel.element[net].setVertNumb(vertNumbP);
				}
			}
		}

		prismModel.region[1].setFirstEl(1);
		prismModel.region[1].setLastEl(K*model.region[1].getLastEl());
		prismModel.region[1].setName(model.region[1].getName());

		for(int i=2;i<=model.numberOfRegions;i++){

			prismModel.region[i].setFirstEl(prismModel.region[i-1].getLastEl()+1);
			prismModel.region[i].setLastEl(prismModel.region[i].getFirstEl()+K*model.region[i].getNumbElements()-1);

			prismModel.region[i].setName(model.region[i].getName());
		}

		prismModel.scaleFactor=	model.scaleFactor;

		String folder=new File(model.meshFilePath).getParentFile().getPath();
		String prismMesh = folder+ "\\prismatic.txt";
		prismModel.writeMesh(prismMesh);

	}

	public Model pileRotate(Model model,int K,double dtt,boolean write){
		
		Model prismModel=new Model();
		
		if(model.elCode>1) return null;
		int nv=3;
		String et="prism";
		if(model.elCode==1) {
			nv=4;
			et="hexahedron";
		}

		int nReg=model.numberOfRegions;
		int nEls=model.numberOfElements*(K);
		int nNodes=model.numberOfNodes*(K+1);

		prismModel.alloc(nReg, nEls, nNodes,et);

		for(int j=0;j<=K;j++){
			Mat R2D=util.rotMat2D(j*dtt);
			Mat R=new Mat(3,3);
			R.el[0][0]=1;
			for(int m=1;m<3;m++)
				for(int n=1;n<3;n++)
					R.el[m][n]=R2D.el[m-1][n-1];

			for(int i=1;i<=model.numberOfNodes;i++){
				Vect P=model.node[i].getCoord();

				Vect Pp=R.mul(P.v3());

				prismModel.node[j*model.numberOfNodes+i].setCoord(Pp);

			}
		}
		int[] vertNumbP=new int[2*nv];
		int nn=model.numberOfNodes;


		int net=0;
		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for(int j=0;j<K;j++){

				boolean touch=(abs((j+1)*dtt-2*PI)<1e-3);
				for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

					net++;

					int[] vertNumb=model.element[i].getVertNumb();

					for(int k=0;k<nv;k++)
						if(!touch)
							vertNumbP[k]=vertNumb[k]+(j+1)*nn;
						else
							vertNumbP[k]=vertNumb[k];
					for(int k=nv;k<2*nv;k++)
						vertNumbP[k]=vertNumb[k-nv]+j*nn;

					prismModel.element[net].setVertNumb(vertNumbP);
				}
			}
		}

		prismModel.region[1].setFirstEl(1);
		prismModel.region[1].setLastEl(K*model.region[1].getLastEl());
		prismModel.region[1].setName(model.region[1].getName());

		for(int i=2;i<=model.numberOfRegions;i++){

			prismModel.region[i].setFirstEl(prismModel.region[i-1].getLastEl()+1);
			prismModel.region[i].setLastEl(prismModel.region[i].getFirstEl()+K*model.region[i].getNumbElements()-1);

			prismModel.region[i].setName(model.region[i].getName());
		}

		prismModel.scaleFactor=	model.scaleFactor;


	
			Mat R=util.rotMat(new Vect(0,0,1), new Vect(1,0,0));	
			
			R=util.rotEuler(new Vect(0,0,1), -PI/2).mul(R);
		
			rotate(prismModel,R,write);
			
			return prismModel;
	}

	
	public void pileRotate(int K,double dtt){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model(bun);

		if(model.elCode>1) return;
		int nv=3;
		String et="prism";
		if(model.elCode==1) {
			nv=4;
			et="hexahedron";
		}



		int nReg=model.numberOfRegions;
		int nEls=model.numberOfElements*(K);
		int nNodes=model.numberOfNodes*(K+1);

		Model prismModel=new Model();
		prismModel.alloc(nReg, nEls, nNodes,et);

		for(int j=0;j<=K;j++){
			Mat R2D=util.rotMat2D(j*dtt);
			Mat R=new Mat(3,3);
			R.el[0][0]=1;
			for(int m=1;m<3;m++)
				for(int n=1;n<3;n++)
					R.el[m][n]=R2D.el[m-1][n-1];

			for(int i=1;i<=model.numberOfNodes;i++){
				Vect P=model.node[i].getCoord();

				Vect Pp=R.mul(P.v3());

				prismModel.node[j*model.numberOfNodes+i].setCoord(Pp);

			}
		}
		int[] vertNumbP=new int[2*nv];
		int nn=model.numberOfNodes;


		int net=0;
		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for(int j=0;j<K;j++){

				boolean touch=(abs((j+1)*dtt-2*PI)<1e-3);
				for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

					net++;

					int[] vertNumb=model.element[i].getVertNumb();

					for(int k=0;k<nv;k++)
						if(!touch)
							vertNumbP[k]=vertNumb[k]+(j+1)*nn;
						else
							vertNumbP[k]=vertNumb[k];
					for(int k=nv;k<2*nv;k++)
						vertNumbP[k]=vertNumb[k-nv]+j*nn;

					prismModel.element[net].setVertNumb(vertNumbP);
				}
			}
		}

		prismModel.region[1].setFirstEl(1);
		prismModel.region[1].setLastEl(K*model.region[1].getLastEl());
		prismModel.region[1].setName(model.region[1].getName());

		for(int i=2;i<=model.numberOfRegions;i++){

			prismModel.region[i].setFirstEl(prismModel.region[i-1].getLastEl()+1);
			prismModel.region[i].setLastEl(prismModel.region[i].getFirstEl()+K*model.region[i].getNumbElements()-1);

			prismModel.region[i].setName(model.region[i].getName());
		}

		prismModel.scaleFactor=	model.scaleFactor;


	
			Mat R=util.rotMat(new Vect(0,0,1), new Vect(1,0,0));	
			
			R=util.rotEuler(new Vect(0,0,1), -PI/2).mul(R);
		
		//	rotate(prismModel,R,true);

			
			//rotate(-90*PI/180);
		String folder=	new File(model.meshFilePath).getParentFile().getPath();

		String prismMesh = folder+ "//prismatic.txt";
		prismModel.writeMesh(prismMesh);

	}


	public void triFiner(double r1, double r2){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);


		String finerMesh = System.getProperty("user.dir") + "//triangElH.txt";


		if(model.elCode==0)
			model.writeMeshTriToTri(finerMesh, r1,  r2);

	}

	public void rescale(double scale){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;


		Model model=new Model();
		model.loadMesh(bun);

		String folder=new File(bun).getParentFile().getPath();
		String sqaledMesh = folder + "//scaled.txt";


		for(int i=1;i<=model.numberOfNodes;i++)
		{
			Vect v=model.node[i].getCoord().times(scale);
		
			//if(util.getAng(v)<.003) v.el[1]-=.00002;

			//	if(v.norm()<.027755) v=v.normalized().times(.02775);
			model.node[i].setCoord(v);
		}
		model.writeMesh(sqaledMesh);



	
		
	}
	
	
	public void rescale(Vect scale){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;


		Model model=new Model();
		model.loadMesh(bun);

		String folder=new File(bun).getParentFile().getPath();
		String sqaledMesh = folder + "//scaled.txt";



		for(int i=1;i<=model.numberOfNodes;i++)
		{
			Vect v=model.node[i].getCoord().times(scale);

			//	if(v.norm()<.027755) v=v.normalized().times(.02775);
			model.node[i].setCoord(v);
		}
		model.writeMesh(sqaledMesh);



	}


	public void quadToTriang(){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);


		String triangMesh = System.getProperty("user.dir") + "//triangEl.txt";

		if(model.elCode==1)
			model.writeMeshQuadToTri(triangMesh);

	}

	public void quadToTriang(String quad,String tri){


		Model model=new Model();
		model.loadMesh(quad);


		String triangMesh = System.getProperty("user.dir") + "//triangEl.txt";

		if(model.elCode==1)
			model.writeMeshQuadToTri(tri);

	}

	public void quadToHexa(){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		String hexMesh = System.getProperty("user.dir") + "//hexaEl.txt";

		if(model.elCode==1)
			model.writeMeshq2h(hexMesh,false);

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
	
	public void pyrToTetra()
	
	{
		String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model();
	model.loadMesh(bun);

	String tetMesh = System.getProperty("user.dir") + "//tetEl.txt";

	if(model.elCode==5)
		writer.writeMeshPyramidToTetra(model,tetMesh);
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
	
public void hexaToTetra()
	
	{
		String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model();
	model.loadMesh(bun);

	String pyrMesh = System.getProperty("user.dir") + "//tetEl.txt";

	if(model.elCode==4)
		writer.writeMeshHexaToTetra(model,pyrMesh);
	}
	
	public void prismToHexa()
	
	{
		String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model();
	model.loadMesh(bun);

	String hexaMesh = System.getProperty("user.dir") + "//hexaEl.txt";

	if(model.elCode==3)
		model.prismToHexa(hexaMesh);

	}

	public void reRegionFullMotor(){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		String[] rgs=new String[1+model.numberOfRegions];
		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(nr==1) rgs[nr]="ironR";

			else if(nr<10) rgs[nr]="PM"+(nr-1);
			else rgs[nr]=model.region[nr-6].getName();

		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(rgs[nr]!=null)
				model.region[nr].setName(rgs[nr]);
		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){

			if(nr==1) continue;

			for(int i=model.region[nr].getFirstEl();i<=model.region[nr].getLastEl();i++){

				if(nr>3) { model.element[i].setRegion(model.element[i].getRegion()+6);}
				else {
					Vect c=model.getElementCenter(i);
					double ang=util.getAng(c);
					int kr=2+(int)(Math.floor(ang/PI*4));
					model.element[i].setRegion(kr);



				}



			}
		}



		reRegionGroupEls(model);

		String file = System.getProperty("user.dir") + "//bunReReg.txt";

		model.writeMesh(file);
	}

	public void reRegionFullMotorLegs(){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		String[] rgs=new String[1+model.numberOfRegions];
		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(nr<28)  rgs[nr]=model.region[nr].getName();

		}
		rgs[28]="legs";
		rgs[29]="airOut";

		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(rgs[nr]!=null)
				model.region[nr].setName(rgs[nr]);
		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){

			if(nr<28) continue;

			for(int i=model.region[nr].getFirstEl();i<=model.region[nr].getLastEl();i++){


				Vect c=model.getElementCenter(i);
				double ang=util.getAng(c);
				if(abs(ang-3*PI/2+PI/6)<PI/19) model.element[i].setRegion(28);
				else if(abs(ang-3*PI/2-PI/6)<PI/19) model.element[i].setRegion(28);
				else model.element[i].setRegion(29);


			}




		}



		reRegionGroupEls(model);

		String file = System.getProperty("user.dir") + "//bunReReg.txt";

		model.writeMesh(file);
	}

	public void reRegionb(){




		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				
/**/
/*				if(ir>7) model.element[i].setRegion(5);
				else if(ir==7) model.element[i].setRegion(4);
				else if(ir<4) model.element[i].setRegion(1);
				else {
					Vect c=model.getElementCenter(i);	
					double tt=util.getAng(c);
					if(tt>PI/2/9 && tt<3*PI/2/9) model.element[i].setRegion(2);
					else if(tt>6*PI/2/9 && tt<8*PI/2/9) model.element[i].setRegion(3);
					else model.element[i].setRegion(1);
				}*/
				//if(ir<4){ model.element[i].setRegion(1);
			//	model.element[i].setRegion(3);
				Vect c=model.getElementCenter(i);	
				
/*				if(ir==1){
				if(c.el[0]<.004&&c.el[1]<.004) model.element[i].setRegion(2);
				else if(c.el[0]<.006&&c.el[1]<.006) model.element[i].setRegion(3);
				else if(c.el[0]>.009||c.el[1]>.009) model.element[i].setRegion(3);
				}else if(c.el[0]<.004&&c.el[1]<.004 && c.el[2]<.01) model.element[i].setRegion(2);
				else				
					 model.element[i].setRegion(3);*/
				
	/*			int[][] regs ={ {24, 26, 31, 25, 27, 10, 8, 7, 9, 6, 3, 2, 4, 5, 1, 39, 41, 40, 43, 37},
						{33, 30, 29, 11, 12, 13, 35, 34},
						{23, 32, 28, 22, 18, 17, 14, 20, 21, 19, 16, 15, 44, 42, 38, 36}};*/
				
				int[][] regs ={ {6},
						{2},
						{1,3,4,5,7}};
				int ix=0;
				
				for(int j=0;j<regs.length;j++)
					for(int k=0;k<regs[j].length;k++){
						if(ir==regs[j][k]){
							ix=j+1;
							break;
						}
					}
				if(ix>0)
				model.element[i].setRegion(ix);
				
				double r=c.v2().norm();
				
				if(r<-.49) model.element[i].setRegion(2);

			//	else model.element[i].setRegion(2);
				double tt=util.getAng(c);
				double rr=c.norm();
				//if(rr<.05) model.element[i].setRegion(1);
				//else
					// model.element[i].setRegion(2);
				//if(rr>1.2 && rr<1.4 && tt>0*PI/18/2 && tt<2*PI/18/2) model.element[i].setRegion(2);
				//else if(rr>1.2 && rr<1.4 && tt>4*PI/18/2 && tt<6*PI/18/2) model.element[i].setRegion(3);
			//	}
				
				//if(model.getElementArea(i)<1e-7) model.element[i].setRegion(2);

			}
		}

		Vect vp=new Vect(70.711,70.711).times(.001);
		Vect vo=new Vect(50,-86.603).times(.001);
		
		for(int ir=1;ir<=0*model.numberOfRegions;ir++){

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


				Vect c=model.getElementCenter(i);		
				Vect v2=c.v2();
				double r=v2.norm();
			
	if(r<.075 ){
		
		if(c.el[2]<.0145) model.element[i].setRegion(2);
		else if( r> .04 && c.el[2]<.0245) model.element[i].setRegion(2);
					
				}
				
			if(r>.083 ){
					if(c.el[2]>.028) model.element[i].setRegion(2);
				}
				else{
				
						if(r>.048 && c.el[2]>.03) model.element[i].setRegion(2);
						else
							if(r>.028 && c.el[2]>.032) model.element[i].setRegion(2);
				}
			
				

			}
		}


			reRegionGroupEls(model);
			
			String folder=new File(bun).getParentFile().getPath();
		String bunFilePath = folder + "//bunReReg.txt";


		model.writeMesh(bunFilePath);


	}
	
	private void reverseHexa(){




		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				
				int[] vn1=model.element[i].getVertNumb();
				int[] vn2=new int[vn1.length];
				
				vn2[0]=vn1[0];
				vn2[1]=vn1[3];
				vn2[2]=vn1[2];
				vn2[3]=vn1[1];
				vn2[4]=vn1[4];
				vn2[5]=vn1[7];
				vn2[6]=vn1[6];
				vn2[7]=vn1[5];
				model.element[i].setVertNumb(vn2);

			}
		}

		
			String folder=new File(bun).getParentFile().getPath();
		String bunFilePath = folder + "//bunReversed.txt";


		model.writeMesh(bunFilePath);


	}

	public void reRegionf(){




		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		
		for(int ir=1;ir<=1*model.numberOfRegions;ir++){

		//	if(ir==4|| ir>7) continue;
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				int[] vn=model.element[i].getVertNumb();
				
				int temp=vn[1];
				vn[1]=vn[3];
				vn[3]=temp;
				model.element[i].setVertNumb(vn);
				Vect c=model.getElementCenter(i);		

				//for(int k=0;k<8;k++)
			
				if(ir==8 || ir==13 || ir==18) model.element[i].setRegion(3);
				else 	if(ir==7 || ir==12 || ir==17 || ir==22) model.element[i].setRegion(4);
				else 	if(ir==23) model.element[i].setRegion(5);
				//if(c.el[2]>.005) 	model.element[i].setRegion(6);

			}
		}
		
		

			reRegionGroupEls(model);
		String bunFilePath = System.getProperty("user.dir") + "//bunReReg.txt";


		model.writeMesh(bunFilePath);


	}
	
	public void reverseFace2D(){


		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		
		for(int ir=1;ir<=1*model.numberOfRegions;ir++){

			if(ir==4|| ir>7) continue;
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				int[] vn=model.element[i].getVertNumb();
				if(model.elCode==0){
					int tmp=vn[1];
					vn[1]=vn[2];
					vn[2]=tmp;
					
				}
				else if(model.elCode==1){
					int tmp1=vn[1];
					vn[1]=vn[3];
					vn[3]=tmp1;
					
				}
				model.element[i].setVertNumb(vn);


			}
		}
		
		
		String folder=new File(bun).getParentFile().getPath();
		String bunFilePath = folder + "//faceReversed.txt";


		model.writeMesh(bunFilePath);


	}

	public void reRegion(){




		String bun=util.getFile();
		//String bun = System.getProperty("user.dir") + "//prismatic.txt";

		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		
		for(int ir=1;ir<=1*model.numberOfRegions;ir++){

			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
	//	;;	if(ir>10) model.element[i].setRegion(ir-10);
				String[] sp1=model.region[ir].getName().split("region");
				int n1=Integer.parseInt(sp1[sp1.length-1]);
				for(int k=1;k<ir;k++){
					String[] sp2=model.region[k].getName().split("region");
					int n2=Integer.parseInt(sp2[sp2.length-1]);
				//	util.pr(sp2[sp2.length-1]);
					//if(ir==18)
					//	util.pr(model.region[ir].getName()+"    "+model.region[ir].getName());
					if(n1==n2){
						util.pr(n1+"  "+n2);

						if(k==18)
						util.pr(model.region[k].getName()+"  ---  "+model.region[ir].getName());

						 model.element[i].setRegion(k);
					}
				}
				
				//if(ir==18) model.element[i].setRegion(1);
				//else if(ir==19) model.element[i].setRegion(4);
				/*if(ir==1 || ir==3)model.element[i].setRegion(1);
					if(ir==2 || ir==4) model.element[i].setRegion(2);*/
				Vect c=model.getElementCenter(i);
				double r=c.norm();
				double tt=util.getAng(c);
				if(tt>6.2){
					
					//model.element[i].setRegion(4);
					}
				/*
				if(ir>8){
					double td=12;
					int v=(int)(tt*180/PI/td);
					model.element[i].setRegion(9+v%2);
					}
				else 	if(ir==5 || ir==6){
					double td=10;
					int v=(int)(tt*180/PI/td);
					model.element[i].setRegion(5+v%2);
					}
				
				else 	if(ir==1 || ir==2){
					double td=60;
					int v=(int)(tt*180/PI/td);
					model.element[i].setRegion(1+v%2);
					}*/
		
					/*if(tt<PI/5) model.element[i].setRegion(1);
					
					else if(tt<2*PI/5) model.element[i].setRegion(2);
					
					else if(tt<3*PI/5) model.element[i].setRegion(1);
					
					else if(tt<4*PI/5) model.element[i].setRegion(2);
					
					else if(tt<5*PI/5) model.element[i].setRegion(1);
					
					else if(tt<6*PI/5) model.element[i].setRegion(2);
					
					else if(tt<7*PI/5) model.element[i].setRegion(1);
					
					else if(tt<8*PI/5) model.element[i].setRegion(2);
					
					else if(tt<9*PI/5) model.element[i].setRegion(1);
					else model.element[i].setRegion(2);*/
					
							
				//}
		/*		if(r<.065 && tt<PI/10) model.element[i].setRegion(1);
				else 
					if(r<.065 && tt>PI/10) model.element[i].setRegion(2);
					else if(r<.066 ) model.element[i].setRegion(3);
					else model.element[i].setRegion(4);
*/
				
/*				if(r<.035 ) model.element[i].setRegion(1);
				else if(r<.04 ) model.element[i].setRegion(2);
				else if(r<.06) model.element[i].setRegion(3);
				else if(r<.065 ) model.element[i].setRegion(4);
				else  model.element[i].setRegion(5);*/
				

			}
		}
		
		double rm=0;
for(int i=1;i<=0*model.numberOfNodes;i++){

			
				Vect c=model.node[i].getCoord();
				double r=c.norm();
				double tt=util.getAng(c);
				
				if(r>rm) rm=r;
				
			/*	if(r>.0351 && r<.041){
					Vect cm=c.normalized();
					//Vect ct=c.times((.037/.04));
					Vect ct=cm.times(r-.002);
					model.node[i].setCoord(ct);
				}*/
				
				if(r>.03651){
					Vect cm=c.normalized();
					Vect ct=c.times(1+.0365/r);
					//Vect ct=cm.times(r-.005);
					model.node[i].setCoord(ct);
				}
/*				if(r<.035 ) model.element[i].setRegion(1);
				else if(r<.04 ) model.element[i].setRegion(2);
				else if(r<.06) model.element[i].setRegion(3);
				else if(r<.065 ) model.element[i].setRegion(4);
				else  model.element[i].setRegion(5);*/
				

			
		}


util.pr(rm);


		reRegionGroupEls(model);
		
		String bunFilePath = System.getProperty("user.dir") + "//bunReReg.txt";


		model.writeMesh(bunFilePath);


	}


	public void reRegionz(){




		String bun=util.getFile();
		//String bun = System.getProperty("user.dir") + "//prismatic.txt";

		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);

		for(int ir=1;ir<=1*model.numberOfRegions;ir++){

		
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				Vect c=model.getElementCenter(i);
				
				if(ir>8) model.element[i].setRegion(ir+7);
				
				else if(c.el[0]>-.115 && c.el[0]<-.105 && c.el[2]<.06)
					model.element[i].setRegion(1);
				else if(c.el[0]>-.115 && c.el[0]<-.105 && c.el[2]>.06)
					model.element[i].setRegion(2);
				else if(c.el[0]>-.105 && c.el[0]<-.055 )
					model.element[i].setRegion(3);
				else if(c.el[0]>-.055 && c.el[0]<-.045 && c.el[2]>.06)
					model.element[i].setRegion(4);
				else if(c.el[0]>-.055 && c.el[0]<-.045 && c.el[2]<.06)
					model.element[i].setRegion(5);
				else if(c.el[0]>-.035 && c.el[0]<-.025 && c.el[2]<.06)
					model.element[i].setRegion(6);
				else if(c.el[0]>-.035 && c.el[0]<-.025 && c.el[2]>.06)
					model.element[i].setRegion(7);
				else if(c.el[0]>-.025 && c.el[0]<.025)
					model.element[i].setRegion(8);
				else if(c.el[0]>.025 && c.el[0]<.035 && c.el[2]>.06)
					model.element[i].setRegion(9);
				else if(c.el[0]>.025 && c.el[0]<.035 && c.el[2]<.06)
					model.element[i].setRegion(10);
				else if(c.el[0]>.045 && c.el[0]<.055 && c.el[2]<.06)
					model.element[i].setRegion(11);
				else if(c.el[0]>.045 && c.el[0]<.055 && c.el[2]>.06)
					model.element[i].setRegion(12);
				else if(c.el[0]>.055 && c.el[0]<.105)
					model.element[i].setRegion(13);
				else if(c.el[0]>.105 && c.el[0]<.115 && c.el[2]>.06)
					model.element[i].setRegion(14);
				else if(c.el[0]>.105 && c.el[0]<.115 && c.el[2]<.06)
					model.element[i].setRegion(15);
				
			
	

			
			}
		}

		for(int i=1;i<=0*model.numberOfNodes;i++){

			Vect v=model.node[i].getCoord();
			if(v.el[2]<-1e-6){
				v.el[2]+=.01;
				model.node[i].setCoord(v);
			}
		
		}



			reRegionGroupEls(model);
		String bunFilePath = System.getProperty("user.dir") + "//bunReReg.txt";


		model.writeMesh(bunFilePath);


	}
	
	public void deform(int ir,double shift){



		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);
		int[] nn=model.getRegNodes(ir);
	
		for(int i=0;i<nn.length;i++){
			int n=nn[i];
			
			Vect v=model.node[n].getCoord();
	
			model.node[n].setCoord(v.add(shift));
		}

		String bunFilePath = System.getProperty("user.dir") + "//deformed.txt";

		model.writeMesh(bunFilePath);


		
	}
	
	public void deform(){


		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);
			
	
		boolean[] nn=new boolean[1+model.numberOfNodes];
		
		for(int ir=1;ir<=model.numberOfRegions;ir++){
			if(ir<=7)
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
				{
				int[] nv=model.element[i].getVertNumb();
				for(int k=0;k<nv.length;k++)
					nn[nv[k]]=true;
				}
		}
		double R=Math.sqrt(3)*.2;
		for(int i=1;i<=1*model.numberOfNodes;i++){

		if(nn[i]){
			Vect v=model.node[i].getCoord();
		
		double r=v.v2().norm();

		double tt=util.getAng(v);
		if(tt>PI) tt-=2*PI;
		
		//double r2=v.el[0];
		//double t2=(v.el[1]/0.08)*9.*PI/180;
	//	Vect v2=new Vect(r2*Math.cos(t2),r2*Math.sin(t2),v.el[2]);
		Vect v2=new Vect(r*Math.cos(22.5/9*tt),r*Math.sin(22.5/9*tt),v.el[2]);
		model.node[i].setCoord(v2);
		
/*		if(v.el[2]>.1 && v.el[2]<.45 && r<.8 ){
			double t2=tt;
			double r2=v.el[0];
			r2=r*(1-.3*Math.random());
			
			Vect v2=new Vect(r2*Math.cos(t2),r2*Math.sin(t2),v.el[2]);
			v2.hshow();
			model.node[i].setCoord(v2);
		}*/
		
/*		if(r>.18 && r<-R){
			v=v.times(R/r);
			model.node[i].setCoord(v);
		}
		if(v.el[2]>-.19){
		//	v=v.times(.25/r);
			v=v.add(new Vect().rand(3, -.001, .001));
			model.node[i].setCoord(v);
		}
		*/
		//double tt=util.getAng(v)/2;
		//Vect v2=new Vect(r*Math.cos(tt),r*Math.sin(tt));
		//model.node[i].setCoord(v2);
		//	model.node[i].setCoord(model.node[i].getCoord().times(.9));
			}
		}
			//list1.add(c.el[2]);
		

		String folder=new File(bun).getParentFile().getPath();
		String bunFilePath = folder+ "//deformed.txt";

		model.writeMesh(bunFilePath);


	}
	
	public void edgeDirec(){



		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model();
		model.loadMesh(bun);
	
	
		
		Vect cc1=model.node[1].getCoord();
		
		List<Double> list1=new ArrayList<Double>();

			
		for(int i=1;i<=0*model.numberOfNodes;i++){

			Vect c=model.node[i].getCoord();
			
			double r=c.norm();
			if(r>.05468){
			double tt=util.getAng(c);
			double td=Math.round(tt/(PI/18/200))*(PI/18/200);
			Vect v=new Vect(r*Math.cos(td),r*Math.sin(td));
			
			model.node[i].setCoord(v);
			}
			//list1.add(c.el[2]);
		}
		int ix=0;
		int[][] nn=new int[model.numberOfNodes][3];
		int[][] vnn=new int[model.numberOfElements+1][3];

		
		for(int i=1;i<=1*model.numberOfElements;i++){

			Vect c=model.getElementCenter(i);
			double tc=util.getAng(c);
			
			Mat R=util.rotMat2D(-tc);
			//if(util.getAng(c)>.003)continue;
			double r=c.norm();
			if(r<.05465 || r >.05475) continue;
			
			//ix++;
			//util.pr(ix);
			
			int[] vn=model.element[i].getVertNumb();
		Vect	v1=model.node[vn[1]].getCoord().sub(model.node[vn[0]].getCoord());
		Vect	v2=model.node[vn[2]].getCoord().sub(model.node[vn[1]].getCoord());
		Vect	v3=model.node[vn[0]].getCoord().sub(model.node[vn[2]].getCoord());
		
		double t1=util.getAng(R.mul(v1));
		double t2=util.getAng(R.mul(v2));
		double t3=util.getAng(R.mul(v3));
		
	/*	util.pr(t1);
		util.pr(t2);
		util.pr(t3);*/

	
		
		if(t1<3 && t1>PI/2+.1 || t1<6 && t1>3*PI/2+.1) { nn[ix][0]=vn[0];nn[ix][1]=vn[1];}
		else if(t2<3 &&t2>PI/2+.1 || t2<6 && t2>3*PI/2+.1) {nn[ix][0]=vn[1];nn[ix][1]=vn[2];}
		else if(t3<3 &&t3>PI/2+.1|| t3<6 && t3>3*PI/2+.1) {nn[ix][0]=vn[2];nn[ix][1]=vn[0];}
		
	
	if(nn[ix][0]>0) {nn[ix][2]=i; ix++;}

		}
		
		
		for(int j=0;j<ix;j++){
			int ie1=nn[j][2];
			int[] vn1=model.element[ie1].getVertNumb();
			int vn13=0;
			for(int u=0;u<3;u++)
				if(vn1[u]!=nn[j][0] && vn1[u]!=nn[j][1]) vn13=vn1[u];
			for(int k=j+1;k<ix;k++)
				if(nn[k][1]==nn[j][0] && nn[k][0]==nn[j][1]){
					int ie2=nn[k][2];
					util.pr(ie1+"  "+ie2);
					int[] vn2=model.element[ie2].getVertNumb();
					
					int vn23=0;
					for(int u=0;u<3;u++)
						if(vn2[u]!=nn[k][0] && vn2[u]!=nn[k][1]) vn23=vn2[u];
					
					 vnn[ie1][0]=vn23;
					 vnn[ie1][1]=nn[k][0];
					 vnn[ie1][2]=vn13;
					 
					 vnn[ie2][0]=vn13;
					 vnn[ie2][1]=nn[j][0];
					 vnn[ie2][2]=vn23;
				}
				
		}
		
			
			for(int j=1;j<=model.numberOfElements;j++)
				if( vnn[j][0]>0) model.element[j].setVertNumb(vnn[j]);
				
/*
			
Vect v=new Vect(c.el[0],c.el[1]);
double rn=v.norm();

if(abs(rn-.025)<-1e-4){
	double r2=.026;
	Vect v1=v.times(r2/rn);
	Vect vx=new Vect(v1.el[0],v1.el[1],c.el[2]);
//	util.pr(rn);

//	vx.hshow();
	model.node[i].setCoord(vx);
}

		}
		
Set<Double> set = new HashSet<Double>(list1);

ArrayList<Double> hh = new ArrayList<Double>(set);

double[] h=new double[hh.size()];
for(int i=0; i<hh.size(); i++){
	h[i]=hh.get(i);
}

double[] dh=new double[hh.size()-1];

Vect v=new Vect(h);
v.bubble();

for(int i=1; i<hh.size(); i++){
	dh[i-1]=v.el[i]-v.el[i-1];
}

for(int i=0; i<dh.length; i++){
	System.out.print(dh[i]+",");
}*/

		

		String bunFilePath = System.getProperty("user.dir") + "//deformed.txt";


		model.writeMesh(bunFilePath);


	}



	public void reRegionIPMS8fold(String bun,String outBun){


		Model model0=new Model(bun);

		int ntr=16;
		Model model=model0.fill(ntr,model0.numberOfElements,model0.numberOfNodes,model0.elType);

		double f0=-PI/36;
		double f1=PI/18;

		String[] rgs=new String[1+model.numberOfRegions];
		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(nr==1) rgs[nr]="ironS";
			else if(nr<3) rgs[nr]="coilU+";
			else if(nr<4) rgs[nr]="coilUh+";
			else if(nr<5) rgs[nr]="coilU-";
			else if(nr<6) rgs[nr]="coilUh-";
			else if(nr<7) rgs[nr]="coilV+";
			else if(nr<8) rgs[nr]="coilVh+";
			else if(nr<9) rgs[nr]="coilV-";
			else if(nr<10) rgs[nr]="coilVh-";
			else if(nr<11) rgs[nr]="coilW+";
			else if(nr<12) rgs[nr]="coilWh+";
			else if(nr<13) rgs[nr]="coilW-";
			else if(nr<14) rgs[nr]="coilWh-";
			else if(nr<15) rgs[nr]="airgapS";
			else if(nr<16) rgs[nr]="frame";
			else if(nr<17) rgs[nr]="airOut";
		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(rgs[nr]!=null)
				model.region[nr].setName(rgs[nr]);
		}

		int[] regNumbs1={2,3,12,6,7,4,10,11,8,2};
		int[] regNumbs2={2,13,12,6,5,4,10,9,8,2};


		for(int nr=1;nr<=model.numberOfRegions;nr++){

			if(nr==1) continue;

			for(int i=model.region[nr].getFirstEl();i<=model.region[nr].getLastEl();i++){

				model.element[i].setRegion(1);


				if(nr>3) { model.element[i].setRegion(nr+10);}


				Vect c=model.getElementCenter(i);
				double ang=util.getAng(c);

				if(nr==2){

					for(int k=0;k<regNumbs1.length;k++){
						if(abs(ang-k*f1+f0)<f1) model.element[i].setRegion(regNumbs1[k]);
					}


				}

				else if(nr==3){

					for(int k=0;k<regNumbs1.length;k++){
						if(abs(ang-k*f1+f0)<f1) model.element[i].setRegion(regNumbs2[k]);
					}
				}

			}
		}

		reRegionGroupEls(model);

		model.writeMesh(outBun);
	}


	public void reRegionIPMSFull(String bun,String outBun){


		Model model0=new Model(bun);

		int ntr=17;
		Model model=model0.fill(ntr,model0.numberOfElements,model0.numberOfNodes,model0.elType);

		double f0=-PI/36;
		double f1=PI/18;

		String[] rgs=new String[1+model.numberOfRegions];
		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(nr==1) rgs[nr]="ironS";
			else if(nr<3) rgs[nr]="coilU+";
			else if(nr<4) rgs[nr]="coilUh+";
			else if(nr<5) rgs[nr]="coilU-";
			else if(nr<6) rgs[nr]="coilUh-";
			else if(nr<7) rgs[nr]="coilV+";
			else if(nr<8) rgs[nr]="coilVh+";
			else if(nr<9) rgs[nr]="coilV-";
			else if(nr<10) rgs[nr]="coilVh-";
			else if(nr<11) rgs[nr]="coilW+";
			else if(nr<12) rgs[nr]="coilWh+";
			else if(nr<13) rgs[nr]="coilW-";
			else if(nr<14) rgs[nr]="coilWh-";
			else if(nr<15) rgs[nr]="airgapS";
			else if(nr<16) rgs[nr]="frame";
			else if(nr<17) rgs[nr]="legs";
			else if(nr<18) rgs[nr]="airOut";
		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(rgs[nr]!=null)
				model.region[nr].setName(rgs[nr]);
		}

		int[] regNumbs1={2,3,12,6,7,4,10,11,8,2};
		int[] regNumbs2={2,13,12,6,5,4,10,9,8,2};


		for(int nr=1;nr<=model.numberOfRegions;nr++){

			if(nr==1) continue;

			for(int i=model.region[nr].getFirstEl();i<=model.region[nr].getLastEl();i++){

				model.element[i].setRegion(1);


				Vect c=model.getElementCenter(i);
				double ang=util.getAng(c);

				if(nr==2){

					int m=(int)(2*ang/PI);
					ang =ang-m*PI/2;

					for(int k=0;k<regNumbs1.length;k++){

						if(abs(ang-k*f1+f0)<f1) model.element[i].setRegion(regNumbs1[k]);
					}


				}

				else if(nr==3){
					int m=(int)(2*ang/PI);
					ang =ang-m*PI/2;
					for(int k=0;k<regNumbs1.length;k++){
						if(abs(ang-k*f1+f0)<f1) model.element[i].setRegion(regNumbs2[k]);
					}
				}

				else if(nr==6){
					if(abs(ang-3*PI/2-PI/8)<PI/24) model.element[i].setRegion(16);
					else if(abs(ang-3*PI/2+PI/8)<PI/24) model.element[i].setRegion(16);
					else model.element[i].setRegion(17);

				}

				else	if(nr>3) 
					model.element[i].setRegion(nr+10);

			}
		}

		reRegionGroupEls(model);

		model.writeMesh(outBun);
	}


	public void reRegionIPMR4th(String bun,String file){

		Model model0=new Model(bun);

		int ntr=7;

		Model model=model0.fill(ntr,model0.numberOfElements,model0.numberOfNodes,model0.elType);


		String[] rgs=new String[1+model.numberOfRegions];
		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(nr==1) rgs[nr]="ironR";

			else if(nr<=3) rgs[nr]=model.region[2].getName()+(nr-1);
			else rgs[nr]=model.region[nr-1].getName();

		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(rgs[nr]!=null)
				model.region[nr].setName(rgs[nr]);
		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){

			if(nr==1) continue;

			for(int i=model.region[nr].getFirstEl();i<=model.region[nr].getLastEl();i++){

				if(nr>2) { model.element[i].setRegion(model.element[i].getRegion()+1);}
				else if(nr==2){
					Vect c=model.getElementCenter(i);
					double ang=util.getAng(c);
					if( ang<PI/4) model.element[i].setRegion(2);
					else if( ang<PI/2) model.element[i].setRegion(3);
					//	else  model.element[i].setRegion(4);


				}



			}
		}

		reRegionGroupEls(model);


		model.writeMesh(file);
	}


	public void reRegionIPMRFull(String bun,String file){


		Model model0=new Model(bun);

		int ntr=13;

		Model model=model0.fill(ntr,model0.numberOfElements,model0.numberOfNodes,model0.elType);


		String[] rgs=new String[1+model.numberOfRegions];
		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(nr==1) rgs[nr]="ironR";

			else if(nr<10) rgs[nr]=model.region[2].getName()+(nr-1);
			else rgs[nr]=model.region[nr-7].getName();

		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(rgs[nr]!=null)
				model.region[nr].setName(rgs[nr]);
		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){

			if(nr==1) continue;

			for(int i=model.region[nr].getFirstEl();i<=model.region[nr].getLastEl();i++){

				if(nr>2) { model.element[i].setRegion(model.element[i].getRegion()+7);}
				else if(nr==2){
					Vect c=model.getElementCenter(i);
					double ang=util.getAng(c);
					if( ang<PI/8 || 2*PI-ang<PI/8  ) model.element[i].setRegion(2);
					else if( ang<3*PI/8) model.element[i].setRegion(3);
					else if( ang<5*PI/8) model.element[i].setRegion(4);
					else if( ang<7*PI/8) model.element[i].setRegion(5);
					else if( ang<9*PI/8) model.element[i].setRegion(6);
					else if( ang<11*PI/8) model.element[i].setRegion(7);
					else if( ang<13*PI/8) model.element[i].setRegion(8);
					else  model.element[i].setRegion(9);


				}



			}
		}

		reRegionGroupEls(model);


		model.writeMesh(file);
	}

	public void reRegionIPMR2fold1deg(String bun,String file){




		Model model0=new Model(bun);

		int ntr=7;

		Model model=model0.fill(ntr,model0.numberOfElements,model0.numberOfNodes,model0.elType);


		String[] rgs=new String[1+model.numberOfRegions];
		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(nr==1) rgs[nr]="ironR";

			else if(nr<4) rgs[nr]=model.region[2].getName()+(nr-1);
			else rgs[nr]=model.region[nr-1].getName();

		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){
			if(rgs[nr]!=null)
				model.region[nr].setName(rgs[nr]);
		}

		for(int nr=1;nr<=model.numberOfRegions;nr++){

			if(nr==1) continue;

			for(int i=model.region[nr].getFirstEl();i<=model.region[nr].getLastEl();i++){

				if(nr>2) { model.element[i].setRegion(model.element[i].getRegion()+1);}
				else if(nr==2){
					Vect c=model.getElementCenter(i);
					double ang=util.getAng(c);
					if( ang<PI/4) model.element[i].setRegion(2);
					else  model.element[i].setRegion(3);


				}



			}
		}

		reRegionGroupEls(model);


		model.writeMesh(file);


	}

	public void reRegionGroupEls(Model model){

		for(int ir=1;ir<=model.numberOfRegions;ir++){
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
				if(model.element[i].getRegion()==0) model.element[i].setRegion(ir);
		}
		int[] mape=new int[model.numberOfElements+1];
		int[][] mapr=new int[model.numberOfRegions+1][2];

		int ix=0;
		int[][] vn=new int[model.numberOfElements+1][model.nElVert];

		for(int ir=1;ir<=model.numberOfRegions;ir++){
			mapr[ir][0]=ix+1;
			for(int i=1;i<=model.numberOfElements;i++)
				if(model.element[i].getRegion()==ir)
				{
					mape[++ix]=i;
					vn[ix]=model.element[i].getVertNumb();
				}
			mapr[ir][1]=ix;
		}


		for(int ir=1;ir<=model.numberOfRegions;ir++){
			model.region[ir].setFirstEl(mapr[ir][0]);
			model.region[ir].setLastEl(mapr[ir][1]);
		}


		for(int i=1;i<=model.numberOfElements;i++){
			model.element[i].setVertNumb(vn[i]);
		}
	}


	public void reRegionIPMrotor8th(){

		String  bun=System.getProperty("user.dir") + "//meshFromDel.txt";

		Model model=new Model(bun);

		double rout=.0544;
		double rshaft=.016;
		double x0=.0455,x1=.0505;
		double y0=-.0145,y1=.0145;
		double gp2=.0001;
		double xg0=x0+gp2,xg1=x1-gp2;
		double yg0=y0+gp2,yg1=y1-gp2;
		double Rs=.040;
		double rsi=.00325;

		Vect bc=new Vect(Rs,0);

		for(int i=1;i<=model.numberOfElements;i++){
			Vect c=model.getElementCenter(i);
			double cn=c.norm();
			double se=model.getElementArea(i);

			if(model.el3angMaxCosAng(i)>.9999 || se<1e-9) model.element[i].setRegion(model.numberOfRegions);
			else if(cn<rshaft){ model.element[i].setRegion(3);	}
			else if(cn>rout){ model.element[i].setRegion(model.numberOfRegions-1);	}

			else if(c.el[0]>x0 && c.el[0]<x1 &&c.el[1]>y0 &&c.el[1]<y1 ){

				if(c.el[0]>xg0 && c.el[0]<xg1 &&c.el[1]>yg0 &&c.el[1]<yg1 )
					model.element[i].setRegion(2);
				else
					model.element[i].setRegion(5);
			}
			else if(c.sub(bc).norm()<rsi){ model.element[i].setRegion(4);	}

		}

		model.region[1].setName("ironR");
		model.region[2].setName("PM");
		model.region[3].setName("shaft");
		model.region[4].setName("bolt");
		model.region[5].setName("gapPM");
		model.region[6].setName("airgap");



		reRegionGroupEls(model);

		String bunFilePath = System.getProperty("user.dir") + "//rotor8th.txt";

		model.writeMesh(bunFilePath);



	}


	public void reRegionIPMSHalfSlot(){


		String bun=System.getProperty("user.dir") + "//meshFromDel.txt";

		Model model=new Model(bun);

		double rin=.0547,gp=.0003,tp=.0007;
		double rsOut=rin+.0223, ric=.00388;
		double rsCenter=rsOut-.00388, rout1=rsOut+.01047, rout2=rout1+.025;
		double rmc=.0675;

		double PI=Math.PI;	
		double f1=5.0/180*PI;
		double f2=3.5/180*PI;
		double f3=2.5/180*PI;
		for(int i=1;i<=model.numberOfElements;i++){
			Vect c=model.getElementCenter(i);
			double ang=util.getAng(c);
			double cn=c.norm();
			double se=model.getElementArea(i);

			if(cn<rin || model.el3angMaxCosAng(i)>.999 || se<1e-9) model.element[i].setRegion(model.numberOfRegions);
			else if(cn>rout1){ model.element[i].setRegion(model.numberOfRegions-1);	}
			else if(cn<rin+gp) model.element[i].setRegion(model.numberOfRegions-2);
			else if(cn<rin+gp+tp && ang>f2 ) model.element[i].setRegion(model.numberOfRegions-2);

			else if(cn<rin+gp+tp ){

				if(ang>f2)model.element[i].setRegion(model.numberOfRegions-2);
			}

			else if( cn<rsOut && ang>f3 ){

				if( cn<rmc) model.element[i].setRegion(2);

				else {
					double d=c.sub(util.rotMat2D(f1).mul(new Vect(rsCenter,0))).norm();
					if(cn<rsCenter || d<ric)
						model.element[i].setRegion(3);
				}
			}

			else
				model.element[i].setRegion(1);

		}



		model.region[1].setName("ironS");
		model.region[2].setName("slotIn");
		model.region[3].setName("slotOut");
		model.region[4].setName("airGap");
		model.region[5].setName("airOut");





	}

	public void reRegionIPMSHalfSlotTip(String bun){



		Model model=new Model(bun);

		double rm=.0547,gp=.0003,tp=.0007,tph=.00033;
		double rsOut=rm+.0223, ric=.00388;
		double rsCenter=rsOut-.00388, rout1=rsOut+.01047, rout2=rout1+.025;
		double rring=rout1+.0045;
		double rmc=.0675;

		double PI=Math.PI;	
		double f1=5.0/180*PI;
		double f2=3.5/180*PI;
		double f3=2.3/180*PI;
		for(int i=1;i<=model.numberOfElements;i++){
			Vect c=model.getElementCenter(i);
			double ang=util.getAng(c);
			double cn=c.norm();
			double se=model.getElementArea(i);

			if(cn<rm || model.el3angMaxCosAng(i)>.999 || se<1e-9) model.element[i].setRegion(model.numberOfRegions);
			else if(cn>rring){ model.element[i].setRegion(model.numberOfRegions-1);}
			else if(cn>rout1){ model.element[i].setRegion(model.numberOfRegions-2);	}
			else if(cn<rm+gp) model.element[i].setRegion(model.numberOfRegions-3);

			else if(cn<rm+gp+tp){
				if(ang>f2 ) model.element[i].setRegion(model.numberOfRegions-3);
				else
					model.element[i].setRegion(1);
			}

			else if(cn<rm+gp+tp+tph   ){

				if( ang>f2)
					model.element[i].setRegion(2);

			}

			else if(cn<rm+gp+tp+tph   ){

				if( 3*c.el[0]+c.el[1]>3*rm+7*tp)
					model.element[i].setRegion(2);

			}
			else if( cn<rsOut && ang>f3 ){

				if( cn<rmc) model.element[i].setRegion(2);

				else {
					double d=c.sub(util.rotMat2D(f1).mul(new Vect(rsCenter,0))).norm();
					if(cn<rsCenter || d<ric)
						model.element[i].setRegion(3);
				}
			}

			else
				model.element[i].setRegion(1);

		}


		model.region[1].setName("ironS");
		model.region[2].setName("slotIn");
		model.region[3].setName("slotOut");
		model.region[4].setName("airGap");
		model.region[5].setName("frame");
		model.region[6].setName("airOut");




		reRegionGroupEls(model);

		String bunFilePath = System.getProperty("user.dir") + "//slotH.txt";

		model.writeMesh(bunFilePath);


	}

	public void reRegionFrame(String bun){

		double r1=.0875;
		double r2=r1+.0045;

		Model model=new Model(bun);


		for(int i=1;i<=model.numberOfElements;i++){
			Vect c=model.getElementCenter(i);
			double cn=c.norm();
			double se=model.getElementArea(i);

			if(cn<r1 || model.el3angMaxCosAng(i)>.999 || se<1e-9) model.element[i].setRegion(model.numberOfRegions);

			else
				model.element[i].setRegion(1);

		}


		model.region[1].setName("steel");

		reRegionGroupEls(model);

		String bunFilePath = System.getProperty("user.dir") + "//frameH.txt";

		model.writeMesh(bunFilePath);


	}

	public void reRegionCap(String bun){

		double r0=.016;

		Model model=new Model(bun);


		for(int i=1;i<=model.numberOfElements;i++){
			Vect c=model.getElementCenter(i);
			double cn=c.norm();
			double se=model.getElementArea(i);

			if(cn<r0 || model.el3angMaxCosAng(i)>.999 || se<1e-9) model.element[i].setRegion(model.numberOfRegions);

			else
				model.element[i].setRegion(1);

		}


		model.region[1].setName("steel");

		reRegionGroupEls(model);

		String bunFilePath = System.getProperty("user.dir") + "//capH.txt";

		model.writeMesh(bunFilePath);


	}

	public void reRegionLegs(String bun){

		double r1=.0875;
		double r2=r1+.0045;
		double h1=5;
		double h2=33;
		double f3=3*PI/2-40*1.0/180*Math.PI;;
		double f4=3*PI/2+40*1.0/180*Math.PI;;
		double x1=.025;

		Model model=new Model(bun);


		for(int i=1;i<=model.numberOfElements;i++){
			Vect c=model.getElementCenter(i);
			double cn=c.norm();

			double se=model.getElementArea(i);
			model.element[i].setRegion(3);

			if(cn<r1 || model.el3angMaxCosAng(i)>.999 || se<1e-9) model.element[i].setRegion(model.numberOfRegions);


			else if(cn<r2)
				model.element[i].setRegion(1);
			else if(abs(util.getAng(c)-3*PI/2)<40*1.0/180*Math.PI)
			{
				if(c.el[1]>-r2-.0045 || abs(c.el[0])>.078-x1 )
					if(c.el[1]>-r2-.0045-.005/* && abs(c.el[0])>.088-x1*/)
						model.element[i].setRegion(2);

			}
			else
				model.element[i].setRegion(3);

			if(c.el[1]<-r2-.0055 && abs(c.el[0])>.095-x1 )
				model.element[i].setRegion(3);
		}


		model.region[1].setName("steel");
		model.region[2].setName("steel2");

		reRegionGroupEls(model);

		String bunFilePath = System.getProperty("user.dir") + "//legs.txt";

		model.writeMesh(bunFilePath);


	}


	public void dropUnusedNodes(){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model(bun);
		util.pr(model.numberOfNodes);

		dropUnusedNodesWrite(model,bun);
		util.pr(model.numberOfNodes);

	}

	public void dropUnusedNodes(String bun){

		Model model=new Model(bun);
		util.pr(model.numberOfNodes);
		dropUnusedNodesWrite(model,bun);
		util.pr(model.numberOfNodes);

	}

	public void dropUnusedNodesWrite(Model model,String bun){


		model.setInUseNodes();

		int ix=0;
		int[] map=new int[1+model.numberOfNodes];
		for(int i=1;i<=model.numberOfNodes;i++)
			if(model.node[i].inUse) {
				map[i]=++ix;
			}

		DecimalFormat formatter;
		if(model.scaleFactor==1)
			formatter= new DecimalFormat("0.0000000");
		else 
			formatter= new DecimalFormat("0.00000");
		
		String folder=new File(bun).getParentFile().getPath();
		String bunFilePath = folder + "//noUnusedNodes.txt";

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(bunFilePath)));		
			pwBun.println(model.elType);
			pwBun.println("//Number_of_Node");
			pwBun.println(ix);

			pwBun.println("//Number_of_Element");
			pwBun.println(model.numberOfElements);

			pwBun.println("//Number_of_Region");
			pwBun.println(model.numberOfRegions);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);


			for(int i=1;i<=model.numberOfElements;i++){
				int[] vertNumb=model.element[i].getVertNumb();
				for(int j=0;j<model.nElVert;j++)
					pwBun.print(map[vertNumb[j]]+",");
				pwBun.println();
			}
			for(int i=1;i<=model.numberOfNodes;i++){ 
				if(map[i]==0) continue;

				Vect xyz=model.node[i].getCoord();

				for(int j=0;j<model.dim;j++){

					pwBun.print(formatter.format(xyz.el[j]*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}

			for(int i=1;i<=model.numberOfRegions;i++){ 

				pwBun.println(model.region[i].getFirstEl()+","+model.region[i].getLastEl()+","+model.region[i].getName());

			}


			model.numberOfNodes=ix;

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

	public Model dropUnusedNodes(Model model){


		model.setInUseNodes();

		int ix=0;
		int[] map=new int[1+model.numberOfNodes];
		for(int i=1;i<=model.numberOfNodes;i++)
			if(model.node[i].inUse) {
				map[i]=++ix;
			}

		Model newModel=new Model(model.numberOfRegions,model.numberOfElements,ix,model.elType);

		for(int i=1;i<=model.numberOfElements;i++){
			int[] vertNumb=model.element[i].getVertNumb();
			for(int j=0;j<model.nElVert;j++)
				newModel.element[i].setVertNumb(j,map[vertNumb[j]]);
		}

		ix=0;
		for(int i=1;i<=model.numberOfNodes;i++){
			if(map[i]>0) 
				newModel.node[map[i]]=model.node[i];
		}

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			newModel.region[ir]=model.region[ir];
		}
		newModel.scaleFactor=model.scaleFactor;
		newModel.commonNodes=model.commonNodes;

		return newModel;

	}


	public void rotate(double rad){
		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model(bun);


		if(rad==0) return;

		rotate(model,rad,true);

	}

	public void rotate(String bun, double rad){

		Model model=new Model(bun);

		if(rad==0) return;

		rotate(model,rad,true);

	}
	public void rotate(Mat R, boolean write){
		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model(bun);
		rotate( model, R,  write);

	}
	public void rotate(Model model,Mat R, boolean write){



		for(int i=1;i<=model.numberOfNodes;i++){
			Vect v=model.node[i].getCoord();
			model.node[i].setCoord(R.mul(v));
		}

		if(write){
			String file = System.getProperty("user.dir") + "//bunRotated.txt";

			model.writeMesh(file);
		}
	}


	public void rotate(Mat R){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model(bun);

		for(int i=1;i<=model.numberOfNodes;i++){
			Vect v=model.node[i].getCoord();
			model.node[i].setCoord(R.mul(v));
		}

		String folder=new File(model.meshFilePath).getParentFile().getPath();
		String file = folder + "//bunRotated.txt";

		model.writeMesh(file);

	}

	public void rotate(Model model,double rad, boolean write){

		int dim=model.dim;
		Mat R=new Mat(dim,dim);
		if(dim==2) R=util.rotMat2D(rad);
		else
			R=util.rotEuler(new Vect(0,0,1), rad);

		for(int i=1;i<=model.numberOfNodes;i++){
			Vect v=model.node[i].getCoord();
			model.node[i].setCoord(R.mul(v));
		}

		if(write){
			String folder=new File(model.meshFilePath).getParentFile().getPath();
			String file = folder + "//bunRotated.txt";

			model.writeMesh(file);
		}
	}


	public void translate(Vect v){
		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model(bun);

		translate(model,v);

	}

	public void translate(Model model,Vect v){


		translate(model,v,true);

	}

	public void translate(Model model,Vect v, boolean write){


		for(int i=1;i<=model.numberOfNodes;i++){
			Vect z=model.node[i].getCoord().add(v);
			model.node[i].setCoord(z);
		}

		if(write){
			String folder=new File(model.meshFilePath).getParentFile().getPath();
			String file = folder+ "\\bunTranslated.txt";

			model.writeMesh(file);
		}
	}


	public void node4Del(){

		String bun=util.getFile(0);
		Model model=new Model(bun);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";
		DecimalFormat formatter= new DecimalFormat("0.0000000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(nodeFile)));		

			pwBun.println(model.numberOfNodes+" "+model.dim+" "+0+" "+0);

			for(int i=1;i<=model.numberOfNodes;i++){ 

				Vect z=model.node[i].getCoord();
				pwBun.print(i+" ");
				for(int j=0;j<model.dim;j++){

					pwBun.print(formatter.format(z.el[j]*model.scaleFactor)+" ");
				}
				pwBun.println();	

			}

			pwBun.close();

			util.pr(" node file for Delauny was written to "+nodeFile);
		}

		catch(IOException e){}



	}

	public void makeIPMmesh4th(){

		meshIPMrotor4th();

		meshIPMstator4th();

		String rtf = System.getProperty("user.dir") + "\\rotorMesh.txt";
		String stf  = System.getProperty("user.dir") + "\\statorMesh.txt";
		Model mt=new Model();
		combineTwoNode(mt, rtf, stf,0,.5);
		String bunm = System.getProperty("user.dir") + "\\motCombined.txt";
		mt.writeMesh(bunm);

	}

	public void makeIPMmeshFull(){

		meshIPMrotorFull();

		meshIPMstatorFull();

		String rtf = System.getProperty("user.dir") + "\\rotorMesh.txt";
		String stf  = System.getProperty("user.dir") + "\\statorMesh.txt";
		Model mt=new Model();
		combineTwoNode(mt, rtf, stf,22.5*PI/180,.5);
		String bunm = System.getProperty("user.dir") + "\\motCombined.txt";
		mt.writeMesh(bunm);

	}

	public void makeIPMmesh1deg(){

		meshIPMrotor1deg();

		meshIPMstator1deg();

		String rtf = System.getProperty("user.dir") + "\\rotorMesh.txt";
		String stf  = System.getProperty("user.dir") + "\\statorMesh.txt";
		Model mt=new Model();
		combineTwoNode(mt, rtf, stf,0,.5);
		String bunm = System.getProperty("user.dir") + "\\motCombined.txt";
		mt.writeMesh(bunm);

	}


	public void meshIPMrotor4th(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMrotor(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(7);


		wait (4*nNodes);

		reRegionIPMrotor8th();

		wait (5*nNodes);
		int[] regList={1,2,3,4,5,6};
		String bun = System.getProperty("user.dir") + "//rotor8th.txt";

		extractReg(bun,regList);


		wait (4*nNodes);

		bun=System.getProperty("user.dir") + "//extReg.txt";

		dropUnusedNodes(bun);

		wait (4*nNodes);



		bun = System.getProperty("user.dir") + "//noUnusedNodes.txt";

		Model mod=new Model(bun);
		extendFlip(mod,1);
		wait (4*nNodes);


		bun = System.getProperty("user.dir") + "//reflected.txt";

		double f1=PI/8;
		rotate(bun,f1);
		wait (4*nNodes);

		bun = System.getProperty("user.dir") + "//bunRotated.txt";
		rotExtendNfold(bun,1);

		bun = System.getProperty("user.dir") + "//extended.txt";
		String outBun = System.getProperty("user.dir") + "//rotorMesh.txt";

		wait (4*nNodes);

		reRegionIPMR4th(bun,outBun);

	}

	public void meshIPMrotorFull(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMrotor(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(7);


		wait (4*nNodes);

		reRegionIPMrotor8th();

		wait (5*nNodes);
		int[] regList={1,2,3,4,5,6};
		String bun = System.getProperty("user.dir") + "//rotor8th.txt";

		extractReg(bun,regList);


		wait (4*nNodes);

		bun=System.getProperty("user.dir") + "//extReg.txt";

		dropUnusedNodes(bun);

		wait (4*nNodes);

		double f1=PI/8;

		bun = System.getProperty("user.dir") + "//noUnusedNodes.txt";
		rotate(bun,-f1);

		wait (4*nNodes);

		bun = System.getProperty("user.dir") + "//bunRotated.txt";
		Model mod=new Model(bun);
		extendFlip(mod,1);
		wait (4*nNodes);

		bun = System.getProperty("user.dir") + "//reflected.txt";

		rotate(bun,f1);
		wait (4*nNodes);

		bun = System.getProperty("user.dir") + "//bunRotated.txt";
		rotExtendNfold(bun,7);

		bun = System.getProperty("user.dir") + "//extended.txt";
		String outBun = System.getProperty("user.dir") + "//rotorMesh.txt";

		wait (4*nNodes);


		reRegionIPMRFull(bun,outBun);


	}


	public void meshIPMrotor1deg(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMrotor1deg(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(7);


		wait (4*nNodes);

		reRegionIPMrotor8th();

		wait (5*nNodes);
		int[] regList={1,2,3,4,5,6};
		String bun = System.getProperty("user.dir") + "//rotor8th.txt";

		extractReg(bun,regList);


		wait (4*nNodes);

		bun=System.getProperty("user.dir") + "//extReg.txt";

		dropUnusedNodes(bun);

		wait (4*nNodes);

		double f1=PI/8;

		bun = System.getProperty("user.dir") + "//noUnusedNodes.txt";
		rotate(bun,f1);

		wait (4*nNodes);

		bun = System.getProperty("user.dir") + "//bunRotated.txt";

		rotExtendNfold(bun,1);

		bun = System.getProperty("user.dir") + "//extended.txt";
		String outBun = System.getProperty("user.dir") + "//rotorMesh.txt";

		wait (4*nNodes);


		reRegionIPMR2fold1deg(bun,outBun);

	}

	public void meshIPMrotor2xx(){

		//==  makes the rotor mesh starting from a full 8th reg segment
		int nNodes=10000;
		int dim=2;
		double scaleFactor=1000;
		double r1=16;
		double rout=54.4;
		double f0=-Math.PI/8;
		double f1=Math.PI/8;

		Vect[] P=new Vect[nNodes];
		int ix=0;

		P[ix++]=new Vect(2);

		double x0=45.5,x1=50.5;
		double y0=-14.5,y1=14.5;

		double eps=1.5;

		double gapHalf=.3;


		double r00=40;
		double r01=3.25;

		r1=10;
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),0, 8, -PI/8, PI/8, 1, 1);

		ix=this.addPointsOnArc(P, ix, new Vect(2),12, -PI/8, PI/8, 2);

		ix=this.addPointsOnArc(P, ix, new Vect(2),16, -PI/8, PI/8, 3);
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),2*r1,3*r1, -PI/8, PI/8,2,4);

		ix=this.addPointsOnArc(P, ix, new Vect(r00,0),1*r01,-PI, PI, 8);


		//  next two lines for having a finer mesh around the bold holes

		//ix=this.addPointsOnArc(P, ix, new Vect(r00,0),1.3*r01,-PI, PI, 8);
		//ix=this.addPointsOnArc(P, ix, new Vect(r00,0),.6*r01,-PI, PI, 8);


		ix=this.addPointsOnPath(P, ix, new Vect(r00,0));

		ix=this.addPointsOnArc(P, ix, new Vect(2),3.5*r1, PI/16, PI/8, 2);
		ix=this.addPointsOnArc(P, ix, new Vect(2),3.5*r1, -PI/8, -PI/18, 2);


		ix=this.addPointsOnArc(P, ix, new Vect(2),4.1*r1, PI/14, PI/8, 3);
		ix=this.addPointsOnArc(P, ix, new Vect(2),4.1*r1, -PI/8, -PI/14, 3);
		ix=this.addPointsOnArc(P, ix, new Vect(2),4.4*r1, PI/20, PI/8, 4);
		ix=this.addPointsOnArc(P, ix, new Vect(2),4.4*r1, -PI/8, -PI/20, 4);




		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),4.6*r1,5.1*r1, PI/9.5, PI/8, 2,1);
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),4.6*r1,5.1*r1, -PI/8, -PI/9.5,2, 1);

		ix=this.addPointsOnArc(P, ix, new Vect(2), rout-1, -PI/8, PI/8,30);
		//ix=this.addPointsOnArc(P, ix, new Vect(2), rout-1.4, 0, PI/10,12);
		//ix=this.addPointsOnArc(P, ix, new Vect(2), rout-1.4, -PI/10, PI/10,12);
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2), rout,rout+gapHalf, -PI/8, PI/8,3,45);
		Vect vx=new Vect(52,0);

		//ix=this.addPointsOnArc(P, ix, util.rotMat2D(-PI/19).mul(vx),util.rotMat2D(PI/19).mul(vx), PI/12,8);



		double eps2=.1;


		ix=addPointsOnArea(P,ix,x0,x0+eps2,y0+eps2,y1-eps2,2,12);
		ix=addPointsOnArea(P,ix,x1-eps2,x1,y0+eps2,y1-eps2,2,12);
		ix=addPointsOnArea(P,ix,x0+eps2,x1-eps2,y0,y0+eps2,4,2);
		ix=addPointsOnArea(P,ix,x0+eps2,x1-eps2,y1-eps2,y1,4,2);

		ix=addPointsOnArea(P,ix,x0,x0+eps2,y0,y0+eps2,1,1);
		ix=addPointsOnArea(P,ix,x0,x0+eps2,y1-eps2,y1,1,1);
		ix=addPointsOnArea(P,ix,x1-eps2,x1,y1-eps2,y1,1,1);
		ix=addPointsOnArea(P,ix,x1-eps2,x1,y0,y0+eps2,1,1);

		ix=addPointsOnArea(P,ix,x0+eps2,x1-eps2,y0+eps2,y1-eps2,4,12);


		nNodes=ix;

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";
		DecimalFormat formatter= new DecimalFormat("0.0000000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(nodeFile)));		

			pwBun.println(nNodes+" "+dim+" "+0+" "+0);

			for(int i=0;i<nNodes;i++){ 

				pwBun.print((i+1)+" ");
				for(int j=0;j<dim;j++){

					pwBun.print(formatter.format(P[i].el[j]/scaleFactor)+" ");
				}
				pwBun.println();	

			}

			pwBun.close();

			util.pr(" node file for Delauny was written to "+nodeFile);
		}


		catch(IOException e){}

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(7);


		wait (4*nNodes);

		reRegionIPMrotor8th();

		wait (4*nNodes);

		String bun = System.getProperty("user.dir") + "//rotor8th.txt";
		dropUnusedNodes(bun);

		wait (4*nNodes);

		f1=PI/8;

		bun = System.getProperty("user.dir") + "//noUnusedNodes.txt";
		rotate(bun,f1);

		wait (4*nNodes);

		bun = System.getProperty("user.dir") + "//bunRotated.txt";

		rotExtendNfold(bun,1);
		bun = System.getProperty("user.dir") + "//extended.txt";
		String outBun = System.getProperty("user.dir") + "//rotorMesh.txt";

		wait (4*nNodes);


		reRegionIPMR4th(bun,outBun);

	}



	private void writeNodesForDelaunay(Vect[] P,int nNodes,double scaleFactor){
		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		int dim=P[0].length;

		DecimalFormat formatter= new DecimalFormat("0.0000000");

		try{
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(nodeFile)));		

			pwBun.println(nNodes+" "+dim+" "+0+" "+0);

			for(int i=0;i<nNodes;i++){ 

				pwBun.print((i+1)+" ");
				for(int j=0;j<dim;j++){

					pwBun.print(formatter.format(P[i].el[j]/scaleFactor)+" ");
				}
				pwBun.println();	

			}

			pwBun.close();

			util.pr(" node file for Delauny was written to "+nodeFile);
		}


		catch(IOException e){}

	}


	public void meshIPMstator4th(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMstator(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";


		wait (5*nNodes);

		Delaunay(nodeFile);

		wait (5*nNodes);

		meshFromDel(7);

		wait (5*nNodes);

		String bun = System.getProperty("user.dir") + "//meshFromDel.txt";		

		reRegionIPMSHalfSlotTip(bun);

		wait (5*nNodes);
		int[] regList={1,2,3,4,5,6};
		bun = System.getProperty("user.dir") + "//slotH.txt";

		extractReg(bun,regList);


		wait (5*nNodes);

		bun=System.getProperty("user.dir") + "//extReg.txt";


		dropUnusedNodes(bun);

		wait (5*nNodes);

		Model mod=new Model(bun);
		extendFlip(mod,1);
		wait(5*nNodes);
		bun = System.getProperty("user.dir") + "//reflected.txt";
		double f1=5*1.0/180*Math.PI;
		rotate(bun,f1);
		wait (5*nNodes);
		bun = System.getProperty("user.dir") + "//bunRotated.txt";
		wait (5*nNodes);

		rotExtendNfold(bun,8);
		bun = System.getProperty("user.dir") + "//extended.txt";
		String outBun = System.getProperty("user.dir") + "//statorMesh.txt";
		wait (5*nNodes);
		reRegionIPMS8fold(bun,outBun);

	}

	public void meshIPMstatorFull(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMstator(P);
		//nNodes=makeNodeIPMstatorRoughQuad(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";


		wait (5*nNodes);

		Delaunay(nodeFile);

		wait (5*nNodes);

		meshFromDel(7);

		wait (5*nNodes);

		String bun = System.getProperty("user.dir") + "//meshFromDel.txt";		

		reRegionIPMSHalfSlotTip(bun);

		wait (5*nNodes);
		int[] regList={1,2,3,4,5,6};
		bun = System.getProperty("user.dir") + "//slotH.txt";




		extractReg(bun,regList);


		wait (5*nNodes);

		bun=System.getProperty("user.dir") + "//extReg.txt";

		dropUnusedNodes(bun);

		wait (5*nNodes);

		Model mod=new Model(bun);
		extendFlip(mod,1);
		wait(5*nNodes);
		bun = System.getProperty("user.dir") + "//reflected.txt";
		double f1=5*1.0/180*Math.PI;
		rotate(bun,f1);
		wait (5*nNodes);
		bun = System.getProperty("user.dir") + "//bunRotated.txt";
		wait (5*nNodes);

		rotExtendNfold(bun,35);
		bun = System.getProperty("user.dir") + "//extended.txt";
		String outBun = System.getProperty("user.dir") + "//statorMesh.txt";
		wait (5*nNodes);
		reRegionIPMSFull(bun,outBun);

	}


	public void meshIPMframe(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMframe(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";



		wait (5*nNodes);

		Delaunay(nodeFile);

		wait (5*nNodes);

		meshFromDel(2);

		wait (5*nNodes);

		String bun = System.getProperty("user.dir") + "//meshFromDel.txt";		

		reRegionFrame(bun);

		bun = System.getProperty("user.dir") + "//frameH.txt";

		wait (5*nNodes);
		int[] regList={1};

		//

		extractReg(bun,regList);


		wait (5*nNodes);

		bun=System.getProperty("user.dir") + "//extReg.txt";


		dropUnusedNodes(bun);

		wait (5*nNodes);

		bun = System.getProperty("user.dir") + "//noUnusedNodes.txt";

		rotExtendNfold(bun,71);


	}

	public void meshIPMCap(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMCap(P);
		//nNodes=makeNodeIPMCapRoughQ(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";



		wait (5*nNodes);

		Delaunay(nodeFile);

		wait (5*nNodes);

		meshFromDel(2);

		wait (5*nNodes);

		String bun = System.getProperty("user.dir") + "//meshFromDel.txt";		

		reRegionCap(bun);

		bun = System.getProperty("user.dir") + "//capH.txt";

		wait (5*nNodes);
		int[] regList={1};

		//

		extractReg(bun,regList);


		wait (5*nNodes);

		bun=System.getProperty("user.dir") + "//extReg.txt";


		dropUnusedNodes(bun);

		wait (5*nNodes);

		bun = System.getProperty("user.dir") + "//noUnusedNodes.txt";


		rotExtendNfold(bun,7);


	}

	public void meshIPMCapQ(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMCap(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";



		wait (5*nNodes);

		Delaunay(nodeFile);

		wait (5*nNodes);

		meshFromDel(2);

		wait (5*nNodes);

		String bun = System.getProperty("user.dir") + "//meshFromDel.txt";		

		reRegionCap(bun);

		bun = System.getProperty("user.dir") + "//capH.txt";

		wait (5*nNodes);
		int[] regList={1};

		//

		extractReg(bun,regList);


		wait (5*nNodes);

		bun=System.getProperty("user.dir") + "//extReg.txt";


		dropUnusedNodes(bun);

		wait (5*nNodes);

		bun = System.getProperty("user.dir") + "//noUnusedNodes.txt";

		rotExtendNfold(bun,7);


	}

	public void meshIPMLegs(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMLegs(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";



		wait (5*nNodes);

		Delaunay(nodeFile);

		wait (5*nNodes);

		meshFromDel(4);

		wait (5*nNodes);

		String bun = System.getProperty("user.dir") + "//meshFromDel.txt";		

		reRegionLegs(bun);


		bun = System.getProperty("user.dir") + "//legs.txt";
		int[] regList={1,2,3};

		//

		extractReg(bun,regList);


		wait (5*nNodes);
		bun = System.getProperty("user.dir") + "//extReg.txt";		
		dropUnusedNodes(bun);

		wait (5*nNodes);





	}


	public void meshIPMstator1deg(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeIPMstator1deg(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";


		wait (5*nNodes);

		Delaunay(nodeFile);

		wait (5*nNodes);

		meshFromDel(6);

		wait (5*nNodes);

		String bun = System.getProperty("user.dir") + "//meshFromDel.txt";		

		reRegionIPMSHalfSlotTip(bun);

		wait (5*nNodes);
		int[] regList={1,2,3,4,5};
		bun = System.getProperty("user.dir") + "//slotH.txt";

		extractReg(bun,regList);


		wait (5*nNodes);

		bun=System.getProperty("user.dir") + "//extReg.txt";


		dropUnusedNodes(bun);

		wait (5*nNodes);

		Model mod=new Model(bun);
		extendFlip(mod,1);
		wait(5*nNodes);
		bun = System.getProperty("user.dir") + "//reflected.txt";
		double f1=5*1.0/180*Math.PI;
		rotate(bun,f1);
		wait (5*nNodes);
		bun = System.getProperty("user.dir") + "//bunRotated.txt";
		wait (5*nNodes);

		rotExtendNfold(bun,8);
		bun = System.getProperty("user.dir") + "//extended.txt";
		String outBun = System.getProperty("user.dir") + "//statorMesh.txt";
		wait (5*nNodes);
		reRegionIPMS8fold(bun,outBun);


	}



	private int  makeNodeIPMstator(Vect[] P){

		double rm=54.4+.3;
		double rin=rm+.3;
		double rout1=175.0/2;
		double rring=rout1+4.5;
		double rout2=rout1+30;
		double rmc=67.298;
		double rs=3.88;
		double Rs=rout1-10.47-3.88;

		double rt=Rs+3.88;


		double f1=5*1.0/180*Math.PI;

		double f2=f1-1.5/55;;

		double f3=f1-2.4/55;
		double f4=f1-rs/Rs;

		int ix=0;

		P[ix++]=new Vect(2);


		Vect[] crn=new Vect[50];



		crn[0]=new Vect(rin,0);
		crn[1]=new Vect(Rs,0);
		crn[2]=new Vect(rout1,0);
		crn[3]=new Vect(rout2,0);

		crn[4]=util.rotMat2D(f1).mul(crn[3]);
		crn[5]=util.rotMat2D(f1).mul(crn[2]);
		crn[6]=util.rotMat2D(f1).mul(new Vect(rt,0));
		crn[7]=util.rotMat2D(f1).mul(crn[1]);
		crn[8]=util.rotMat2D(f1).mul(crn[0]);
		crn[9]=util.rotMat2D(f2).mul(new Vect(rin+.7,0));
		crn[10]=util.rotMat2D(f2).mul(new Vect(rin,0));
		crn[11]=util.rotMat2D(f3).mul(new Vect(rin+1.03,0));
		crn[12]=util.rotMat2D(f4).mul(crn[1]);
		//crn[12]=util.rotMat2D(f2).mul(crn[9]);

		/*		for(int i=0;i<crn.length;i++){
				if(crn[i]!=null)
				ix=addPointsOnPath(P,ix,crn[i]);
			}*/


		ix=addPointsOnPath(P,ix,crn[9],crn[10],2);
		ix=addPointsOnPath(P,ix,crn[9],crn[11],2);


		ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,0,f3,5);
		ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,f3,.5*(f2+f3),1);
		ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,f1-.666*(f1-f2),f1,2);


		ix=addPoint(P,ix,util.rotMat2D(1.01*f2).mul(new Vect(rm+.2,0)));



		ix=addPointsOnArc(P,ix,new Vect(2), rin,0,f3,5);
		ix=addPointsOnArc(P,ix,new Vect(2), rin,f3,f2,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rin,f2,f1,3);





		ix=addPointsOnAreaArc(P,ix,new Vect(2), rm,rm+.1,0,f1,1,10);



		ix=addPointsOnArc(P,ix,crn[12],crn[6],.6*PI,3);
		ix=addPointsOnArc(P,ix,util.rotMat2D(-f3/3).mul(crn[12]),crn[6].times(1.03),.6*PI,3);

		ix=addPointsOnArc(P,ix,new Vect(2),rout1-6,0,f1,3);
		ix=addPointsOnArc(P,ix,new Vect(2),rout1-3,0,f1,3);
		//	ix=addPointsOnArc(P,ix,new Vect(2),rout1-9,0,f1,2);
		//ix=addPointsOnArc(P,ix,crn[1].times(1.08),crn[6].times(1.09),.6*PI,2);

		ix=addPoint(P,ix,new Vect(rout1-8,0));
		ix=addPoint(P,ix,new Vect(rout1-12,0));

		//))))))))))) fine tooth
		/*ix=addPointsOnArc(P,ix,new Vect(2),rin+.2,0,f3,6);
			ix=addPointsOnArc(P,ix,new Vect(2),rin+.2,f3,f2,3);
			ix=addPointsOnArc(P,ix,new Vect(2),rin+.5,0,f3,6);
			ix=addPointsOnArc(P,ix,new Vect(2),rin+.5,f3,f2,3);*/
		//)))))))
		ix=addPointsOnArc(P,ix,new Vect(2),rin+.35,0,f3,4);
		ix=addPointsOnArc(P,ix,new Vect(2),rin+.35,f3,f2,2);
		//ix=addPointsOnArc(P,ix,new Vect(2),rin+.35,f2,f1,3);

		ix=addPointsOnArc(P,ix,new Vect(2), rin+.7,f2,f1,3);

		ix=addPointsOnArc(P,ix,new Vect(2), rin+1.03,0,f3,3);
		//	ix=addPointsOnArc(P,ix,new Vect(2), rin+1.03,f3,f1,2);





		double D=crn[12].sub(crn[11]).norm();
		double rx;

		rx=2.2-1.03;
		double fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));
		ix=addPointsOnArc(P,ix,new Vect(2), rin+2.2,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rin+2.2,fx,f1,2);




		rx=4-1.03;
		fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));

		ix=addPointsOnArc(P,ix,new Vect(2), rin+4,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rin+4,fx,f1,2);

		rx=8-1.03;
		fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));

		ix=addPointsOnArc(P,ix,new Vect(2), rin+8,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rin+8,fx,f1,2);

		rx=16-1.03;
		fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));
		ix=addPointsOnArc(P,ix,new Vect(2), rin+16,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rin+16,fx,f1,2);

		rx=rmc-crn[11].el[0];
		fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));
		ix=addPointsOnArc(P,ix,new Vect(2), rmc,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rmc,fx,f1,2);



		ix=addPointsOnArc(P,ix,new Vect(2), rout1,0,f1,3);

		ix=addPointsOnAreaArc(P,ix,new Vect(2), rout1,rring,0,f1,3,3);
		ix=addPointsOnArc(P,ix,new Vect(2), .8*rring+.2*rout2,0,f1,2);
		ix=addPointsOnArc(P,ix,new Vect(2), .4*rout1+.6*rout2,0,f1,2);
		//ix=addPointsOnArc(P,ix,new Vect(2), rring,0,f1,3);
		//	ix=addPointsOnArc(P,ix,new Vect(2), .7*rout1+.3*rout2,0,f1,3);
		//ix=addPointsOnArc(P,ix,new Vect(2), .5*rout1+.5*rout2,0,f1,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rout2,0,f1,2);


		return ix;
	}


	private int  makeNodeIPMframe(Vect[] P){


		double r1=87.5;
		double r2=r1+4.5;


		double f1=5*1.0/180*Math.PI;


		int ix=0;

		P[ix++]=new Vect(2);

		ix=addPointsOnAreaArc(P,ix,new Vect(2), r1,r2,0,f1,3,3);
		//	ix=addPointsOnAreaArc(P,ix,new Vect(2), r2,r2+10,f3,f4,3,3);


		return ix;
	}



	private int  makeNodeIPMCap(Vect[] P){


		double r0=16.0;

		double r1=87.5;
		double r2=r1+4.5;


		double f1=45*1.0/180*Math.PI;


		int ix=0;

		P[ix++]=new Vect(2);

		ix=addPointsOnAreaArc(P,ix,new Vect(2), r1,r2,0,f1,3,27);
		ix=addPointsOnArc(P,ix,new Vect(2),.95*r1,0,f1,15);
		//	ix=addPointsOnArc(P,ix,new Vect(2),.9*r1,0,f1,6);
		ix=addPointsOnArc(P,ix,new Vect(2),.85*r1,0,f1,6);
		ix=addPointsOnArc(P,ix,new Vect(2),.7*r1,0,f1,4);
		ix=addPointsOnArc(P,ix,new Vect(2),.5*r1,0,f1,3);
		ix=addPointsOnArc(P,ix,new Vect(2),.35*r1,0,f1,3);
		//ix=addPointsOnArc(P,ix,new Vect(2),.3*r1,0,f1,2);
		ix=addPointsOnArc(P,ix,new Vect(2),r0,0,f1,2);


		//	ix=addPointsOnAreaArc(P,ix,new Vect(2), r2,r2+10,f3,f4,3,3);


		return ix;
	}


	private int  makeNodeIPMCapRoughQ(Vect[] P){


		double r0=16.0;

		double r1=87.5;
		double r2=r1+4.5;


		double f1=2.5*1.0/180*Math.PI;


		int ix=0;

		P[ix++]=new Vect(2);

		ix=addPointsOnAreaArc(P,ix,new Vect(2), r0,r2,0,f1,12,1);



		return ix;
	}

	private int  makeNodeIPMLegs(Vect[] P){


		//double r1=87.5;
		double r2=87.5;
		double h1=5;
		double h2=33;
		double x1=25;


		double f1=5*1.0/180*Math.PI;
		double f3=3*PI/2-40*1.0/180*Math.PI;;
		double f4=3*PI/2+40*1.0/180*Math.PI;;


		int ix=0;

		P[ix++]=new Vect(2);

		Vect[] z=new Vect[50];

		z[0]=util.rotMat2D(f3).mul(new Vect(r2,0));
		z[1]=util.rotMat2D(f3).mul(new Vect(r2+h2,0));
		z[2]=z[1].add(new Vect(0,-10));
		z[3]=z[2].add(new Vect(x1,0));
		z[4]=z[3].add(new Vect(0,5));
		//double d=-r2*Math.cos(f3)*(r2+h1)/r2-15;

		z[9]=util.rotMat2D(f4).mul(new Vect(r2,0));
		z[8]=util.rotMat2D(f4).mul(new Vect(r2+h2,0));
		z[7]=z[8].add(new Vect(0,-10));
		z[6]=z[7].add(new Vect(-x1,0));
		z[5]=z[6].add(new Vect(0,5));

		/*	z[5]=z[4].add(new Vect(2*d,0));
			z[6]=z[5].add(new Vect(0,-5));
			z[7]=z[6].add(new Vect(15,0));
			z[8]=z[7].add(new Vect(0,10));
			z[9]=util.rotMat2D(f4).mul(new Vect(r2,0));*/

		ix=addPointsOnAreaArc(P,ix,new Vect(2), r2,r2+h1,0,2*PI,3,3*72);
		int kk=0;
		for(int i=0;i<z.length;i++)
			if(z[i]!=null) kk++;
		Vect[] z2=new Vect[kk];

		for(int i=0;i<z2.length;i++)
			z2[i]=z[i].deepCopy();

		int[] N=new int[kk];

		for(int i=0;i<z2.length;i++){
			N[i]=1;
			if(i==4) N[i]=22;
		}
		N[0]=1;
		N[1]=1;
		N[2]=2;
		N[6]=2;
		N[7]=1;
		N[8]=1;
		ix=addPointsOnPath(P,ix,z2,N);

		ix=addPointsOnArc(P,ix,new Vect(2),r2+2.2*h1,f3,f3+PI/9,7);
		ix=addPointsOnArc(P,ix,new Vect(2),r2+4*h1,f3,f3+PI/22,2);

		ix=addPointsOnArc(P,ix,new Vect(2),r2+2.2*h1,f4,f4-PI/9,7);
		ix=addPointsOnArc(P,ix,new Vect(2),r2+4*h1,f4,f4-PI/22,2);


		//	ix=addPointsOnArc(P,ix,new Vect(2),r2+2.1*h1,0,PI+50*PI/180,72);
		ix=addPointsOnArc(P,ix,new Vect(2),r2+2.5*h1,-45*PI/180,PI+45*PI/180,90);
		ix=addPointsOnArc(P,ix,new Vect(2),r2+5*h1,-45*PI/180,PI+45*PI/180,50);
		ix=addPointsOnArc(P,ix,new Vect(2),r2+5*h1,1.5*PI-30*PI/180,1.5*PI+30*PI/180,10);
		//ix=addPointsOnArc(P,ix,new Vect(2),z[2].norm(),-60*PI/180,PI+60*PI/180,50);
		/*	ix=addPointsOnArc(P,ix,new Vect(2),r2+4*h1,0,1.5*PI-PI/4,60);
			ix=addPointsOnArc(P,ix,new Vect(2),r2+4*h1,1.5*PI+PI/4,2*PI,60);*/
		ix=addPointsOnArc(P,ix,new Vect(2),z[2].norm(),-50*PI/180,PI+50*PI/180,60);
		ix=addPointsOnArc(P,ix,new Vect(2),z[2].norm(),1.5*PI-30*PI/180,1.5*PI+30*PI/180,10);

		ix=addPointsOnPath(P,ix,z[1],z[4],2);
		ix=addPointsOnPath(P,ix,z[5],z[8],2);

		return ix;
	}

	private int  makeNodeIPMstatorRoughQuad(Vect[] P){



		double rm=54.4+.3;
		double rin=rm+.3;
		double rout1=175.0/2;
		double rring=rout1+4.5;
		double rout2=92;
		double rmc=67.298;
		double rs=3.88;
		double Rs=rout1-10.47-3.88;

		double rt=Rs+3.88;


		double f1=5*1.0/180*Math.PI;

		double f2=f1-1.5/55;;

		double f3=f1-2.4/55;
		double f4=f1-rs/Rs;

		int ix=0;

		P[ix++]=new Vect(2);


		Vect[] crn=new Vect[50];



		crn[0]=new Vect(rin,0);
		crn[1]=new Vect(Rs,0);
		crn[2]=new Vect(rout1,0);
		crn[3]=new Vect(rout2,0);

		crn[4]=util.rotMat2D(f1).mul(crn[3]);
		crn[5]=util.rotMat2D(f1).mul(crn[2]);
		crn[6]=util.rotMat2D(f1).mul(new Vect(rt,0));
		crn[7]=util.rotMat2D(f1).mul(crn[1]);
		crn[8]=util.rotMat2D(f1).mul(crn[0]);
		crn[9]=util.rotMat2D(f2).mul(new Vect(rin+.7,0));
		crn[10]=util.rotMat2D(f2).mul(new Vect(rin,0));
		crn[11]=util.rotMat2D(f3).mul(new Vect(rin+1.03,0));
		crn[12]=util.rotMat2D(f4).mul(crn[1]);
		//crn[12]=util.rotMat2D(f2).mul(crn[9]);

		for(int i=0;i<crn.length;i++){
			if(crn[i]!=null)
				ix=addPointsOnPath(P,ix,crn[i]);
		}


		ix=addPointsOnPath(P,ix,crn[9],crn[10],1);
		ix=addPointsOnPath(P,ix,crn[9],crn[11],1);


		ix=addPointsOnPath(P,ix,crn[11],crn[12],5);
		//ix=addPointsOnPath(P,ix,new Vect(rin+1.03,crn[12].el[1]/2),new Vect(rt,crn[12].el[1]/2),5);
		ix=addPointsOnPath(P,ix,new Vect(rin+1.03,0),new Vect(Rs,0),5);
		//	ix=addPointsOnPath(P,ix,util.rotMat2D(-f3/2).mul(crn[11]),util.rotMat2D(-f3/2).mul(crn[12]),5);
		//ix=addPointsOnPath(P,ix,util.rotMat2D(-f3).mul(crn[11]),util.rotMat2D(-f3/2).mul(crn[12]),5);
		//	ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,0,f3,2);
		/*ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,f3,.5*(f2+f3),1);
			ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,f1-.666*(f1-f2),f1,1);*/





		ix=addPointsOnArc(P,ix,new Vect(2), rin,0,f3,1);
		//	ix=addPointsOnAreaArc(P,ix,new Vect(2),crn[11].norm(),Rs-5,0,f3,8,2);
		ix=addPointsOnArc(P,ix,new Vect(2),Rs+6,0,f1,2);
		ix=addPointsOnArc(P,ix,new Vect(2),Rs+10,0,f1,2);
		ix=addPointsOnArc(P,ix,new Vect(2),rout1,0,f1,2);
		//ix=addPointsOnAreaArc(P,ix,new Vect(2),Rs+8.5,rout1,0,f1,1,2);
		ix=addPointsOnAreaArc(P,ix,new Vect(2),rout1,rring,0,f1,1,2);


		//ix=addPointsOnPath(P,ix,util.rotMat2D(.4*f3).mul(new Vect(rt,0)));
		ix=addPointsOnPath(P,ix,new Vect(rt,0));
		//	ix=addPointsOnArc(P,ix,new Vect(2), rin,f2,f1,3);





		//ix=addPointsOnAreaArc(P,ix,new Vect(2), rm,rm+.1,0,f1,1,10);



		ix=addPointsOnArc(P,ix,crn[12],crn[6],.6*PI,3);
		//	ix=addPointsOnArc(P,ix,util.rotMat2D(-f3/3).mul(crn[12]),crn[6].times(1.03),.6*PI,3);

		//	ix=addPointsOnArc(P,ix,new Vect(2),rout1-6,0,f1,3);
		//	ix=addPointsOnArc(P,ix,new Vect(2),rout1-3,0,f1,3);
		//	ix=addPointsOnArc(P,ix,new Vect(2),rout1-9,0,f1,2);
		//ix=addPointsOnArc(P,ix,crn[1].times(1.08),crn[6].times(1.09),.6*PI,2);

		//	ix=addPoint(P,ix,new Vect(rout1-8,0));
		//ix=addPoint(P,ix,new Vect(rout1-12,0));

		//))))))))))) fine tooth
		/*ix=addPointsOnArc(P,ix,new Vect(2),rin+.2,0,f3,6);
			ix=addPointsOnArc(P,ix,new Vect(2),rin+.2,f3,f2,3);
			ix=addPointsOnArc(P,ix,new Vect(2),rin+.5,0,f3,6);
			ix=addPointsOnArc(P,ix,new Vect(2),rin+.5,f3,f2,3);*/
		//)))))))
		/*	ix=addPointsOnArc(P,ix,new Vect(2),rin+.35,0,f3,4);
			ix=addPointsOnArc(P,ix,new Vect(2),rin+.35,f3,f2,2);
			//ix=addPointsOnArc(P,ix,new Vect(2),rin+.35,f2,f1,3);

			ix=addPointsOnArc(P,ix,new Vect(2), rin+.7,f2,f1,3);

			ix=addPointsOnArc(P,ix,new Vect(2), rin+1.03,0,f3,3);
			//	ix=addPointsOnArc(P,ix,new Vect(2), rin+1.03,f3,f1,2);





			double D=crn[12].sub(crn[11]).norm();
			double rx;

			rx=2.2-1.03;
			double fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));
			ix=addPointsOnArc(P,ix,new Vect(2), rin+2.2,0,fx,2);
			ix=addPointsOnArc(P,ix,new Vect(2), rin+2.2,fx,f1,2);




			rx=4-1.03;
			fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));

			ix=addPointsOnArc(P,ix,new Vect(2), rin+4,0,fx,2);
			ix=addPointsOnArc(P,ix,new Vect(2), rin+4,fx,f1,2);

			rx=8-1.03;
			fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));

			ix=addPointsOnArc(P,ix,new Vect(2), rin+8,0,fx,2);
			ix=addPointsOnArc(P,ix,new Vect(2), rin+8,fx,f1,2);

			rx=16-1.03;
			fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));
			ix=addPointsOnArc(P,ix,new Vect(2), rin+16,0,fx,2);
			ix=addPointsOnArc(P,ix,new Vect(2), rin+16,fx,f1,2);

			rx=rmc-crn[11].el[0];
			fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));
			ix=addPointsOnArc(P,ix,new Vect(2), rmc,0,fx,2);
			ix=addPointsOnArc(P,ix,new Vect(2), rmc,fx,f1,2);



			ix=addPointsOnArc(P,ix,new Vect(2), rout1,0,f1,3);

			ix=addPointsOnAreaArc(P,ix,new Vect(2), rout1,rring,0,f1,3,3);
			ix=addPointsOnArc(P,ix,new Vect(2), .8*rring+.2*rout2,0,f1,2);
			ix=addPointsOnArc(P,ix,new Vect(2), .4*rout1+.6*rout2,0,f1,2);
			//ix=addPointsOnArc(P,ix,new Vect(2), rring,0,f1,3);
			//	ix=addPointsOnArc(P,ix,new Vect(2), .7*rout1+.3*rout2,0,f1,3);
			//ix=addPointsOnArc(P,ix,new Vect(2), .5*rout1+.5*rout2,0,f1,2);
			ix=addPointsOnArc(P,ix,new Vect(2), rout2,0,f1,2);
		 */

		return ix;

	}

	private int  makeNodeIPMstator1deg(Vect[] P){

		double rm=54.4+.3;
		double rin=rm+.3;
		double rout1=175.0/2;
		double rout2=rout1+30;
		double rmc=67.298;
		double rs=3.88;
		double Rs=rout1-10.47-3.88;

		double rt=Rs+3.88;


		double f1=5*1.0/180*Math.PI;

		double f2=f1-1.5/55;;

		double f3=f1-2.4/55;
		double f4=f1-rs/Rs;

		int ix=0;

		P[ix++]=new Vect(2);


		Vect[] crn=new Vect[50];



		crn[0]=new Vect(rin,0);
		crn[1]=new Vect(Rs,0);
		crn[2]=new Vect(rout1,0);
		crn[3]=new Vect(rout2,0);

		crn[4]=util.rotMat2D(f1).mul(crn[3]);
		crn[5]=util.rotMat2D(f1).mul(crn[2]);
		crn[6]=util.rotMat2D(f1).mul(new Vect(rt,0));
		crn[7]=util.rotMat2D(f1).mul(crn[1]);
		crn[8]=util.rotMat2D(f1).mul(crn[0]);
		crn[9]=util.rotMat2D(f2).mul(new Vect(rin+.7,0));
		crn[10]=util.rotMat2D(f2).mul(new Vect(rin,0));
		crn[11]=util.rotMat2D(f3).mul(new Vect(rin+1.03,0));
		crn[12]=util.rotMat2D(f4).mul(crn[1]);
		//crn[12]=util.rotMat2D(f2).mul(crn[9]);

		for(int i=0;i<crn.length;i++){
			if(crn[i]!=null)
				ix=addPointsOnPath(P,ix,crn[i]);
		}


		ix=addPointsOnPath(P,ix,crn[9],crn[10],2);
		ix=addPointsOnPath(P,ix,crn[9],crn[11],2);





		/*ix=addPointsOnArc(P,ix,new Vect(2), rin,0,f3,3);
			ix=addPointsOnArc(P,ix,new Vect(2), rin,f3,f2,2);
			ix=addPointsOnArc(P,ix,new Vect(2), rin,f2,f1,2);*/

		ix=addPointsOnArc(P,ix,new Vect(2),  rm+.2,0,1.2*f3,3);

		ix=addPointsOnArc(P,ix,new Vect(2),  rin,0,1.2*f3,3);

		ix=addPointsOnAreaArc(P,ix,new Vect(2), rm,rm+.1,0,f1,1,5);


		//ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,0,f1,5);

		/*	ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,0,f3,3);
			ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,f3,.5*(f2+f3),1);
			ix=addPointsOnArc(P,ix,new Vect(2), rm+.2,f1-.666*(f1-f2),f1,2);
		 */

		//ix=addPoint(P,ix,util.rotMat2D(1.01*f2).mul(new Vect(rm+.2,0)));




		ix=addPointsOnArc(P,ix,crn[12],crn[6],.6*PI,3);
		ix=addPointsOnArc(P,ix,util.rotMat2D(-f3/2).mul(crn[12]),crn[6].times(1.03),.6*PI,3);


		ix=addPointsOnArc(P,ix,crn[1].times(1.08),crn[6].times(1.09),.6*PI,2);


		//ix=addPointsOnArc(P,ix,new Vect(2),rin+.35,0,f2,3);
		//	ix=addPointsOnArc(P,ix,new Vect(2),rin+.35,f3,f2,2);
		//ix=addPointsOnArc(P,ix,new Vect(2),rin+.35,f2,f1,3);

		ix=addPointsOnArc(P,ix,new Vect(2), rin+.7,f2,f1,2);

		ix=addPointsOnArc(P,ix,new Vect(2), rin+1.03,0,f3,2);
		//	ix=addPointsOnArc(P,ix,new Vect(2), rin+1.03,f3,f1,2);





		double D=crn[12].sub(crn[11]).norm();
		double rx;

		rx=2.2-1.03;
		double fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));
		ix=addPointsOnArc(P,ix,new Vect(2), rin+2.2,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rin+2.2,fx,f1,2);




		rx=4-1.03;
		fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));

		ix=addPointsOnArc(P,ix,new Vect(2), rin+4,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rin+4,fx,f1,2);

		rx=8-1.03;
		fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));

		ix=addPointsOnArc(P,ix,new Vect(2), rin+8,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rin+8,fx,f1,2);

		rx=16-1.03;
		fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));
		ix=addPointsOnArc(P,ix,new Vect(2), rin+16,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rin+16,fx,f1,2);

		rx=rmc-crn[11].el[0];
		fx=util.getAng(crn[11].times(1-rx/D).add(crn[12].times(rx/D)));
		ix=addPointsOnArc(P,ix,new Vect(2), rmc,0,fx,2);
		ix=addPointsOnArc(P,ix,new Vect(2), rmc,fx,f1,2);



		ix=addPointsOnArc(P,ix,new Vect(2), rout1,0,f1,2);
		ix=addPointsOnArc(P,ix,new Vect(2), .85*rout1+.15*rout2,0,f1,2);
		ix=addPointsOnArc(P,ix,new Vect(2), .7*rout1+.3*rout2,0,f1,2);
		ix=addPointsOnArc(P,ix,new Vect(2), .4*rout1+.6*rout2,0,f1,1);


		return ix;
	}


	private int  makeNodeIPMrotor(Vect[] P){


		double r1=16;
		double rout=54.4;


		int ix=0;

		P[ix++]=new Vect(2);

		double x0=45.5,x1=50.5;
		double y0=0,y1=14.5;

		double gapHalf=.3;


		double r00=40;
		double r01=3.25;

		r1=10;
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),0, 8, 0, PI/8, 1, 1);

		//ix=this.addPointsOnArc(P, ix, new Vect(2),12, 0, PI/8, 1);

		ix=this.addPointsOnArc(P, ix, new Vect(2),16,0, PI/8, 2);
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),2.2*r1,3.0*r1,0, PI/8,1,2);

		ix=this.addPointsOnArc(P, ix, new Vect(r00,0),1*r01,0, PI, 4);
		ix=this.addPointsOnArc(P, ix, new Vect(r00,0),2*r01,PI/3, PI, 4);


		//  next two lines for having a finer mesh around the bold holes

		//ix=this.addPointsOnArc(P, ix, new Vect(r00,0),1.3*r01,-PI, PI, 8);
		//ix=this.addPointsOnArc(P, ix, new Vect(r00,0),.6*r01,-PI, PI, 8);


		ix=this.addPointsOnPath(P, ix, new Vect(r00,0));

		ix=this.addPointsOnArc(P, ix, new Vect(2),3.5*r1, PI/16, PI/8, 2);

		ix=this.addPointsOnArc(P, ix, new Vect(2),4.0*r1, PI/14, PI/8, 2);
		ix=this.addPointsOnArc(P, ix, new Vect(2),4.4*r1, PI/16, PI/8, 5);




		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),4.6*r1,5.3*r1, PI/9.5, PI/8, 4,2);

		ix=this.addPointsOnArc(P, ix, new Vect(2), rout-.5, 0, PI/8,22);
		//ix=this.addPointsOnArc(P, ix, new Vect(2), rout-1.2, 0, PI/8,12);
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2), rout,rout+gapHalf,0, PI/8,3,45);
		//"?
		ix=this.addPointsOnArc(P, ix, new Vect(52,0),util.rotMat2D(PI/19).mul(new Vect(52,0)), PI/12,4);



		double eps2=.1;


		ix=addPointsOnArea(P,ix,x0,x0+eps2,y0,y1-eps2,2,8);
		ix=addPointsOnArea(P,ix,x1-eps2,x1,y0,y1-eps2,2,8);
		ix=addPointsOnArea(P,ix,x0+eps2,x1-eps2,y1-eps2,y1,4,2);

		ix=addPointsOnArea(P,ix,x0,x0+eps2,y1-eps2,y1,1,1);
		ix=addPointsOnArea(P,ix,x1-eps2,x1,y1-eps2,y1,1,1);

		ix=addPointsOnArea(P,ix,x0+eps2,x1-eps2,y0,y1-eps2,4,8);

		eps2=.3;
		ix=addPointsOnArc(P,ix,new Vect(x1+eps2,y1+eps2),new Vect(x0-eps2,y1+eps2),PI/4,5);



		return ix;

	}


	private int  makeNodeIPMrotor1deg(Vect[] P){



		double r1=16;
		double rout=54.4;
		int ix=0;

		P[ix++]=new Vect(2);

		double x0=45.5,x1=50.5;
		double y0=-14.5,y1=14.5;

		double eps=1.5;

		double gapHalf=.3;


		double r00=40;
		double r01=3.25;

		r1=10;
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),0, 8, -PI/8, PI/8, 1, 1);

		ix=this.addPointsOnArc(P, ix, new Vect(2),12, -PI/8, PI/8, 2);

		ix=this.addPointsOnArc(P, ix, new Vect(2),16, -PI/8, PI/8, 3);
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),2*r1,3*r1, -PI/8, PI/8,2,4);

		ix=this.addPointsOnArc(P, ix, new Vect(r00,0),1*r01,-PI, PI, 8);


		//  next two lines for having a finer mesh around the bold holes

		//ix=this.addPointsOnArc(P, ix, new Vect(r00,0),1.3*r01,-PI, PI, 8);
		//ix=this.addPointsOnArc(P, ix, new Vect(r00,0),.6*r01,-PI, PI, 8);


		ix=this.addPointsOnPath(P, ix, new Vect(r00,0));

		ix=this.addPointsOnArc(P, ix, new Vect(2),3.5*r1, PI/16, PI/8, 2);
		ix=this.addPointsOnArc(P, ix, new Vect(2),3.5*r1, -PI/8, -PI/18, 2);


		ix=this.addPointsOnArc(P, ix, new Vect(2),4.1*r1, PI/14, PI/8, 3);
		ix=this.addPointsOnArc(P, ix, new Vect(2),4.1*r1, -PI/8, -PI/14, 3);
		ix=this.addPointsOnArc(P, ix, new Vect(2),4.4*r1, PI/20, PI/8, 4);
		ix=this.addPointsOnArc(P, ix, new Vect(2),4.4*r1, -PI/8, -PI/20, 4);




		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),4.6*r1,5.1*r1, PI/9.5, PI/8, 2,1);
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2),4.6*r1,5.1*r1, -PI/8, -PI/9.5,2, 1);

		ix=this.addPointsOnArc(P, ix, new Vect(2), rout-1, -PI/8, PI/8,30);
		//ix=this.addPointsOnArc(P, ix, new Vect(2), rout-1.4, 0, PI/10,12);
		//ix=this.addPointsOnArc(P, ix, new Vect(2), rout-1.4, -PI/10, PI/10,12);
		ix=this.addPointsOnAreaArc(P, ix, new Vect(2), rout,rout+gapHalf, -PI/8, PI/8,3,45);
		Vect vx=new Vect(52,0);

		//ix=this.addPointsOnArc(P, ix, util.rotMat2D(-PI/19).mul(vx),util.rotMat2D(PI/19).mul(vx), PI/12,8);



		double eps2=.1;


		ix=addPointsOnArea(P,ix,x0,x0+eps2,y0+eps2,y1-eps2,2,12);
		ix=addPointsOnArea(P,ix,x1-eps2,x1,y0+eps2,y1-eps2,2,12);
		ix=addPointsOnArea(P,ix,x0+eps2,x1-eps2,y0,y0+eps2,4,2);
		ix=addPointsOnArea(P,ix,x0+eps2,x1-eps2,y1-eps2,y1,4,2);

		ix=addPointsOnArea(P,ix,x0,x0+eps2,y0,y0+eps2,1,1);
		ix=addPointsOnArea(P,ix,x0,x0+eps2,y1-eps2,y1,1,1);
		ix=addPointsOnArea(P,ix,x1-eps2,x1,y1-eps2,y1,1,1);
		ix=addPointsOnArea(P,ix,x1-eps2,x1,y0,y0+eps2,1,1);

		ix=addPointsOnArea(P,ix,x0+eps2,x1-eps2,y0+eps2,y1-eps2,4,12);

		return ix;

	}




	private void wait(int secs){
		try {
			new Thread().sleep(secs);
		} catch (InterruptedException exception) {
			// TODO Auto-generated catch-block stub.
			exception.printStackTrace();
		}
	}


	private int addPointsOnArea(Vect[] P, int ix,double x0,double x1,double y0,double y1,int Nx,int Ny){

		double x,y;
		double dx=(x1-x0)/Nx;
		double dy=(y1-y0)/Ny;


		for(int kx=0;kx<=Nx;kx++){

			for(int ky=0;ky<=Ny;ky++){

				x=x0+kx*dx;
				y=y0+ky*dy;

				Vect z=new Vect(x,y);
				P[ix++]=z;
			}
		}

		return ix;

	}

	private int addPointsOnPath(Vect[] P,int ix, Vect z1,Vect z2,int N){





		Vect dv=z2.sub(z1).times(1.0/N);
		for(int k=0;k<=N;k++){
			P[ix++]=z1.add(dv.times(k));
		}

		return ix;


	}

	private int addPointsOnPath(Vect[] P,int ix, Vect z){
		P[ix++]=z.deepCopy();
		return ix;
	}

	private int addPointsOnPath(Vect[] P,int ix, Vect[] z){
		int[] N=new int[z.length];
		for(int i=0;i<N.length;i++)
			N[i]=1;

		if(z.length==0) return ix;
		if(z.length==1){
			P[ix++]=z[0];
			return ix;
		}
		for(int i=0;i<z.length-1;i++){

			Vect dv=z[i+1].sub(z[i]).times(1.0/N[i]);
			for(int k=0;k<=N[i];k++){
				P[ix++]=z[i].add(dv.times(k));
			}
		}

		return ix;

	}


	private int addPointsOnPath(Vect[] P,int ix, Vect[] z,int[] N){

		if(z.length==0) return ix;
		if(z.length==1){
			P[ix++]=z[0];
			return ix;
		}
		for(int i=0;i<z.length-1;i++){

			Vect dv=z[i+1].sub(z[i]).times(1.0/N[i]);
			for(int k=0;k<=N[i];k++){
				P[ix++]=z[i].add(dv.times(k));
			}
		}

		return ix;

	}

	private int addPointsOnAreaArc(Vect[] P,int ix, Vect center,double r0,double r1,double phi0,double phi1,int Nr,int Nphi){



		double x,y;

		double dr=(r1-r0)/Nr;
		double df=(phi1-phi0)/Nphi;

		for(int kr=0;kr<=Nr;kr++){
			double r=r0+kr*dr;

			for(int kf=0;kf<=Nphi;kf++){
				double phi=phi0+kf*df;
				x=r*Math.cos(phi);
				y=r*Math.sin(phi);
				Vect z=new Vect(x,y).add(center);

				P[ix++]=z;
			}
		}


		return ix;

	}


	private int addPoint(Vect[] P,int ix, Vect z){


		P[ix++]=z.deepCopy();
		return ix;

	}

	private int addPointsOnArc(Vect[] P,int ix, Vect center,double r,double phi0,double phi1,int N){

		double df=(phi1-phi0)/N;

		for(int i=0;i<=N;i++){
			double phi=phi0+i*df;
			P[ix++]=new Vect(r*Math.cos(phi),r*Math.sin(phi)).add(center);
		}

		return ix;

	}

	private int addPointsOnArc(Vect[] P,int ix, Vect z1,Vect z2,  double theta,int N){

		double d=z2.sub(z1).norm()/2;
		double r=d/Math.sin(theta/2);
		double h=abs(r*Math.cos(theta/2));
		Vect c1=z1.add(z2).times(.5);
		Vect dir=util.rotMat2D(PI/2).mul(z2.sub(c1)).normalized();
		Vect center=c1.add(dir.times(h));

		double phi0=util.getAng(z1.sub(center));

		//double phi1=util.getAng(z2.sub(center));

		double df=theta/N;

		for(int i=0;i<=N;i++){
			double phi=phi0+i*df;
			P[ix++]=new Vect(r*Math.cos(phi),r*Math.sin(phi)).add(center);
		}

		return ix;

	}


	public void Delaunay(String nodeFile){



		String line;
		try {
			Process p = Runtime.getRuntime().exec("triangle "+ nodeFile);
			BufferedReader bri = new BufferedReader
					(new InputStreamReader(p.getInputStream()));
			BufferedReader bre = new BufferedReader
					(new InputStreamReader(p.getErrorStream()));
			while ((line = bri.readLine()) != null) {
				System.out.println(line);
			}
			bri.close();
			while ((line = bre.readLine()) != null) {
				System.out.println(line);
			}
			bre.close();
			p.waitFor();
			System.out.println("Done.");
		}
		catch (Exception err) {
			err.printStackTrace();
		}

	}

	public void Delaunay(){

		String nodeFile=util.getFile(0);


		String line;
		try {
			Process p = Runtime.getRuntime().exec("triangle "+ nodeFile);
			BufferedReader bri = new BufferedReader
					(new InputStreamReader(p.getInputStream()));
			BufferedReader bre = new BufferedReader
					(new InputStreamReader(p.getErrorStream()));
			while ((line = bri.readLine()) != null) {
				System.out.println(line);
			}
			bri.close();
			while ((line = bre.readLine()) != null) {
				System.out.println(line);
			}
			bre.close();
			p.waitFor();
			System.out.println("Done.");
		}
		catch (Exception err) {
			err.printStackTrace();
		}

	}
	public void meshFromDel(){
		meshFromDel(1);
	}

	public void meshFromDel(int nReg) {
		String nodeFile=System.getProperty("user.dir") + "//node4Del.1.node";
		String elementFile=System.getProperty("user.dir") + "//node4Del.1.ele";

		Model model=new Model();

		try{
			BufferedReader br1 = new BufferedReader(new FileReader(nodeFile));
			BufferedReader br2 = new BufferedReader(new FileReader(elementFile));
			String[] s1,s2;
			String line1,line2;
			line1=br1.readLine();
			s1=line1.split(regex);	
			line2=br2.readLine();
			s2=line2.split(regex);	
			int nNodes=Integer.parseInt(s1[0]);
			int nElements=Integer.parseInt(s2[0]);
			model.alloc(nReg, nElements, nNodes, "triangle");
			model.scaleFactor=1;
			for(int i=1;i<=model.numberOfElements;i++){

				line2=br2.readLine();
				s2=line2.split(regex);	
				int m=1;
				if(s2[0].equals("")) m=2;
				for(int k=0;k<model.nElVert;k++){
					int nn=Integer.parseInt(s2[k+m]);

					model.element[i].setVertNumb(k, nn);
				}

			}

			for(int i=1;i<=model.numberOfNodes;i++){

				line1=br1.readLine();
				s1=line1.split(regex);	
				int m=1;
				if(s1[0].equals("")) m=2;
				Vect v=new Vect(model.dim);
				for(int k=0;k<model.dim;k++){
					v.el[k]=Double.parseDouble(s1[m+k]);

				}
				model.node[i].setCoord(v);
			}
			model.region[1].setFirstEl(1);
			model.region[1].setLastEl(nElements);
			model.region[1].setName("iron");

			for(int ir=2;ir<=nReg;ir++){
				model.region[ir].setFirstEl(1);
				model.region[ir].setLastEl(0);
				model.region[ir].setName("reg"+ir);
			}

			br1.close();
			br2.close();
			String bunFilePath = System.getProperty("user.dir") + "//meshFromDel.txt";

			model.writeMesh(bunFilePath);


		}

		catch(IOException e){}



	}

	private double[] getTabbedData(String line){
		String[] sp=line.split(regex);	
		int L=sp.length;
		double[] v=new double[L];
		for( int p=0;p<L;p++)
			v[p]=Double.parseDouble(sp[p]);

		return v;
	}

	public void extractReg(int[] regList)
	{


		String bun=util.getFile(0);

		extractReg(bun,regList);

	}

	public void extractReg(String bun,int[] regList)
	{


		Model model=new Model(bun);

		int nEls=0;

		for(int ir=0;ir<regList.length;ir++)
			nEls+=model.region[regList[ir]].getNumbElements();

		int nReg=regList.length;
		String folder=new File(bun).getParentFile().getPath();
		String bunFilePath = folder + "//extReg.txt";


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
			pwBun.println(nEls);

			pwBun.println("//Number_of_Region");
			pwBun.println(nReg);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			double r0=model.node[1].getCoord().norm();

			for(int k=0;k<regList.length;k++){
				int ir=regList[k];
				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
					int[] vertNumb=model.element[i].getVertNumb();
					for(int j=0;j<model.nElVert;j++)
						pwBun.print(vertNumb[j]+",");
					pwBun.println();
				}
			}

			for(int i=1;i<=model.numberOfNodes;i++){ 
				Vect xyz;;

				xyz=model.node[i].getCoord();
				for(int j=0;j<model.dim;j++){
					pwBun.print(formatter.format((xyz.el[j])*model.scaleFactor)+" ,");
				}
				pwBun.println();	

			}


			int npr=1;
			for(int k=0;k<regList.length;k++){
				int ir=regList[k];
				pwBun.println(npr+","+(npr-1+model.region[ir].getNumbElements())+","+model.region[ir].getName());
				npr+=model.region[ir].getNumbElements();
			}


			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);


			pwBun.close();
		}
		catch(IOException e){}

	}


	public void extractReg(int ir){

		String bun=util.getFile(0);

		Model model=new Model(bun);

		int nEls=0;

		nEls=model.region[ir].getNumbElements();

		String folder=new File(bun).getParentFile().getPath();

		String bunFilePath = folder + "//extReg.txt";



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
			pwBun.println(nEls);

			pwBun.println("//Number_of_Region");
			pwBun.println(1);
			pwBun.println("//Factor");
			pwBun.println(model.scaleFactor);

			double r0=model.node[1].getCoord().norm();

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				//if(i>model.region[ir].getFirstEl()) continue;
				int[] vertNumb=model.element[i].getVertNumb();
				for(int j=0;j<model.nElVert;j++)
					pwBun.print(vertNumb[j]+",");
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



			pwBun.println(1+","+model.region[ir].getNumbElements()+","+model.region[ir].getName());


			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+bunFilePath);


			pwBun.close();
		}
		catch(IOException e){}



	}

	public void extractReg( double r1, double r2,double t1,double t2,double z1, double z2){

		String bun=util.getFile(0);
		extractReg(bun,r1,r2,t1,t2,z1,z2,null);
	}


	public void extractReg( double r1, double r2,double t1,double t2,int[] regList){

		String bun=util.getFile(0);
		extractReg(bun,r1,r2,t1,t2,-1e10,1e10,regList);
	}

	public void extractReg( double r1, double r2,double t1,double t2){

		String bun=util.getFile(0);
		extractReg(bun,r1,r2,t1,t2,-1e10,1e10,null);
	}

	public void extractReg(String bun, double r1, double r2,double t1,double t2,double z1,double z2,int[] regList){

		if(t1<0) t1+=2*PI;

		Model model1=new Model(bun);

		int nReg=model1.numberOfRegions+1;
		int nEls=model1.numberOfElements;
		int nNodes=model1.numberOfNodes;
		String type=model1.elType;
		
		 model1=null;

		 System.gc();
		
		Model model2=new Model();
		model2.alloc(nReg,nEls,nNodes,type);
		model2.loadMesh(bun);

		double z=0;
		int[] list;
		if(regList==null){
			list=new int[model2.numberOfRegions];
			for(int k=0;k<list.length;k++) list[k]=k+1;

		}
		else{
			list=regList;
		}

		for(int k=0;k<list.length;k++){

			int ir=list[k];
			for(int i=model2.region[ir].getFirstEl();i<=model2.region[ir].getLastEl();i++){

			

				Vect c=model2.getElementCenter(i);
				if(model2.dim==3)
					z=c.el[2];
				c=c.v2();		
				double tt=util.getAng(c);
				//model.element[i].setRegion(ir);

				double r=c.norm();
				//	if(c.norm()>.0876)model.element[i].setRegion(2);
				if((tt<t1 || tt>t2 )|| (r<r1 || r>r2 ) || (model2.dim==3 &&(z<z1 || z>z2 )))
				{

					model2.element[i].setRegion(model2.numberOfRegions);
				}


			}
		}


		reRegionGroupEls(model2);
		String folder=new File(bun).getParentFile().getPath();
		String bunFilePath =folder+ "//bunReReg.txt";


		model2.writeMesh(bunFilePath);


	}
	
	public void elimitateZeroAreaElements(double err)
	{


		String bun=util.getFile(0);

		elimitateZeroAreaElements(bun,err);

	}
	
	public void elimitateZeroAreaElements(String bun,double err){


		Model model1=new Model(bun);
		Model model2=new Model(model1.numberOfRegions+1,model1.numberOfElements,model1.numberOfNodes,model1.elType);
		model2.loadMesh(bun);

		double z=0;
		int[] list;
			list=new int[model2.numberOfRegions];
			for(int k=0;k<list.length;k++) list[k]=k+1;




			for(int k=0;k<list.length;k++){

			int ir=list[k];
			for(int i=model2.region[ir].getFirstEl();i<=model2.region[ir].getLastEl();i++){

				
				double vol=0;
				if(model2.dim==2)
					vol=model2.getElementArea(i);
				else
					vol=model2.getElementVolume(i);
				
				if(vol<0)	
					System.out.println(" volume of element "+i+" is minus ...:"+vol);
				
				if(vol<err) model2.element[i].setRegion(model2.numberOfRegions);

			}
		}


		reRegionGroupEls(model2);
		String folder=new File(bun).getParentFile().getPath();
		String bunFilePath =folder+ "//bunReReg.txt";


		model2.writeMesh(bunFilePath);
	}
	




	public void cut2D(double h1, double h2)
	{
		
		String bun=util.getFile(0);

		if(bun==null) return;
		
			
		//double dh=0.001;
	
		Model model=new Model(bun);


		String e2type="quadrangle";
		
		if(model.elType.equals("prism")) e2type="triangle";
		
		int ne2=0;
		for(int ir=1;ir<=model.numberOfRegions;ir++)	{	
		//	double dhx=dh;
		//if(ir<6) dhx=dh2;
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)	{

			double ze=model.getElementCenter(i).el[2];
			if(ze>h2 ||ze<h1 ) continue;
			ne2++;
			
			}
		}
		
		Model m2d=new Model(model.numberOfRegions,ne2,model.numberOfNodes,e2type);

		m2d.scaleFactor=model.scaleFactor;
		
		for(int i=1;i<=model.numberOfNodes;i++)		
			m2d.node[i].setCoord(model.node[i].getCoord().v2());


		int ix=0;
		int nex=0;

		for(int ir=1;ir<=model.numberOfRegions;ir++)		{
			m2d.region[ir].setName(model.region[ir].getName());
			
			//double dhx=dh;
		//	if(ir<6) dhx=dh2;
			if(ir==1)
			m2d.region[ir].setFirstEl(1);
			else
			m2d.region[ir].setFirstEl(m2d.region[ir-1].getLastEl()+1);
			
			nex=0;
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)	{

			double ze=model.getElementCenter(i).el[2];
			if(ze>h2 ||ze<h1) continue;
			ix++;
			nex++;
			m2d.element[ix].setRegion(model.element[i].getRegion());
			int[] e2vn=new int[model.nElVert/2];
			int[] vn=model.element[i].getVertNumb();
			for(int j=0;j<model.nElVert/2;j++)
				e2vn[j]=vn[j+model.nElVert/2];
			m2d.element[ix].setVertNumb(e2vn);
		}
			//util.pr(nex);
			m2d.region[ir].setLastEl(m2d.region[ir].getFirstEl()+nex-1);

		}

		String bunFilePath = System.getProperty("user.dir") + "//ext2Dfrom3D.txt";

		m2d.writeMesh(bunFilePath);
	}


	public Model extendRotationalSym(){

		String bun=util.getFile();
		Model forth=new Model(bun);

		boolean[] nodeIn=new boolean[forth.numberOfNodes+1];
		for(int i=1;i<=forth.numberOfElements;i++){
			int[] vertNum=forth.element[i].getVertNumb();
			for(int j=0;j<forth.nElVert;j++)
				nodeIn[vertNum[j]]=true;
		}


		double tmax=-10,tmin=10;
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(!nodeIn[i]) continue;
			Vect z=forth.node[i].getCoord();
			double t;
			if(z.norm()==0) t=0;
			else if(z.el[0]>=0) t=Math.atan(z.el[1]/z.el[0]);
			else t=Math.atan(z.el[1]/z.el[0])+Math.PI;

			if(t>tmax) tmax=t;
			if(t<tmin) tmin=t;

		}
		double alpha1=tmin;
		double alpha2=tmax;


		int nRegions=forth.numberOfRegions;

		int nElements=2*forth.numberOfElements;

		int ix=0;
		boolean[] onBound1=new boolean[forth.numberOfNodes+1];
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(!nodeIn[i]) continue;
			Vect z=forth.node[i].getCoord();
			double t;
			if(z.norm()==0) t=alpha1;
			else if(z.el[0]>=0) t=Math.atan(z.el[1]/z.el[0]);
			else t=Math.atan(z.el[1]/z.el[0])+Math.PI;
			if( Math.abs(t-alpha1)<1e-6){
				ix++;
				onBound1[i]=true;
			}
		}

		Vect xx=new Vect(ix);
		int[] nx=new int[ix];
		ix=0;
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(onBound1[i]) {
				nx[ix]=i;
				xx.el[ix++]=forth.node[i].getCoord().norm();
			}
		}


		Vect yy=new Vect(ix);
		int[] ny=new int[ix];
		ix=0;
		boolean[] onBound2=new boolean[forth.numberOfNodes+1];
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(!nodeIn[i]) continue;
			Vect z=forth.node[i].getCoord();
			double t;
			if(z.norm()==0) t=alpha2;
			else if(z.el[0]>=0 ) t=Math.atan(z.el[1]/z.el[0]);
			else t=Math.atan(z.el[1]/z.el[0])+Math.PI;
			if( Math.abs(t-alpha2)<1e-6){
				ix++;
				onBound2[i]=true;
			}
		}
		ix=0;
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(onBound2[i]) {
				ny[ix]=i;
				yy.el[ix++]=forth.node[i].getCoord().norm();
			}
		}


		int[] indX=xx.bubble();
		int[] indY=yy.bubble();


		int[] map1=new int[forth.numberOfNodes+1];
		int[] map2=new int[forth.numberOfNodes+1];
		for(int i=0;i<ix;i++){
			map1[nx[indX[i]]]=ny[indY[i]];
			map2[ny[indY[i]]]=nx[indX[i]];
		}


		int nNodes=2*forth.numberOfNodes-xx.length;


		Model half=new Model(nRegions,nElements,nNodes,forth.elType);
		half.numberOfElements=nElements;

		int n1=1,N;
		for(int i=1;i<=forth.numberOfRegions;i++){

			if(i>1)
				n1=half.region[i-1].getLastEl()+1;
			N=forth.region[i].getNumbElements();
			half.region[i].setFirstEl(n1);
			half.region[i].setLastEl(n1+2*N-1);
			half.region[i].setMaterial(forth.region[i].getMaterial());
			half.region[i].setMur(forth.region[i].getMur());
			half.region[i].setM(forth.region[i].getM());
			half.region[i].setJ(forth.region[i].getJ());
			half.region[i].setSigma(forth.region[i].getSigma());
			half.region[i].setPois(forth.region[i].getPois());
			half.region[i].setYng(forth.region[i].getYng());

		}



		double theta=alpha2-alpha1;

		Mat R;

		for(int i=1;i<=forth.numberOfNodes;i++)	
			half.node[i].setCoord(forth.node[i].getCoord());


		R=util.rotMat2D(theta);
		for(int i=1;i<=forth.numberOfNodes;i++)	

			half.node[i+forth.numberOfNodes].setCoord(R.mul(forth.node[i].getCoord()));

		ix=0;

		for(int ir=1;ir<=forth.numberOfRegions;ir++)		{

			int ne=forth.region[ir].getNumbElements();
			ix=0;
			for(int i=forth.region[ir].getFirstEl();i<=forth.region[ir].getLastEl();i++)	{

				ix++;
				int p=1;	
				int[] vertNumb=new int[3];
				p=half.region[ir].getFirstEl();
				vertNumb=forth.element[i].getVertNumb();

				for(int j=0;j<2;j++){
					int ie=p-1+j*ne+ix;
					if(j==0)
						half.element[ie].setVertNumb(vertNumb);	
					else{
						for(int k=0;k<forth.nElVert;k++){
							if(!onBound1[vertNumb[k]] && !onBound2[vertNumb[k]])
								half.element[ie].setVertNumb(k,vertNumb[k]+j*forth.numberOfNodes);
							else if(onBound1[vertNumb[k]])
								half.element[ie].setVertNumb(k,map1[vertNumb[k]]);
							else{

								Vect z=half.node[vertNumb[k]+j*forth.numberOfNodes].getCoord();
								double t;
								if(z.norm()==0) t=alpha1;
								else t=util.getAng(z);
								if(Math.abs(t-alpha1)<1e-6 || Math.abs(t-2*Math.PI-alpha1)<1e-6){
									half.element[ie].setVertNumb(k,map2[vertNumb[k]]);
								}
								else {
									half.element[ie].setVertNumb(k,vertNumb[k]+j*forth.numberOfNodes);
								}


							}
						}


					}
				}
			}

		}

		String mesh = System.getProperty("user.dir") + "//extended.txt";

		half.scaleFactor=forth.scaleFactor;

		half.writeMesh(mesh);



		return half;
	}

	public Model extendOneFold(){

		String bun=util.getFile();
		Model forth=new Model(bun);

		boolean[] nodeIn=new boolean[forth.numberOfNodes+1];
		for(int i=1;i<=forth.numberOfElements;i++){
			int[] vertNum=forth.element[i].getVertNumb();
			for(int j=0;j<forth.nElVert;j++)
				nodeIn[vertNum[j]]=true;
		}


		int m=2;
		double tmax=-10,tmin=10;
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(!nodeIn[i]) continue;
			Vect z=forth.node[i].getCoord();
			double t;
			if(z.norm()==0) t=0;
			else if(z.el[0]>=0) t=Math.atan(z.el[1]/z.el[0]);
			else t=Math.atan(z.el[1]/z.el[0])+Math.PI;

			if(t>tmax) tmax=t;
			if(t<tmin) tmin=t;

		}
		double alpha1=tmin;
		double alpha2=tmax;


		int nRegions=forth.numberOfRegions;

		int nElements=m*forth.numberOfElements;

		int ix=0;
		boolean[] onBound1=new boolean[forth.numberOfNodes+1];
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(!nodeIn[i]) continue;
			Vect z=forth.node[i].getCoord();
			double t;
			if(z.norm()==0) t=alpha1;
			else if(z.el[0]>=0) t=Math.atan(z.el[1]/z.el[0]);
			else t=Math.atan(z.el[1]/z.el[0])+Math.PI;
			if( Math.abs(t-alpha1)<1e-6){
				ix++;
				onBound1[i]=true;
			}
		}

		Vect xx=new Vect(ix);
		int[] nx=new int[ix];
		ix=0;
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(onBound1[i]) {
				nx[ix]=i;
				xx.el[ix++]=forth.node[i].getCoord().norm();
			}
		}


		Vect yy=new Vect(ix);
		int[] ny=new int[ix];
		ix=0;
		boolean[] onBound2=new boolean[forth.numberOfNodes+1];
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(!nodeIn[i]) continue;
			Vect z=forth.node[i].getCoord();
			double t;
			if(z.norm()==0) t=alpha2;
			else if(z.el[0]>=0 ) t=Math.atan(z.el[1]/z.el[0]);
			else t=Math.atan(z.el[1]/z.el[0])+Math.PI;
			if( Math.abs(t-alpha2)<1e-6){
				ix++;
				onBound2[i]=true;
			}
		}
		ix=0;
		for(int i=1;i<=forth.numberOfNodes;i++)
		{
			if(onBound2[i]) {
				ny[ix]=i;
				yy.el[ix++]=forth.node[i].getCoord().norm();
			}
		}


		int[] indX=xx.bubble();
		int[] indY=yy.bubble();


		int[] map1=new int[forth.numberOfNodes+1];
		int[] map2=new int[forth.numberOfNodes+1];
		for(int i=0;i<ix;i++){
			map1[nx[indX[i]]]=ny[indY[i]];
			map2[ny[indY[i]]]=nx[indX[i]];
		}


		int nNodes=m*forth.numberOfNodes;


		Model half=new Model(nRegions,nElements,nNodes,forth.elType);
		half.numberOfElements=nElements;

		int n1=1,N;
		for(int i=1;i<=forth.numberOfRegions;i++){

			if(i>1)
				n1=half.region[i-1].getLastEl()+1;
			N=forth.region[i].getNumbElements();
			half.region[i].setFirstEl(n1);
			half.region[i].setLastEl(n1+2*N-1);
			half.region[i].setMaterial(forth.region[i].getMaterial());
			half.region[i].setMur(forth.region[i].getMur());
			half.region[i].setM(forth.region[i].getM());
			half.region[i].setJ(forth.region[i].getJ());
			half.region[i].setSigma(forth.region[i].getSigma());
			half.region[i].setPois(forth.region[i].getPois());
			half.region[i].setYng(forth.region[i].getYng());

		}




		double theta=alpha2-alpha1;

		Mat R;

		for(int i=1;i<=forth.numberOfNodes;i++)	
			half.node[i].setCoord(forth.node[i].getCoord());

		for(int j=1;j<m;j++)	{
			R=util.rotMat2D(j*theta);
			for(int i=1;i<=forth.numberOfNodes;i++)	

				half.node[i+j*forth.numberOfNodes].setCoord(R.mul(forth.node[i].getCoord()));
		}

		ix=0;

		for(int ir=1;ir<=forth.numberOfRegions;ir++)		{

			int ne=forth.region[ir].getNumbElements();
			ix=0;
			for(int i=forth.region[ir].getFirstEl();i<=forth.region[ir].getLastEl();i++)	{

				ix++;
				int p=1;	
				int[] vertNumb=new int[3];
				p=half.region[ir].getFirstEl();
				vertNumb=forth.element[i].getVertNumb();

				for(int j=0;j<m;j++){
					int ie=p-1+j*ne+ix;
					if(j==0)
						half.element[ie].setVertNumb(vertNumb);	
					else{
						for(int k=0;k<forth.nElVert;k++){
							if(!onBound1[vertNumb[k]] && !onBound2[vertNumb[k]])
								half.element[ie].setVertNumb(k,vertNumb[k]+j*forth.numberOfNodes);
							else if(onBound1[vertNumb[k]])
								half.element[ie].setVertNumb(k,map1[vertNumb[k]]);
							else{
								Vect z=half.node[vertNumb[k]+j*forth.numberOfNodes].getCoord();
								double t;
								if(z.norm()==0) t=alpha1;
								else t=util.getAng(z);

								if(Math.abs(t-alpha1)<1e-6)
									half.element[ie].setVertNumb(k,map2[vertNumb[k]]);
								else {
									half.element[ie].setVertNumb(k,vertNumb[k]+j*forth.numberOfNodes);
								}


							}
						}


					}
				}
			}

		}

		String mesh = System.getProperty("user.dir") + "//extended.txt";

		half.scaleFactor=forth.scaleFactor;

		half.writeMesh(mesh);



		return half;
	}

	public Model rotExtendNfold( int ns){
		String bun=util.getFile();

		Model slice=new Model(bun);
		util.pr(slice.numberOfElements);
		return rotExtendNfold(slice, ns);

	}
	

	public Model rotExtendNfold(Model slice, int ns){
		int nT=ns+1;

		setSliceBounds(slice);

		setNodeOnMotorBound(slice);


		double alpha1=slice.alpha1;
		double alpha2=slice.alpha2;

		
		double theta=alpha2-alpha1;

		if(nT*theta>2*PI+1e-6) {nT=(int)(2*PI/theta)+1; ns=nT-1;}

		int[][] mapping=getBoundaryNodesMap(slice,2,3);

		int L=mapping.length;
		int nRegions=slice.numberOfRegions;

		int nElements=nT*slice.numberOfElements;
		int nNodes=nT*slice.numberOfNodes;

		Model extenSlice=new Model();
		extenSlice.alloc(nRegions,nElements,nNodes,slice.elType);

		for(int i=1;i<=slice.numberOfRegions;i++){


			extenSlice.region[i].setMaterial(slice.region[i].getMaterial());
			extenSlice.region[i].setName(slice.region[i].getName());
			extenSlice.region[i].setMur(slice.region[i].getMur());
			extenSlice.region[i].setM(slice.region[i].getM());
			extenSlice.region[i].setJ(slice.region[i].getJ());
			extenSlice.region[i].setSigma(slice.region[i].getSigma());
			extenSlice.region[i].setPois(slice.region[i].getPois());
			extenSlice.region[i].setYng(slice.region[i].getYng());

		}


		int[] map=new int[slice.numberOfNodes+1];
		int[] mapEnd=new int[slice.numberOfNodes+1];
		int pivot=0;
		for(int i=0;i<L;i++){
	/*		map[mapping[i][0]]=mapping[i][1];
			mapEnd[mapping[i][1]]=mapping[i][0];
			if(mapping[i][0]==mapping[i][1])
				pivot=mapping[i][0];*/
		}


		Mat R;
		int nNodes1=slice.numberOfNodes;
		int nElem=slice.numberOfElements;

		int[] map2=new int[1+extenSlice.numberOfNodes];
		for(int j=0;j<=ns;j++){

			if(slice.dim==2)
				R=util.rotMat2D(j*theta);
			else {
				Mat R2D=util.rotMat2D(j*theta);
				R=new Mat(slice.dim,slice.dim);
				for(int m=0;m<2;m++)
					for(int n=0;n<2;n++)
						R.el[m][n]=R2D.el[m][n];

				R.el[2][2]=1;

			}

			for(int i=1;i<=slice.numberOfNodes;i++)	{
				int nn=i+j*nNodes1;
				extenSlice.node[nn].setCoord(R.mul(slice.node[i].getCoord()));
				if(slice.node[i].hasF()){
					extenSlice.node[nn].F=R.mul(slice.node[i].F);
				}
				if (j>0 && map[i]>0 )
					if(i==pivot)
						map2[nn]=pivot;
					else
						map2[nn]=map[i]+(j-1)*nNodes1;
			}
		}


		for(int j=0;j<=ns;j++){
			boolean touch=(abs((j+1)*theta-2*PI)<1e-4);
			

			for(int ir=1;ir<=extenSlice.numberOfRegions;ir++){

				for(int i=slice.region[ir].getFirstEl();i<=slice.region[ir].getLastEl();i++){	

					int[] vertNumb=slice.element[i].getVertNumb();
					for(int k=0;k<slice.nElVert;k++){
						int n0=vertNumb[k];
						int nn=n0+j*nNodes1;
						int ne=i+j*nElem;

						if(j==0 || (map[vertNumb[k]]==0 && mapEnd[vertNumb[k]]==0)) 
							extenSlice.element[ne].setVertNumb(k,nn);
						else if(n0==pivot)
							extenSlice.element[ne].setVertNumb(k,pivot);
						else{

							if(map[n0]>0)
								extenSlice.element[ne].setVertNumb(k,map[n0]+(j-1)*nNodes1);
							else if(touch &&mapEnd[n0]>0)
								extenSlice.element[ne].setVertNumb(k,mapEnd[n0]);
							else	 
								extenSlice.element[ne].setVertNumb(k,nn);
						}	
						


					}

					extenSlice.element[i+j*nElem].setRegion(ir);

					



				}
			}
		}

		int[][] mapr=new int[extenSlice.numberOfRegions+1][2];

		int ix=0;
		int[][] vn=new int[extenSlice.numberOfElements+1][extenSlice.nElVert];

		for(int ir=1;ir<=extenSlice.numberOfRegions;ir++){
			mapr[ir][0]=ix+1;
			for(int i=1;i<=extenSlice.numberOfElements;i++)
				if(extenSlice.element[i].getRegion()==ir)
				{
					ix++;
					vn[ix]=extenSlice.element[i].getVertNumb();
				}
			mapr[ir][1]=ix;
		}

		for(int ir=1;ir<=extenSlice.numberOfRegions;ir++){
			extenSlice.region[ir].setFirstEl(mapr[ir][0]);
			extenSlice.region[ir].setLastEl(mapr[ir][1]);
		}


		for(int i=1;i<=extenSlice.numberOfElements;i++){
			extenSlice.element[i].setVertNumb(vn[i]);
		}

		extenSlice.scaleFactor=slice.scaleFactor;



		Model model=dropUnusedNodes(extenSlice);

		model.defMode=slice.defMode;
		model.deform=slice.deform;

		String folder=new File(slice.meshFilePath).getParentFile().getPath();
		String fout=folder+"\\rotExtended.txt";

		model.writeMesh(fout);


		return model;
	}
	
	public Model revolveLine(Vect lineElems,int[] regs, int seg,double dtt){


		int nRegions=0;
		for(int ir:regs)
			if(ir>nRegions) nRegions=ir;
			

		int L=lineElems.length;

		int nElements=seg*(L-1);
		int nNodes=(seg+1)*L;

		Model revolved=new Model(nRegions,nElements,nNodes,"quadrangle");

		for(int i=1;i<=nRegions;i++){


			revolved.region[i].setMaterial("mat"+i);
			revolved.region[i].setName("reg"+i);

		}




		Mat R;


	
		for(int j=0;j<=seg;j++){
			
			if(j*dtt-2*PI>1e-6) break;

				R=util.rotMat2D(j*dtt);
				

			for(int i=0;i<L;i++)	
			{
				
				int nn=i+1+j*L;
								
				revolved.node[nn].setCoord(R.mul(new Vect(lineElems.el[i],0)));
			}
		}
		

			
		int[] vertnumb=new int[4];
		for(int j=0;j<seg;j++){
			
			if(j*dtt-2*PI>1e-6) break;
			
			boolean repeat=false;
			if(abs(j*dtt-2*PI)<1e-6)
				repeat=true;
			
			for(int i=0;i<L-1;i++)	{
				
				int ne=i+1+j*(L-1);
				
				vertnumb[0]=i+1+j*L;
				vertnumb[1]=i+2+j*L;
				
			
				if(!repeat){
				vertnumb[2]=i+2+(j+1)*L;
				vertnumb[3]=i+1+(j+1)*L;
				}
				else{
					vertnumb[2]=i+2;
					vertnumb[3]=i+1;
					}
				
				revolved.element[ne].setVertNumb(vertnumb);
				
				revolved.element[ne].setRegion(regs[i]);
				
		
		}
		}
		
		int[][] mapr=new int[revolved.numberOfRegions+1][2];
		int ix=0;
		for(int ir=1;ir<=revolved.numberOfRegions;ir++){
			mapr[ir][0]=ix+1;
			for(int i=1;i<=revolved.numberOfElements;i++)
				if(revolved.element[i].getRegion()==ir)
				{
					++ix;
				}
			mapr[ir][1]=ix;
		}


		for(int ir=1;ir<=revolved.numberOfRegions;ir++){
			revolved.region[ir].setFirstEl(mapr[ir][0]);
			revolved.region[ir].setLastEl(mapr[ir][1]);
		}
		
		reRegionGroupEls(revolved);

		
		revolved.scaleFactor=1;

	
		String bunFilePath = System.getProperty("user.dir") + "//revolved.txt";


		revolved.writeMesh(bunFilePath);

		return revolved;
	}
	
	
	
	
	public void reorderRot( ){
		String bun=util.getFile();

		Model m3d=new Model(bun);

	
			Model m2d= new Model(System.getProperty("user.dir") + "\\mot4thNewest.txt");
			
			
		

			int ir1=8;
			int ir2=1;
			
			int[][] mapEls2D3D=new int[1][1];

			int[] ind2D3D=new int[1];
			
			Mat[] R=new Mat[4];
			for(int k=0;k<4;k++){
				R[k]=util.rotMat2D(k*PI/2);
			}
			
			mapEls2D3D=new int[m2d.region[ir1].getNumbElements()][100];
			ind2D3D=new int[m2d.region[ir1].getNumbElements()];
			int p=0;
			for(int i=m2d.region[ir1].getFirstEl();i<=m2d.region[ir1].getLastEl();i++){
				Vect v1=m2d.getElementCenter(i);
				for(int j=m3d.region[ir2].getFirstEl();j<=m3d.region[ir2].getLastEl();j++){
					Vect v2=m3d.getElementCenter(j).v2();


					for(int k=0;k<4;k++){
				
					if(R[k].mul(v1).sub(v2).norm()<1e-4) mapEls2D3D[p][ind2D3D[p]++]=j;
				}
				}
				p++;
							
		}
			
			int[][] vns=new int[m3d.region[ir2].getNumbElements()][m2d.nElVert];
		
			int L=ind2D3D[0];

			int[][] map=new int[L][ind2D3D.length];
		
			
				for(int k=0;k<ind2D3D.length;k++)
						for(int j=0;j<L;j++)
							map[j][k]=mapEls2D3D[k][j];
				
				
				int ix=0;
				
				
			

			int[] map2=new int[ind2D3D.length*L];
			for(int j=0;j<L;j++)
				for(int k=0;k<ind2D3D.length;k++){
					vns[ix]=Arrays.copyOf(m3d.element[map[j][k]].getVertNumb(),m3d.nElVert);
					map2[ix++]=map[j][k];
				}
			
			for(int i=0;i<map2.length;i++)
					m3d.element[map2[i]].setVertNumb(vns[i]);
			
			
			double tm=0;
			for(int j=m3d.region[ir2].getFirstEl();j<=m3d.region[ir2].getLastEl()/4;j++){
				Vect v2=m3d.getElementCenter(j).v2();
				double t=util.getAng(v2);
				if(t>tm) tm=t;
			}

			util.pr(tm);
	
			
			m3d.writeMesh((System.getProperty("user.dir") + "\\reordered.txt"));
			
		}
	
	



	public Model extendFlip(Model md0){

		int hv=0;

		boolean[] nodeIn=new boolean[md0.numberOfNodes+1];
		for(int i=1;i<=md0.numberOfElements;i++){
			int[] vertNum=md0.element[i].getVertNumb();
			for(int j=0;j<md0.nElVert;j++)
				nodeIn[vertNum[j]]=true;
		}



		int nRegions=md0.numberOfRegions;

		int nElements=2*md0.numberOfElements;

		int ix=0;
		boolean[] onAx=new boolean[md0.numberOfNodes+1];
		for(int i=1;i<=md0.numberOfNodes;i++)
		{
			if(!nodeIn[i]) continue;
			Vect z=md0.node[i].getCoord();
			if(Math.abs(z.el[hv])<1e-10) onAx[i]=true;

		}

		int nNodes=2*md0.numberOfNodes;


		Model md1=new Model(nRegions,nElements,nNodes,md0.elType);

		int n1=1,N;
		for(int i=1;i<=md0.numberOfRegions;i++){

			if(i>1)
				n1=md1.region[i-1].getLastEl()+1;
			N=md0.region[i].getNumbElements();
			md1.region[i].setFirstEl(n1);
			md1.region[i].setLastEl(n1+2*N-1);
			md1.region[i].setMaterial(md0.region[i].getMaterial());
			md1.region[i].setMur(md0.region[i].getMur());
			md1.region[i].setM(md0.region[i].getM());
			md1.region[i].setJ(md0.region[i].getJ());
			md1.region[i].setSigma(md0.region[i].getSigma());
			md1.region[i].setPois(md0.region[i].getPois());
			md1.region[i].setYng(md0.region[i].getYng());

		}





		for(int i=1;i<=md0.numberOfNodes;i++)	
			md1.node[i].setCoord(md0.node[i].getCoord());

		int[] map=new int[md0.numberOfNodes+1];
		for(int i=1;i<=md0.numberOfNodes;i++)	
		{

			if(onAx[i]) continue;

			Vect v=md0.node[i].getCoord();

			v.el[hv]*=-1;
			md1.node[i+md0.numberOfNodes].setCoord(v);
			map[i]=i+md0.numberOfNodes;
		}
		ix=0;

		for(int ir=1;ir<=md0.numberOfRegions;ir++)		{

			int ne=md0.region[ir].getNumbElements();
			ix=0;
			for(int i=md0.region[ir].getFirstEl();i<=md0.region[ir].getLastEl();i++)	{

				ix++;
				int p=1;	
				p=md1.region[ir].getFirstEl();
				int[] vertNumb=md0.element[i].getVertNumb();

				for(int j=0;j<2;j++){
					int ie=p-1+j*ne+ix;
					if(j==0)
						md1.element[ie].setVertNumb(vertNumb);	
					else{

						int[] vert1=new int[md0.nElVert];
						for(int k=0;k<md0.nElVert;k++)
							if(!onAx[vertNumb[k]] )
								vert1[k]=map[vertNumb[k]];
							else
								vert1[k]=vertNumb[k];

						int[] vert2=new int[md0.nElVert];
						for(int k=0;k<md0.nElVert;k++)				
							vert2[k]=vert1[md0.nElVert-1-k];			

						for(int k=0;k<md0.nElVert;k++){
							md1.element[ie].setVertNumb(k,vert2[k]);

						}


					}
				}
			}

		}


		md1.scaleFactor=md0.scaleFactor;
		String file = System.getProperty("user.dir") + "//reflected.txt";
		md1.writeMesh(file);

		return md1;
	}

	
	public void flip(int mode)
	{
		
		String bun=util.getFile();
		if(bun==null || bun.equals("") )  throw new NullPointerException("file not found.");
		Model model=new Model(bun);
		flip(model,mode,true);
		
	}
	public void flip(Model model,int mode,boolean write){

		for(int i=1;i<=model.numberOfNodes;i++)
		{
			
			Vect z=model.node[i].getCoord();
			if(mode==0) z.el[0]*=-1;
			else
				if(mode==1) z.el[1]*=-1;
				else if(mode==2) z.el[2]*=-1;
			
			model.node[i].setCoord(z);
			
			
		}
		
		int L=model.nElVert;
		for(int i=1;i<=model.numberOfElements;i++){
			int[] vn=model.element[i].getVertNumb();
			
			for(int j=0;j<L;j++)
				model.element[i].setVertNumb(j,vn[L-1-j]);
		}
		
			


		if(write){
		String file = System.getProperty("user.dir") + "//flipped.txt";
		model.writeMesh(file);		
		}

	}


	public void combineTwoNode(Model motor,String rotf, String statf,double ang,double step){
		Combiner cb=new Combiner();
		cb.combineTwoNode(motor,rotf,statf,ang,step,false);
	}

	public void combineOneNode(Model motor,String rotf, String statf,double ang,double step){
		Combiner cb=new Combiner();
		cb.combineOneNode(motor,rotf,statf,ang,step);
	}

/*	public void combineFull(Model motor,String rotf, String statf,double ang){
		Combiner cb=new Combiner();
		cb.combineFull(motor,rotf,statf,ang);
	}*/

	public void putCap(Model model,String body, String cap,int mode){

		Model mbody=new Model(body);

		Model mcap=new Model(cap);



		setSliceBounds(mbody);

		double eps=1e-4;
		double hm;
		if(mode==0) hm=mbody.h1;
		else hm=mbody.h2;

		int[] map=new int[mcap.numberOfNodes+1];

		for(int i=1;i<=mbody.numberOfNodes;i++){

			Vect v=mbody.node[i].getCoord();
			if(abs(v.el[2]-hm)>eps) continue;
			for(int j=1;j<=mcap.numberOfNodes;j++){

				Vect v2=mcap.node[j].getCoord();
				//v2.hshow();

				if(abs(v2.el[2]-hm)>eps) continue;
				if(v2.sub(v).norm()<eps) {	
					map[j]=i;
				}

			}
		}

		int nc=0;
		for(int i=1;i<=mcap.numberOfNodes;i++)
			if(map[i]>0) nc++;

		util.pr(nc);

		int nNodes1=mbody.numberOfNodes;;
		int nNodes2=mcap.numberOfNodes;

		int nRegions1=mbody.numberOfRegions;
		int nRegions2=mcap.numberOfRegions;
		int nElements1=mbody.numberOfElements;
		int nElements2=mcap.numberOfElements;

		int nRegions=nRegions1+nRegions2;
		int nElements=nElements1+nElements2;
		int nNodes=nNodes1+nNodes2;

		model.alloc(nRegions,nElements,nNodes,mbody.elType);

		for(int i=1;i<=nRegions1;i++){
			model.region[i]=mbody.region[i].deepCopy();

		}
		for(int i=1;i<=nRegions2;i++)	{

			model.region[i+nRegions1]=mcap.region[i].deepCopy();
			model.region[i+nRegions1].setFirstEl(mcap.region[i].getFirstEl()+nElements1);
			model.region[i+nRegions1].setLastEl(mcap.region[i].getLastEl()+nElements1);

		}

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
				model.element[i].setRegion(ir);


		for(int i=1;i<=nElements1;i++)	
			model.element[i]=mbody.element[i];


		for(int i=1;i<=nElements2;i++)	{
			int[] vertNumb=mcap.element[i].getVertNumb();
			model.element[i+nElements1]=new Element(mcap.elType);

			int[] vn2=new int[mcap.nElVert];
			for(int j=0;j<vn2.length;j++)
				if(map[vertNumb[j]]==0)
					vn2[j]=vertNumb[j]+nNodes1;
				else
					vn2[j]=map[vertNumb[j]];
			model.element[i+nElements1].setVertNumb(vn2);			

		}

		for(int i=1;i<=nNodes1;i++)
			model.node[i]=mbody.node[i];

		for(int i=1;i<=nNodes2;i++)	
			model.node[i+nNodes1]=mcap.node[i];

		model.scaleFactor=mbody.scaleFactor;
		model.setElType(mbody.elType);



		model.writeMesh(System.getProperty("user.dir") + "//capped.txt");





	}



	public void putBase(Model model,String body, String cap,int mode){

		Model mbody=new Model(body);

		Model mcap=new Model(cap);


		setSliceBounds(mbody);

		double hm;
		if(mode==0) hm=mbody.h1;
		else hm=mbody.h2;

		int[] map=new int[mcap.numberOfNodes+1];

		for(int i=1;i<=mbody.numberOfNodes;i++){

			Vect v=mbody.node[i].getCoord();
			//	if(abs(v.el[2]-hm)>1e-6) continue;
			for(int j=1;j<=mcap.numberOfNodes;j++){

				Vect v2=mcap.node[j].getCoord();

				//if(abs(v2.el[2]-hm)>1e-6) continue;
				if(v2.sub(v).norm()<1e-4) {	
					util.pr(i);
					map[j]=i;
				}

			}
		}

		int nc=0;
		for(int i=1;i<=mcap.numberOfNodes;i++)
			if(map[i]>0) nc++;

		util.pr(nc);

		int nNodes1=mbody.numberOfNodes;;
		int nNodes2=mcap.numberOfNodes;

		int nRegions1=mbody.numberOfRegions;
		int nRegions2=mcap.numberOfRegions;
		int nElements1=mbody.numberOfElements;
		int nElements2=mcap.numberOfElements;

		int nRegions=nRegions1+nRegions2;
		int nElements=nElements1+nElements2;
		int nNodes=nNodes1+nNodes2;

		model.alloc(nRegions,nElements,nNodes,mbody.elType);

		for(int i=1;i<=nRegions1;i++){
			model.region[i]=mbody.region[i].deepCopy();

		}
		for(int i=1;i<=nRegions2;i++)	{

			model.region[i+nRegions1]=mcap.region[i].deepCopy();
			model.region[i+nRegions1].setFirstEl(mcap.region[i].getFirstEl()+nElements1);
			model.region[i+nRegions1].setLastEl(mcap.region[i].getLastEl()+nElements1);

		}

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
				model.element[i].setRegion(ir);


		for(int i=1;i<=nElements1;i++)	
			model.element[i]=mbody.element[i];


		for(int i=1;i<=nElements2;i++)	{
			int[] vertNumb=mcap.element[i].getVertNumb();
			model.element[i+nElements1]=new Element(mcap.elType);

			int[] vn2=new int[mcap.nElVert];
			for(int j=0;j<vn2.length;j++)
				if(map[vertNumb[j]]==0)
					vn2[j]=vertNumb[j]+nNodes1;
				else
					vn2[j]=map[vertNumb[j]];
			model.element[i+nElements1].setVertNumb(vn2);			

		}

		for(int i=1;i<=nNodes1;i++)
			model.node[i]=mbody.node[i];

		for(int i=1;i<=nNodes2;i++)	
			model.node[i+nNodes1]=mcap.node[i];

		model.scaleFactor=mbody.scaleFactor;
		model.setElType(mbody.elType);



		model.writeMesh(System.getProperty("user.dir") + "//capped.txt");





	}

	public void connect(String bun1file, String bun2file){
		double[] bound= {-1e10,1e10,-1e10,1e10,-1e10,1e10};
		connect(bun1file,bun2file,bound);
		
	}
	
	public void connect(String bun1file, String bun2file,double[] bound){

		double eps=1e-5;
		Model bun1=new Model(bun1file);

		Model bun2=new Model(bun2file);
		Model model=new Model();

		int[] map=new int[bun2.numberOfNodes+1];

		for(int i=1;i<=bun1.numberOfNodes;i++){

			Vect v=bun1.node[i].getCoord();

			if(v.el[0]<bound[0] || v.el[0]>bound[1] ) continue;
			if(v.el[1]<bound[2] || v.el[1]>bound[3] ) continue;
			if(bun1.dim==3)
			if(v.el[2]<bound[4] || v.el[2]>bound[5] ) continue;
			
			for(int j=1;j<=bun2.numberOfNodes;j++){

				Vect v2=bun2.node[j].getCoord();
				
				if(v2.el[0]<bound[0] || v2.el[0]>bound[1] ) continue;
				if(v2.el[1]<bound[2] || v2.el[1]>bound[3] ) continue;
				if(bun1.dim==3)
				if(v2.el[2]<bound[4] || v2.el[2]>bound[5] ) continue;


				double d=v2.sub(v).norm();
				
				if(d<eps) {	
				
					map[j]=i;
				}

			}
		}

		int nc=0;
		for(int i=1;i<=bun2.numberOfNodes;i++)

			if(map[i]>0) {nc++;

		}

		
		util.pr("Number of Matched nodes: "+nc);

		int nNodes1=bun1.numberOfNodes;;
		int nNodes2=bun2.numberOfNodes;

		int nRegions1=bun1.numberOfRegions;
		int nRegions2=bun2.numberOfRegions;
		int nElements1=bun1.numberOfElements;
		int nElements2=bun2.numberOfElements;

		int nRegions=nRegions1+nRegions2;
		int nElements=nElements1+nElements2;
		int nNodes=nNodes1+nNodes2;

		model.alloc(nRegions,nElements,nNodes,bun1.elType);

		for(int i=1;i<=nRegions1;i++){
			model.region[i]=bun1.region[i].deepCopy();

		}
		for(int i=1;i<=nRegions2;i++)	{

			model.region[i+nRegions1]=bun2.region[i].deepCopy();
			model.region[i+nRegions1].setFirstEl(bun2.region[i].getFirstEl()+nElements1);
			model.region[i+nRegions1].setLastEl(bun2.region[i].getLastEl()+nElements1);

		}

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
				model.element[i].setRegion(ir);


		for(int i=1;i<=nElements1;i++)	
			model.element[i]=bun1.element[i];


		for(int i=1;i<=nElements2;i++)	{
			int[] vertNumb=bun2.element[i].getVertNumb();
			model.element[i+nElements1]=new Element(bun2.elType);

			int[] vn2=new int[bun2.nElVert];
			for(int j=0;j<vn2.length;j++)
				if(map[vertNumb[j]]==0)
					vn2[j]=vertNumb[j]+nNodes1;
				else
					vn2[j]=map[vertNumb[j]];
			model.element[i+nElements1].setVertNumb(vn2);			

		}

		for(int i=1;i<=nNodes1;i++)
			model.node[i]=bun1.node[i];

		for(int i=1;i<=nNodes2;i++)	
			model.node[i+nNodes1]=bun2.node[i];

		model.scaleFactor=bun1.scaleFactor;
		model.setElType(bun1.elType);



		model.writeMesh(System.getProperty("user.dir") + "//connected.txt");





	}
	
	public void connectCyl(String bun1file, String bun2file,double[] bound,double eps){

		Model bun1=new Model(bun1file);

		Model bun2=new Model(bun2file);
		Model model=new Model();

		int[] map=new int[bun2.numberOfNodes+1];

		for(int i=1;i<=bun1.numberOfNodes;i++){
			

			Vect v0=bun1.node[i].getCoord();
			Vect v=new Vect(model.dim);
			Vect v2d=v0.v2();
			v.el[0]=v2d.norm();
			v.el[1]=util.getAng(v2d);
			if(model.dim==3)
			v.el[2]=v0.el[2];

			if(v.el[0]<bound[0] || v.el[0]>bound[1] ) continue;
			if(v.el[1]<bound[2] || v.el[1]>bound[3] ) continue;
			if(bun1.dim==3)
			if(v.el[2]<bound[4] || v.el[2]>bound[5] ) continue;
			
			//if(Math.abs(v.el[2])<.08 ) continue;
			
			for(int j=1;j<=bun2.numberOfNodes;j++){

				Vect v00=bun2.node[j].getCoord();
				Vect v2=new Vect(model.dim);
				Vect v22d=v00.v2();
				v2.el[0]=v22d.norm();
				v2.el[1]=util.getAng(v22d);
				if(model.dim==3)
				v2.el[2]=v00.el[2];
				
				//if(Math.abs(v2.el[2])<.08 ) continue;

			
				if(v2.el[0]<bound[0] || v2.el[0]>bound[1] ) continue;
				if(v2.el[1]<bound[2] || v2.el[1]>bound[3] ) continue;
				if(bun1.dim==3)
				if(v2.el[2]<bound[4] || v2.el[2]>bound[5] ) continue;

				double d=v2.sub(v).norm();
				
				if(d<eps) {	
				
					map[j]=i;
				}

			}
		}

		int nc=0;
		for(int i=1;i<=bun2.numberOfNodes;i++)

			if(map[i]>0) {nc++;

		}

		
		util.pr("Number of Matched nodes: "+nc);

		int nNodes1=bun1.numberOfNodes;;
		int nNodes2=bun2.numberOfNodes;

		int nRegions1=bun1.numberOfRegions;
		int nRegions2=bun2.numberOfRegions;
		int nElements1=bun1.numberOfElements;
		int nElements2=bun2.numberOfElements;

		int nRegions=nRegions1+nRegions2;
		int nElements=nElements1+nElements2;
		int nNodes=nNodes1+nNodes2;

		model.alloc(nRegions,nElements,nNodes,bun1.elType);

		for(int i=1;i<=nRegions1;i++){
			model.region[i]=bun1.region[i].deepCopy();

		}
		for(int i=1;i<=nRegions2;i++)	{

			model.region[i+nRegions1]=bun2.region[i].deepCopy();
			model.region[i+nRegions1].setFirstEl(bun2.region[i].getFirstEl()+nElements1);
			model.region[i+nRegions1].setLastEl(bun2.region[i].getLastEl()+nElements1);

		}

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
				model.element[i].setRegion(ir);


		for(int i=1;i<=nElements1;i++)	
			model.element[i]=bun1.element[i];


		for(int i=1;i<=nElements2;i++)	{
			int[] vertNumb=bun2.element[i].getVertNumb();
			model.element[i+nElements1]=new Element(bun2.elType);

			int[] vn2=new int[bun2.nElVert];
			for(int j=0;j<vn2.length;j++)
				if(map[vertNumb[j]]==0)
					vn2[j]=vertNumb[j]+nNodes1;
				else
					vn2[j]=map[vertNumb[j]];
			model.element[i+nElements1].setVertNumb(vn2);			

		}

		for(int i=1;i<=nNodes1;i++)
			model.node[i]=bun1.node[i];

		for(int i=1;i<=nNodes2;i++)	
			model.node[i+nNodes1]=bun2.node[i];

		model.scaleFactor=bun1.scaleFactor;
		model.setElType(bun1.elType);



		model.writeMesh(System.getProperty("user.dir") + "//connected.txt");





	}
	
	public void connectivity(double eps){

		String bun=util.getFile();
		if(bun==null || bun.equals("") )  throw new NullPointerException("file not found.");
		Model model=new Model(bun);
		int[] map=new int[1+model.numberOfNodes];
		int ix=0;
		boolean[] nnc=new boolean[1+model.numberOfNodes];
		
		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
		{
			int[] vn=model.element[i].getVertNumb();
			
			
			for(int j=0;j<model.nElVert;j++){
				Vect v=model.node[vn[j]].getCoord();
				if(abs(v.el[0])>1e-4)continue;
					nnc[vn[j]]=true;
			}
		}
		}
		
		for(int i=1;i<=model.numberOfNodes;i++){
			
			if(!nnc[i]) {
			
				continue;
			}
			    
			Vect v1=model.node[i].getCoord();
			for(int j=i+1;j<=model.numberOfNodes;j++){
				
				if(!nnc[j]) continue;
				
				Vect v2=model.node[j].getCoord();
				double dist=v1.sub(v2).norm();

				if(map[j]==0 && dist<eps){

					map[j]=i;
					ix++;
				}
				}
			}


		for(int i=1;i<=1*model.numberOfElements;i++){
			int[] vn=model.element[i].getVertNumb();
			int[] vn2=new int[vn.length];
			for(int j=0;j<model.nElVert;j++){
		
				
				if(map[vn[j]]>0){
				 vn2[j]=map[vn[j]];
				
				}
				else{
					vn2[j]=vn[j];
				}
		}
		
			model.element[i].setVertNumb(vn2);
		}
		

		util.pr(ix+"/"+model.numberOfNodes);
		String folder=new File(bun).getParentFile().getPath();

		String sewedMesh = folder + "//sewed.txt";
	
		model.writeMesh(sewedMesh);

	}

	public void putRotor(Model model,String rot, String body){

		Model mbody=new Model(body);

		Model mcap=new Model(rot);


		setSliceBounds(mbody);

		double h1=mbody.h1;
		double h2=mbody.h2;

		double t=.0045;
		double hm1=h1+t/2;
		double hm2=h2-t/2;
		util.pr(hm1);
		double eps=1e-4;
		double et=t;

		int[] map=new int[mbody.numberOfNodes+1];

		for(int i=1;i<=mcap.numberOfNodes;i++){

			Vect v=mcap.node[i].getCoord();
			if(abs(v.el[2]-hm1)>et && abs(v.el[2]-hm2)>et) continue;
			for(int j=1;j<=mbody.numberOfNodes;j++){

				Vect v2=mbody.node[j].getCoord();

				if(abs(v2.el[2]-hm1)>et && abs(v2.el[2]-hm2)>et) continue;

				if(v2.sub(v).norm()<eps) {	
					map[j]=i;
				}

			}
		}

		int nc=0;
		for(int i=1;i<=mbody.numberOfNodes;i++)
			if(map[i]>0) nc++;



		int nNodes2=mbody.numberOfNodes;;
		int nNodes1=mcap.numberOfNodes;

		int nRegions2=mbody.numberOfRegions;
		int nRegions1=mcap.numberOfRegions;
		int nElements2=mbody.numberOfElements;
		int nElements1=mcap.numberOfElements;

		int nRegions=nRegions1+nRegions2;
		int nElements=nElements1+nElements2;
		int nNodes=nNodes1+nNodes2;

		model.alloc(nRegions,nElements,nNodes,mbody.elType);

		for(int i=1;i<=nRegions1;i++){
			model.region[i]=mcap.region[i].deepCopy();

		}
		for(int i=1;i<=nRegions2;i++)	{

			model.region[i+nRegions1]=mbody.region[i].deepCopy();
			model.region[i+nRegions1].setFirstEl(mbody.region[i].getFirstEl()+nElements1);
			model.region[i+nRegions1].setLastEl(mbody.region[i].getLastEl()+nElements1);

		}

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
				model.element[i].setRegion(ir);


		for(int i=1;i<=nElements1;i++)	
			model.element[i]=mcap.element[i];


		for(int i=1;i<=nElements2;i++)	{
			int[] vertNumb=mbody.element[i].getVertNumb();
			model.element[i+nElements1]=new Element(mbody.elType);

			int[] vn2=new int[mbody.nElVert];
			for(int j=0;j<vn2.length;j++)
				if(map[vertNumb[j]]==0)
					vn2[j]=vertNumb[j]+nNodes1;
				else
					vn2[j]=map[vertNumb[j]];
			model.element[i+nElements1].setVertNumb(vn2);			

		}

		for(int i=1;i<=nNodes1;i++)
			model.node[i]=mcap.node[i];

		for(int i=1;i<=nNodes2;i++)	
			model.node[i+nNodes1]=mbody.node[i];

		model.scaleFactor=mbody.scaleFactor;
		model.setElType(mcap.elType);



		model.writeMesh(System.getProperty("user.dir") + "//rotorPut.txt");





	}
	public void assemble(){

		String body=util.getFile();
		if(body==null || body.equals("") )return;
		
		String cap=util.getFile();
		if(cap==null || cap.equals("") )return;
	
	assemble(cap,body);

	}

	public void assemble(String rot, String body){
		
		assemble(rot,body,true);
	}
	
	public Model assemble(String rot, String body, boolean write){

		Model mbody=new Model(body);

		Model mcap=new Model(rot);
		
		return assemble(mcap,mbody,write);
		
	}
		
		public Model assemble(Model mcap, Model mbody, boolean write){

	
		Model model=new Model();
		
		int nNodes2=mbody.numberOfNodes;;
		int nNodes1=mcap.numberOfNodes;

		int nRegions2=mbody.numberOfRegions;
		int nRegions1=mcap.numberOfRegions;
		int nElements2=mbody.numberOfElements;
		int nElements1=mcap.numberOfElements;

		int nRegions=nRegions1+nRegions2;
		int nElements=nElements1+nElements2;
		int nNodes=nNodes1+nNodes2;

		model.alloc(nRegions,nElements,nNodes,mbody.elType);

		for(int i=1;i<=nRegions1;i++){
			model.region[i]=mcap.region[i].deepCopy();

		}
		for(int i=1;i<=nRegions2;i++)	{

			model.region[i+nRegions1]=mbody.region[i].deepCopy();
			model.region[i+nRegions1].setFirstEl(mbody.region[i].getFirstEl()+nElements1);
			model.region[i+nRegions1].setLastEl(mbody.region[i].getLastEl()+nElements1);

		}

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
				model.element[i].setRegion(ir);


		for(int i=1;i<=nElements1;i++)	
			model.element[i]=mcap.element[i];


		for(int i=1;i<=nElements2;i++)	{
			int[] vertNumb=mbody.element[i].getVertNumb();
			model.element[i+nElements1]=new Element(mbody.elType);

			int[] vn2=new int[mbody.nElVert];
			for(int j=0;j<vn2.length;j++)
				vn2[j]=vertNumb[j]+nNodes1;

			model.element[i+nElements1].setVertNumb(vn2);			
			

		}

		for(int i=1;i<=nNodes1;i++)
			model.node[i]=mcap.node[i];

		for(int i=1;i<=nNodes2;i++)	
			model.node[i+nNodes1]=mbody.node[i];

		model.scaleFactor=mbody.scaleFactor;
		model.setElType(mcap.elType);



		if(write){
			String folder=new File(mbody.meshFilePath).getParentFile().getPath();
		model.writeMesh(folder + "//assembled.txt");
		}

		return model;

		
	}


		public Model meshPlunger(double x,boolean write){
			
			x=x*100;
			

			
				//double[][] bb={{-2,6.5,0,3},{0,4,0,2},{.4,3.6,0,1.6},{.5,3.5,1.0,1.5},{.4,1.7,0,.825},{3.6,4.0,.825,.925},{2.7-x,6-x,0,.825}};

			// ANSYS
			//	double[][] bb={{0,32,0,12},{0,24,0,12},{4,2,0,8},{4,20,4.2,8},{4,6,0,42},{20,24,4,4.2},{8-x,28-x,0,4}};
		
			double[][] bb={{0,3.2,0,1.2},{0,2.4,0,1.2},{.4,2,0,.8},{.4,2,.42,.8},{.4,.6,0,.42},{2,2.4,.4,.42},{.8-x,2.8-x,0,.4}};
			
				double scale=10;

				for(int j=0;j<bb.length;j++)
					for(int k=0;k<bb[0].length;k++){
						bb[j][k]*=scale;

						if(k==2&& bb[j][k]==0) 	bb[j][k]=.2;
					}

				Geometry mg=new Geometry(bb);


				mg.blockName[0]="air";
				mg.blockName[1]="iron";
				mg.blockName[2]="air";
				mg.blockName[3]="coil";
				mg.blockName[4]="iron";
				mg.blockName[5]="air";
				mg.blockName[6]="plunger";
			
				for(int j=0;j<bb.length;j++){
					if(j==0){
						for(int k=0;k<bb[0].length;k++){

							mg.discLeft[j][k]=false;
							mg.discRight[j][k]=false;
						}
					}



					for(int k=0;k<bb[0].length;k++){

						mg.minMeshRight[j][k]=5;
						mg.minMeshLeft[j][k]=5;

					}


					if(j==6){
						for(int k=0;k<bb[0].length;k++){
							if(k==0){

								mg.minMeshLeft[j][k]=5;
								mg.baseLeft[j][k]=1;
							}
							else if(k==2){

								mg.minMeshRight[j][k]=2;

							}
							else if(k==3){

									mg.minMeshRight[j][k]=1;
									mg.baseRight[j][k]=1;

								}
						}
					}


				}
				//mg.setDisc();
			

				Model model=getOrthogMesh(mg,"",false);
		
				model.rotate(PI/2);
				this.flip(model,0,false);
				//model=pileRotate(model,90, PI/90,false);
				
				if(write){
				String bun=System.getProperty("user.dir") + "\\plungers\\plungerX.txt";
				model.writeMesh(bun);
				}
				

				model.setBounds();

				model.setFemCalc();

			return model;

		}
		
	public Model meshEM(double x,boolean write){
			
		x=.00;
			x=x*100;
			
			x=.1;
	
			//double[][] bb={{0,6,0,8},{0,2.5,0,4},{.5,2,.5,4},{0,2.5,5-x,5.5-x},{.7,1.8,.7,3.5}};
			
			double[][] bb={{0,8,0,8},{0,2.5,0,5},{.5,1.5,1,5},{0,2.5,5+x,5.5+x}};
			
				double scale=10;

				for(int j=0;j<bb.length;j++)
					for(int k=0;k<bb[0].length;k++){
						bb[j][k]*=scale;

						if(k==0&& bb[j][k]==0) 	bb[j][k]=.2;
					}

				Geometry mg=new Geometry(bb);


				mg.blockName[0]="xair";
				mg.blockName[1]="iron";
				mg.blockName[2]="coil";
				mg.blockName[3]="moving";
			//	mg.blockName[4]="coil";
			
			
				for(int j=0;j<bb.length;j++){
					if(j==0){
						for(int k=0;k<bb[0].length;k++){

							mg.discLeft[j][k]=false;
							mg.discRight[j][k]=false;
						}
					}



					for(int k=0;k<bb[0].length;k++){

						mg.minMeshRight[j][k]=.5;
						mg.minMeshLeft[j][k]=.5;

					}


					if(j==3){
						for(int k=0;k<bb[0].length;k++){
							if(k==2){

								mg.minMeshLeft[j][k]=.2;
								mg.baseLeft[j][k]=1;
							}
									}
					}


				}
			

				Model model=getOrthogMesh(mg,"",false);
		
			//	model.rotate(-PI/2);
				
				//this.flip(model,1,false);
			//	model=pileRotate(model,1, PI/180,false);
				
				if(write){
				String bun=System.getProperty("user.dir") + "\\plungers\\plungerX.txt";
				model.writeMesh(bun);
				}
				

				model.setBounds();

				model.setFemCalc();

			return model;

		}
		
		public Model meshPlunger97(double x,boolean write){
			
			x=x*100;
			

	
			double[][] bb={{0,4.3,0,2},{.3,4,1,1.9},{.3,4.3,0,.8},{1.3,2.5,0.8,1},{1.1,3.8,0,.7},{1.1,2.4,0,.4},{3.2,3.8,0,.4}};
			
				double scale=10;

				for(int j=0;j<bb.length;j++)
					for(int k=0;k<bb[0].length;k++){
						bb[j][k]*=scale;

						if(k==2&& bb[j][k]==0) 	bb[j][k]=1;
					}

				Geometry mg=new Geometry(bb);


				mg.blockName[0]="iron";
				mg.blockName[1]="coil";
				mg.blockName[2]="air";
				mg.blockName[3]="air";
				mg.blockName[4]="plunger";
				mg.blockName[5]="air";
				mg.blockName[6]="air";
			
				for(int j=0;j<bb.length;j++){
					if(j==0){
						for(int k=0;k<bb[0].length;k++){

							mg.discLeft[j][k]=false;
							mg.discRight[j][k]=false;
						}
					}



					for(int k=0;k<bb[0].length;k++){

						mg.minMeshRight[j][k]=.8;
						mg.minMeshLeft[j][k]=.8;

					}


					if(j==4){
						for(int k=0;k<bb[0].length;k++){
						
						 if(k==3){

									mg.minMeshRight[j][k]=.2;
									mg.baseRight[j][k]=1;
									mg.minMeshLeft[j][k]=.4;
									//mg.baseLeft[j][k]=1;

								}
						}
					}
					
					if(j==5){
						for(int k=0;k<bb[0].length;k++){
						
						 if(k==3){

								mg.minMeshLeft[j][k]=.6;
								mg.baseLeft[j][k]=1;

								}
						}
					}
					
					if(j==1){
						for(int k=0;k<bb[0].length;k++){
							if(k==0){

								mg.minMeshLeft[j][k]=.4;
								mg.baseLeft[j][k]=1;
							}
							else if(k==1){

								mg.minMeshRight[j][k]=.4;
								mg.baseRight[j][k]=1;

							}
							else if(k==3){

									mg.minMeshRight[j][k]=.25;
									mg.baseRight[j][k]=1;

								}
						}
					}


				}
			

				Model model=getOrthogMesh(mg,"",false);
		
				model.rotate(-PI/2);
				this.translate(model,new Vect(0,.043));
			//	this.flip(model,0,false);
				//model=pileRotate(model,90, PI/90,false);
				
				double y1=.02446666;
				double y2=.03533333;
				double x1=.004;
				double x2=.00725;
				
				
				double y3=.029;
				double y4=.042;
				double x3=.007;
				double x4=.01;
			
		
				
				for(int i=1;i<=0*model.numberOfNodes;i++){
					Vect z=model.node[i].getCoord();
					 x=z.el[0];
					double y=z.el[1];
					
					if(x>x1 &&  x<x2 && y>y1 && y<y2){
						Vect z2=z.deepCopy();
						z2.el[0]=x-(y-y1)/(y2-y1)*(x-x1)/(x2-x1)*(x2-x1)/2;
						
						model.node[i].setCoord(z2);
					}
					
					else if(x>x3 &&  x<x4 && y>y3 && y<y4){
						Vect z2=z.deepCopy();
						z2.el[0]=x-(y-y3)/(y4-y3)*(x-x3)/(x4-x3)*.002;
						//z2.el[0]=x-.002;
						
						model.node[i].setCoord(z2);
					}
						
					
				}
				if(write){
				String bun=System.getProperty("user.dir") + "\\plungers\\plungerX.txt";
				model.writeMesh(bun);
				}
				

				model.setBounds();

				model.setFemCalc();

			return model;

		}

	public Model meshPlungerOld(double x,boolean write){
	
		x=x*100;
		

			double[][] bb={{-10,25,0,10},{0,12,0,4.4},{.5,11.5,0,4},{.5,4,0,1},{1.25,10.75,1.5,3.5},{5-x,17-x,0,1},
					{15-x,15.3-x,1,1.6},{12,15-x,1.2,1.5},{11.05,12,1.0,1.03},{4,5-x,0,1}};



			double scale=10;

			for(int j=0;j<bb.length;j++)
				for(int k=0;k<bb[0].length;k++){
					bb[j][k]*=scale;

					if(k==2&& bb[j][k]==0) 	bb[j][k]=.4;
				}

			Geometry mg=new Geometry(bb);


			mg.blockName[0]="air";
			mg.blockName[2]="air";
			mg.blockName[1]="iron";
			mg.blockName[3]="iron";
			mg.blockName[4]="coil";
			mg.blockName[5]="plunger";
			mg.blockName[6]="plunger";
			mg.blockName[7]="spring";
			mg.blockName[9]="air";
			mg.blockName[8]="air";

			for(int j=0;j<bb.length;j++){
				if(j==0){
					for(int k=0;k<bb[0].length;k++){

						mg.discLeft[j][k]=false;
						mg.discRight[j][k]=false;
					}
				}



				for(int k=0;k<bb[0].length;k++){

					mg.minMeshRight[j][k]=.8;
					mg.minMeshLeft[j][k]=.8;

				}

				if(j==8){
					for(int k=0;k<bb[0].length;k++){
						if(k==2){

							mg.minMeshRight[j][k]=.08;
							mg.baseRight[j][k]=1;

						}
						else
							if(k==3){

								mg.minMeshLeft[j][k]=.08;
								mg.baseLeft[j][k]=1;

							}
					}
				}

				if(j==9){
					for(int k=0;k<bb[0].length;k++){
						if(k==0){

							mg.minMeshRight[j][k]=1.5*(1-x);
							mg.baseRight[j][k]=1;
						}
						else
							if(k==1){

								mg.minMeshLeft[j][k]=1.5*(1-x);
								mg.baseLeft[j][k]=1;
							}
					}
				}


			}
			//mg.setDisc();
		

			Model model=getOrthogMesh(mg,"",false);
	
			model.rotate(PI/2);
			this.flip(model,0,false);
		//	model=pileRotate(model,1, 1*PI/90,false);
			
			if(write){
			String bun=System.getProperty("user.dir") + "\\plungers\\plungerX.txt";
			model.writeMesh(bun);
			}
			

			model.setBounds();

			model.setFemCalc();

		return model;

	}
	
	public void meshTrans3(){
	
			double[][] bb={{-6,6,0,6},{-2.5,2.5,0,4},{-1.5,-.5,0,3},{.5,1.5,0,3},{-2.9,-2.6,0,2},{-1.4,-1.1,0,2},{-1,-.7,0,2},{.6,.9,0,2},{1.1,1.4,0,2},{2.6,2.9,0,2}};
			
		//	double[][] bb={{-6,6,0,6},{-2.5,2.5,0,4},{-1.5,1.5,0,3},{-3.6,-2.6,0,2},{-1.4,-.4,0,2}};


	//	//	double[][] bb={{-10,25,0,10},{0,12,0,4.4},{.95,11.05,1,3.45},{1.25,10.75,1.2,3.3},{6,6.1,0,1}};


			double scale=10;

			for(int j=0;j<bb.length;j++)
				for(int k=0;k<bb[0].length;k++){
					bb[j][k]*=scale;

					if(k==2&& bb[j][k]==0) 	bb[j][k]=.4;
				}

			Geometry mg=new Geometry(bb);
/*			mg.blockName[0]="air";
			mg.blockName[2]="air";

					mg.blockName[4]="air";*/

			mg.blockName[0]="air";
			mg.blockName[1]="iron";
			mg.blockName[2]="air";
			mg.blockName[3]="air";
			mg.blockName[4]="coilU+";
			mg.blockName[5]="coilU-";
			mg.blockName[6]="coilV+";
			mg.blockName[7]="coilV-";
			mg.blockName[8]="coilW+";
			mg.blockName[9]="coilW-";
			
		

			for(int j=0;j<bb.length;j++){
				
				for(int k=0;k<bb[0].length;k++){
				mg.minMeshRight[j][k]=1;
				mg.minMeshLeft[j][k]=1;
				mg.baseLeft[j][k]=2;
				mg.baseRight[j][k]=2;
				
				if(j==0){
					mg.minMeshRight[j][0]=10;
					mg.minMeshLeft[j][1]=10;
					mg.minMeshLeft[j][3]=10;

				}
				}
			}



			Model model=getOrthogMesh2D(mg,"",false);
	
			
			//model=pileRotate(model,1, 1*PI/90,false);
			
			String bun=System.getProperty("user.dir") + "\\model2D.txt";
			model.writeMesh(bun);

		

	}



	public Model getOrthogMesh2D( Geometry mg, String file, boolean b){

		Model model=new Model();
		util.show(mg.blockBoundary);

		model.nBlocks=mg.nBlocks;
		model.nElVert=4;
		model.dim=2;;
		model.scaleFactor=mg.scaleFactor;
		model.elType="quadrangle";
		//**************************
		Discretizer discretizer=new Discretizer(mg);

		int I,J;
		double[] X=discretizer.X;
		double[] Y=discretizer.Y;


		I=X.length; J=Y.length;
		int nNodes=I*J;
		int nElements=(I-1)*(J-1); 
		Vect[] coords=new Vect[nNodes];
		int[][] elVerts=new int[nElements][model.nElVert];


		int ix=0;


		for(int j=0;j<J;j++)
			for(int i=0;i<I;i++)
			{

				coords[ix++]=new Vect(X[i],Y[j]);	

			}

		int elementNumber,Nx;


		Nx=I;
		elementNumber=0;

		for(int j=0;j<J-1;j++)
			for(int i=0;i<I-1;i++){

				elVerts[elementNumber][0]=(j+1)*Nx+i+2;
				elVerts[elementNumber][1]=(j+1)*Nx+i+1;
				elVerts[elementNumber][2]=j*Nx+i+1;
				elVerts[elementNumber][3]=j*Nx+i+2;
				elementNumber++;


			}

		Vect centerOfMass;

		int[] elementBlock=new int[nElements];
		for(int i=0;i<nElements;i++){
			centerOfMass=new Vect(2);
			for(int v=0;v<4; v++)
				centerOfMass=centerOfMass.add(coords[elVerts[i][v]-1]);
			centerOfMass=centerOfMass.times(.25);

			for(int ir=0;ir<model.nBlocks;ir++)
				if(discretizer.blockBoundary[ir][0]<centerOfMass.el[0])
					if(centerOfMass.el[0]<discretizer.blockBoundary[ir][1])
						if(discretizer.blockBoundary[ir][2]<centerOfMass.el[1])
							if(centerOfMass.el[1]<discretizer.blockBoundary[ir][3])
							{
								elementBlock[i]=ir+1;



							}
		} 




		ix=0;
		int iy=0;
		boolean[] nc=new boolean[nNodes+1];


		for(int i=0;i<nElements;i++){
			if(elementBlock[i]>0) {

				++ix;
				for(int j=0;j<model.nElVert;j++)
					if(!nc[elVerts[i][j]]) 
					{
						iy++;
						nc[elVerts[i][j]]=true;
					}
			}

		}

		model.numberOfElements=ix;
		model.numberOfNodes=iy;



		model.node=new Node[model.numberOfNodes+1];

		int[] mapNd=new int[nNodes+1];
		ix=0;
		for(int i=1;i<=nNodes;i++){

			if(nc[i]){
				ix++;
				model.node[ix]=new Node(ix,model.dim);		
				model.node[ix].setCoord(coords[i-1]);	
				mapNd[i]=ix;
			}
		}





		int nRegX=model.nBlocks;
		int[] nRegEls=new int[nRegX+1];
		Region[] regionx=new Region[nRegX+1];




		for(int ir=1;ir<=nRegX;ir++){
			regionx[ir]=new Region(model.dim);
		}

		regionx[1].setFirstEl(1);


		int[] mapEl=new int[nElements+1];
		int ire=0;
		for(int ir=1;ir<=nRegX;ir++){
			if(ir>1) regionx[ir].setFirstEl(regionx[ir-1].getLastEl()+1);
			for(int i=0;i<elementBlock.length;i++){
				if(elementBlock[i]==ir){
					ire++;
					mapEl[ire]=i;
				}

			}

			regionx[ir].setLastEl(ire);
			regionx[ir].setName(mg.blockName[ir-1]);	
		}


		model.element=new Element[model.numberOfElements+1];

		for(int i=1;i<=model.numberOfElements;i++){
			model.element[i]=new Element("hexahedron");
			for(int j=0;j<model.nElVert;j++)
				model.element[i].setVertNumb(j,mapNd[elVerts[mapEl[i]][j]]);


		}




		List<String> list1=new ArrayList<String>();
		for(int i=1;i<=nRegX;i++){
			list1.add(regionx[i].getName());

		}


		Set<String> set = new HashSet<String>(list1);

		ArrayList<String> regName = new ArrayList<String>(set);

		int nair=0;

		model.numberOfRegions=regName.size();


		boolean joinRegs=false;
		if(model.numberOfRegions<mg.nBlocks)

			joinRegs=true;



		if(!joinRegs){


			model.numberOfRegions=nRegX;
			model.region=new Region[model.numberOfRegions+1];

			int elNumb=0;
			for(int ir=1; ir<=model.numberOfRegions; ir++){
				model.region[ir]=new Region(model.dim);

				model.region[ir].setFirstEl(regionx[ir].getFirstEl());

				model.region[ir].setLastEl(regionx[ir].getLastEl());


				model.region[ir].setName(regionx[ir].getName());



			}
		}

		else{


			for(int i=0; i<model.numberOfRegions; i++){
				if(regName.get(i).startsWith("air"))
					nair++;
			}

			String[] sr1=new String[model.numberOfRegions-nair];
			String[] sr2=new String[nair];

			int i1=0,i2=0;
			for(int i=0; i<model.numberOfRegions; i++)
				if(!regName.get(i).startsWith("air"))
					sr1[i1++]=regName.get(i);
				else
					sr2[i2++]=regName.get(i);

			Arrays.sort(sr1);
			Arrays.sort(sr2);



			ix=0;

			String[] sortedRegions=new String[model.numberOfRegions];
			for(int i=0; i<sr1.length; i++)
				if(sr1[i]!=null)
					sortedRegions[ix++]=sr1[i];
			for(int i=0; i<sr2.length; i++)
				if(sr2[i]!=null)
					sortedRegions[ix++]=sr2[i];


			int[][] br=new int[model.numberOfRegions+1][nRegX];

			for(int ir=1; ir<=model.numberOfRegions; ir++){
				ix=0;
				for(int ib=1; ib<=nRegX;ib++){
					if(list1.get(ib-1).equals(sortedRegions[ir-1])){

						br[ir][ix++]=ib;

					}
				}

			}



			model.region=new Region[model.numberOfRegions+1];
			int[][] elVertsP=new int[model.numberOfElements+1][model.nElVert];
			for(int i=1; i<=model.numberOfElements; i++)
				elVertsP[i]=model.element[i].getVertNumb();

			int elNumb=0;
			for(int ir=1; ir<=model.numberOfRegions; ir++){
				model.region[ir]=new Region(model.dim);

				model.region[ir].setFirstEl(elNumb+1);
				for(int ib=0; (ib<nRegX && br[ir][ib]!=0) ;ib++){
					for(int i=regionx[br[ir][ib]].getFirstEl();i<=regionx[br[ir][ib]].getLastEl();i++)
					{		elNumb++;

					model.element[elNumb].setVertNumb(elVertsP[i]);

					}	
				}

				model.region[ir].setLastEl(elNumb);
			}



			for(int ir=1; ir<=model.numberOfRegions; ir++)
			{

				int i=br[ir][0];
				if(i==nRegX) i=0;

				model.region[ir].setName(sortedRegions[ir-1]);




			}

		}


		if(b)
			model.writeMesh(file);


		return model;



	}






	public void getOrthogMesh( Geometry mg){

		String bun = System.getProperty("user.dir") + "\\orBun.txt";
		getOrthogMesh(mg,bun,false);
	}

	public Model getOrthogMesh( Geometry mg, String file, boolean b){




		Model model=new Model();

		model.nBlocks=mg.nBlocks;

		model.scaleFactor=mg.scaleFactor;
		
		String elType="hexahedron";
		

		if(mg.nBoundary==6){ 

			model.nElVert=8;
			model.dim=3;

		}
		else{
			model.nElVert=4;
			model.dim=2;
			elType="quadrangle";
		}
		
		model.setElType(elType);

		

		Discretizer discretizer=new Discretizer(mg);

		int I,J,K=0;
		double[] X=discretizer.X;
		double[] Y=discretizer.Y;

		I=X.length; J=Y.length;

		int nNodes=I*J;
		int nElements=(I-1)*(J-1); 

		double[] Z=null;

		if(model.dim==3){
			Z=discretizer.Z;	
			K=Z.length;

			nNodes*=K;
			nElements*=(K-1); 
		}




		Vect[] coords=new Vect[nNodes];
		int[][] elVerts=new int[nElements][model.nElVert];

		int[] elementBlock=new int[nElements];

		int ix=0;

		if(model.dim==3){

			for(int k=0;k<K;k++)
				for(int j=0;j<J;j++)
					for(int i=0;i<I;i++)
					{

						coords[ix++]=new Vect(X[i],Y[j],Z[k]);	

					}

			int elementNumber,Nx,Nxy;


			Nx=I;
			Nxy=I*J;
			elementNumber=0;

			for(int k=0;k<K-1;k++)
				for(int j=0;j<J-1;j++)
					for(int i=0;i<I-1;i++){

						elVerts[elementNumber][0]=(k+1)*Nxy+(j+1)*Nx+i+2;
						elVerts[elementNumber][1]=(k+1)*Nxy+(j+1)*Nx+i+1;
						elVerts[elementNumber][2]=(k+1)*Nxy+j*Nx+i+1;
						elVerts[elementNumber][3]=(k+1)*Nxy+j*Nx+i+2;
						elVerts[elementNumber][4]=k*Nxy+(j+1)*Nx+i+2;
						elVerts[elementNumber][5]=k*Nxy+(j+1)*Nx+i+1;
						elVerts[elementNumber][6]=k*Nxy+j*Nx+i+1;
						elVerts[elementNumber][7]=k*Nxy+j*Nx+i+2;

						elementNumber++;


					}

			Vect centerOfMass;

		
			for(int i=0;i<nElements;i++){
				centerOfMass=new Vect(model.dim);
				for(int v=0;v<model.nElVert; v++)
					centerOfMass=centerOfMass.add(coords[elVerts[i][v]-1]);
				centerOfMass=centerOfMass.times(.125);
				for(int ir=0;ir<model.nBlocks;ir++)
					if(discretizer.blockBoundary[ir][0]<centerOfMass.el[0])
						if(centerOfMass.el[0]<discretizer.blockBoundary[ir][1])
							if(discretizer.blockBoundary[ir][2]<centerOfMass.el[1])
								if(centerOfMass.el[1]<discretizer.blockBoundary[ir][3])
									if(discretizer.blockBoundary[ir][4]<centerOfMass.el[2])
										if(centerOfMass.el[2]<discretizer.blockBoundary[ir][5]){
											elementBlock[i]=ir+1;



										}
			} 



		}

		else{	
			ix=0;


			for(int j=0;j<J;j++)
				for(int i=0;i<I;i++)
				{

					coords[ix++]=new Vect(X[i],Y[j]);	

				}

			int elementNumber,Nx;


			Nx=I;
			elementNumber=0;

			for(int j=0;j<J-1;j++)
				for(int i=0;i<I-1;i++){

					elVerts[elementNumber][0]=(j+1)*Nx+i+2;
					elVerts[elementNumber][1]=(j+1)*Nx+i+1;
					elVerts[elementNumber][2]=j*Nx+i+1;
					elVerts[elementNumber][3]=j*Nx+i+2;
					elementNumber++;


				}

			Vect centerOfMass;

			for(int i=0;i<nElements;i++){
				centerOfMass=new Vect(2);
				for(int v=0;v<4; v++)
					centerOfMass=centerOfMass.add(coords[elVerts[i][v]-1]);
				centerOfMass=centerOfMass.times(.25);

				for(int ir=0;ir<model.nBlocks;ir++)
					if(discretizer.blockBoundary[ir][0]<centerOfMass.el[0])
						if(centerOfMass.el[0]<discretizer.blockBoundary[ir][1])
							if(discretizer.blockBoundary[ir][2]<centerOfMass.el[1])
								if(centerOfMass.el[1]<discretizer.blockBoundary[ir][3])
								{
									elementBlock[i]=ir+1;



								}
			} 
		}


		ix=0;
		int iy=0;
		boolean[] nc=new boolean[nNodes+1];


		for(int i=0;i<nElements;i++){
			if(elementBlock[i]>0) {

				++ix;
				for(int j=0;j<model.nElVert;j++)
					if(!nc[elVerts[i][j]]) 
					{
						iy++;
						nc[elVerts[i][j]]=true;
					}
			}

		}

		model.numberOfElements=ix;
		model.numberOfNodes=iy;



		model.node=new Node[model.numberOfNodes+1];

		int[] mapNd=new int[nNodes+1];
		ix=0;
		for(int i=1;i<=nNodes;i++){

			if(nc[i]){
				ix++;
				model.node[ix]=new Node(ix,model.dim);		
				model.node[ix].setCoord(coords[i-1]);	
				mapNd[i]=ix;
			}
		}





		int nRegX=model.nBlocks;
		int[] nRegEls=new int[nRegX+1];
		Region[] regionx=new Region[nRegX+1];



		for(int ir=1;ir<=nRegX;ir++){
			regionx[ir]=new Region(model.dim);
		}

		regionx[1].setFirstEl(1);


		int[] mapEl=new int[nElements+1];
		int ire=0;
		for(int ir=1;ir<=nRegX;ir++){
			if(ir>1) regionx[ir].setFirstEl(regionx[ir-1].getLastEl()+1);
			for(int i=0;i<elementBlock.length;i++){
				if(elementBlock[i]==ir){
					ire++;
					mapEl[ire]=i;
				}

			}

			regionx[ir].setLastEl(ire);
			regionx[ir].setName(mg.blockName[ir-1]);	
		}


		model.element=new Element[model.numberOfElements+1];

		for(int i=1;i<=model.numberOfElements;i++){
			model.element[i]=new Element(model.elType);
			for(int j=0;j<model.nElVert;j++)
				model.element[i].setVertNumb(j,mapNd[elVerts[mapEl[i]][j]]);


		}




		List<String> list1=new ArrayList<String>();
		for(int i=1;i<=nRegX;i++){
			list1.add(regionx[i].getName());

		}


		Set<String> set = new HashSet<String>(list1);

		ArrayList<String> regName = new ArrayList<String>(set);

		int nair=0;

		model.numberOfRegions=regName.size();


		boolean joinRegs=false;
		if(model.numberOfRegions<mg.nBlocks)

			joinRegs=true;



		if(!joinRegs){


			model.numberOfRegions=nRegX;
			model.region=new Region[model.numberOfRegions+1];

			int elNumb=0;
			for(int ir=1; ir<=model.numberOfRegions; ir++){
				model.region[ir]=new Region(model.dim);

				model.region[ir].setFirstEl(regionx[ir].getFirstEl());

				model.region[ir].setLastEl(regionx[ir].getLastEl());


				model.region[ir].setName(regionx[ir].getName());



			}
		}

		else{


			for(int i=0; i<model.numberOfRegions; i++){
				if(regName.get(i).startsWith("air"))
					nair++;
			}

			String[] sr1=new String[model.numberOfRegions-nair];
			String[] sr2=new String[nair];

			int i1=0,i2=0;
			for(int i=0; i<model.numberOfRegions; i++)
				if(!regName.get(i).startsWith("air"))
					sr1[i1++]=regName.get(i);
				else
					sr2[i2++]=regName.get(i);

			Arrays.sort(sr1);
			Arrays.sort(sr2);



			ix=0;

			String[] sortedRegions=new String[model.numberOfRegions];
			for(int i=0; i<sr1.length; i++)
				if(sr1[i]!=null)
					sortedRegions[ix++]=sr1[i];
			for(int i=0; i<sr2.length; i++)
				if(sr2[i]!=null)
					sortedRegions[ix++]=sr2[i];


			int[][] br=new int[model.numberOfRegions+1][nRegX];

			for(int ir=1; ir<=model.numberOfRegions; ir++){
				ix=0;
				for(int ib=1; ib<=nRegX;ib++){
					if(list1.get(ib-1).equals(sortedRegions[ir-1])){

						br[ir][ix++]=ib;

					}
				}

			}



			model.region=new Region[model.numberOfRegions+1];
			int[][] elVertsP=new int[model.numberOfElements+1][model.nElVert];
			for(int i=1; i<=model.numberOfElements; i++)
				elVertsP[i]=model.element[i].getVertNumb();

			int elNumb=0;
			for(int ir=1; ir<=model.numberOfRegions; ir++){
				model.region[ir]=new Region(model.dim);

				model.region[ir].setFirstEl(elNumb+1);
				for(int ib=0; (ib<nRegX && br[ir][ib]!=0) ;ib++){
					for(int i=regionx[br[ir][ib]].getFirstEl();i<=regionx[br[ir][ib]].getLastEl();i++)
					{		elNumb++;

					model.element[elNumb].setVertNumb(elVertsP[i]);

					}	
				}

				model.region[ir].setLastEl(elNumb);
			}



			for(int ir=1; ir<=model.numberOfRegions; ir++)
			{

				int i=br[ir][0];
				if(i==nRegX) i=0;

				model.region[ir].setName(sortedRegions[ir-1]);



			}

		}



		if(b)
			model.writeMesh(file);

		return model;

	}
	
	public void meshCapCorner(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeCapCorner(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(1);



	}
	
	private int  makeNodeCapCorner(Vect[] P){

		double r=184./2;
		double r2=190./2;

		int ix=0;

		P[ix++]=new Vect(2);


		Vect[] crn=new Vect[50];



		ix=addPointsOnArc(P,ix,new Vect(2), r,0,PI/4,216/8);
		
		ix=addPointsOnArc(P,ix,new Vect(70.711,70.711), 5,-3*PI/4,PI/4,4);
		
		//ix=addPointsOnArc(P,ix,new Vect(70.711,70.711), 7.5,-3*PI/4,PI/4,3);
		
		double tx=.00*PI;
	//	ix=addPointsOnArc(P,ix,new Vect(2), r+2.5,tx,PI/4-.08,15);
		tx=.5;
		ix=addPointsOnArc(P,ix,new Vect(2), r+5,.35,PI/4-.1,5);
		ix=addPointsOnArc(P,ix,new Vect(2), r+10,.5,PI/4-.2,2);
		
	/*	tx=.45*PI;
		ix=addPointsOnArc(P,ix,new Vect(2), r+20,tx,PI/2-tx,10);*/
		ix=addPointsOnPath(P,ix,new Vect(r2,0), new Vect(r2,r2/3),5);
		
		ix=addPoint(P,ix,new Vect(70.711,70.711));
		ix=addPointsOnPath(P,ix,new Vect(r2,r2/3), new Vect(r2,r2-3),6);
		
		ix=addPointsOnPath(P,ix,new Vect(r2-8,.65*r2), new Vect(r2-8,r2-15),3);
		
		ix=addPointsOnPath(P,ix,new Vect(r2/Math.sqrt(2),r2/Math.sqrt(2)).times(1.15),  new Vect(r2,r2).times(1-1.5/r2),4);

		return ix;
	}


	public void meshCapCornerLow(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		nNodes=makeNodeCapCornerLow(P);

		double scaleFactor=1000;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(1);



	}
	
	private int  makeNodeCapCornerLow(Vect[] P){

		double r=184./2;
		double r2=190./2;

		int ix=0;

		P[ix++]=new Vect(2);





		ix=addPointsOnArc(P,ix,new Vect(2), r,-PI/2,0,216/4);
		
		ix=addPointsOnArc(P,ix,new Vect(50,-86.603), 5,0,2*PI,8);
		ix=addPoint(P,ix,new Vect(50,-86.603));

		ix=addPointsOnPath(P,ix,new Vect(r2,0), new Vect(r2,-79),12);
		ix=addPointsOnPath(P,ix,new Vect(60,-82), new Vect(92,-82),4);
		ix=addPointsOnPath(P,ix,new Vect(60,-82), new Vect(60,-92),2);
		ix=addPointsOnPath(P,ix,new Vect(60,-92),new Vect(57,-95),1);
		ix=addPointsOnPath(P,ix,new Vect(0,-95),new Vect(57,-95),8);
		
		ix=addPoint(P,ix,new Vect(90,-77));

		
		double tx=.01*PI;
		ix=addPointsOnArc(P,ix,new Vect(2), r+5,-PI/3.3,-PI/10,10);
		
		ix=addPointsOnArc(P,ix,new Vect(2), r+5,-PI/3.3,-PI/10,10);

		ix=addPointsOnArc(P,ix,new Vect(2), r+15,-PI/4,-PI/5,2);
		
		ix=addPointsOnArc(P,ix,new Vect(2), r+5,-PI/2.4,-PI/2.8,4);
		return ix;
	}
	
	
	public void meshHub(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];



		double r1=28;
		double r2=86./2;
		double r3=166./2;
		double r4=177./2;
		double r5=184./2;
		double rx=150./2;
		double rb1=80./2;
		double rb2=50./2;
	//	double rin=12.5;
		

		int ix=0;

		P[ix++]=new Vect(2);

		int ndiv=216;
		

		ix=addPointsOnArc(P,ix,new Vect(2), r5,.0,2*PI,ndiv);


		ix=addPointsOnArc(P,ix,new Vect(2), r4,.0,2*PI,ndiv);
		
		ix=addPointsOnArc(P,ix,new Vect(2), r3,0.0,2*PI,12*8);
		
		ix=addPointsOnArc(P,ix,new Vect(2), rx,0.0,2*PI,8*8);
		ix=addPointsOnArc(P,ix,new Vect(2), r3-20,0.0,2*PI,6*8);
		ix=addPointsOnArc(P,ix,new Vect(2), r2,0.0,2*PI,6*8);
		ix=addPointsOnArc(P,ix,new Vect(2), rb1,0.0,2*PI,4*8);
		ix=addPointsOnArc(P,ix,new Vect(2), rb2,0.0,2*PI,4*8);
		ix=addPointsOnArc(P,ix,new Vect(2), r1,0,2*PI,4*8);

		double scaleFactor=1000;

		nNodes=ix;
		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(1);



	}
	

	public void meshHubFront(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		double r2=150./2;
		double r3=166./2;
		double r4=175./2;
		double r5=184./2;
		
		double rb1=86./2;
		double rb2=62./2;
		
		

		int ix=0;

		P[ix++]=new Vect(2);

		int nd=27;

		ix=addPointsOnArc(P,ix,new Vect(2), r5,.0,2*PI,nd*8);


		ix=addPointsOnArc(P,ix,new Vect(2), r4,.0,2*PI,nd*8);
		
		ix=addPointsOnArc(P,ix,new Vect(2), r3,0.0,2*PI,12*8);
		
		ix=addPointsOnArc(P,ix,new Vect(2), r2,0.0,2*PI,6*8);
		ix=addPointsOnArc(P,ix,new Vect(2), r2-15,0.0,2*PI,3*8);
		ix=addPointsOnArc(P,ix,new Vect(2), rb1,0.0,2*PI,2*8);
		ix=addPointsOnArc(P,ix,new Vect(2), rb2,0.0,2*PI,2*8);


		double scaleFactor=1000;

		nNodes=ix;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(1);



	}
	
	
	public void meshAirgap(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		double r1=54.4;
		double r2=55;
		double rm=(r1+r2)/2;
		

		int ix=0;

		P[ix++]=new Vect(2);

		double f1=PI/18;

	/*	ix=addPointsOnArc(P,ix,new Vect(2), r1,.0,f1,20);
		ix=addPointsOnArc(P,ix,new Vect(2), r1+.1,.0,f1,40);
		ix=addPointsOnArc(P,ix,new Vect(2), r1+.2,.0,f1,60);
		ix=addPointsOnArc(P,ix,new Vect(2), r1+.25,.0,f1,100);
		ix=addPointsOnArc(P,ix,new Vect(2), rm,.0,f1,100);
		ix=addPointsOnArc(P,ix,new Vect(2), r2-.25,.0,f1,100);
		ix=addPointsOnArc(P,ix,new Vect(2), r2-.2,.0,f1,60);
		ix=addPointsOnArc(P,ix,new Vect(2), r2-.1,.0,f1,40);
		ix=addPointsOnArc(P,ix,new Vect(2), r2,.0,f1,20);*/
		
		ix=addPointsOnArc(P,ix,new Vect(2), r1,.0,f1,20);
		ix=addPointsOnArc(P,ix,new Vect(2), r1+.1,.0,f1,50);
		ix=addPointsOnArc(P,ix,new Vect(2), r1+.2,.0,f1,100);
		ix=addPointsOnArc(P,ix,new Vect(2), r1+.25,.0,f1,200);
		ix=addPointsOnArc(P,ix,new Vect(2), rm,.0,f1,200);
		ix=addPointsOnArc(P,ix,new Vect(2), r2-.25,.0,f1,200);
		ix=addPointsOnArc(P,ix,new Vect(2), r2-.2,.0,f1,100);
		ix=addPointsOnArc(P,ix,new Vect(2), r2-.1,.0,f1,50);
		ix=addPointsOnArc(P,ix,new Vect(2), r2,.0,f1,20);

		double scaleFactor=1000;

		nNodes=ix;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(1);



	}
	
	public void meshRotAirgapEnd(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];


		double r1=54.65;
		double r2=54.7;
		double rm=(r1+r2)/2;
		

		int ix=0;

		P[ix++]=new Vect(2);

		double f1=PI/18;

		ix=addPointsOnArc(P,ix,new Vect(2), r1,.0,f1,100);
		ix=addPointsOnArc(P,ix,new Vect(2), rm,.0,f1,200);
		ix=addPointsOnArc(P,ix,new Vect(2), r2,.0,f1,200);

		double scaleFactor=1000;

		nNodes=ix;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(1);



	}
	
	public void meshAirgap2(){

		int nNodes=10000;
		Vect[] P=new Vect[nNodes];

		String bun = System.getProperty("user.dir") + "//assembled.txt";

		Model model=new Model(bun);
		
		double r1=54.4;
		double r2=55;
		double rm=(r1+r2)/2;
		

		int ix=0;

		P[ix++]=new Vect(2);

		for(int i=1;i<=model.numberOfNodes;i++){
			Vect v=model.node[i].getCoord();
			double r=v.norm();
			if(util.getAng(v)<PI/17.99&& r>.0548999 && r<.055001) P[ix++]=v.deepCopy();
		}
		
		ix=addPointsOnArc(P,ix,new Vect(2), .05495,.0,PI/18,10);
		ix=addPointsOnArc(P,ix,new Vect(2), .054,.0,PI/18,30);
		ix=addPointsOnArc(P,ix,new Vect(2), .052,.0,PI/18,20);
		ix=addPointsOnArc(P,ix,new Vect(2), .050,.0,PI/18,20);
		ix=addPointsOnArc(P,ix,new Vect(2), .04,.0,PI/18,20);
		ix=addPointsOnArc(P,ix,new Vect(2), .02,.0,PI/18,20);




		double scaleFactor=1000;

		nNodes=ix;

		writeNodesForDelaunay(P,nNodes,scaleFactor);

		String nodeFile = System.getProperty("user.dir") + "//node4Del.node";

		wait(4*nNodes);


		Delaunay(nodeFile);

		wait(4*nNodes);

		meshFromDel(1);



	}
	
	public void RCM(){
		String bun=util.getFile();
		if(bun==null || bun.equals("") )return;
		Model model=new Model(bun);


		RCM(model);
		
	}
	
	public void RCM(Model model){
		
		SpMat Rs=this.adjMat(model);
		int[] rr=Rs.RCMorder();
		
		int[] map=new int[rr.length];
		for(int i=0;i<rr.length;i++){
			map[rr[i]]=i;
			
		
		}
		

		
					Model rmc=model.deepCopy();
					
				for(int i=1;i<=rmc.numberOfElements;i++){						
					int[] vertNumb=model.element[i].getVertNumb();
						for(int j=0;j<rmc.nElVert;j++){
							
						int nn=map[vertNumb[j]-1]+1;
							rmc.element[i].setVertNumb(j, nn);
						}
				}
			
				for(int i=1;i<=rmc.numberOfNodes;i++){
					int nn=map[i-1]+1;
					rmc.node[nn].setCoord(model.node[i].getCoord());
				}
				
				String bun=System.getProperty("user.dir") + "\\reNumb.txt";
				rmc.writeMesh(bun);

							

				}
	
	
				

	public SpMat adjMat(){
		return adjMat(this.getModel());
		//String bun=System.getProperty("user.dir") + "\\cube3.txt";
		//return adjMat(new Model(bun));
	}
	
	public SpMat adjMat(Model model){

			System.out.println(" Calculating adjacency matrix ...");
		
			int m,row,col,ext=4;
			int[] nz=new int[model.numberOfNodes];
			SpMat Rs=new SpMat(model.numberOfNodes,model.numberOfNodes,model.nNodNod);


			for(int i=1;i<=model.numberOfElements;i++){

				int[] vertNumb=model.element[i].getVertNumb();

				for(int j=0;j<model.nElVert;j++){

					row=vertNumb[j]-1;

					for(int k=0;k<model.nElVert;k++){

						col=vertNumb[k]-1;

						m=util.search(Rs.row[row].index,nz[row]-1,col);

						if(m<0)
						{	

							Rs.row[row].index[nz[row]]=col;															
				
							nz[row]++;

							if(nz[row]==Rs.row[row].nzLength-1){
								Rs.row[row].extend(ext);
							}

						}

			
					}

				}
			}


			Rs.sortAndTrim(nz);


			//Rs.band().plot();
return Rs;


		}
	
	
	public void meshStatQ(){
		
		
	//	double[][] bb={{-6,6,0,6},{-2.5,2.5,0,4},{-1.5,1.5,0,3},{-3.6,-2.6,0,2},{-1.4,-.4,0,2}};


		double[][] bb={{55,55.7,4.3,10.9},{55.7,77,5.1,10.1},{77,87.5,0,15.2},{87.5,92,0,15.2}};


		double scale=1;

		for(int j=0;j<bb.length;j++)
			for(int k=0;k<bb[0].length;k++){
				bb[j][k]*=scale;
				if(k>1)
				bb[j][k]-=7.6;

				//if(k==2&& bb[j][k]==0) 	bb[j][k]=.4;
			}

		Geometry mg=new Geometry(bb);
/*			mg.blockName[0]="air";
		mg.blockName[2]="air";

				mg.blockName[4]="air";*/

		mg.blockName[0]="stat";
		mg.blockName[1]="stat";
		mg.blockName[2]="stat";
		mg.blockName[3]="frame";

		
	

		for(int j=0;j<bb.length;j++){
			
			for(int k=0;k<bb[0].length;k++){
			mg.baseLeft[j][k]=1.2;
			mg.baseLeft[j][k]=1.2;
			mg.baseLeft[j][k]=1.2;
			mg.baseRight[j][k]=1.2;
			
			if((j==3 && k<1) || (j==2 && k==1)){
				mg.minMeshRight[j][k]=.05;
				mg.minMeshLeft[j][k]=.05;
	

			}
			else{
				mg.minMeshRight[j][k]=.7;

				mg.minMeshLeft[j][k]=.7;
			}
			}
		}



		Model model=getOrthogMesh2D(mg,"",false);
		
		for(int i=1;i<=1*model.numberOfNodes;i++)
		{
			if(model.node[i].getCoord(0)<-55.1e-3 || model.node[i].getCoord(0)>76.8e-3){
			double rh=model.node[i].getCoord(1)*1000/7.6;
			 
			double tt=PI/36*rh;
			util.pr(tt);
			double r=model.node[i].getCoord(0);
			Vect v=new Vect(r*Math.cos(tt),r*Math.sin(tt));
			model.node[i].setCoord(v);
			}


		}

		
		
		String bun=System.getProperty("user.dir") + "\\model2D.txt";
		model.writeMesh(bun);

	

}
	
	
	public void meshQ(){
		
		
			//double[][] bb={{471,492.32,0,12.5},{478.43,481.23,.55,6.25},{478.43,481.23,.55,6.25},{478.43,481.23,.55,6.25},{478.43,481.23,.55,6.25}};
		//	double[][] bb={{0,1000,-1000,1000},{-1000,1000,-1000,-300},{-400,-300,0,100},{300,400,0,100}};
			double[][] bb={{-300,300,-300,300},{-60,60,-60,60},{-40,40,-40,40},{40,60,-5,5}, {-35,0,-35,35},{-100,-65,-35,35}};
			
			/*
			double[] bair={470,700,0,12.5};
			double[] balum={0,2.8,0,5.7};
			 double[][] bb1=new double[4][4];
			 
			 for(int j=0;j<2;j++)
			 for(int k=0;k<2;k++)
					for(int m=0;m<4;m++)
						if(m<2)
						bb1[j*2+k][m]=balum[m]+j*3.66;
						else
							bb1[j*2+k][m]=balum[m]+k*(5.7+.06)+.55;
			 
			 int nc=11;
			 
			 double[][] bb12=new double[nc*4+1][4];
		
						
			 
			 for(int k=0;k<4;k++)
				 bb12[0][k]=bair[k];
			 
				for(int j=0;j<nc;j++)
					for(int k=0;k<4;k++)
						for(int m=0;m<4;m++){
							if(m<2)
						bb12[j*4+k+1][m]=bb1[k][m]+j*(14+3.66)+478;
							else
								bb12[j*4+k+1][m]=bb1[k][m];
					}*/
			 
/*			double[][] bb={{-100,100,-100,100},{-80,80,-80,80},
					{-50,-15,-50,50},{15,50,-50,50}
					, {-15,15,-1,1},{-45,-20,-40,40},{20,45,-40,40}};*/


			double scale=1;
		//	scale=1000;

			for(int j=0;j<bb.length;j++)
				for(int k=0;k<bb[0].length;k++){
					bb[j][k]*=scale;
				/*if(k<2)
					bb12[j][k]+=-400;*/
					//if(k==2&& bb[j][k]==0) 	bb[j][k]=.4;
				}

			Geometry mg=new Geometry(bb);
	/*			mg.blockName[0]="air";
			mg.blockName[2]="air";

					mg.blockName[4]="air";*/

			mg.blockName[0]="air";
		//	for(int j=1;j<bb12.length;j++)
				//mg.blockName[j]="coil";
	

			mg.blockName[1]="core";
			mg.blockName[2]="air";
			mg.blockName[3]="air";
			mg.blockName[4]="coil1";
			mg.blockName[5]="coil2";

			
for(int j=0;j<bb.length;j++){
				
				for(int k=0;k<bb[0].length;k++){
				
			
				mg.baseLeft[j][k]=1.3;
				mg.baseRight[j][k]=1.3;
					
			
			
		
				if(j<2){
				mg.minMeshRight[j][k]=1;
				mg.minMeshLeft[j][k]=1;
				}
				else{
					mg.minMeshRight[j][k]=1;
					mg.minMeshLeft[j][k]=1;
					
				}
				}
			
				}			

		

/*			for(int j=0;j<bb12.length;j++){
				
				for(int k=0;k<bb12[0].length;k++){
				
				
				mg.baseLeft[j][k]=1.3;
				mg.baseRight[j][k]=1.3;
				
				if(j>0){
		
				mg.minMeshRight[j][k]=1.;
				mg.minMeshLeft[j][k]=1.;
				
	
					
				if(k<2){
					mg.minMeshRight[j][k]=.9;
					mg.minMeshLeft[j][k]=.9;
				}
					}
			
			
				
				}
			
				}			
*/

			


			Model model=getOrthogMesh2D(mg,"",false);
			


			
			
			String bun=System.getProperty("user.dir") + "\\model2D.txt";
			model.writeMesh(bun);

		

	}
	

	
	public void meshQx(){
		
			//double[][] bb={{50,65,0,5},{65,67,0,5}};
		double[][] bb={{88,89,0,5}};
			//double[][] bb={{67,87,0,5},{87,89,0,5}};
			//double[][] bb={{89,99,0,5}};
		
		

		double scale=1;
	//	scale=1000;

		for(int j=0;j<bb.length;j++)
			for(int k=0;k<bb[0].length;k++){
				bb[j][k]*=scale;
			
			}

		Geometry mg=new Geometry(bb);

		mg.blockName[0]="outPM";
		//mg.blockName[1]="inAir2";
	

		
for(int j=0;j<bb.length;j++){
			
			for(int k=0;k<bb[0].length/2;k++){
			
		
			mg.baseLeft[j][k]=1;
			mg.baseRight[j][k]=1;
	
			if(j<1){
				mg.minMeshRight[j][k]=.24;
				mg.minMeshLeft[j][k]=.26;
			} else {
				mg.minMeshRight[j][k]=.3;
				mg.minMeshLeft[j][k]=.3;
			}
		
		
			}	
}


		


		Model model=getOrthogMesh2D(mg,"",false);

		double tt=1.0*Math.PI/180;
		
		for(int i=1; i<=model.numberOfNodes;i++){
			if(model.node[i].getCoord(1)>1e-4){
				
			double r=model.node[i].getCoord(0);
				model.node[i].setCoord(0,r*Math.cos(tt));
				model.node[i].setCoord(1,r*Math.sin(tt));
				
			}
		}
/*
		model.motor=true;
		
		Model model2=this.rotExtendNfold(model,1);*/
		
		
		String bun=System.getProperty("user.dir") + "\\model2D.txt";
		model.writeMesh(bun);

	

}
	
	public void meshQz(){
		
		  int nr=20;
		  int nt=361;
		   int nNodes=nt*nr;
		   double r0=.020;
		   double r1=.080;
		   double dr=(r1-r0)/nr;
		   double dtt=2*PI/(nt-1);
		   int ne=(nt-1)*(nr-1);
		   
		   Model model=new Model();
		   model.alloc(1, ne, nNodes, "quadrangle");
		   int nx=0;
		   for(int i=0;i<nt;i++){
			   double tt=i*dtt;
			   for(int j=0;j<nr;j++){
				   nx++;
				   double r=r0+j*dr;
				   Vect v=new Vect(r*Math.cos(tt),r*Math.sin(tt));
				   model.node[nx].setCoord(v);
				   if(i<nt-1 && j<nr-1){
					   {
						   int kx=i*(nr-1)+j+1;
						   model.element[kx].setVertNumb(0, i*nr+j+1);  
						   model.element[kx].setVertNumb(1, i*nr+j+1+1); 
						   model.element[kx].setVertNumb(2, (i+1)*nr+j+1+1); 
						   model.element[kx].setVertNumb(3, (i+1)*nr+j+1); 
					   }
				  
			   }
			   
		   }
	}
		   
		   model.region[1].setFirstEl(1);
		   model.region[1].setLastEl(ne);

		   

		   model.scaleFactor=1000;

	
	String bun=System.getProperty("user.dir") + "\\model2D.txt";
	model.writeMesh(bun);



}

	public void getNeuMeshQ(int mode){
		String regex="[ ,\\t]+";
		String s=util.getFile();
		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];

	if(mode==0)
	while(!br.readLine().startsWith("<<< End Solid Transmit <<<")){		}

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
	
	
	
	public void getNeuMeshHexa(int mode){
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


			int[][] vernumb=new int[10*nnMax+1][8];
			int[] nReg=new int[1000000+1];
			
			for(int i=0;i<nReg.length;i++)
				nReg[i]=-1;
				
			int ix=0;
		
				line=br.readLine();
			//	line=br.readLine();
	for(int i=0;i<7;i++)
		line=br.readLine();

			sp=new String[13];
			for(int i=1;i<=vernumb.length;i++)
			{
				
				linep=br.readLine();
				line=br.readLine();
	
				
				sp=line.split(regex);
	
				if(sp.length<5) break;
				
				if(sp[6].equals("0")) {
					for(int j=0;j<5;j++)
						line=br.readLine();
					continue;
				}
				

				
				sp=linep.split(regex);
				ix++;
				nReg[ix]=Integer.parseInt(sp[2]);
			
			
				sp=line.split(regex);


				for(int j=0;j<8;j++){

					vernumb[ix][j]=map[Integer.parseInt(sp[j])];
				}
				
			

				for(int j=0;j<8;j++){
					if(vernumb[ix][j]==0) {
						if(j>0) vernumb[ix][j]=vernumb[ix][j-1];
						ix--;
						break;
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


			int nNodes=nnMax;
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
			for(int i=1;i<=nnMax;i++){ 
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
	
	
	public void twistNeuMeshToMatchHelical(double radius,double pitchFactor){
		String regex="[ ,\\t]+";
		String bun=util.getFile();
		try{
			File f=new File(bun);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
			

			String folder=new File(bun).getParentFile().getPath();
			String fout=folder+"\\twisted.neu";


			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));	
			
			line=br.readLine();
			pwBun.println(line);
			
			boolean atNodeData=false;
			
			while(true)

			{
		
				line=br.readLine();
				
			if(line==null || line.equals("")) {
				//pwBun.println(line);
				break;
			}
				
			//	ix++;
				//if(ix==38972) util.pr(line);
				
				//util.pr(ix);
				sp=line.split(regex);
			//if(!atNodeData)

				if(sp.length==15) {
			//	util.pr

				Vect v3d=new Vect(Double.parseDouble(sp[11]),Double.parseDouble(sp[12]),Double.parseDouble(sp[13]));
				Vect v2d=v3d.v2();
			//	double ang=util.getAng(v2d);

				if(abs(v2d.norm()-radius)<1e-3){
					//double tt=10*v3d.el[2]*pitchFactor/180;
					int nl=(int)Math.floor(v3d.el[2]*1.01/.05);
				if(nl%2==1){
					double tt=2*pitchFactor/180;
				
					Mat R2D=util.rotMat2D(tt);
					
					Vect v2dRotated=R2D.mul(v2d);
					
					v2dRotated.sub(v2d).hshow();
					sp[11]=Double.toString(v2dRotated.el[0]);
					sp[12]=Double.toString(v2dRotated.el[1]);
					
					line=sp[0]+",";
					for(int j=1;j<sp.length;j++)
						line=line+sp[j]+",";
					util.pr(line);
					System.out.println();
					System.out.println(" Node was rotated by "+tt+" radians");
				}
					}
				
			
				}
			
				
				pwBun.println(line);

				
			}
				pwBun.flush();
				
			System.out.println();
			System.out.println(" Bun data was written to:");
			System.out.println("    "+fout);
			System.out.println();
		
		
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
	
	
	public void getPostMeshHex(){
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
					if(vernumb[ix][j]==0) {

						numbAddedNodes++;
				
						 vernumb[ix][j]=nnMax+numbAddedNodes;//vernumb[ix][j-1];
							if(j>=0){
						 coord1[numbAddedNodes]=coord1[vernumb[ix][j-1]].deepCopy().times(0);
							}
						//ix--;
						//break;
					}
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
	
	public void getPostMeshHexAtlas(){

		String s=util.getFile();
		
		String line;
		int max=1000000;
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
		
			model.scaleFactor=1;
			
			for(int i=1;i<=nEls;i++)
			model.element[i].setRegion(regMap[regNumb[i]]);
			
			reRegionGroupEls(model);
			
			String fout=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";
			
			model.writeMesh(fout);



	}
	
	public void getPostMeshHexAtlasOrig(){

		String s=util.getFile();
		
		String line;
		int max=1000000;
		Vect[] coord1=new Vect[max];
		int[][] vertNumb=new int[max][8];
		int[] nodeNumb=new int[max];
		int[] elNumb=new int[max];
		int[] regNumb=new int[max];


		int nNodes=0,nEls=0;

	
		for(int i=0;i<regNumb.length;i++)
			regNumb[i]=-1;


		try{
			File f=new File(s);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);

	
			while(true){
				line=br.readLine();
				
				if(line==null) break;
				if(line.startsWith("GRID")){
					nNodes=readAtlasNodesOrig(br,nodeNumb,coord1,nNodes);


				}
				

				
				if(line.startsWith("CONC")){	
		
					nEls=readAtlasElementsOrig(br,elNumb,regNumb,vertNumb,nEls);
				}

			}
		
		}
			catch(Exception e){System.err.println("error");	e.printStackTrace(); }

			int nNodeMax=0;
			for(int i=0;i<nNodes;i++)
				if(nodeNumb[i]>nNodeMax) nNodeMax=nodeNumb[i];
			
			int[] vertNumb0=new int[8];
			int jx=0;
			for(int i=0;i<nNodes;i++){
				if(nodeNumb[i]>0) vertNumb0[jx++]=nodeNumb[i];
				if(jx==8) break;
			}

			int nElMax=0;
			for(int i=0;i<nEls;i++)
				if(elNumb[i]>nElMax) nElMax=elNumb[i];
			
			int nRegMax=0;
			for(int i=0;i<regNumb.length;i++)
				if(regNumb[i]>nRegMax) nRegMax=regNumb[i];
			
			
			int[] nRegEls=new int[nRegMax+1];
			for(int i=0;i<regNumb.length;i++)
				if(regNumb[i]>0) nRegEls[regNumb[i]]++;
						
			int[][]regEnd=new int[nRegMax+1][2];
			regEnd[0][0]=1;

					
			for(int i=0;i<nRegEls.length;i++){
		
				if(i>0)
				regEnd[i][0]=regEnd[i-1][1]+1;
				
					regEnd[i][1]=regEnd[i][0]+nRegEls[i]-1;
				}		
			

			
			Model model=new Model(nRegMax,nElMax,nNodeMax,"hexahedron");
			
			for(int i=1;i<=nElMax;i++){
				if(vertNumb[i][0]>0){
					model.element[i].setVertNumb(vertNumb[i]);
				}
				else{
					model.element[i].setVertNumb(vertNumb0);
				}
			}
				
				
			for(int ir=1;ir<=nRegMax;ir++)
			{
				model.region[ir].setFirstEl(regEnd[ir][0]);;
				model.region[ir].setLastEl(regEnd[ir][1]);;
				model.region[ir].setName("region"+ir);
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
	/*			for(int j=0;j<8;j++){
					int mpn=vertNumb[i][j];

					model.element[i].setVertNumb(j,mpn);
					
				}*/
				model.element[i].setRegion(ir);

			}
			
		
			}
		
			util.pr(nRegMax);
			
			for(int i=1;i<=nNodeMax;i++){
				if(coord1[i]==null)
					model.node[i].setCoord(new Vect(3));
				else
				model.node[i].setCoord(coord1[i]);
				}
		
			model.scaleFactor=1;
			
			reRegionGroupEls(model);
			

			
			String fout=System.getProperty("user.dir")+"\\EMSol\\HexaOrig.txt";
			
			model.writeMesh(fout);



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
			if(sp1.length<3) {return nEl;}

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
			if(sp1.length<3) return nEl;

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
	
	
	public void getEMSolFlux(String bbf,int dim, int numb){

		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
		
			Mat BB=new Mat(100*1000,dim);
			
			for(int tt=0;tt<numb;tt++){

				line="";
				while(!line.startsWith("STEP")){
					line=br.readLine();
					
					}
				int[] nx=new int[dim];;

				String fout=System.getProperty("user.dir")+"\\EMSol\\flux"+tt+".txt";
				
				PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));	
				
	while(!line.startsWith("BMAG")){
	line=br.readLine();
			}
	
	for(int k=0;k<6;k++)
		line=br.readLine();
	
	int k=0;

	
			for(int i=1;i<180000;i++){
				
	
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
	
	
	
	
	
	public void modifyEMSolFlux(int dim){

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
			
			
			Mat BB=new Mat(100*1000,dim);
			
			int[] nx=new int[dim];;
			
			
	while(!line.startsWith("BMAG")){
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
	
	



public void pileHelic(double totalTwistDeg, double height,int K){
	
	double dtt=totalTwistDeg/K/180*PI;
	
	double dheight=height/K;

	String bun=util.getFile();
	if(bun==null || bun.equals("") )return;
	Model model=new Model(bun);

	if(model.elCode>1) return;
	int nv=3;
	String et="prism";
	if(model.elCode==1) {
		nv=4;
		et="hexahedron";
	}

	

	int nReg=model.numberOfRegions;
	int nEls=model.numberOfElements*(K);
	int nNodes=model.numberOfNodes*(K+1);

	Model prismModel=new Model();
	prismModel.alloc(nReg, nEls, nNodes,et);

	for(int j=0;j<=K;j++){
		
		double tt=j*dtt;
		Mat R2D=util.rotMat2D(tt);
		Mat R=new Mat(3,3);
		R.el[2][2]=1;

		for(int m=0;m<2;m++)
			for(int n=0;n<2;n++)
				R.el[m][n]=R2D.el[m][n];
		

		for(int i=1;i<=model.numberOfNodes;i++){
			Vect P=model.node[i].getCoord();
			

			Vect Pp=R.mul(P.v3());
			
		
			double hh=j*dheight;//pitch*(tt);
			
		
			Pp.el[2]=Pp.el[2]+hh;

			prismModel.node[j*model.numberOfNodes+i].setCoord(Pp);

		}
	}
	int[] vertNumbP=new int[2*nv];
	int nn=model.numberOfNodes;


	int net=0;
	for(int ir=1;ir<=model.numberOfRegions;ir++){

		for(int j=0;j<K;j++){

			boolean touch=(abs((j+1)*dtt-2*PI)<-1e-3);
			for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				net++;

				int[] vertNumb=model.element[i].getVertNumb();

				for(int k=0;k<nv;k++)
					if(!touch)
						vertNumbP[k]=vertNumb[k]+(j+1)*nn;
					else
						vertNumbP[k]=vertNumb[k];
				for(int k=nv;k<2*nv;k++)
					vertNumbP[k]=vertNumb[k-nv]+j*nn;

				prismModel.element[net].setVertNumb(vertNumbP);
			}
		}
	}

	prismModel.region[1].setFirstEl(1);
	prismModel.region[1].setLastEl(K*model.region[1].getLastEl());
	prismModel.region[1].setName(model.region[1].getName());

	for(int i=2;i<=model.numberOfRegions;i++){

		prismModel.region[i].setFirstEl(prismModel.region[i-1].getLastEl()+1);
		prismModel.region[i].setLastEl(prismModel.region[i].getFirstEl()+K*model.region[i].getNumbElements()-1);

		prismModel.region[i].setName(model.region[i].getName());
	}

	prismModel.scaleFactor=	model.scaleFactor;



		//Mat R=util.rotMat(new Vect(0,0,1), new Vect(1,0,0));	
		
		//R=util.rotEuler(new Vect(0,0,1), -PI/2).mul(R);
	
	//	rotate(prismModel,R,true);

		
		//rotate(-90*PI/180);
	
	String folder=new File(bun).getParentFile().getPath();

	String prismMesh = folder + "//prismatic.txt";
	prismModel.writeMesh(prismMesh);

}


public void renumberNeuMeshElNodes(){

	
	boolean writeFile=true;
	
	PrintWriter pw=null;
	
	String regex="[ ,\\t]+";
	String bun=util.getFile();
	
	String[] lines=new String[7];

	String line="";
	String[] sp=new String[15];
	
	int nElems=0;
	int maxElemNumb=0;
	
	
	int nNodes=0;
	int maxNodeNumb=0;
	

	
	int[] nodeMap=new int[10000*10000];
	
	try{
		File f=new File(bun);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);

		
		

		

		String folder=new File(bun).getParentFile().getPath();
		String fout=folder+"\\renumbered.neu";
		if(writeFile)
		pw = new PrintWriter(new BufferedWriter(new FileWriter(fout)));	
		
		line=br.readLine();
		if(writeFile)
		pw.println(line);


		while(!line.startsWith("   403")){

			line=br.readLine();
			if(writeFile)
				pw.println(line);
			if(line==null) break;

		}
		

		
		while(true){
			line=br.readLine();
			if(line==null) break;
			sp=line.split(regex);
		//	util.pr(line);
		//	util.pr(sp.length);

			//if(sp.length!=15) break;

			sp=line.split(regex);

			
			if(!util.first(line).equals("-1")) {
				int ib=0;
				if(sp[0].equals("")) ib=1;
				int nn=Integer.parseInt(sp[ib]);
				nNodes++;
				nodeMap[nn]=nNodes;
				
				sp[ib]=Integer.toString(nNodes);
				
				
				line=sp[0]+",";
				for(int j=1;j<sp.length;j++)
					line=line+sp[j]+",";

				if(writeFile)
					pw.println(line);
	
				if(nn>maxNodeNumb) maxNodeNumb=nn;
				
			} 
			else break;
			}
		
		if(writeFile)
			pw.println(line);
		
		System.out.println("no_nodes ----  nnodeMax    "+nNodes+"    "+maxNodeNumb);

		while(!line.startsWith("   404")){

			line=br.readLine();
			
			if(writeFile)
				pw.println(line);
			
			if(line==null) break;

		}


		int[] nVerts=new int[8];
		
		while(true){

			for(int i=0;i<7;i++){
				lines[i]=br.readLine();

			}

			if(util.first(lines[0]).equals("-1")) {
				break;
			}
			
			sp=lines[0].split(regex);
			int ib=0;
			if(sp[0].equals("")) ib=1;
			
			int ne=Integer.parseInt(sp[ib]);
			
			nElems++;
			
			sp[ib]=Integer.toString(nElems);
			lines[0]=sp[ib]+",";
			for(int j=1;j<sp.length;j++)
				lines[0]=lines[0]+sp[j]+",";
			
			if(ne>maxElemNumb) maxElemNumb=ne;
			
			
			sp=lines[1].split(regex);

			ib=0;
			if(sp[0].equals("")) ib=1;
			for(int i=0;i<8;i++){
				nVerts[i]=nodeMap[Integer.parseInt(sp[i+ib])];
				sp[i+ib]=Integer.toString(nVerts[i]);
			}
			lines[1]=sp[ib]+",";
			for(int j=1;j<sp.length;j++)
				lines[1]=lines[1]+sp[j]+",";


			if(writeFile)
				for(int i=0;i<7;i++)
					pw.println(lines[i]);		


		}
		
		while(true){

			line=br.readLine();
			if(line==null) break;
	
			if(writeFile)
				pw.println(line);
			

		}


		pw.flush();

		System.out.println("nElems ----  maxElemNumb    "+nElems+"    "+maxElemNumb);

		System.out.println();
		System.out.println(" Bun data was written to:");
		System.out.println("    "+fout);
		System.out.println();


		br.close();
		fr.close();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	


}



public void renumberNeuMeshRegions(){

	boolean twoD=true;
	boolean byCoord=true;

	//int[][] hexVert=new int[1000000][8];
	
	int[] NodeIndex=new int[2000*10000];
	Vect[] coords=new Vect[2000*100000];
	
	int[] hexVert=new int[8];
//	Vect angs=new Vect().linspace(0, 2*PI, 7);
	Vect angs=new Vect().linspace(0, 2*PI, 6);
	//	angs.hshow();
	//	angs.hshow();
	int nElems=0;

	int nNodes=0;
	int nnMax=0;
	String file=util.getFile();
	String folder=new File(file).getParentFile().getPath();
	String reRegFile=folder+"\\reRegion.neu";

	String line="";


	String[] sp;


	try{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reRegFile)));
		
		File f=new File(file);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		while(!line.startsWith("   403")){

			line=br.readLine();
			pw.println(line);
			if(line==null) break;

		}
		

		
		while(true){
			line=br.readLine();
			pw.println(line);
			sp=line.split(regex);
			
			

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

			
			sp=lines[1].split(regex);

			int ib=0;
			if(sp[0].equals("")) ib=1;
			for(int i=0;i<8;i++)
				hexVert[i]=Integer.parseInt(sp[i+ib]);

		
			
			
			sp=lines[0].split(regex);
			ib=0;
			if(sp[0].equals("")) ib=1;

	
			int nReg=Integer.parseInt(sp[ib+2]);
			int nRegNew=nReg;

			if(byCoord)	{
				
				Vect center=new Vect(3);
				int nVerts=0;
				for(int i=0;i<8;i++){
					int nn=hexVert[i];
					if(nn>0){
						center=center.add(coords[NodeIndex[nn]]);
						nVerts++;
					}
				}
	
			center=center.times(1./nVerts);
			
			double ang=util.getAng(center.v2());
			int angBlock=0;
			for(int k=0;k<angs.length-1;k++){

				if(ang>=angs.el[k] && ang<angs.el[k+1]){
					angBlock=k;
					break;
				}
			}
			

			
			if(nRegNew>=111 & nRegNew<=120){
				//int n=nRegNew/100;
				//nRegNew-=n*100;
				//nRegNew=1000+n*100+nRegNew+angBlock*1000;
				
				nRegNew-=100;
				nRegNew=100+nRegNew+angBlock*100;
				util.pr(" DDD ------------->>  "+nRegNew);
			}
			
			if(nReg==1){
				int ix=(int)((center.el[0]+.5)/1.);
				int iy=(int)((center.el[1]+.5)/1.);
				nRegNew=(iy+1)*100+ix+1;
				util.pr(" SSS ------------->>  "+nRegNew);
			}
			else
				nRegNew=nReg;
				
			
		//	if(nReg>1009 && nReg<1117){
	/*			if(nRegNew>=40010 & nRegNew<=40013){
				nRegNew-=40000;
				nRegNew=nRegNew+angBlock*10;
				util.pr(" AAA ------------->>  "+nRegNew);
			}	else if(nReg>=1010000&& nReg<=1012000){
				nRegNew-=1000;
				nRegNew=nRegNew+1000+angBlock*10;
				util.pr(" BBB ------------->>  "+nRegNew);
			}
			else if(nReg==1120){
				nRegNew=nRegNew+angBlock*100;
				util.pr(" CCC ------------->>  "+nRegNew);
			}
			if(nRegNew>=111 & nRegNew<=120){
				nRegNew-=100;
				nRegNew=100+nRegNew+angBlock*100;
				util.pr(" DDD ------------->>  "+nRegNew);
			}
			}	else if(nReg>=1010000&& nReg<=1012000){*/
			}else{
				if(twoD){
				if(nRegNew<=5) nRegNew=98+nRegNew;
			//	else if(nRegNew>=111 && nRegNew<=120)
				else{
					
					if(nReg==14)nRegNew=120;
					else if(nReg==23)nRegNew=220;
					else if(nReg==32)nRegNew=320;
					else if(nReg==41)nRegNew=420;
					else if(nReg==50)nRegNew=520;
					else if(nReg==59)nRegNew=620;
					
					else if(nReg>=6&& nReg<=14)
					nRegNew=100+11+nRegNew-6;
					else if(nReg>=15&& nReg<=23)
						nRegNew=200+11+nRegNew-15;
					else if(nReg>=24&& nReg<=32)
						nRegNew=300+11+nRegNew-24;
					else if(nReg>=33&& nReg<=41)
						nRegNew=400+11+nRegNew-33;
					else if(nReg>=42&& nReg<=50)
						nRegNew=500+11+nRegNew-42;
					else if(nReg>=51&& nReg<=59)
						nRegNew=600+11+nRegNew-51;
					util.pr(" EEE ------------->>  "+nRegNew);
					
				}
/*				else if(nRegNew==7) nRegNew=120;
				else if(nRegNew==8) nRegNew=218;
				else if(nRegNew==9) nRegNew=220;
				else if(nRegNew==10) nRegNew=318;
				else if(nRegNew==11) nRegNew=320;
				else if(nRegNew==12) nRegNew=418;
				else if(nRegNew==13) nRegNew=420;
				else if(nRegNew==14) nRegNew=518;
				else if(nRegNew==15) nRegNew=520;
				else if(nRegNew==16) nRegNew=618;
				else if(nRegNew==17) nRegNew=620;
				else if(nRegNew<=26)nRegNew=111+nRegNew-20;
				else if(nRegNew<=34)nRegNew=211+nRegNew-27;
				else if(nRegNew<=42)nRegNew=311+nRegNew-35;
				else if(nRegNew<=50)nRegNew=411+nRegNew-43;
				else if(nRegNew<=58)nRegNew=511+nRegNew-51;
				else if(nRegNew<=66)nRegNew=611+nRegNew-59;*/
				}else{
					if(nRegNew<=5) nRegNew=998+nRegNew;
					else if(nRegNew<=13) nRegNew=1110+nRegNew-5;
					else if(nRegNew==14) nRegNew=1120;
					else if(nRegNew<=22) nRegNew=1210+nRegNew-14;
					else if(nRegNew==23) nRegNew=1220;
					else if(nRegNew<=31) nRegNew=1310+nRegNew-23;
					else if(nRegNew==32) nRegNew=1320;
					else if(nRegNew<=40) nRegNew=1410+nRegNew-32;
					else if(nRegNew==41) nRegNew=1420;
					else if(nRegNew<=49) nRegNew=1510+nRegNew-41;
					else if(nRegNew==50) nRegNew=1520;
					else if(nRegNew<=58) nRegNew=1610+nRegNew-50;
					else if(nRegNew==59) nRegNew=1620;
					else if(nRegNew<=26)nRegNew=1111+nRegNew-20;
					else if(nRegNew<=34)nRegNew=1211+nRegNew-27;
					else if(nRegNew<=42)nRegNew=1311+nRegNew-35;
					else if(nRegNew<=50)nRegNew=1411+nRegNew-43;
					else if(nRegNew<=58)nRegNew=1511+nRegNew-51;
					else if(nRegNew<=66)nRegNew=1611+nRegNew-59;
					util.pr(nRegNew);
				}
				
				
			}
			


			sp[ib+2]=Integer.toString(nRegNew);
			
			lines[0]=sp[ib]+",";
			for(int j=1;j<sp.length;j++)
				lines[0]=lines[0]+sp[j]+",";
	

			for(int i=0;i<7;i++)
				pw.println(lines[i]);		


		}

		pw.flush();

		br.close();
		fr.close();
		
		System.out.println();
		System.out.println(" Bun data was written to:");
		System.out.println("    "+reRegFile);
		System.out.println();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
	



}




public void renumber6CoilsByCoord(int ir, Vect axis, int nSegment){


	
	int[] NodeIndex=new int[2000*10000];
	Vect[] coords=new Vect[2000*100000];
	
	int[] hexVert=new int[8];
	Vect angs=new Vect().linspace(0, 2*PI, nSegment+1);
	angs.hshow();
	int nElems=0;

	int nNodes=0;
	int nnMax=0;
	String file=util.getFile();
	String folder=new File(file).getParentFile().getPath();
	String reRegFile=folder+"\\reRegion.neu";

	String line="";


	String[] sp;


	try{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reRegFile)));
		
		File f=new File(file);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		while(!line.startsWith("   403")){

			line=br.readLine();
			pw.println(line);
			if(line==null) break;

		}
		

		
		while(true){
			line=br.readLine();
			pw.println(line);
			sp=line.split(regex);
			
			

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

			
			sp=lines[1].split(regex);

			int ib=0;
			if(sp[0].equals("")) ib=1;
			for(int i=0;i<8;i++)
				hexVert[i]=Integer.parseInt(sp[i+ib]);

		
			
			
			sp=lines[0].split(regex);
			ib=0;
			if(sp[0].equals("")) ib=1;

	
			int nReg=Integer.parseInt(sp[ib+2]);
			int nRegNew=nReg;

			if(nReg==ir)	{
				
				Vect center=new Vect(3);
				int nVerts=0;
				for(int i=0;i<8;i++){
					int nn=hexVert[i];
					if(nn>0){
						center=center.add(coords[NodeIndex[nn]]);
						nVerts++;
					}
				}
	
			center=center.times(1./nVerts);
			
			Vect center_shifted2D=center.v2().sub(axis.v2());
			double ang=util.getAng(center_shifted2D);
			//util.pr(ang);
			int angBlock=0;
			for(int k=0;k<angs.length-1;k++){

				if(ang>=angs.el[k] && ang<angs.el[k+1]){
					angBlock=k;

					break;
					
				}
			}

				nRegNew=nRegNew+angBlock*100;
				util.pr(nRegNew);
	
			}


			sp[ib+2]=Integer.toString(nRegNew);
			
			lines[0]=sp[ib]+",";
			for(int j=1;j<sp.length;j++)
				lines[0]=lines[0]+sp[j]+",";
	

			for(int i=0;i<7;i++)
				pw.println(lines[i]);		


		}

		pw.flush();

		br.close();
		fr.close();
		
		System.out.println();
		System.out.println(" Bun data was written to:");
		System.out.println("    "+reRegFile);
		System.out.println();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
	



}




public void renumberRegionAcualModel(){


	
	int[] NodeIndex=new int[2000*10000];
	Vect[] coords=new Vect[2000*100000];
	
	int[] hexVert=new int[8];
	
	Vect corner1=new Vect( 147.2 ,0.,-77.8 );
	Vect x1s=new Vect(6);
	x1s.el[0]=corner1.el[0];
	for(int i=1;i<x1s.length;i++)
		x1s.el[i]=x1s.el[i-1]+1.12;
	
	x1s.hshow();
	
	Vect z1s=new Vect(6);
	
	z1s.el[0]=corner1.el[2];
	for(int i=1;i<z1s.length;i++)
		z1s.el[i]=z1s.el[i-1]+1.12;
	z1s.hshow();
	
/*	Vect corner2=new Vect( 144.4 ,0.,72.2 );
	Vect x2s=new Vect(11);
	x2s.el[0]=corner2.el[0];
	for(int i=1;i<x2s.length;i++)
		x2s.el[i]=x2s.el[i-1]+1.12;
	
	x2s.hshow();
	Vect z2s=new Vect(6);
	
	z2s.el[0]=corner2.el[2];
	for(int i=1;i<z2s.length;i++)
		z2s.el[i]=z2s.el[i-1]+1.12;*/
	int nElems=0;
	int nRegNew=0;
	int nBlock=0;
	int nNodes=0;
	int nnMax=0;
	String file=util.getFile();
	String folder=new File(file).getParentFile().getPath();
	String reRegFile=folder+"\\reRegion.neu";

	String line="";


	String[] sp;


	try{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reRegFile)));
		
		File f=new File(file);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		while(!line.startsWith("   403")){

			line=br.readLine();
			pw.println(line);
			if(line==null) break;

		}
		

		
		while(true){
			line=br.readLine();
			pw.println(line);
			sp=line.split(regex);
			
			

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

			
			sp=lines[1].split(regex);

			int ib=0;
			if(sp[0].equals("")) ib=1;
			for(int i=0;i<8;i++)
				hexVert[i]=Integer.parseInt(sp[i+ib]);

		
			
			
			sp=lines[0].split(regex);
			ib=0;
			if(sp[0].equals("")) ib=1;

	
			int nReg=Integer.parseInt(sp[ib+2]);
			if(nReg>1000000) nReg-=1000000;

			nRegNew=nReg;
			
			nBlock=0;
			
			if((nReg>=111 && nReg<=116) || (nReg>=211 && nReg<=216)|| (nReg>=311 && nReg<=316)|| (nReg>=411 && nReg<=416)|| (nReg>=511 && nReg<=516)|| (nReg>=611 && nReg<=616))	{
				
				Vect center=new Vect(3);
				int nVerts=0;
				for(int i=0;i<8;i++){
					int nn=hexVert[i];
					if(nn>0){
						center=center.add(coords[NodeIndex[nn]]);
						nVerts++;
					}
				}
	
			center=center.times(1./nVerts);
			
			int nb=0;
			boolean found=false;
			for(int ix=0;ix<x1s.length-1;ix++){
				for(int iz=0;iz<z1s.length-1;iz++){

					nb=ix*(x1s.length-1)+iz+1;
				if(center.el[0]>x1s.el[ix] && center.el[0]<x1s.el[ix+1] && center.el[2]>z1s.el[iz] && center.el[2]<z1s.el[iz+1]){
					nBlock=nb;
					found=true;
					//util.pr(nBlock);
					
					break;
				
				}
				}
				if(found) break;
				
			}
			


			if(found){
				nRegNew=nRegNew+nBlock*1000;
				//util.pr(nRegNew);
			}
	
			}
	


			sp[ib+2]=Integer.toString(nRegNew);
			
			lines[0]=sp[ib]+",";
			for(int j=1;j<sp.length;j++)
				lines[0]=lines[0]+sp[j]+",";
	
			//if(nBlock>=1000) util.pr(nRegNew+"       "+lines[0]);

			for(int i=0;i<7;i++)
				pw.println(lines[i]);		


		}

		pw.flush();

		br.close();
		fr.close();
		
		System.out.println();
		System.out.println(" Bun data was written to:");
		System.out.println("    "+reRegFile);
		System.out.println();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
	



}


public void twist3DNeu(Vect ax,double pitch, double range){

	double factor=2*PI/pitch;
	DecimalFormat formatter= new DecimalFormat("0.0000000");
	Mat R2D=new Mat(2,2);

	

	String file=util.getFile();
	String folder=new File(file).getParentFile().getPath();
	String twistedFile=folder+"\\twisted.neu";

	String line="";


	String[] sp;


	try{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(twistedFile)));
		
		File f=new File(file);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		while(!line.startsWith("   403")){

			line=br.readLine();
			pw.println(line);
			if(line==null) break;

		}
		
	
		while(true){
			line=br.readLine();
			sp=line.split(regex);
			
			if(!util.first(line).equals("-1")) {
				int ib=0;
				if(sp[0].equals("")) ib=1;
				//int nn=Integer.parseInt(sp[ib]);
				Vect v2d=new Vect(Double.parseDouble(sp[ib+11]),Double.parseDouble(sp[ib+12]));
				Vect v2d_shifted=v2d.sub(ax.v2());

				double r=v2d_shifted.norm();
				
				if(r<range){
				
		
	
				double z=Double.parseDouble(sp[ib+13]);

				double angle=z*factor;
				R2D=util.rotMat2D(angle);
				Vect v2d_shifted_rotated=R2D.mul(v2d_shifted);
				Vect v2d_rotated=v2d_shifted_rotated.add(ax.v2());
				sp[ib+11]=formatter.format(v2d_rotated.el[0]);
				sp[ib+12]=formatter.format(v2d_rotated.el[1]);
				
				line=sp[ib]+",";
				for(int j=1;j<sp.length;j++)
					line=line+sp[j]+",";
			
				}
				
				pw.println(line);

				
			} 
			else{
				pw.println(line);

				break;
			}
				
		}
		
		
		while(true){
			line=br.readLine();
			pw.println(line);	
			if(line==null) break;

		}

		pw.flush();

		pw.close();
		br.close();
		fr.close();
		
		System.out.println();
		System.out.println(" Bun data was written to:");
		System.out.println("    "+twistedFile);
		System.out.println();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
}
	

public void wind(Vect ax,double pitch){

	double factor=2*PI/pitch;
	DecimalFormat formatter= new DecimalFormat("0.0000000");
	Mat R3D=new Mat(3,3);

	double theta=0*PI/4;

	String file=util.getFile();
	String folder=new File(file).getParentFile().getPath();
	String twistedFile=folder+"\\twisted.neu";

	String line="";


	String[] sp;


	try{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(twistedFile)));
		
		File f=new File(file);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		while(!line.startsWith("   403")){

			line=br.readLine();
			pw.println(line);
			if(line==null) break;

		}
		
	
		while(true){
			line=br.readLine();
			sp=line.split(regex);
			
			if(!util.first(line).equals("-1")) {
				int ib=0;
				if(sp[0].equals("")) ib=1;
				//int nn=Integer.parseInt(sp[ib]);
				Vect v3d=new Vect(Double.parseDouble(sp[ib+11]),Double.parseDouble(sp[ib+12]),Double.parseDouble(sp[ib+13]));
				
				Vect v2d=new Vect(Double.parseDouble(sp[ib+11]),Double.parseDouble(sp[ib+12]));
				Vect v2d_shifted=v2d.sub(ax.v2());

				double r=v2d_shifted.norm();
				
						
		
	
				double z=Double.parseDouble(sp[ib+13]);
				

				double angle=z*factor;
				
				Vect newAxis=new Vect(Math.sin(theta)*Math.cos(angle),Math.sin(theta)*Math.sin(angle),Math.cos(theta));
				R3D=util.rotMat(newAxis,new Vect(0,0,1));
				//Vect v2d_shifted_rotated=R2D.mul(v2d_shifted);
				Vect v3d_rotated=R3D.mul(v3d);;//=v2d_shifted_rotated.add(ax.v2());
				sp[ib+11]=formatter.format(v3d_rotated.el[0]);
				sp[ib+12]=formatter.format(v3d_rotated.el[1]);
				sp[ib+13]=formatter.format(v3d_rotated.el[2]);
				
				line=sp[ib]+",";
				for(int j=1;j<sp.length;j++)
					line=line+sp[j]+",";
			
				
				
				pw.println(line);

				
			} 
			else{
				pw.println(line);

				break;
			}
				
		}
		
		
		while(true){
			line=br.readLine();
			pw.println(line);	
			if(line==null) break;

		}

		pw.flush();

		pw.close();
		br.close();
		fr.close();
		
		System.out.println();
		System.out.println(" Bun data was written to:");
		System.out.println("    "+twistedFile);
		System.out.println();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
}
	


public void multipleTwist3DNeu(double pitch, double range1){
	double scaleFactor=1e0;
	//double pitch=15;


	double factor=2*PI/pitch;
	DecimalFormat formatter= new DecimalFormat("0.0000000");
	Mat R2D=new Mat(2,2);

	//double range1=.158*scaleFactor;
	//double range1=.1985*scaleFactor;
	//double range1=.9*.158*scaleFactor;
	
	// range1=.1598*scaleFactor;
	
	//Vect[] axs=new Vect[6];
	int nSegments=6;
	Vect[] axs=new Vect[nSegments];
	Vect ranges1=new Vect(axs.length);
	//axs[0]=new Vect(0.294627 , 0.170103 ).times(scaleFactor);
	//axs[0]=new Vect(  0.282009  , 0.204925  ).times(scaleFactor);
	axs[0]=new Vect(  0.484403  , 0.279682  ).times(scaleFactor);
//	axs[0]=new Vect(  1.424699  ,  0.82261 ).times(scaleFactor);
	 
	ranges1.el[0]=range1;
	for(int k=1;k<axs.length;k++){
		R2D=util.rotMat2D(2*PI/nSegments);
		axs[k]=R2D.mul(axs[k-1].v2()).v3();
		ranges1.el[k]=range1;
		
	}
	
	for(int k=0;k<axs.length;k++){
		axs[k].hshow();
		
	}
	String file=util.getFile();
	String folder=new File(file).getParentFile().getPath();
	String twistedFile=folder+"\\twisted.neu";

	String line="";


	String[] sp;


	try{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(twistedFile)));
		
		File f=new File(file);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		while(!line.startsWith("   403")){

			line=br.readLine();
			pw.println(line);
			if(line==null) break;

		}
		
	
		while(true){
			line=br.readLine();
			sp=line.split(regex);
			
			if(!util.first(line).equals("-1")) {
				int ib=0;
				if(sp[0].equals("")) ib=1;
				//int nn=Integer.parseInt(sp[ib]);
				Vect v2d=new Vect(Double.parseDouble(sp[ib+11]),Double.parseDouble(sp[ib+12]));
				
				int kx=-1;
				double rk=0;
				//double rout=v2d.sub(axOut.v2()).norm();
				
				//if(rout<range2){
					
				for(int k=0;k<axs.length;k++){
					rk=v2d.sub(axs[k].v2()).norm();
					//util.pr(k+"  "+rk+"  "+ranges1.el[k]);

					if(rk<ranges1.el[k]){
						kx=k;
						break;
					}
			
				}
				
				Vect axis=null;
				
				if(kx>=0){ axis=axs[kx];
				//else axis=axOut;
				//util.pr(kx);

				Vect v2d_shifted=v2d.sub(axis.v2());
	
				double z=Double.parseDouble(sp[ib+13]);

				double angle=z*factor;
			//	util.pr(angle);
				//if(r>range1){
					//angle*=(range2-r)/(range2-range1);
			//	}
				
				R2D=util.rotMat2D(angle);
				Vect v2d_shifted_rotated=R2D.mul(v2d_shifted);
				Vect v2d_rotated=v2d_shifted_rotated.add(axis.v2());
				sp[ib+11]=formatter.format(v2d_rotated.el[0]);
				sp[ib+12]=formatter.format(v2d_rotated.el[1]);
				
				line=sp[ib]+",";
				for(int j=1;j<sp.length;j++)
					line=line+sp[j]+",";
			
				}
				//}
				pw.println(line);

				
			} 
			else{
				pw.println(line);

				break;
			}
				
		}
		
		
		while(true){
			line=br.readLine();
			pw.println(line);	
			if(line==null) break;

		}

		pw.flush();

		pw.close();
		br.close();
		fr.close();
		
		System.out.println();
		System.out.println(" Bun data was written to:");
		System.out.println("    "+twistedFile);
		System.out.println();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
}




public void copyTranslateNeu(int ns, Vect transV){

	DecimalFormat df= new DecimalFormat("0.0000000000");
	
	
	PrintWriter pw=null;
	
	String regex="[ ,\\t]+";
	String bun=util.getFile();
	
	String[] lines=new String[7];

	String line="";
	String[] sp=new String[15];
	

	int nElems=0;
	int nNodes=0;
	int maxNodeNumb=0;
	

	
	int[] nodeMap=new int[1000*10000];
	
//	Vect[] coords=new Vect[nodeMap.length];
	
	try{
		File f=new File(bun);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);

		String folder=new File(bun).getParentFile().getPath();
		String fout=folder+"\\copy_translated.neu";
		pw = new PrintWriter(new BufferedWriter(new FileWriter(fout)));	
		
		line=br.readLine();
		pw.println(line);


		while(!line.startsWith("   403")){

			line=br.readLine();
				pw.println(line);
			if(line==null) break;

		}
		
		while(true){
			line=br.readLine();
			if(line==null) break;
			sp=line.split(regex);
			
			if(!util.first(line).equals("-1")) {
				int ib=0;
				if(sp[0].equals("")) ib=1;
				int nn=Integer.parseInt(sp[ib]);
				
				nNodes++;
				nodeMap[nn]=nNodes;
				Vect coord=new Vect(Double.parseDouble(sp[ib+11]),Double.parseDouble(sp[ib+12]),Double.parseDouble(sp[ib+13]));

				sp[ib]=Integer.toString(nNodes);

				pw.println(line);
	
				int nnNew=ns+nn;
				Vect nodeCopy=coord.add(transV);
				String linex=nnNew+",0,0,1,46,0,0,0,0,0,0,"+df.format(nodeCopy.el[0])+","+df.format(nodeCopy.el[1])+","+df.format(nodeCopy.el[2])+",0,";
				pw.println(linex);
			} 
			else {
				
			//	for(int i=1;i<=nNodes;i++){
					//Vect nodeCopy=coords[i].add(transV);
				//	int nnNew=ns+i;
					//String linex=nnNew+",0,0,1,46,0,0,0,0,0,0,"+df.format(nodeCopy.el[0])+","+df.format(nodeCopy.el[1])+","+df.format(nodeCopy.el[2])+",0,";
				//	if(writeFile)
					//	pw.println(linex);
			//	}
				break;
			}
		}
		
			pw.println(line);
		

		while(!line.startsWith("   404")){

			line=br.readLine();
			
				pw.println(line);
			
			if(line==null) break;

		}


				
		int[] nVerts=new int[8];
		
		while(true){

			for(int i=0;i<7;i++){
				lines[i]=br.readLine();

			}

			if(util.first(lines[0]).equals("-1")) {

				break;
			}
			
			
			for(int i=0;i<7;i++)
				pw.println(lines[i]);	
			
			
			nElems++;

			
				for(int i=0;i<7;i++)
					pw.println(lines[i]);	
				
					sp=lines[0].split(regex);
					int ib=0;
					if(sp[0].equals("")) ib=1;
					int ne=Integer.parseInt(sp[ib]);
					int neNew=ne+ns;
					sp[0]=Integer.toString(neNew);
					lines[0]=sp[0]+",";
					for(int j=1;j<sp.length;j++)
						lines[0]=lines[0]+sp[j]+",";
					
					sp=lines[1].split(regex);

					ib=0;
					if(sp[0].equals("")) ib=1;
					for(int i=0;i<8;i++)
						nVerts[i]=Integer.parseInt(sp[i+ib]);
					
					for(int i=0;i<8;i++){
						int nv=	nVerts[i]+ns;
						sp[i]=Integer.toString(nv);
					}
					
					lines[1]=sp[0]+",";
					for(int j=1;j<sp.length;j++)
						lines[1]=lines[1]+sp[j]+",";

					for(int i=0;i<7;i++)
						pw.println(lines[i]);	
					
				}
					

		while(true){

			line=br.readLine();
			if(line==null) break;

				pw.println(line);
			

		}


		pw.flush();


		System.out.println();
		System.out.println(" Bun data was written to:");
		System.out.println("    "+fout);
		System.out.println();


		br.close();
		fr.close();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	


}


public void cutCorners(){
	
	
	Vect[] coords=new Vect[10000000];
	PrintWriter pw=null;
	
	String regex="[ ,\\t]+";
	String bun=util.getFile();
	
	String[] lines=new String[7];

	String line="";
	String[] sp=new String[15];
	
	int nElems=0;

		
	try{
		File f=new File(bun);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);

		
		

		

		String folder=new File(bun).getParentFile().getPath();
		String fout=folder+"\\cornersCut.neu";

		pw = new PrintWriter(new BufferedWriter(new FileWriter(fout)));	
		
		line=br.readLine();
		pw.println(line);


		while(!line.startsWith("   403")){

			line=br.readLine();
				pw.println(line);
			if(line==null) break;

		}
		

		
		while(true){
			line=br.readLine();
			if(line==null) break;
			sp=line.split(regex);

			

				if(!util.first(line).equals("-1")) {
					int ib=0;
					if(sp[0].equals("")) ib=1;
					int nn=Integer.parseInt(sp[ib]);
					coords[nn]=new Vect(Double.parseDouble(sp[ib+11]),Double.parseDouble(sp[ib+12]),Double.parseDouble(sp[ib+13]));

				
			} 
			else break;
			}
		
			pw.println(line);
		

		while(!line.startsWith("   404")){

			line=br.readLine();
			
				pw.println(line);
			
			if(line==null) break;

		}


		int[] nVerts=new int[8];
		
		while(true){

			for(int i=0;i<7;i++){
				lines[i]=br.readLine();

			}

			if(util.first(lines[0]).equals("-1")) {
				break;
			}
			
			sp=lines[0].split(regex);
			int ib=0;
			if(sp[0].equals("")) ib=1;
			

				
			
			sp=lines[1].split(regex);

			ib=0;
			if(sp[0].equals("")) ib=1;
			for(int i=0;i<8;i++){
				nVerts[i]=Integer.parseInt(sp[i+ib]);
			}
			
			
			Vect center=new Vect(3);
			int nV=0;
			for(int i=0;i<8;i++){
				int nn=nVerts[i];
				if(nn>0){
					center=center.add(coords[nn]);
					nV++;
				}
			}

		center=center.times(1./nV);
			if(center.el[2]>2.5){
			
			nElems++;
			
				for(int i=0;i<7;i++)
					pw.println(lines[i]);		
			}


		}
		
		while(true){

			line=br.readLine();
			if(line==null) break;
	
			pw.println(line);
			

		}


		pw.flush();

		System.out.println(" Number of elements:"+ nElems);
		System.out.println();
		System.out.println(" Bun data was written to:");
		System.out.println("    "+fout);
		System.out.println();


		br.close();
		fr.close();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
	
}

public void cutCornersbb(){

	

	//int[][] hexVert=new int[1000000][8];
	
	int[] NodeIndex=new int[10000000];
	Vect[] coords=new Vect[10000000];

	int nElems=0;

	int nNodes=0;
	int nnMax=0;
	String file=util.getFile();
	String folder=new File(file).getParentFile().getPath();
	String tetMidFile=folder+"\\cornerCut.neu";

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



		br.close();
		fr.close();

	}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }


	int[] hexVert=new int[8];
	
	int ix=0;
	try{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(tetMidFile)));
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

			
			pw.println(line);
			if(line==null) break;

		}

		String[] lines=new String[7];
		
		ix=-1;
	
		while(true){

			for(int i=0;i<7;i++){
				lines[i]=br.readLine();

			}

			if(util.first(lines[0]).equals("-1")) {
				break;
			}
			//if(!sp1[3].equals("25")) continue;
			//	sp1[1]=Integer.toString(elNum);
			//sp1[3]="26";
			//sp1[4]="10";
			//lines[0]=sp1[ib];
		//	for(int i=1+ib;i<sp1.length;i++)
			//	lines[0]=lines[0]+","+sp1[ib+i];



			
			sp=lines[1].split(regex);

			int ib=0;
			if(sp[0].equals("")) ib=1;
			for(int i=0;i<8;i++)
				hexVert[i]=Integer.parseInt(sp[i+ib]);
			

			int n1=hexVert[0];
			int n2=hexVert[1];
			int n3=hexVert[2];
			int n4=hexVert[4];

			Vect mid1=coords[NodeIndex[n1]].add(coords[NodeIndex[n2]]).times(.5);
			Vect mid2=coords[NodeIndex[n3]].add(coords[NodeIndex[n4]]).times(.5);

			Vect mid=mid1.add(mid2).times(.5);
		//if(mid.el[0]<-.3 || mid.el[0]>.539 || mid.el[1]<-.3 ||  mid.el[1]>.539 ) continue;
		if(mid.el[2]<15. ) continue;

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
