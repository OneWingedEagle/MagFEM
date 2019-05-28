package fem;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;

import femSolver.FEMsolver;
import femSolver.StaticElectricSolver;
import femSolver.StaticNonlinearMagSolver;
import io.Loader;
import io.Writer;
import main.Main;
import materialData.BHCurve;
import materialData.CurrentWaveForm;
import materialData.LamBCurve;
import math.*;
import meshFactory.MeshManipulator;
import static java.lang.Math.*;


public class Model{

	public MagMatAssembler magMat;
	public int dim=3,iterMax=10000,nonLinIterMax=30,cpb=1,nRotReg;
	public double cpm=PI/2,Rg=1e5;
	public int numberOfRegions, numberOfNodes,numberOfElements,numberOfEdges,nXYedges,
	nBlocks,nBH,nLam;
	public int nodeNodeConnectionMax=27,nElVert=8,nBoundary=0;
	public int nElEdge=12;
	public double[] spaceBoundary;
	public int[] BCtype;
	public String[] blockMaterial;
	public Vect unifB;
	public int unifBTimeId;
	public Region[] region;
	public PhiCoil[] phiCoils;
	public Node[] node;
	public Element[] element;
	public Edge[] edge;
	public Loader loader=new Loader();
	public Writer writer=new Writer();
	public BoundarySet bcSet=new BoundarySet();
	public Calculator femCalc;
	public Force forceCalc;
	public int[] edgeUnknownIndex,unknownEdgeNumber,nodeUnknownIndex,unknownNodeNumber,
	U_unknownIndex,unknownUnumber,A_unknownIndex,unknownAnumber,T_unknownIndex,unknownTnumber;
	public int[] knownEdgeNumber,nodeVarIndex,varNodeNumber,unCurRegNumb;
	public double[] knownEdgeValue;
	public double scaleFactor=1,Bmax1=0,Bmin1=0,Bmax,Bmin,stressMax=0,nodalStressMax=0,
			stressMin,nodalStressMin=0,Jmin,Jmax,Jemin,Jemax,maxDim,minEdgeLength,maxEdgeLength,
			FmsMin,FmsMax=0,FreluctMax=0,FreluctMin=0,uMax=0,AuMax,defScale;
	public int numberOfUnknownEdges,numberOfMappedEdges,numberOfKnownEdges,numberOfVarNodes
	,numberOfKnownPhis,numberOfUnknowns,analysisMode,stressViewCode=-1
	,numberOfUnknownU,numberOfUnknownUcomp,numberOfCurrents,numberOfUnknownCurrents,nNeutral,defMode=1,numberOfUnknownT;
	public boolean deform,thermal,hasJ,hasM,forceLoaded,fluxLoaded,potentialLoaded,
	stressLoaded,forceCalcLoaded,coupled,cpvms,calCurve,nonLin,hasPBC,hasMechPBC=true,rotIndex;
	public int nEdEd=34,nNodNod=27,nEdNod=18,nNodEd=54,nEdgeHasJ;
	static int[][] facetVert={{6,2,5},{7,4,3},{6,7,2},{0,4,1},{4,7,5},{0,1,3}};

	public byte elCode=4;
	public double nu0=1e7/(4*PI);
	public double freq=1,dt,errCGmax=1,
			errNRmax=1e-6,errCG_NR=1e-1,errFluxMax=1e-3;
	public int nTsteps,nBegin,nEnd,nInc,nRotorElements,currentTimeStep;
	public int coordCode=0,timeIntegMode=0,eddyTimeIntegMode;
	public LamBCurve[] lamB;
	public BHCurve[] BH;
	public CurrentWaveForm ia,ib,ic,va,vb,vc;
	public String elType="hexahedron";
	public SpMat Hs,Ms,Ks,Ls,Cs,Ss,Ps,Qs,Fs;
	public SpMat Rs,Bs,BtBs;
	public Mat eigVects,bigMat;
	public Vect lams,RHS,bU,bT,HpAp,HkAk,RHS_boundary;
	public boolean AC,writeFiles,Tmethod,
	circuit,stranded,wavePC,loadFlux,loadPotentioal,loadPrevMag,loadForce,saveFlux,saveJe,saveForce,
	magAnalysis=true,axiSym;
	public int forceCalcMode=1,dataType;
	public double alpha1,alpha2,r1,r2,h1,h2,rm,TrqZ,height=1,gw=1e4,vNeutral,tet,tetp,tetpp,spwmLevel,rayAlpha, rayBeta;
	public int[] PBCpair;
	public int[][][] commonNodes;
	public SpMatSolver solver=new SpMatSolver();
	public StaticNonlinearMagSolver nonlinearSolver=null;
	FEMsolver femSolver=new FEMsolver();;

	public Vect up,upp,xp,Ci,ud,udd,hp;
	public String eddyFolderxx,filePath,meshFilePath,dataFilePath,eddyFilePathxx,fluxFilePathxx,resultFolder,fluxFolderIn;
	public double[][] H2,H3;
	public double[] C,Cj2d;
	public Vect[] Cj;

	public String[] forceFile;
	public String forceFolder;
	public SpVect[] lastRows,lastRowsAll;
	public Model m2d;
	public Vect[][] forceLamin;
	public int[] mapnr;
	public boolean hasBunif,open_vps;
	public int nCLNstages;
	public Network network;
	public TimeFunction[] timeFunctions;
	public StaticElectricSolver phiSolver=null;
	public Main main;

	public Model(){}


	public Model(int nRegions, int nElements,int nNodes, String elType){

		this.alloc(nRegions, nElements, nNodes, elType);
		/*this.numberOfRegions=nRegions;
		this.numberOfElements=nElements;

		this.numberOfNodes=nNodes;

		this.setElType(elType);

		region=new Region[this.numberOfRegions+1];
		for(int i=1;i<=this.numberOfRegions;i++)
			region[i]=new Region(dim);

		element=new Element[this.numberOfElements+1];
		for(int i=1;i<=this.numberOfElements;i++)
			element[i]=new Element(elType);

		node=new Node[this.numberOfNodes+1];
		for(int i=1;i<=this.numberOfNodes;i++)
			node[i]=new Node(dim,i);*/
	}


	public Model(String bun){
		new Model();
		loadMesh(bun);
	}


	public void setFemCalc(){
		femCalc=new Calculator(this);

	}
	public void setForceCalc(){
		forceCalc=new Force(this);

	}

	public void alloc(int nRegions, int nElements,int nNodes, String elType){

		this.numberOfRegions=nRegions;
		this.numberOfElements=nElements;

		this.numberOfNodes=nNodes;

		this.setElType(elType);


		region=new Region[this.numberOfRegions+1];
		for(int i=1;i<=this.numberOfRegions;i++)
			region[i]=new Region(dim);

		element=new Element[this.numberOfElements+1];
		for(int i=1;i<=this.numberOfElements;i++)
			element[i]=new Element(elType);

		node=new Node[this.numberOfNodes+1];
		for(int i=1;i<=this.numberOfNodes;i++)
			node[i]=new Node(i,dim);
	}
	public void loadMesh(String bunFilePath){

		loader.loadMesh(this, bunFilePath);

	}

	public void loadData(String dataFilePath)
	{
		//if(dim==3)
			loader.loadData(this,dataFilePath);
	//	else
		//	loader.loadData2D(this,dataFilePath);
		this.femCalc=new Calculator(this);
		this.forceCalc=new Force(this);	
		setMagMatAssembler();

	}

	public Model deepCopy(){

		Model copy=new Model(this.numberOfRegions,this.numberOfElements,this.numberOfNodes,this.elType);

		for(int i=1;i<=this.numberOfRegions;i++){
			copy.region[i].setFirstEl(this.region[i].getFirstEl());
			copy.region[i].setLastEl(this.region[i].getLastEl());
		}


		for(int i=1;i<=this.numberOfElements;i++){
			copy.element[i].setVertNumb(this.element[i].getVertNumb());
			copy.element[i].setEdgeNumb(this.element[i].getEdgeNumb());
			copy.element[i].setRegion(this.element[i].getRegion());
		}


		for(int i=1;i<=this.numberOfNodes;i++){
			copy.node[i].setCoord(this.node[i].getCoord());
		}

		if(this.edge!=null){
			copy.edge=new Edge[this.numberOfEdges+1];
			for(int i=1;i<=this.numberOfEdges;i++){
				//copy.edge[i]=new Edge(this.edge[i].endNodeNumber[0],this.edge[i].endNodeNumber[1]);
				copy.edge[i]=new Edge(this.edge[i].node[0],this.edge[i].node[1]);
				copy.edge[i].A=this.edge[i].A;
				copy.edge[i].Ap=this.edge[i].Ap;
			}
		}

		copy.elCode=this.elCode;
		copy.scaleFactor=this.scaleFactor;

		return copy;
	}

	public Model fill(int nRegs,int nEls,int nNodes, String elType){

		Model copy=new Model(nRegs,nEls,nNodes,elType);
		for(int i=1;i<=this.numberOfRegions;i++){
			copy.region[i].setFirstEl(this.region[i].getFirstEl());
			copy.region[i].setLastEl(this.region[i].getLastEl());
			copy.region[i].setName(this.region[i].getName());
		}

		for(int i=this.numberOfRegions+1;i<=nRegs;i++){
			copy.region[i].setFirstEl(1);
			copy.region[i].setLastEl(0);
		}

		for(int i=1;i<=this.numberOfElements;i++){
			copy.element[i].setVertNumb(this.element[i].getVertNumb());
			copy.element[i].setRegion(this.element[i].getRegion());
		}


		for(int i=1;i<=this.numberOfNodes;i++){
			copy.node[i].setCoord(this.node[i].getCoord());
		}


		copy.scaleFactor=this.scaleFactor;
		return copy;
	}


	public double edgeLength(int i){
	//	double length=node[edge[i].endNodeNumber[1]].getCoord().sub(node[edge[i].endNodeNumber[0]].getCoord()).norm();
		double length=edge[i].node[0].getCoord().sub(edge[i].node[1].getCoord()).norm();

		return length;
	}

	public void setMinEdge(){
		double minEdge=1e40;
		for(int i=1;i<=this.numberOfEdges;i++){
			if(edge[i].length<minEdge) minEdge=edge[i].length;
		}
		this.minEdgeLength=minEdge;
	}
	public void setMaxEdge(){
		double maxEdge=0;
		for(int i=1;i<=this.numberOfEdges;i++)
			if(edge[i].length>maxEdge) maxEdge=edge[i].length;

		this.maxEdgeLength=maxEdge;
	}

	public  void setMaxDim(){

		maxDim=getRmax();
	}

	public void setHasJ(){

		for(int ir=1;ir<=numberOfRegions;ir++){
			if(region[ir].hasJ){
			hasJ=true;
			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

					element[i].setHasJ(true);
			}
			}
			else for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setHasJ(false);
		}
			}
	}
	public void setHasM(){
		for(int i=1;i<=numberOfRegions;i++)
			if(region[i].hasM){
				hasM=true;
				break;
			}
	}

	public void setNonLin(boolean nonlin){

		nonLin=nonlin;
		
		if(nonlin) nonlinearSolver=new StaticNonlinearMagSolver();

		for(int ir=1;ir<=numberOfRegions;ir++){
			if(region[ir].isNonLinear)
				for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
					element[i].setNonlin(true);
		}	

	}
	public void setNonLinToElememts(){

		if(nonLin)
			for(int ir=1;ir<=numberOfRegions;ir++){
				if(region[ir].BHnumber>0){
					region[ir].isNonLinear=true;
					for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
						element[i].setNonlin(true);
				}
			}	

	}


	public void setHasThermal(){		

		for(int ir=1;ir<=numberOfRegions;ir++){
			if(region[ir].thermal)
				for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

					element[i].setHasThermal(true);

				}

		}	
	}



	public void setMagMatAssembler(){
		magMat=new MagMatAssembler(this);
	}
	public void setMagMat(){
		magMat.setMagMat(this);
	}

	public void femSolver(){
		this.femSolver=new FEMsolver(this);

	}


	public void writeMesh(String bunFilePath){
		writer.writeMesh(this,bunFilePath);
	}

	public void writeMesh(String bunFilePath, boolean deform){
		writer.writeMesh(this,bunFilePath,deform);
	}

	public void writeMeshOneNode(String bunFilePath){
		writer.writeMeshOneNode(this,bunFilePath);
	}

	public void writeModeShape(String bunFilePath , double a){
		writer.writeModeShape(this,bunFilePath,a);

	}

	public void writeMeshTriToQuad2(String bunFilePath , boolean deformed){
		writer.writeMeshTriToQuad2(this,bunFilePath ,deformed);
	}

	public void writeMeshq2h(String bunFilePath , boolean deformed){
		writer.writeMeshq2h(this,bunFilePath ,deformed);
	}

	public void hexaToPyramid(String bunFilePath){
		writer.writeMeshHexaToPyramid(this,bunFilePath);
	}

	public void pyramidToTetra(String bunFilePath){
		writer.writeMeshPyramidToTetra(this,bunFilePath);
	}

	public void writeMeshTriToTri(String bunFilePath,double r1,double r2){
		//	writer.writeMesh323(this,bunFilePath,r1,r2);
		writer.writeMeshTriToTriW(this,bunFilePath);

	}

	public void writeMeshTriToQuad(String bunFilePath){
		writer.writeMeshTriToQuad(this,bunFilePath);
	}

	public void writeMeshQuadToTri(String bunFilePath){
		writer.writeMeshQuadToTri(this,bunFilePath);

	}

	public void writeMeshHexaToPrism(String bunFilePath){
		writer.writeMeshHexaToPrism(this,bunFilePath );
	}
	public void writeMeshHexaToPrism(String bunFilePath,int dir){
		writer.writeMeshHexaToPrism(this,bunFilePath,dir );
	}

	public void prismToHexa(String bunFilePath){
		writer.writeMeshPrismToHexa(this,bunFilePath );
	}

	public void writeData(String dataFilePathOut) {
		try{
			this.writer.copyFile(this.dataFilePath,dataFilePathOut);
			System.out.println(" Data file was written to "+dataFilePathOut);

		}
		catch(IOException e){System.out.println("writing data file failed.");}

	}



	public Node[] elementNodes(int i){
		Node[] elementNode=new Node[nElVert];
		int[] vertNumb=element[i].getVertNumb();
		for(int j=0;j<nElVert;j++){

			elementNode[j]=node[vertNumb[j]];
		}

		return elementNode;
	}
	public Edge[] elementEdges(int i){
		Edge[] elementEdge=new Edge[nElEdge];
		int[] edgeNumb=element[i].getEdgeNumb();
		for(int j=0;j<nElEdge;j++)
			elementEdge[j]=edge[edgeNumb[j]];
		return elementEdge;
	}

	public int[] getContainingElement(Vect P){
		if(elCode==0) return getContaining3angElement(P);
		else if(elCode==1) return getContainingQuadElement(P);
		int[] elResult;
		int[] result=new int[7];

		for(int ir=1;ir<=this.numberOfRegions;ir++){

			//	if(ir!=1 && ir!=2) continue;
			//if(ir!=1 ) continue;


			for(int i=this.region[ir].getFirstEl();i<=this.region[ir].getLastEl();i++){

				elResult=elContains(i,P);
				result[0]=0;
				for(int j=0;j<6;j++)
					if(elResult[j]<0){
						result[0]=-1;
						break;
					}

				if(result[0]==-1) continue;		


				result[0]=i;
				for(int j=0;j<6;j++)
					result[j+1]=elResult[j];
				break;

			}
		}



		return result;
	}

	public int[] getContaining3angElement(Vect P){
		int[] ne=new int[1];


		double[] S;
		double S0,St;
		for(int i=1;i<=numberOfElements;i++){
			S0=el3angArea(i);
			S=subS3ang(i,P);
			St=S[0]+S[1]+S[2];

			if(abs(1-St/S0)<1e-6){ ne[0]= i; break;}

		}

		return ne;
	}



	public int[] elContains(int i,Vect P){
		int n;
		int[] result=new int [nBoundary];
		for(int ns=0;ns<nBoundary;ns++){
			n=pointOnFace(i,ns,P);
			result[ns]=n;

		}

		return result;
	}

	public int pointOnFace(int i, int ns,Vect P){
		Vect v0,v1,v2,vP;

		v0=node[element[i].getVertNumb(facetVert[ns][0])].getCoord();
		v1=node[element[i].getVertNumb(facetVert[ns][1])].getCoord();
		v2=node[element[i].getVertNumb(facetVert[ns][2])].getCoord();
		vP=P.sub(v0);

		if(vP.norm()==0)
			return 0;
		double crossdot=v1.sub(v0).cross(v2.sub(v0)).dot(P.sub(v0));
		if(crossdot==0)
			return 0;
		else if(crossdot>1e-10)
			return -1;
		else
			return 1;

	}

	public int[] getContainingQuadElement(Vect P){

		int[] ne=new int[1];

		Vect[] v=new Vect[4];

		for(int i=1;i<=numberOfElements;i++){
			int[] vertNumb=element[i].getVertNumb();
			for(int j=0;j<4;j++)
				v[j]=node[vertNumb[j]].getCoord();

			double S0=v[1].sub(v[0]).cross(v[3].sub(v[0])).norm()+
					v[1].sub(v[2]).cross(v[3].sub(v[2])).norm();
			double S=0;
			for(int j=0;j<4;j++)
			{
				S=S+v[j].sub(P).cross(v[(j+1)%4].sub(P)).norm();
			}
			if(abs(1-S/S0)<1e-6) {ne[0]=i; break;}

		}

		return ne;

	}

	public double getRegionArea(int ir)

	{

		if(this.dim>2) throw new NullPointerException(" Region is not 2D. ");

		double S=0;


		for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
			S+=getElementArea( i);

		return S;

	}

	public double getRegionAreaXY(int ir)

	{

		if(this.dim!=3) throw new NullPointerException(" Region is not 3D. ");

		double S=0;


		S=3.08E-4;


		/*for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
				S+=getElementArea( i);
		 */
		return S;

	}

	public double getRegionVolume(int ir)

	{

		if(this.dim!=3) throw new NullPointerException(" Region is not 3D. ");

		double V=0;


		for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++)
			V+=getElementVolume( i);

		return V;

	}

	public double getElementArea(int i){
		if(elCode==0) return el3angArea(i);
		else 	if(elCode==1) return elQuadArea(i);
		else throw new NullPointerException(" Element is not 2D. ");

	}

	public double el3angArea(int i){
		Node[] vertexNode=elementNodes(i);
		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=v1.cross(v2).norm()/2;
		return S;
	}

	public double el3angMaxCosAng(int i){
		Node[] vertexNode=elementNodes(i);
		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double dd1= abs(v1.normalized().dot(v2.normalized()));
		v1=v1.times(-1);
		v2=vertexNode[2].getCoord().sub(vertexNode[1].getCoord());
		double dd2= abs(v1.normalized().dot(v2.normalized()));
		v2=v2.times(-1);
		v1=vertexNode[0].getCoord().sub(vertexNode[2].getCoord());
		double dd3= abs(v1.normalized().dot(v2.normalized()));
		double dm=max(dd3,max(dd1,dd2));

		return dm;
	}

	public double elQuadArea(int i){
		Node[] vertexNode=elementNodes(i);

		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[3].getCoord().sub(vertexNode[0].getCoord());
		Vect v3=vertexNode[1].getCoord().sub(vertexNode[2].getCoord());
		Vect v4=vertexNode[3].getCoord().sub(vertexNode[2].getCoord());
		double S=(v1.cross(v2).norm()+v4.cross(v3).norm())/2;
		return S;
	}

	public void setSolvedAL(Vect x){

		for(int i=1;i<=numberOfUnknownEdges;i++){


			edge[unknownEdgeNumber[i]].setSolvedAL(x.el[i-1]);	

		}


		if(this.hasPBC ){

			for(int i=1;i<=numberOfEdges;i++){
				if(edge[i].map>0 && !edge[i].edgeKnown)
				{
					if(edge[i].aPBC){

						edge[i].setSolvedAL(this.cpb*edge[edge[i].map].A);	

					}


					else{

						edge[i].setSolvedAL(edge[edge[i].map].A);	


					}
				}

			}	
		}




	}




	public void deformMesh(){

		for(int i=1;i<=this.numberOfNodes;i++){
			if(node[i].isDeformable() && node[i].u.norm()!=0)
				node[i].setCoord(node[i].getCoord().add(node[i].getU()));

		}

	}


	private void setElementB(int i){

		if(elCode==0&& !this.axiSym)set3angElementB(i);
		else if(elCode==0) set3angElementAxiB(i);
		else if(elCode==1&& !this.axiSym)setQuadElementB(i);
		else if(elCode==1) setQuadElementAxiB(i);
		//else if(elCode==2) setTetraElementB(i);
	//	else if(elCode==3) setPrismElementB(i);
		else{
			boolean[] edgeDir=element[i].getEdgeReverse();

			Node[] vertexNode=elementNodes(i);
			
			Vect zero=new Vect(3);

			Mat jac=femCalc.jacobian(vertexNode,zero);
			Vect B;
			Vect[] rotNe=femCalc.rotNe(jac,zero,edgeDir);
	
			B=getElementB(i,rotNe);

			element[i].setB(B);


		}
	}	



	private void setQuadElementB(int i){

		Node[] vertexNode=elementNodes(i);
		Vect zero=new Vect(2);
		Mat jac=new Mat();
		Vect[] rotNe;

		jac=femCalc.jacobian(vertexNode,zero);
		rotNe=femCalc.rotNeQuad(jac,zero);

		Vect B=getElementB(i,rotNe);

		element[i].setB(B);

	}


	private void setQuadElementAxiB(int i){

		Node[] vertexNode=elementNodes(i);
		Vect zero=new Vect(2);
		Mat jac=new Mat();
		Vect[] rotNe;

		jac=femCalc.jacobian(vertexNode,zero);
		rotNe=femCalc.rotNeQuad(jac,zero);


		double rr=0;
		for(int j=0;j<4;j++){
			double r1=vertexNode[j].getCoord(0)+1e-5;
			rr+=.25/r1;
		}

		//	rr=1;


		Vect B=getElementB(i,rotNe).times(rr);

		//	if(this.getElementCenter(i).el[1]<.052) B=new Vect(0,1);

		element[i].setB(B);

	}

	private void set3angElementAxiB(int i){

		Edge[] elEdge=elementEdges(i);
		double[] A=new double[3];
		for(int j=0;j<3;j++)
			A[j]=elEdge[j].A;


		Node[] vertexNode=elementNodes(i);


		double rr=0;
		for(int j=0;j<3;j++){
			double r1=vertexNode[j].getCoord(0)+1e-6;
			rr+=.33333333/r1;
		}

		Vect[] rotNe=femCalc.rotNe3ang(this.elementNodes(i));
		Vect B=new Vect(2);		
		for(int j=0;j<3;j++)
			B=B.add(rotNe[j].times(A[j]*rr));


		element[i].setB(B);

	}



	private void set3angElementB(int i){

		Edge[] elEdge=elementEdges(i);
		double[] A=new double[3];
		for(int j=0;j<3;j++)
			A[j]=elEdge[j].A;

		Vect[] rotNe=femCalc.rotNe3ang(this.elementNodes(i));
		Vect B=new Vect(2);		
		for(int j=0;j<3;j++)
			B=B.add(rotNe[j].times(A[j]));

		//	B=new Vect(this.getElement3angA(i),0);

		element[i].setB(B);
	}


	public void setElType(String type){
		elType=type;
		if(type.equals("triangle") ){
			elCode=0;
			nElVert=3;
			nElEdge=3;
			this.dim=2;
		}
		else if(type.equals("quadrangle") ){
			elCode=1;
			nElVert=4;
			nElEdge=4;
			this.dim=2;
		}
		else if(type.equals("tetrahedron") ){
			elCode=2;
			nElVert=4;
			nElEdge=6;
			dim=3;
		}
		else if(type.equals("prism") ){
			elCode=3;
			nElVert=6;
			nElEdge=9;
			dim=3;
		}
		else if(type.equals("hexahedron") ){
			elCode=4;
			nElVert=8;
			nElEdge=12;
			dim=3;
		}

		else if(type.equals("pyramid") ){
			elCode=5;
			nElVert=5;
			nElEdge=8;
			dim=3;
		}
		nBoundary=2*dim;
	}


	private Vect getVeclocity(int i, Vect[] gradN){

		int[] vert=element[i].getVertNumb();
		Vect Kt=element[i].getSigma();
		Vect B=new Vect(dim);
		for(int j=0;j<nElVert;j++)		{
			int nn=vert[j];
			B=B.add(gradN[j].times(node[nn].T));
		}

		return Kt.times(B);

	}


	private Vect getElementB(int i, Vect[] rotNe){

		Edge[] edge=elementEdges(i);
		Vect B=new Vect(dim);
		for(int j=0;j<nElEdge;j++)		{
		//	util.pr(edge[j].A);
		//	rotNe[j].hshow();
			
		B=B.add(rotNe[j].times(edge[j].A));
		}
	

		return B;

	}

	public double getElementQuadA(int i){
		Edge[] elEdge=elementEdges(i);
		double[] Ae=new double[4];
		for(int j=0;j<4;j++)
			Ae[j]=elEdge[j].A;
		double A=0;
		Vect zero=new Vect(2);

		double[] Ne=femCalc.NeQuad(zero);

		for(int j=0;j<nElEdge;j++)	{		
			A= A+Ne[j]*Ae[j];
		}
		return  A;	

	}



	public void setReluctForce(){
		forceCalc.setReluctForce(this);
	}

	public void resetReluctForce(){
		for(int i=1;i<=this.numberOfNodes;i++)
			if(this.node[i].F!=null) this.node[i].F=new Vect(dim);
	}


	public double getElementVolume(int i){

		if( this.elCode==2) return getElementVolumeTetra( i);
		if( this.elCode==3) return getElementVolumePrism( i);

		if(dim==2) return 0;
		double vol=0;
		Node[] vertexNode=elementNodes(i);
		Mat jac;
		double detJac,ws;
		Vect localCo=new Vect(dim);
		ws=8;
		jac=femCalc.jacobian(vertexNode,localCo);
		//detJac=abs(jac.determinant());
		detJac=jac.determinant();

		vol=detJac*ws;

		return vol;

	}

	public double getElementVolumePrism(int i){

		Node[] vertexNode=elementNodes(i);
		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=abs(v1.cross(v2).norm())/2;


		double h=vertexNode[3].getCoord().sub(vertexNode[0].getCoord()).norm();
		double vol=S*h;



		return vol;

	}


	public double getElementVolumeTetra(int i){

		double vol=0;
		Node[] vertexNode=elementNodes(i);
		Mat jac;
		double detJac,ws;
		Vect localCo=new Vect(4);
		jac=femCalc.jacobian(vertexNode,localCo);
		detJac=abs(jac.determinant())/6;

		vol=detJac;

		return vol;

	}


	public void setMagBC(){
		this.bcSet.setMagBC(this);
	}


	public void setEdge(){

		EdgeSet edgeSet=new EdgeSet();
		edgeSet.setEdge(this);

	}


	public void setSliceBounds(){
		this.bcSet.setSliceBounds(this);
	}

	public void setBounds(){
		if(this.coordCode==1)
			this.bcSet.setSliceBounds(this);
		else if(coordCode==0){

			double[] spb=new double[2*this.dim];


			for(int i=0;i<this.dim;i++){
				spb[2*i]=1e10;
				spb[2*i+1]=-1e10;
			}

			
			boolean[] nc=new boolean[1+this.numberOfNodes];
			int ix=0;
			for(int i=1;i<=this.numberOfElements;i++){		
				int[] vertNumb=this.element[i].getVertNumb();
				for(int j=0;j<nElVert;j++){
					int nx=vertNumb[j];
					if(!nc[nx]){

						nc[nx]=true;
					}
				}
			}

			for(int i=1;i<=this.numberOfNodes;i++){
				if(!nc[i]) continue;

				if(this.node[i].getCoord(0)<spb[0]) spb[0]=this.node[i].getCoord(0);
				else if(this.node[i].getCoord(0)>spb[1]) spb[1]=this.node[i].getCoord(0);

				if(this.node[i].getCoord(1)<spb[2]) spb[2]=this.node[i].getCoord(1);
				else if(this.node[i].getCoord(1)>spb[3]) spb[3]=this.node[i].getCoord(1);
				if(this.dim==3){
					if(this.node[i].getCoord(2)<spb[4]) spb[4]=this.node[i].getCoord(2);
					else if(this.node[i].getCoord(2)>spb[5]) spb[5]=this.node[i].getCoord(2);
				}


			}
			this.spaceBoundary=spb;

		}
	}


	public void setNodeOnBound(){
		this.bcSet.setNodeOnBound(this);
	}

	public void mapPBC(){
		this.bcSet.mapPBC(this);
	}


	public void setInUseNodes(){
		//=====================	identifying the used nodes
		for(int i=1;i<=this.numberOfElements;i++){

			int[] vert=this.element[i].getVertNumb();
			for(int j=0;j<this.nElVert;j++)	{	
				if(vert[j]>0)
					this.node[vert[j]].inUse=true;
			}
		}
	}



	public void setJ0(double dt){



		double Jn2=0,Jmax2=0,Jmin2=0;


		Vect regJ=new Vect(3);
		for(int ir=1;ir<=numberOfRegions;ir++){
			if(!region[ir].hasJ) continue;
			regJ=region[ir].getJ();
			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				if(this.coordCode==1 && this.dim==3)
				{

					Mat R=new Mat();


					Mat R2D=util.rotMat2D(util.getAng(this.getElementCenter(i).v2()));
					R=new Mat(dim,dim);
					for(int m=0;m<2;m++)
						for(int n=0;n<2;n++)
							R.el[m][n]=R2D.el[m][n];

					R.el[2][2]=1;

					element[i].setJ(R.mul(regJ));


				}

				else{
		
					element[i].setJ(regJ);
				}

				Jn2=regJ.dot(regJ);
				if(Jn2>Jmax2)
					Jmax2=Jn2;
				if(Jn2<Jmin2)
					Jmin2=Jn2;

			}
		}

		Jmax=sqrt(Jmax2);
		Jmin=sqrt(Jmin2);


	}

	public void setM(){

		Vect MReg=new Vect(dim);

		for(int ir=1;ir<=numberOfRegions;ir++){
			if(!region[ir].hasM) continue;


			MReg=region[ir].getM().times(region[ir].getNu());

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setM(MReg);
			}
		}

	}

	public void setFreq(double f){
		this.freq=f;

	}
	public void setDt(double dt){

		this.dt=dt;
	}
	public double getCurrentTime(){

		double time =this.currentTimeStep*dt;
		return time;
	}


	public void setnTsteps(int N){

		this.nTsteps=N;
	}

	public void setDefMode(int defMode){
		this.defMode=defMode;

	}


	public void writeNodalScalar(String file){
		writer.writeNodalScalar(this,file);
	}
	
	public void writePhi(String file){
		writer.writePhi(this,file);
	}


	public void setTorque(double r){
		setTorque(new Vect(dim),0,r,1);

	}

	public void setTorque(double r1, double r2){
		setTorque(new Vect(dim),r1,r2,1);
	}

	public void setTorque(double r1, double r2, int mode){

		setTorque(new Vect(dim),r1,r2,mode);

	}


	public void setTorque(Vect origin, double r1,double r2, int mode){

		int nx=0;

		double trq=0;
		Vect F=new Vect();
		for(int i=1;i<=this.numberOfNodes;i++)
		{

			//if(node[i].getMap()>0) {continue;}	
			
			
			Vect r2d=this.node[i].getCoord().v2();
	

			double rn=r2d.norm();
			if(rn<r1 ||rn>r2) continue;

			if(mode>=0) F=this.node[i].getNodalVect(mode);


			if(F!=null){
				nx++;
				double nodeTrq=0;
				if(this.coordCode==1){
					nodeTrq=rn*F.el[1];

				}
				else if(this.coordCode==0){
					nodeTrq=F.v2().cross(r2d).el[2];
				}
				trq=trq+nodeTrq;
			}

		}


		if(this.dim==2)
			this.TrqZ=trq*this.height;
		else this.TrqZ=trq;


		if(this.hasPBC)
			if(this.alpha2-this.alpha1>0){
				this.TrqZ*=2*PI/(alpha2-alpha1);
			}


	}

	public void setTorque(int ir){

		int mode=1;

		int[] nn=this.getRegNodes(ir);
		double trq=0;
		Vect F=new Vect();
		for(int i=0;i<nn.length;i++)
		{
			//if(node[i].getMap()>0) {continue;}		
			int nx=nn[i];
			double rn=this.node[nx].getCoord().v2().norm();

			//if(rn<r1 ||rn>r2) continue;

			if(mode>=0) F=this.node[nx].getNodalVect(mode);


			if(F!=null){

				trq=trq+rn*F.el[1];
			}

		}

		if(this.dim==2)
			this.TrqZ=trq*this.height;
		else this.TrqZ=trq;


		if(this.hasPBC)
			if(this.alpha2-this.alpha1>0)
				this.TrqZ*=2*PI/(alpha2-alpha1);


	}

	public void setJe(){
		double Jn2=0,Jmax2=0,Jmin2=0;
		for(int i=1;i<=numberOfElements;i++){

			if(!element[i].isConductor()) continue;

			Vect center=new Vect(dim);
			Vect Je=getElementJe(i, center);

			element[i].setJe(Je);
		
		
				Jn2=Je.dot(Je);
				if(Jn2>Jmax2)
					Jmax2=Jn2;
				if(Jn2<Jmin2)
					Jmin2=Jn2;

		}

		Jmax=sqrt(Jmax2);
		Jmin=sqrt(Jmin2);
	}
	
	public void setJPhiCoil(){
		
		double Jn2=0,Jmax2=0,Jmin2=0;
		
		int[] coilIndices=new int[numberOfRegions+1];

		for(int ir=1;ir<=numberOfRegions;ir++)
			coilIndices[ir]=-1;

		for(int ic=0;ic<phiCoils.length;ic++){
			int nr=phiCoils[ic].regNo;
			coilIndices[nr]=ic;

		}
		
		Vect lc=this.centerLocalCo();

		
		for(int ir=1;ir<=numberOfRegions;ir++){

			if(coilIndices[ir]<0) continue;

			double conductivity=1;//model.phiCoils[coilIndices[ir]].conductivity;
			if(coilIndices[ir]>0){
				double conductivity1=phiCoils[0].conductivity;
				double conductivity2=phiCoils[coilIndices[ir]].conductivity;
			if(conductivity1>0) conductivity=conductivity2/conductivity1;
			}

			region[ir].hasJ=true;
			
			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

			Vect J=getElementJPhiCoil(i,lc,conductivity);
					
			element[i].setJ(J);

				Jn2=J.dot(J);
				if(Jn2>Jmax2)
					Jmax2=Jn2;
				if(Jn2<Jmin2)
					Jmin2=Jn2;


		}
		}
		Jmax=sqrt(Jmax2);
		Jmin=sqrt(Jmin2);
	}

	
	public void setJStatic(){
		double Jn2=0,Jmax2=0,Jmin2=0;
		for(int i=1;i<=numberOfElements;i++){

			if(!element[i].isConductor()) continue;

			setElementJStatic(i);

			Vect Je=element[i].getJe();
			//Je.hshow();
			if(element[i].isConductor()){
				Jn2=Je.dot(Je);
				if(Jn2>Jmax2)
					Jmax2=Jn2;
				if(Jn2<Jmin2)
					Jmin2=Jn2;

			}

		}

		Jmax=sqrt(Jmax2);
		Jmin=sqrt(Jmin2);
	}

	
	public Vect getElementJe(int i, Vect lc){
		if(dim==2) {return getElement2DJe(i,lc); }
		if(elCode==2) {return getElementJeTetra(i,lc); }
		Vect sigma=element[i].getSigma();
		
		Node[] vertex=elementNodes(i);
		Edge[] edge=elementEdges(i);
		double[] dAe=new double[nElEdge];
		Vect dA=new Vect(dim);
		for(int j=0;j<nElEdge;j++)
			dAe[j]=edge[j].getDiffA();

		dA=getElementdA(vertex,dAe,i,lc);

		double rdt=1.0/dt;
		
		if(this.AC) rdt=this.freq;
		
		Vect Je=new Vect(3);
			
			Je=dA.times(rdt);

			if(analysisMode==2){
			double[] nodePhi=new double[nElVert];
			Vect gradPhi=new Vect(dim);

			for(int j=0;j<nElVert;j++)
				nodePhi[j]=vertex[j].getPhi();
			
			gradPhi=femCalc.gradPhi(vertex,nodePhi,lc);
			
			
			Je=Je.add(gradPhi);
		}
		
			Je=Je.times(sigma.times(-1.));

	return Je;

	}

	private void setElementJStatic(int i){

		if(this.dim==2){
			setElementJStaticJe2D(i);
			return;
		}
		Node[] vertexNode=elementNodes(i);
		Edge[] elemEdges=elementEdges(i);
		boolean[] edgeDir=element[i].getEdgeReverse();

		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	
			A[j]=elemEdges[j].getA();

			double[] nodePhi=new double[nElVert];
			Vect gradPhi=new Vect(dim);

			Vect localCo=new Vect(this.dim);

			for(int j=0;j<nElVert;j++){
				nodePhi[j]=vertexNode[j].getPhi();
			gradPhi=femCalc.gradPhi(vertexNode,nodePhi,localCo);
			}
			
			
			Mat jac=femCalc.jacobian(vertexNode,localCo);

			
			Vect[] Ne=femCalc.Ne(jac,localCo,edgeDir);
			
			Vect Adot=new Vect(dim);
			
			for(int j=0;j<this.nElEdge;j++){
				Adot=Adot.add(Ne[j].times(A[j]));
			
			}
			
			element[i].setJe(gradPhi.add(Adot).times(element[i].getSigma()).times(-1));
			
}
	private Vect getElementJPhiCoil(int i, Vect lc,double sigma){

		if(dim==2){
			return getElementJPhiCoil2D(i,lc,sigma);

		}

		Node[] vertexNode=elementNodes(i);

			double[] nodePhi=new double[nElVert];
			Vect gradPhi=new Vect(dim);


			for(int j=0;j<nElVert;j++){
				nodePhi[j]=vertexNode[j].getPhi();

			}
			if(elCode==4)
				gradPhi=femCalc.gradPhi(vertexNode,nodePhi,lc);	
			else if(elCode==2)
				gradPhi=femCalc.gradPhiTet(vertexNode,nodePhi);

			Vect J=gradPhi.times(sigma).times(-1);
			

			return J;
		//	element[i].setJ(J);


}
	

	private Vect getElementJPhiCoil2D(int i,Vect lc,double sigma){

			

		Vect J=new Vect(3);
	
				double gradPhiz=0;
				if(elCode==0){
					gradPhiz=getElement3angPhi(i,lc);
				}
				else if(elCode==1){
					gradPhiz=getElementQuadPhi(i,lc);
					
				}
				J=new Vect(0,0,-(gradPhiz)*sigma);

	
		return J;

}

	private Vect getElement2DJe(int i, Vect lc){

		double rdt=1.0/dt;
		
		if(this.AC) rdt=this.freq;

		Edge[] edge=elementEdges(i);
		double[] dAe=new double[nElEdge];
		double dA=0;
		for(int j=0;j<nElEdge;j++){
			dAe[j]=edge[j].getDiffA();
		}

		if(elCode==0){
			dA=getElement3angdA(i,dAe,lc);
		}
		else if(elCode==1)
			dA=getElementQuaddA(dAe,new Vect(2));

		
		Vect Je=new Vect(0,0,-dA*rdt);


		 if(analysisMode==2){
			double gradPhiz=0;
			if(elCode==0) gradPhiz=getElement3angPhi(i,lc);
			else if(elCode==1) gradPhiz=getElementQuadPhi(i,lc);
			
			Je=Je.add(new Vect(0,0,-gradPhiz));
		}
		 Je=Je.times(element[i].getSigmaZ());
		 
		 return Je;
	}
	
	private Vect getElementJeTetra(int i, Vect lc){

		double rdt=1.0/dt;
		if(this.AC) rdt=this.freq;

		Edge[] edge=elementEdges(i);
		double[] dAe=new double[nElEdge];
		
		for(int j=0;j<nElEdge;j++){
			dAe[j]=edge[j].getDiffA();
		}
		if(lc.norm()==0)
		 lc=new Vect().ones(4).times(.25);
		
		Node[] vertexNode=this.elementNodes(i);
		Mat jac=femCalc.jacobianTetra(vertexNode,lc);
		boolean[] edgeReverse=element[i].getEdgeReverse();

		Vect dA=new Vect(3);

		Vect[] Ne=femCalc.NeTetra(jac, lc, edgeReverse);
		
		for(int j=0;j<nElEdge;j++)	{		
			dA=dA.add(Ne[j].times(dAe[j]));
		}

		Vect Je=dA.times(rdt);


		if(analysisMode==2){
		double[] nodePhi=new double[nElVert];
		Vect gradPhi=new Vect(dim);

		for(int j=0;j<nElVert;j++)
			nodePhi[j]=vertexNode[j].getPhi();
		
		gradPhi=femCalc.gradPhiTet(vertexNode,nodePhi);
		
		
		Je=Je.add(gradPhi);
	}

		Vect sigma=element[i].getSigma();
		
		Je=Je.times(sigma.times(-1.));

		return Je;

	}
	
	private void setElementJStaticJe2D(int i){

		double rdt=1.0/this.height;
	
		Vect lc=new Vect(this.dim);
		if(this.elCode==0) lc=new Vect(1./3,1./3,1./3);

		double A=0;

		if(elCode==0)
			A=getElementA3ang(i,lc);
		else if(elCode==1)
			A=getElementAQuad(i,lc);
	
		
			double gradPhiz=0;
			if(elCode==0) gradPhiz=getElement3angPhi(i,lc);
			else if(elCode==1) gradPhiz=getElementQuadPhi(i,lc);
			element[i].setJe(new Vect(0,0,-(A+gradPhiz)*element[i].getSigmaZ()*rdt));
			
		}	

	public Vect getElementA(int ie){

		Vect localCo=centerLocalCo();
		
		return getElementA(ie,localCo);
	}
	
	public Vect centerLocalCo(){
		
		Vect lc=null;
		if(this.elCode==0) lc=new Vect(1,1,1).times(1./3);
		else if(this.elCode==1) lc=new Vect(2);
		else if(this.elCode==2) lc=new Vect().ones(4).times(1./4);
		else if(this.elCode==4) lc=new Vect(3);
		
		return lc;
	}

	public Vect getElementA(int ie,Vect lc){

		boolean[] edgeDir=element[ie].getEdgeReverse();

		Vect A=new Vect(dim);
		Node[] vertexNode=this.elementNodes(ie);
		Edge[] edge=this.elementEdges(ie);
		Mat jac=femCalc.jacobian(vertexNode,lc);
		Vect[] Ne=femCalc.Ne(jac,lc,edgeDir);

		for(int j=0;j<nElEdge;j++)	{		
			//util.ph(edge[j].A);
			//Ne[j].hshow();
			A= A.add(Ne[j].times(edge[j].A));
		}
	//	A.show();
		return  A;	

	}

	public double getElementPhi(int ie){

		double phi=0;
		Vect zero=new Vect(dim);
		Node[] vertexNode=this.elementNodes(ie);
		double[] N=femCalc.N(zero);

		for(int j=0;j<nElVert;j++)	{		
			phi+=N[j]*vertexNode[j].getPhi();
		}
		return  phi;	

	}



	public double getElementT(int ie){

		double T=0;
		Vect zero=new Vect(dim);
		Node[] vertexNode=this.elementNodes(ie);
		double[] N=femCalc.N(zero);

		for(int j=0;j<nElVert;j++)	{		
			T+=N[j]*vertexNode[j].T;
		}
		return  T;	

	}

	public Vect getElementdA(Node[] vertexNode,double[] Ae,int i,Vect lc){

		boolean[] edgeDir=element[i].getEdgeReverse();

		Vect dA=new Vect(dim);

		Mat jac=femCalc.jacobian(vertexNode,lc);
		
		Vect[] Ne=femCalc.Ne(jac,lc,edgeDir);

		for(int j=0;j<nElEdge;j++)	{		
			dA= dA.add(Ne[j].times(Ae[j]));
		}
		return  dA;	

	}
	
	public double getElementQuadA(double[] dAe,Vect lc){

		double dA=0;
		double[] Ne=femCalc.NeQuad(lc);
		for(int j=0;j<nElEdge;j++)	{		
			dA+=Ne[j]*dAe[j];
		}
		return  dA;	

	}


	public double getElementQuadPhi(int ie,Vect lc){

		double phi=0;

		Node[] vertexNode=this.elementNodes(ie);
		double[] N=femCalc.NQuad(lc);

		for(int j=0;j<nElVert;j++)	{		
			phi+=N[j]*vertexNode[j].getPhi();
		}
		return  phi;	

	}

	public double getElementQuaddA(double[] dAe, Vect lc){

		double dA=0;

		double[] Ne=femCalc.NeQuad(lc);
		for(int j=0;j<nElEdge;j++)	{		
			dA+=Ne[j]*dAe[j];
		}
		return  dA;	

	}

	public double getElement2DA(int ie,Vect lc){
		double A=0;
		if(this.elCode==0)
			A= getElementA3ang(ie,lc);
		else if(this.elCode==1)
			A= getElementAQuad(ie,lc);


		return  A;	

	}
	


public double getElementA3ang(int ie,Vect lc){

		int[] en=this.element[ie].getEdgeNumb();


		double[] localNe=lc.el;
		double A=0;
		for(int j=0;j<3;j++)	{		
			double Aj=0;
			if(this.edge[en[j]].map==0)
				Aj=this.edge[en[j]].A;
			else
				Aj=this.edge[this.edge[en[j]].map].A;

			A+=localNe[j]*Aj;
		}
		return  A;	

	}

	public double getElementAQuad(int ie,Vect lc){
	
	
	
		double Az=0;

		double[] N=femCalc.NeQuad(lc);
	
		for(int j=0;j<N.length;j++)	{		

			Az+= N[j]*edge[j].A;
		}
	
		return  Az;	
	
	}

	public double getElement3angPhi(int ie,Vect lc){

		double phi=0;
		Node[] vertexNode=this.elementNodes(ie);
	//	double a=1.0/3;
		double[] localNe=lc.el;
		for(int j=0;j<nElVert;j++)	{		
			phi+=localNe[j]*vertexNode[j].getPhi();
		}
		return  phi;	

	}

	public double getElement3angdA(int ie,double[] dAe, Vect lc){

		
		double[] localNe=lc.el;
		double dA=0;
		for(int j=0;j<3;j++)	{		
			dA+=localNe[j]*dAe[j];
		}
		return  dA;	

	}

	public void setElementsParam(){



		int nx=0,ny=0;
		for(int ir=1;ir<=this.numberOfRegions;ir++){
			if(this.region[ir].isNonLinear){				
				nx++;
				this.region[ir].BHnumber=nx;
			}else{
				this.region[ir].BHnumber=0;
				}

			if(this.region[ir].MS){				
				ny++;
				this.region[ir].lamBNumber=ny;
			}

		}
		this.nBH=nx;
		this.nLam=ny;


			this.BH=new BHCurve[this.nBH+1];
			this.lamB=new LamBCurve[this.nLam+1];
	

		try {

			for(int ir=1;ir<=this.numberOfRegions;ir++){

				if(this.region[ir].lamBNumber>0)		
					{
					String file=this.resultFolder+"\\"+region[ir].getMaterial()+".txt";
						this.lamB[this.region[ir].lamBNumber]=new LamBCurve(file) ;
					}


				if(this.region[ir].BHnumber>0)	{
			
					String file=this.resultFolder+"\\"+region[ir].getMaterial()+".txt";
					
						this.BH[this.region[ir].BHnumber]=new BHCurve(file) ;
						double nux=this.BH[this.region[ir].BHnumber].getNu(1.0);							
						region[ir].setNu(nux);
				}
					


			}


		} catch (Exception e) {
			System.err.println("BH data file was not found!");
		}


		setM();

		nEdgeHasJ=0;


		for(int ir=1;ir<=numberOfRegions;ir++){

			boolean regCond=region[ir].isConductor;

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

					element[i].setNu(region[ir].getNu());
				element[i].setSigma(region[ir].getSigma());

				if(regCond)
					element[i].setSigma(region[ir].getSigma());

				element[i].setRegion(ir);



			}
		}
	

		if(elType.equals("triangle")) elCode=0;
		else if(elType.equals("quadrangle")) elCode=1;
		else if(elType.equals("tetrahedron")) elCode=2;
		else if(elType.equals("prism")) elCode=3;
		else if(elType.equals("hexahedron")) elCode=4;

	}


	public void setElementsParamMech(){

		for(int ir=1;ir<=numberOfRegions;ir++){

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){


				element[i].setRegion(ir);

			}
		}


		if(elType.equals("triangle")) elCode=0;
		else if(elType.equals("quadrangle")) elCode=1;
		else if(elType.equals("tetrahedron")) elCode=2;
		else if(elType.equals("prism")) elCode=3;
		else if(elType.equals("hexahedron")) elCode=4;

	}

	public void setElementsParamSeep(){

		for(int ir=1;ir<=numberOfRegions;ir++){

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setSigma(region[ir].getSigma());

			}
		}


		if(elType.equals("triangle")) elCode=0;
		else if(elType.equals("quadrangle")) elCode=1;
		else if(elType.equals("tetrahedron")) elCode=2;
		else if(elType.equals("prism")) elCode=3;
		else if(elType.equals("hexahedron")) elCode=4;

	}


	public void setElementsParamPC(){

		for(int ir=1;ir<=numberOfRegions;ir++){

			for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

				element[i].setSigma(region[ir].getEr());
				element[i].setNu(region[ir].getMur());

			}
		}


		if(elType.equals("triangle")) elCode=0;
		else if(elType.equals("quadrangle")) elCode=1;
		else if(elType.equals("tetrahedron")) elCode=2;
		else if(elType.equals("prism")) elCode=3;
		else if(elType.equals("hexahedron")) elCode=4;

	}



	public void scaleKnownEdgeAL(double c){
		for(int i=1;i<=numberOfKnownEdges;i++)
			if(knownEdgeValue[i]!=0)
				edge[knownEdgeNumber[i]].setSolvedAL(knownEdgeValue[i]*c);	
	}


	public void saveAp(){
		for(int i=1;i<=this.numberOfEdges;i++){
			edge[i].saveAp();	


		}

	}



	public Vect solveMagLin(int step,Vect x_init){

		return femSolver.solveMagLin(this,step,x_init);
	}

	public Vect solveNonLinear(Vect x, boolean b,int step){
		
		return femSolver.solveMagNonlin(this,x, b,step);
	}



	public void setNodePhi(Vect x){
		int nodeNumber;

		for(int i=1;i<=numberOfVarNodes;i++){
			nodeNumber=varNodeNumber[i];	
			if(nodeNumber>0){
				node[varNodeNumber[i]].setPhi(x.el[i+numberOfUnknownEdges-1]);	
			}
		}

	}

	public void setSolution(Vect x){

		setSolvedAL(x);

		if(analysisMode==2)
			setNodePhi(x);
		

		this.setB();	

	}

	public Vect getUnknownA(){
		Vect x=new Vect(this.numberOfUnknownEdges);
		for(int i=0;i<x.length;i++){
			x.el[i]=edge[unknownEdgeNumber[i+1]].A;	
		}

		return x;
	}

	public Vect getUnknownAp(){
		Vect x=new Vect(this.numberOfUnknownEdges);
		for(int i=0;i<x.length;i++)
			x.el[i]=edge[unknownEdgeNumber[i+1]].Ap;	

		return x;
	}

	public int[] getRegNodes(int ir){

		boolean[] nc=new boolean[1+this.numberOfNodes];
		int[] nn=new int[this.numberOfNodes];
		int ix=0;
		for(int i=this.region[ir].getFirstEl();i<=this.region[ir].getLastEl();i++){		
			int[] vertNumb=this.element[i].getVertNumb();
			for(int j=0;j<nElVert;j++){
				int nx=vertNumb[j];
				if(!nc[nx]){

					nc[nx]=true;
					nn[ix]=nx;
					ix++;
				}
			}
		}

		int[] regNodes=new int[ix];
		for(int i=0;i<ix;i++)
			regNodes[i]=nn[i];

		return regNodes;


	}

	public void setB(){

		double Bn2,Bmax2=0,Bmin2=0;

		for(int i=1;i<=numberOfElements;i++){
			
			setElementB(i);
			//element[i].getB().hshow();
			Bn2=element[i].getB().dot(element[i].getB());
			if(Bn2>Bmax2)
				Bmax2=Bn2;
			if(Bn2<Bmin2)
				Bmin2=Bn2;}

		Bmax=sqrt(Bmax2);
		Bmin=sqrt(Bmin2);

	}


	public void setHead(Vect h){
		for(int i=1;i<=this.numberOfNodes;i++){
			if(this.T_unknownIndex[i]!=0)
			{
				this.node[i].T=h.el[this.T_unknownIndex[i]-1];

			}
		}
	}
	public void setVelocity(){

		double Bn2,Bmax2=0,Bmin2=0;

		for(int i=1;i<=numberOfElements;i++){

			if (!element[i].isConductor()) continue;

			setElementVelocity(i);

			Bn2=element[i].getB().dot(element[i].getB());
			if(Bn2>Bmax2)
				Bmax2=Bn2;
			if(Bn2<Bmin2)
				Bmin2=Bn2;}

		Bmax=sqrt(Bmax2);
		Bmin=sqrt(Bmin2);


	}


	public void setElementVelocity(int i){
		Node[] vertexNode=elementNodes(i);
		Vect zero=new Vect(2);
		Mat jac=new Mat();
		Vect[] gradN;
		Vect B=new Vect(2);
		Vect Kt=element[i].getSigma();

		jac=femCalc.jacobian(vertexNode,zero);
		gradN=femCalc.gradN(jac,zero);

		int[] vn=element[i].getVertNumb();

		for(int j=0;j<this.nElVert;j++){
			double T=node[vn[j]].T;
			B=B.add(gradN[j].times(Kt).times(-T));
		}

		element[i].setB(B);

	}

	public Vect getBAt(Vect P){
		return femCalc.getBAt(this, P);
	}

	public double[] getApAnAt(Vect P){

		return femCalc.getApAnAt(this, P);
	}



	public Vect[] getAllB(){	

		Vect[] B=new Vect[numberOfElements];

		for(int ie=0;ie<numberOfElements;ie++)
			B[ie]=element[ie+1].getB();
		return B;
	}

	public double getDiffMax(Vect[] u, Vect[] v){
		double diff;
		double diffMax=0;
		for(int i=0;i<u.length;i++){


			diff=u[i].sub(v[i]).norm();
			if(diff>diffMax)
				diffMax=diff;
		}

		return diffMax;
	}

	public double getBmax(){	

		double Bmax=0;

		for(int i=1;i<=numberOfElements;i++){
			double Bn=element[i].getB().norm();
			if(Bn>Bmax) Bmax=Bn;

		}
		return Bmax;
	}
	
	public double getFluxErrSquared(Vect[] u, Vect[] v){	

		double err_sqared=getErrorSquared(u,v)/getSumB_Squared(u,v);
		
		return err_sqared;
	}
	
	public double getErrorSquared(Vect[] u, Vect[] v){	

		double sum=0;

		for(int i=0;i<u.length;i++){


			sum+=u[i].sub(v[i]).norm2();
		
		}
		return sum;
	}

	public double getSumB_Squared(Vect[] u, Vect[] v){	

		double sum=0;
		for(int i=0;i<u.length;i++){

			double Bn2=u[i].norm2()+v[i].norm2();
			sum+=Bn2;

		}
		return sum;
	}

	public double getSumB_Squared(){	

		double sum=0;

		for(int i=1;i<=numberOfElements;i++){
			double Bn2=element[i].getB().norm2();
			sum+=Bn2;

		}
		return sum;
	}

	public double getRmax(){	

		double rmax=0;

		for(int i=1;i<=numberOfNodes;i++){

			double rn=node[i].getCoord().v2().norm();
			if(rn>rmax) rmax=rn;

		}
		return rmax;
	}

	public void setAuMax(){	

		AuMax=0;

		for(int i=1;i<=numberOfEdges;i++){
			double aum=abs(edge[i].A);
			if(aum>AuMax) AuMax=aum;

		}

	}

	public void setuMax(){	

		uMax=0;

		for(int i=1;i<=numberOfNodes;i++){
			if(node[i].u==null) continue;
			double a=node[i].u.norm();
			if(a>uMax) uMax=a;

		}

	}

	public Vect[] getAllU(){	

		Vect[] u=new Vect[this.numberOfUnknownU];

		for(int i=0;i<this.numberOfUnknownU;i++){

			u[i]=node[unknownUnumber[i+1]].u.deepCopy();
		}
		return u;
	}


	public Vect getFOf(Vect globalCo){

		int m[]=getContainingElement(globalCo);

		if(m[0]<=0) throw new NullPointerException("given point outside the space ");

		if(element[m[0]].getNu().norm()==nu0)
			return new Vect(dim);
		Vect F=new Vect(dim);
		for(int j=0;j<nElVert;j++)
			F=F.add(node[element[m[0]].getVertNumb(j)].F);

		return  F.times(1.0/nElVert);	

	}

	public Vect getJeAt(Vect globalCo){

		int[] m=getContainingElement(globalCo);
		if(m[0]<=0) throw new NullPointerException("given point outside the space ");
		else if(!this.element[m[0]].isConductor()) return new Vect(3);
		return this.element[m[0]].getJe();

	}

	public Vect getJeOf(int m){
		if(!element[m].isConductor())
			return new Vect(3);

		Vect Jn=new Vect(3);
		if(m>0)
			Jn=element[m].getJe();

		else if(m<-1)
			Jn=element[-m-1].getJe();

		return Jn;
	}

	public void writeB(String file){
		writer.writeB(this,file);

	}


	public void writeA(String file){
		writer.writeA(this, file);

	}
	public void writeNodalA2D(String file){
		

		for(int n=1;n<=this.numberOfEdges;n++)
		{	
			this.edge[n].node[0].T=this.edge[n].A;
		}		
		
		writer.writeNodalScalar(this, file);

	}


	public void writeNodalField(String nodalForceFile,int mode){

		writer.writeNodalField(this, nodalForceFile, mode);
	}


	public void writeJe(String eddyFile){
		writer.writeJe(this, eddyFile);

	}
	

	public double obtainLoss(int ir){

			
		double loss=0;

	for(int i=region[ir].getFirstEl();i<=region[ir].getLastEl();i++){

		double elementLoss=femCalc.obtainElementLoss(this,i);
		
		loss+=elementLoss;
			}
		
		
			return loss;
	}

	public void writeJ0(String j0filr){
		writer.writeJ0(this, j0filr);

	}



	public void rotate(double rad){
		MeshManipulator mf=new MeshManipulator();
		mf.rotate(this,rad,false);	
	}



	private double[] subS3ang(int ie, Vect P){

		Node[] vertexNode=elementNodes(ie);
		Vect v1=vertexNode[0].getCoord().sub(P);
		Vect v2=vertexNode[1].getCoord().sub(P);
		Vect v3=vertexNode[2].getCoord().sub(P);

		double[] N=new double[3];

		N[0]=v2.cross(v3).norm()/2;
		N[1]=v3.cross(v1).norm()/2;
		N[2]=v1.cross(v2).norm()/2;
		return N;
	}


	public Vect getElementCenter(int ie){


		Vect center=new Vect(dim);
		int[] vertNumb=element[ie].getVertNumb();
		for(int j=0;j<nElVert;j++)
			center=center.add(node[vertNumb[j]].getCoord());


		return center.times(1.0/nElVert);
	}

	public void resetMotor(){


		for(int i=1;i<=numberOfNodes;i++){
			if(node[i].getMap()>0){
				node[i].common=false;
				node[i].sPBC=false;
				node[i].aPBC=false;
				node[i].setMap(0);
			}
		}

		for(int i=1;i<=numberOfEdges;i++){
			if(edge[i].map>0){
				edge[i].common=false;
				edge[i].sPBC=false;
				edge[i].aPBC=false;
				edge[i].map=0;
			}
		}

		int L=0;
		if(commonNodes!=null){
			L=commonNodes[0][0].length;

		for(int iz=0;iz<commonNodes.length;iz++)
		for(int i=0;i<L;i++){
			int nn=commonNodes[iz][0][i];
			node[nn].setMap(commonNodes[iz][1][i]);
			node[nn].common=true;
		}
		}

		if(this.hasPBC)
			this.mapPBC();
	}


	public void resultAt(Vect P){
		System.out.println();
		Vect Bn=this.getBAt(P);
		if(Bn.el[0]==1e10 && Bn.el[1]==-1e10) util.pr(" >>>>  Given point is located outside the analyzed space!");
		{
			System.out.print("At : ");	P.hshow();
			System.out.print("B(muT) : ");	Bn.times(1e6).hshow();

			System.out.println();

			if(this.analysisMode>0){

				Vect Jen=this.getJeAt(P);
				System.out.print("J : ");	Jen.hshow();
				System.out.println();
				System.out.println();
			}
		}
	}

	public double getSliceAngle(){
		return this.alpha2-this.alpha1;
	}

public void solveCoils(){
	
	if(phiCoils==null) return;
						

	phiSolver= new StaticElectricSolver();

	phiSolver.setBoundaryCondition(this);
	phiSolver.setMatrix(this);
	phiSolver.setRHS0(this);

	Vect fi=phiSolver.solve(this);

	phiSolver.setSolution(this,fi);
	this.setJPhiCoil();
	this.writeJ0(this.resultFolder+"\\J"+0+".txt");
	//writePhi(this.resultFolder+"\\phi"+0+".txt");

	for(int i=1;i<=this.numberOfNodes;i++){
	this.node[i].setPhiVar(false);
	this.node[i].setPhiKnown(false);
	}
	
	
	
}
}




