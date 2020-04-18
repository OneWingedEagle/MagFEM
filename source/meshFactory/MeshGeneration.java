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



public class MeshGeneration {

	String regex="[ : ,=\\t]+";
	Writer writer=new Writer();

	public static void main(String[] args){
		
		MeshGeneration mf=new MeshGeneration();

		mf.meshQ();
		//mf.mesh16WiresTwisted();
		//mf.meshHexa();
		
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
			model.element[i]=new Element("quadrangle");
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
	

	
	public void meshQ(){
		
		
			//double[][] bb={{471,492.32,0,12.5},{478.43,481.23,.55,6.25},{478.43,481.23,.55,6.25},{478.43,481.23,.55,6.25},{478.43,481.23,.55,6.25}};
		//	double[][] bb={{0,1000,-1000,1000},{-1000,1000,-1000,-300},{-400,-300,0,100},{300,400,0,100}};
			//double[][] bb={{-300,300,-300,300},{-60,60,-60,60},{-40,40,-40,40},{40,60,-5,5}, {-35,0,-35,35},{-100,-65,-35,35}};
			//double[][] bb={{-100,100,-100,100}};
			
			double[][] bb={{1,10,0,10}};
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
	

			//mg.blockName[1]="core";
			//mg.blockName[2]="air";
		//	mg.blockName[3]="air";
			//mg.blockName[4]="coil1";
			//mg.blockName[5]="coil2";

			
for(int j=0;j<bb.length;j++){
				
				for(int k=0;k<bb[0].length;k++){
				
			
				mg.baseLeft[j][k]=1.;
				mg.baseRight[j][k]=1.;
					
			
			
		
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
	

	
	public void meshHexa(){

	double[][] bb={{-250,250,-250,250,0,600},{-50,50,-50,50,0,600}};
	//	double[][] bb={{-100,100,-100,100,-100,100},{-30,30,-30,30,-5,5}};


		double scale=1;
	//	scale=1000;

		for(int j=0;j<bb.length;j++)
			for(int k=0;k<bb[j].length;k++){
				bb[j][k]*=scale;
			
			}

		Geometry mg=new Geometry(bb);

		mg.blockName[0]="outPM";
		mg.blockName[1]="hhh";
	

		
for(int j=0;j<bb.length;j++){
			
			for(int k=0;k<bb[j].length;k++){
			
		
			mg.baseLeft[j][k]=1;
			mg.baseRight[j][k]=1;
	
			if(j!=-1){
				mg.minMeshRight[j][k]=99;
				mg.minMeshLeft[j][k]=99;
			} else {
				//mg.minMeshRight[j][k]=.3;
				//mg.minMeshLeft[j][k]=.3;
			}
		
		
			}	
}


		


		Model model=getOrthogMesh(mg,"",false);

		String bun=System.getProperty("user.dir") + "\\model.txt";
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
	

	public void mesh16WiresTwisted(){

	double[][] bb={{0,5,100,105,0,5},{0,5,110,115,0,5},{0,5,100,105,10,15},{0,5,110,115,10,15}};
	//	double[][] bb={{-100,100,-100,100,-100,100},{-30,30,-30,30,-5,5}};


		double scale=1;
		//scale=1000;

		for(int j=0;j<bb.length;j++)
			for(int k=0;k<bb[0].length;k++){
				bb[j][k]*=scale;
			
			}

		Geometry mg=new Geometry(bb);

		//mg.blockName[0]="outPM";
	//	mg.blockName[1]="hhh";
	

		
for(int j=0;j<bb.length;j++){
			
			for(int k=0;k<bb[0].length;k++){
			
		
			//mg.baseLeft[j][k]=1;
			//mg.baseRight[j][k]=1;
	
			//if(j==1){
				mg.minMeshRight[j][k]=5;
				mg.minMeshLeft[j][k]=5;
		//	} else {
				//mg.minMeshRight[j][k]=.3;
				//mg.minMeshLeft[j][k]=.3;
			//}
		
		
			}	
}


		


		Model model=getOrthogMesh(mg,"",false);

		String bun=System.getProperty("user.dir") + "\\model.txt";
		model.writeMesh(bun);

	

}


}
