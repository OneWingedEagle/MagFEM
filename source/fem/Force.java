package fem;

import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import math.Eigen;
import math.IntVect;
import math.Mat;
import math.Vect;
import math.util;

public class Force {

	public int dim,elCode,nElVert,nElEdge,numberOfElements,numberOfNodes,numberOfRegions;
	private Calculator femCalc;
	private boolean coupled,centerMST;
	private double[][] PW,PW3ang,PWNL;
	Eigen eg=new Eigen();

	public Force()
	{	}

	public Force(Model model)
	{
		this.nElVert=model.nElVert;
		this.nElEdge=model.nElEdge;
		this.numberOfElements=model.numberOfElements;
		this.numberOfNodes=model.numberOfNodes;
		this.numberOfRegions=model.numberOfRegions;
		this.dim=model.dim;
		this.elCode=model.elCode;
		this.femCalc=model.femCalc;
		this.coupled=model.coupled;
		this.PW=this.femCalc.gaussInteg(2);
		this.PWNL=this.femCalc.gaussInteg(3);
		this.PW3ang=this.femCalc.gaussInteg3(4);

		this.centerMST=true;
	}

	public void setReluctForce(Model model){
			boolean sf=false;
			//setMagForceDensity(model,0);

			int ext=4;
		IntVect[] nodeElement=new 	IntVect[this.numberOfNodes+1];
		
		for(int i=1;i<=this.numberOfNodes;i++)
			nodeElement[i]=new IntVect(10);
		
		int[] indx=new int[this.numberOfNodes+1];
		
		for(int ir=1;ir<=this.numberOfRegions;ir++)
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){	
			int[] vertNumb=model.element[i].getVertNumb();
			for(int j=0;j<this.nElVert;j++){
				int nn=vertNumb[j];
				if(model.node[nn].common) continue;
				nodeElement[nn].el[indx[nn]++]=i;
				if(indx[nn]==nodeElement[nn].length-1)
					nodeElement[nn].extend(ext);
			}
		}
		
		boolean[] elementInAir= new boolean[this.numberOfElements+1];
		boolean[] elementHasF= new boolean[this.numberOfElements+1];
		
		Vect nu0=new Vect().ones(this.dim).times(model.nu0);
		double err=1e-6*model.nu0;
		boolean nonLin,conduct,hasM;
		for(int ir=1;ir<=this.numberOfRegions;ir++)
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
			elementInAir[i]=true; 
			int[] vertNumb=model.element[i].getVertNumb();

	
			for(int j=0;j<this.nElVert;j++){
				int nn=vertNumb[j];				
				for(int k=0;k<indx[nn];k++){
					int ne=nodeElement[nn].el[k];
					double dn=model.element[ne].getNu().sub(nu0).norm();
					nonLin=model.element[ne].isNonlin();
					conduct=model.element[ne].hasJ();
					hasM=model.element[ne].hasM();
					

					if(nonLin || conduct || hasM || dn>err) {
						elementInAir[i]=false; 
						elementHasF[ne]=true;
					}
				}
			}
		}
		
		for(int ir=1;ir<=this.numberOfRegions;ir++)
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				if(elementHasF[i])
				{
				int[] vertNumb=model.element[i].getVertNumb();
				for(int j=0;j<this.nElVert;j++){
					int nn=vertNumb[j];
					//if(model.node[nn].common && !model.node[nn].PBC) continue;					
					model.node[nn].setHasF(true);
			
				}
					
				}
			}
				
			for(int i=1;i<=this.numberOfNodes;i++)
				if(model.node[i].hasF())
					model.node[i].setF(new Vect(this.dim));
			Vect[] nodalForce;
			int nBH;
			for(int ir=1;ir<=this.numberOfRegions;ir++){
				nBH=model.region[ir].BHnumber;
				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
					if(elementInAir[i]) continue;
					if(sf){

						Vect[] fs=surfForce(model,0,nBH,i);

						int[] edgeNumb=model.element[i].edgeXYNumb;

						for(int j=0;j<this.nElEdge;j++){

/*							int n0=model.edge[edgeNumb[j]].endNodeNumber[0];
							int n1=model.edge[edgeNumb[j]].endNodeNumber[1];*/
							
							int n0=model.edge[edgeNumb[j]].node[0].id;
							int n1=model.edge[edgeNumb[j]].node[1].id;


							if(model.node[n0].hasF())
								model.node[n0].F=model.node[n0].F.add(fs[j].times(0.5));

							if(model.node[n1].hasF())
								model.node[n1].F=model.node[n1].F.add(fs[j].times(0.5));
						}
					}

					//===========
						else {
							nodalForce=nodalForce(model,nBH,i);


							int[] vertNumb=model.element[i].getVertNumb();
										

							for(int j=0;j<this.nElVert;j++){

								int nn=vertNumb[j];

								if(model.node[nn].hasF() )
								{					
								model.node[nn].F= model.node[nn].F.add(nodalForce[j]);
								}
							
							}

						}
				}
			}
			
			


			
			boolean hasNeumann=false;

			for(int j=0;j<2*this.dim;j++)
				if(model.BCtype[j]==0 || model.PBCpair[j]==-2 ) {hasNeumann=true; break;}
			
			if(hasNeumann){

				for(int i=1;i<=this.numberOfNodes;i++){

					for(int j=0;j<this.dim;j++)
						if(model.BCtype[2*j]==0 || model.PBCpair[2*j]==-2) {
							if(model.node[i].onBound[2*j] || model.node[i].onBound[2*j+1])
								if(model.node[i].hasF())
									for(int k=0;k<this.dim;k++){
										if(k==j)
											model.node[i].F.el[k]=0;

									}

						}
				}
			}
			//================================

			//======== for periodic boundary condition
			

			if(model.hasPBC && model.coordCode==1 ){
				Mat R1=new Mat();
			
				if(model.dim==2)
					 R1=util.rotMat2D(model.alpha2-model.alpha1);
					 else
						 R1=util.rotEuler(new Vect(0,0,1),model.alpha2-model.alpha1);
				

				Mat	R2=R1.transp();
				
				
				for(int i=1;i<=this.numberOfNodes;i++){
					if (! model.node[i].hasF()) continue;

					
					if(model.node[i].hasPBC()){
						int nmap=model.node[i].getMap();
						if (! model.node[nmap].hasF()) continue;
						Vect temp=model.node[i].F.deepCopy();
						model.node[i].F=model.node[i].F.add(R1.mul(model.node[nmap].F)).times(.5);
						model.node[nmap].F=model.node[nmap].F.add(R2.mul(temp)).times(.5);
					}
				}
			}
			
			
			
			if(model.hasPBC && model.coordCode==0 ){
				
				
				for(int i=1;i<=this.numberOfNodes;i++){
					if (! model.node[i].hasF()) continue;

					
					if(model.node[i].hasPBC()){
						int nmap=model.node[i].getMap();
						if (! model.node[nmap].hasF()) continue;
						Vect temp=model.node[i].F.deepCopy();
						model.node[i].F=model.node[i].F.add(model.node[nmap].F).times(.5);
						model.node[nmap].F=model.node[nmap].F.add(temp).times(.5);
					}
				}
			}

				
			
			//=========== converting force components to radial/ tangential components

			if(model.coordCode==1){
				Mat R=new Mat();

				for(int i=1;i<=this.numberOfNodes;i++){
					if(model.dim==2)
					 R=util.rotMat2D(-util.getAng(model.node[i].getCoord()));
					else 
					 R=util.rotEuler(new Vect(0,0,1),-util.getAng(model.node[i].getCoord().v2().v3()));
			
				
					if(model.node[i].hasF()){
						model.node[i].F=R.mul(model.node[i].F);
					}
					}
			
			
					
			}
			//===================
			
			
		
			double frmax=0;
			int im=0;

			for(int i=1;i<=this.numberOfNodes;i++){

				if(model.node[i].hasF()){

					double frn=model.node[i].F.norm();
					if(frn>frmax) 
					{
						frmax=frn;
						im=i;
					}
				}
			}

			
			
			model.FreluctMax=frmax;	

			util.pr("Fmax "+frmax +"at node "+im+" with coordinates" );
			if(im>0){
				model.node[im].getCoord().hshow();
				util.pr("Force vector :" );
				model.node[im].F.hshow();
			}
		}

		



		private Vect[] nodalForce(Model model,int nBH,int i){

			if(this.elCode==0) return nodalForce3ang(model,nBH,i);
			if(this.elCode==3) return nodalForcePrism(model,nBH,i);

			double[][] PW;
		
			if(model.element[i].isNonlin())
				PW=this.PWNL;
			else
				PW=this.PW;

			Node[] vertexNode=model.elementNodes(i);
			Vect[] F=new Vect[this.nElVert];
			for(int j=0;j<this.nElVert;j++)
				F[j]=new Vect(this.dim);
			

			Vect[] gradN=new Vect[this.nElVert];

			Mat T=new Mat(this.dim,this.dim);

			if(this.centerMST)
				T=maxwellTensor(model,nBH,i);

			Mat jac=new Mat(this.dim,this.dim);;
			double ws,detJ,wsJ;
			Vect localCo=new Vect(this.dim);

			int n=PW[0].length; 

			if(this.dim==2){
				for(int p=0;p<n;p++)
					for(int q=0;q<n;q++){

						localCo.el[0]=PW[0][p];
						localCo.el[1]=PW[0][q];

						jac=this.femCalc.jacobian(vertexNode,localCo);
						if(!this.centerMST)
							T=maxwellTensor(model,nBH,i,localCo);
						detJ=abs(jac.determinant());
						gradN=this.femCalc.gradN(jac,localCo);

						ws=1;
						if(n!=2)
							ws=PW[1][p]*PW[1][q];

						wsJ=-ws*detJ;
						for(int j=0;j<this.nElVert;j++){
							F[j]=F[j].add(T.mul(gradN[j]).times(wsJ));
						}
					}

			}
			else
			{

				for(int p=0;p<n;p++)
					for(int q=0;q<n;q++)
						for(int r=0;r<n;r++){

							localCo.el[0]=PW[0][p];
							localCo.el[1]=PW[0][q];
							localCo.el[2]=PW[0][r];

							jac=this.femCalc.jacobian(vertexNode,localCo);
							if(!this.centerMST)
								T=maxwellTensor(model,nBH,i,localCo);
							detJ=abs(jac.determinant());
							gradN=this.femCalc.gradN(jac,localCo);

							ws=1;
							if(n!=2)
								ws=PW[1][p]*PW[1][q]*PW[1][r];

							wsJ=-ws*detJ;

							for(int j=0;j<this.nElVert;j++){
								F[j]=F[j].add(T.mul(gradN[j]).times(wsJ));
							}
						}



			}

			return F;

		}

		private Vect[] nodalForcePrismLLL(Model model,int nBH,int ie){
			
			
			Node[] vertexNode=model.elementNodes(ie);
			
			Vect[] F=new Vect[6];
			for(int j=0;j<6;j++)
				F[j]=new Vect(this.dim);

			Mat T=maxwellTensor(model,nBH,ie);

			Vect[] gradN;

			int nr=this.PW[0].length;
			int n=this.PW3ang.length;
			Vect lc=new Vect(3);
			Mat jac;
			double detJac,wsJ;
			for(int r=0;r<nr;r++)
				for(int p=0;p<n;p++)
				{

					lc.el[0]=this.PW3ang[p][0];
					lc.el[1]=this.PW3ang[p][1];				
					lc.el[2]=this.PW[0][r];

					jac=femCalc.jacobianPrism(vertexNode,lc);
					
					detJac=abs(jac.determinant());
					wsJ=this.PW3ang[p][2]*this.PW[1][r]*detJac;

				
					gradN=femCalc.gradN(jac,lc);
		
					for(int j=0;j<this.nElVert;j++){
						F[j]=F[j].add(T.mul(gradN[j]).times(-wsJ));
					}
				}




			return F;

		}

		private Vect[] nodalForcePrism(Model model,int nBH,int i){

			
			double[][] PW;
		
			
			if(model.element[i].isNonlin())
				PW=this.PWNL;
			else
				PW=this.PW;
			
			int n3ang=this.PW3ang.length; 	
			int nGauss=PW[0].length; 

			Node[] vertexNode=model.elementNodes(i);
			Vect[] F=new Vect[this.nElVert];
		
			for(int j=0;j<this.nElVert;j++)
				F[j]=new Vect(this.dim);
			

			Vect[] gradN=new Vect[this.nElVert];

			Mat T=new Mat(this.dim,this.dim);

			if(this.centerMST)
				T=maxwellTensor(model,nBH,i);

			Mat jac=new Mat(this.dim,this.dim);;
			double ws=0,detJ,wsJ;
			Vect localCo=new Vect(this.dim);
			for(int p=0;p<n3ang;p++)
				for(int q=0;q<nGauss;q++)
				{
		
							localCo.el[0]=this.PW3ang[p][0];
							localCo.el[1]=this.PW3ang[p][1];
							localCo.el[2]=PW[0][q];
							
							jac=this.femCalc.jacobianPrism(vertexNode,localCo);
							
							detJ=abs(jac.determinant());
							
							gradN=this.femCalc.gradN(jac,localCo);

							
							if(!this.centerMST)
								T=maxwellTensor(model,nBH,i,localCo);

							if(nGauss!=2)
								ws=PW[1][q]*this.PW3ang[p][2];
							else
								ws=this.PW3ang[p][2];
							
							wsJ=-ws*detJ;
							
							for(int j=0;j<this.nElVert;j++){
								F[j]=F[j].add(T.mul(gradN[j]).times(wsJ));
							}
					
				}

	
			return F;

		}
		
		private Vect[] nodalForce3ang(Model model,int nBH,int i){
			Vect[] F=new Vect[3];
			for(int j=0;j<3;j++)
				F[j]=new Vect(this.dim);
			Vect[] gradN=this.femCalc.gradN3ang(model,i);

			Mat T=maxwellTensor(model,nBH,i);

			double S=model.getElementArea(i);


			for(int j=0;j<this.nElVert;j++)
				F[j]=T.mul(gradN[j]).times(-S);


			return F;

		}

		private Vect[] surfForce(Model model,int mode, int nCurve,int i){

			if(model.elCode==0)
				return surfForce3ang(model, mode,nCurve,i);
			else {
				Vect[] F=new Vect[model.nElVert];
				for(int j=0;j<model.nElVert;j++)
					F[j]=new Vect(this.dim);
				return F;
			}
		}

		private Vect[] surfForce3ang(Model model, int mode, int nCurve,int i){
			Vect[] F=new Vect[3];
			for(int j=0;j<3;j++)
				F[j]=new Vect(this.dim);


			Vect ec=model.getElementCenter(i);

			Mat T=new Mat(dim,dim);
	
				T=maxwellTensor(model,nCurve,i);


			int[] ednumb=model.element[i].edgeXYNumb;

			for(int j=0;j<3;j++){
				
				/*int n0=model.edge[edgeNumb[j]].endNodeNumber[0];
				int n1=model.edge[edgeNumb[j]].endNodeNumber[1];*/
				
				int n1=model.edge[ednumb[j]].node[0].id;
				int n2=model.edge[ednumb[j]].node[1].id;
				double length=model.edge[ednumb[j]].length;

				Vect v1=model.node[n1].getCoord();
				Vect v2=model.node[n2].getCoord();

				Vect v3=v2.sub(v1);
				v3=v3.normalized();
				Vect v4=v2.sub(ec);
				double proj=v4.dot(v3);
				Vect vn=v4.sub(v3.times(proj)).normalized();
				length=1;
				F[j]=T.mul(vn.times(-length));
			}

			return F;

		}



		

		
		
		

		public Mat maxwellTensor(Model model,int nBH,int ie){
			return 	maxwellTensor(model,nBH,ie,new Vect(this.dim));
		}

		public Mat maxwellTensor(Model model,int nBH,int ie, Vect lc){

			Vect B;
			if(lc.abs().max()==0)
				B=model.element[ie].getB();
			else
				B=this.femCalc.getBAt(model,ie,lc);
			

			double Bn=B.norm();
			Vect H=new Vect();
			double pem;

			Mat T=new Mat(this.dim,this.dim);

			Vect nu=new Vect(this.dim);

			if( model.element[ie].isNonlin()){
				double nux;
				
					nux=model.BH[nBH].getNu(Bn);
				
				for(int k=0;k<this.dim;k++)
					nu.el[k]=nux;

				H=B.times(nu);

				
					pem=model.BH[nBH].getPem(0, H.norm());

	
				T= makeTm(B,H,pem);
				

			}


			else{
				nu=model.element[ie].getNu();

				H=B.times(nu);
				pem=0.5*B.dot(H);
				T= makeTm(B,H,pem);
			}
			
	
			return T;
		}
		
		public Mat makeTm(Vect B, Vect H, double pem){
			Mat T=new Mat(this.dim,this.dim);
			
			for(int i=0;i<this.dim;i++)
				for(int j=0;j<this.dim;j++)
					if(i==j)
						T.el[i][j]=.5*(B.el[i]*H.el[j]+B.el[j]*H.el[i])-pem;

					else
						T.el[i][j]=.5*(B.el[i]*H.el[j]+B.el[j]*H.el[i]);	
	
		
			return T;
		}



		
		
		



	}
