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
	

		public void setMSForce(Model model){
			boolean sf=false;
			Vect[] nodalForce;




			for(int ir=1;ir<=this.numberOfRegions;ir++){
				if(!model.region[ir].MS && !model.region[ir].thermal) continue;
		

				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){	
					if(!model.element[i].hasMS() && !model.element[i].isThermal()) continue;
					int[] vertNumb=model.element[i].getVertNumb();

					for(int j=0;j<this.nElVert;j++){
						model.node[vertNumb[j]].setHasFms(true);
						model.node[vertNumb[j]].Fms=new Vect(this.dim);		
					}
				}
			}

			int nLam;
			for(int ir=1;ir<=this.numberOfRegions;ir++){
	
				if(!model.region[ir].MS && !model.region[ir].thermal) continue;
				nLam=model.region[ir].lamBNumber;
				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){	
					if(!model.element[i].hasMS() && !model.element[i].isThermal()) continue;

					int[] vertNumb=model.element[i].getVertNumb();

					if(sf){/*

						Vect[] fs=surfForce(model,1,nLam,i);

						int[] edgeNumb=model.element[i].edgeXYNumb;

						for(int j=0;j<this.nElEdge;j++){

							int n0=model.edge[edgeNumb[j]].endNodeNumber[0];
							int n1=model.edge[edgeNumb[j]].endNodeNumber[1];

							if(model.node[n0].hasFms())
								model.node[n0].Fms=model.node[n0].Fms.add(fs[j].times(0.5));

							if(model.node[n1].hasFms())
								model.node[n1].Fms=model.node[n1].Fms.add(fs[j].times(0.5));

						}


					*/}

					else {
						nodalForce=nodalMSForce(model,nLam,i);
						

						for(int j=0;j<this.nElVert;j++){

							int nn=vertNumb[j];
							if(model.node[nn].hasFms())
							{	
								model.node[nn].Fms= model.node[nn].Fms.add(nodalForce[j]);

							}

						}
					}

				}
			}
			//======== for Neumann boundary condition
			boolean hasNeumann=false;

			for(int j=0;j<2*this.dim;j++)
				if(model.BCtype[j]==0) {hasNeumann=true; break;}
			if(hasNeumann){

				for(int i=1;i<=this.numberOfNodes;i++){
					Vect v=model.node[i].getCoord();

					for(int j=0;j<this.dim;j++){

						if(model.BCtype[2*j]==0) {
							if(v.el[j]==model.spaceBoundary[2*j] || v.el[j]==model.spaceBoundary[2*j+1])
								if(model.node[i].Fms!=null) 
									for(int k=0;k<this.dim;k++){
										if(k==j)
											model.node[i].Fms.el[k]=0;

									}

						}
					}
				}
			}
			//================================

			//======== for periodic boundary condition
			if(model.hasPBC &&  model.coordCode==1){

				Mat R1;
				if(model.dim==2)
				 R1=util.rotMat2D(model.alpha2-model.alpha1);
				else  R1=util.rotEuler(new Vect(0,0,1),model.alpha2-model.alpha1);
				
				Mat R2=R1.transp();
				for(int i=1;i<=this.numberOfNodes;i++){
					if (! model.node[i].hasFms()) continue;
					if(model.node[i].hasPBC()){
						int nmap=model.node[i].getMap();

						Vect temp=model.node[i].Fms.deepCopy();
						model.node[i].Fms=model.node[i].Fms.add(R1.mul(model.node[nmap].Fms)).times(.5);
						model.node[nmap].Fms=model.node[nmap].Fms.add(R2.mul(temp)).times(.5);
					}
				}
				}
			
			if(model.hasPBC &&  model.coordCode==0){
	
				for(int i=1;i<=this.numberOfNodes;i++){
					if (! model.node[i].hasFms()) continue;
					if(model.node[i].hasPBC()){
						int nmap=model.node[i].getMap();
						Vect temp=model.node[i].Fms.deepCopy();
						model.node[i].Fms=model.node[i].Fms.add(model.node[nmap].Fms);
						model.node[nmap].Fms=model.node[nmap].Fms.add(temp);
					}
				}
				}
			
			//===========

			if(model.coordCode==1){
				Mat R=new Mat();
				
				for(int i=1;i<=this.numberOfNodes;i++){
					if(model.dim==2)
					 R=util.rotMat2D(-util.getAng(model.node[i].getCoord()));
					else R=util.rotEuler(new Vect(0,0,1),-util.getAng(model.node[i].getCoord()));
				
				
					if(model.node[i].hasFms()){
						model.node[i].Fms=R.mul(model.node[i].Fms);
					}
					}
			}
			//===================

			double frmax=0;
			int im=0;
			for(int i=1;i<=this.numberOfNodes;i++){
				if(model.node[i].hasFms()){
					double frn=model.node[i].Fms.norm();
					if(frn>frmax) 
					{
						frmax=frn;
						im=i;
					}
				}
			}

			model.FmsMax=frmax;	

			util.pr("FMsmax "+frmax +"at node "+im+" with coordinates" );
			if(im>0){
				model.node[im].getCoord().hshow();
				util.pr("Force vector :" );
				model.node[im].Fms.hshow();
			}

		}
		
		

		public void setThermalForce(Model model){

			Vect[] nodalForce;
			for(int ir=1;ir<=this.numberOfRegions;ir++){
				if(!model.region[ir].thermal) continue;
				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){	
					if(!model.element[i].isThermal()) continue;
					int[] vertNumb=model.element[i].getVertNumb();

					for(int j=0;j<this.nElVert;j++){
						model.node[vertNumb[j]].setHasFms(true);
						model.node[vertNumb[j]].Fms=new Vect(this.dim);		
					}
				}
			}

			for(int ir=1;ir<=this.numberOfRegions;ir++){
				if(!model.region[ir].thermal) continue;
				for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){	
					if(!model.element[i].isThermal()) continue;

					int[] vertNumb=model.element[i].getVertNumb();

						nodalForce=nodalThermalForce(model,i);


						for(int j=0;j<this.nElVert;j++){

							int nn=vertNumb[j];
							if(model.node[nn].hasFms())
							{	
								model.node[nn].Fms= model.node[nn].Fms.add(nodalForce[j]);

							}

											}

				}
			}
		
			//================================

			//======== for periodic boundary condition
			if(model.hasPBC && model.alpha2-model.alpha1<4){
				if(model.dim==2){
				Mat R1=util.rotMat2D(model.alpha2-model.alpha1);
			
				Mat R2=R1.transp();
				for(int i=1;i<=this.numberOfNodes;i++){
					if (! model.node[i].hasFms()) continue;
					if(model.node[i].hasPBC()){
						int nmap=model.node[i].getMap();
						Vect temp=model.node[i].Fms.deepCopy();
						model.node[i].Fms=model.node[i].Fms.add(R1.mul(model.node[nmap].Fms));
						model.node[nmap].Fms=model.node[nmap].Fms.add(R2.mul(temp));
					}
				}
				}
				else{
					Mat R2D=util.rotMat2D(model.alpha2-model.alpha1);
					Mat R1=new Mat(dim,dim);
					for(int m=0;m<2;m++)
						for(int n=0;n<2;n++)
							R1.el[m][n]=R2D.el[m][n];
					
					R1.el[2][2]=1;
					Mat R2=R1.transp();
					
				
				for(int i=1;i<=this.numberOfNodes;i++){
					if (! model.node[i].hasFms()) continue;
					if(model.node[i].hasPBC()){
						int nmap=model.node[i].getMap();
						Vect temp=model.node[i].Fms.deepCopy();
					
						model.node[i].Fms=model.node[i].Fms.add(R1.mul(model.node[nmap].Fms));
						model.node[nmap].Fms=model.node[nmap].Fms.add(R2.mul(temp));
					}
				}
				}
			}
			
			//===========

			if(model.coordCode==1){
				Mat R=new Mat();
				
				for(int i=1;i<=this.numberOfNodes;i++){
					if(model.dim==2)
					 R=util.rotMat2D(-util.getAng(model.node[i].getCoord()));
					 else {
						Mat R2D=util.rotMat2D(-util.getAng(model.node[i].getCoord().v2()));
						R=new Mat(dim,dim);
						for(int m=0;m<2;m++)
							for(int n=0;n<2;n++)
								R.el[m][n]=R2D.el[m][n];
						
						R.el[2][2]=1;
						 
						//R.show();
					 }
				
					if(model.node[i].hasFms()){
						model.node[i].Fms=R.mul(model.node[i].Fms);
					}
					}
			}
			//===================

			double frmax=0;
			int im=0;
			for(int i=1;i<=this.numberOfNodes;i++){
				if(model.node[i].hasFms()){
					double frn=model.node[i].Fms.norm();
					if(frn>frmax) 
					{
						frmax=frn;
						im=i;
					}
				}
			}

			model.FmsMax=frmax;	

			util.pr("FMsmax "+frmax +"at node "+im+" with coordinates" );
			if(im>0){
				model.node[im].getCoord().hshow();
				util.pr("Force vector :" );
				model.node[im].Fms.hshow();
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
			if(mode==0)
				T=maxwellTensor(model,nCurve,i);
			else
				T=Tms(model,nCurve,i);

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

		private Vect[] nodalMSForce(Model model,int nLam,int i){

			if(this.elCode==0) return nodalMSForce3ang(model,nLam,i);
			if(this.elCode==3) return nodalMSForcePrism(model,nLam,i);
			Mat Tms=new Mat(this.dim,this.dim);
			double[][] PW=this.femCalc.gaussInteg(2);
			Node[] vertexNode=model.elementNodes(i);
			Vect[] F=new Vect[this.nElVert];
			for(int j=0;j<this.nElVert;j++)
				F[j]=new Vect(this.dim);

			Vect[] gradN=new Vect[this.nElVert];

			Mat T=new Mat(this.dim,this.dim);


			Mat jac=new Mat(this.dim,this.dim);;
			double ws,wsJ,detJ;
			Vect localCo=new Vect(this.dim);

			int n=PW[0].length; 
			if(this.dim==2){
				for(int p=0;p<n;p++)
					for(int q=0;q<n;q++){

						localCo.el[0]=PW[0][p];
						localCo.el[1]=PW[0][q];

						jac=this.femCalc.jacobian(vertexNode,localCo);
						T=Tms(model,nLam,i,localCo);
						detJ=abs(jac.determinant());
						gradN=this.femCalc.gradN(jac,localCo);

						ws=1;
						if(n!=2) 
							ws=PW[1][p]*PW[1][q];

						wsJ=-ws*detJ;

						Tms=Tms.add(T.times(ws));

						for(int j=0;j<this.nElVert;j++){
							F[j]=F[j].add(T.mul(gradN[j]).times(wsJ));
						}
					}
			}
			else{

				for(int p=0;p<n;p++)
					for(int q=0;q<n;q++)
						for(int r=0;r<n;r++){

							localCo.el[0]=PW[0][p];
							localCo.el[1]=PW[0][q];
							localCo.el[2]=PW[0][r];

							jac=this.femCalc.jacobian(vertexNode,localCo);
							T=Tms(model,nLam,i,localCo);
							detJ=abs(jac.determinant());
							gradN=this.femCalc.gradN(jac,localCo);

							ws=1;
							if(n!=2) 
								ws=PW[1][p]*PW[1][q]*PW[1][r];

							wsJ=-ws*detJ;

							Tms=Tms.add(T.times(ws));

							for(int j=0;j<this.nElVert;j++){
								F[j]=F[j].add(T.mul(gradN[j]).times(wsJ));
							}
						}

			}

			double rv=1.0/pow(2,this.dim);
			Vect MSS=util.vectorize(Tms.times(rv));

			model.element[i].setStress(MSS.times(1e-6));


			return F;

		}


		private Vect[] nodalThermalForce(Model model,int i){

			if(this.elCode==0) return nodalThermalForce3ang(model,i);
			if(this.elCode==3) return nodalThermalForcePrism(model,i);
			Mat Tth=new Mat(this.dim,this.dim);
			double[][] PW=this.femCalc.gaussInteg(2);
			Node[] vertexNode=model.elementNodes(i);
			Vect[] F=new Vect[this.nElVert];
			for(int j=0;j<this.nElVert;j++)
				F[j]=new Vect(this.dim);

			Vect[] gradN=new Vect[this.nElVert];

			Mat T=new Mat(this.dim,this.dim);


			Mat jac=new Mat(this.dim,this.dim);;
			double ws,wsJ,detJ;
			Vect localCo=new Vect(this.dim);

			int n=PW[0].length; 
			if(this.dim==2){
				for(int p=0;p<n;p++)
					for(int q=0;q<n;q++){

						localCo.el[0]=PW[0][p];
						localCo.el[1]=PW[0][q];

						jac=this.femCalc.jacobian(vertexNode,localCo);
						T=Tth(model,i);
						detJ=abs(jac.determinant());
						gradN=this.femCalc.gradN(jac,localCo);

						ws=1;
						if(n!=2) 
							ws=PW[1][p]*PW[1][q];

						wsJ=-ws*detJ;

						Tth=Tth.add(T.times(ws));

						for(int j=0;j<this.nElVert;j++){
							F[j]=F[j].add(T.mul(gradN[j]).times(wsJ));
						}
					}
			}
			else{

				for(int p=0;p<n;p++)
					for(int q=0;q<n;q++)
						for(int r=0;r<n;r++){

							localCo.el[0]=PW[0][p];
							localCo.el[1]=PW[0][q];
							localCo.el[2]=PW[0][r];

							jac=this.femCalc.jacobian(vertexNode,localCo);
							T=Tth(model,i);
							detJ=abs(jac.determinant());
							gradN=this.femCalc.gradN(jac,localCo);

							ws=1;
							if(n!=2) 
								ws=PW[1][p]*PW[1][q]*PW[1][r];

							wsJ=-ws*detJ;

							Tth=Tth.add(T.times(ws));

							for(int j=0;j<this.nElVert;j++){
								F[j]=F[j].add(T.mul(gradN[j]).times(wsJ));
							}
						}

			}

			double rv=1.0/pow(2,this.dim);
			Vect MSS=util.vectorize(Tth.times(rv));

			model.element[i].setStress(MSS.times(1e-6));


			return F;

		}
		

		private Vect[] nodalMSForce3ang(Model model,int nLam,int i){
			Vect[] F=new Vect[3];
			for(int j=0;j<3;j++)
				F[j]=new Vect(this.dim);
			Vect[] gradN=this.femCalc.gradN3ang(model,i);

			Mat T=Tms(model,nLam,i);

			double S=this.femCalc.el3angArea(model,i);

			for(int j=0;j<this.nElVert;j++)
				F[j]=T.mul(gradN[j]).times(-S);

			Vect MSS=util.vectorize(T);

			
			model.element[i].setStress(MSS.times(1e-6));
		

			return F;

		}
		
		private Vect[] nodalThermalForce3ang(Model model,int i){
			Vect[] F=new Vect[3];
			for(int j=0;j<3;j++)
				F[j]=new Vect(this.dim);
			Vect[] gradN=this.femCalc.gradN3ang(model,i);

			Mat T=Tth(model,i);

			double S=this.femCalc.el3angArea(model,i);

			for(int j=0;j<this.nElVert;j++)
				F[j]=T.mul(gradN[j]).times(-S);

			Vect MSS=util.vectorize(T);

			
			model.element[i].setStress(MSS.times(1e-6));


			return F;

		}
		
		private Vect[] nodalMSForcePrism(Model model,int nLam,int ie){
			
			
			Node[] vertexNode=model.elementNodes(ie);
			
			Vect[] F=new Vect[6];
			for(int j=0;j<6;j++)
				F[j]=new Vect(this.dim);

			Mat T=Tms(model,nLam,ie);

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


			Vect MSS=util.vectorize(T);

			
			model.element[ie].setStress(MSS.times(1e-6));


			return F;

		}
		
private Vect[] nodalThermalForcePrism(Model model,int ie){
			

			
			Node[] vertexNode=model.elementNodes(ie);
			
			Vect[] F=new Vect[6];
			for(int j=0;j<6;j++)
				F[j]=new Vect(this.dim);

			Mat T=Tth(model,ie);

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


			Vect MSS=util.vectorize(T);

			
			model.element[ie].setStress(MSS.times(1e-6));


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

		public Mat Tms(Model model,int nLam,int ie){
			return Tms(model,nLam,ie,new Vect(this.dim));

		}

		public Mat Tms(Model model,int nLam,int ie,Vect lc)
		{
		


			double seq=0;
			
			Vect B;
			if(lc.abs().max()==0)
				B=model.element[ie].getB();
			else
				B=this.femCalc.getBAt(model,ie,lc);
			
		//	model.m2d
		//	Mat S1=model.element[ie].getStressTensor();
			if(1>2 && model.element[ie].getRegion()==8){
		
				Mat S1=model.m2d.element[ie-model.region[8].getFirstEl()+1].getStressTensor();
			
				if(S1!=null){
			
				Mat S2=S1.deepCopy();
			
			double p=0;
			for(int i=0;i<this.dim;i++)
				p+=S2.el[i][i];

			S2.addToDiag(new Vect().ones(dim).times(-p/3));
			
			Vect h=B.normalized();
			
			 seq=1.5*S2.mul(h).dot(h);	
			 
			 
					}
			}
			
	
			double Bn2=B.dot(B);
			
			Mat T=new Mat(this.dim,this.dim);
			
			if(model.element[ie].hasMS()){
			
			if(Bn2<1e-3) return T;
			
			
			double Bn=sqrt(Bn2);
			double E=model.element[ie].getYng().el[0];
			double v=model.element[ie].getPois().el[0];
			double G=E/(1+v);

			Mat S=new Mat(this.dim,this.dim);
			double kms=.5;

			if(model.dim==3){

			for(int i=0;i<this.dim;i++)
				for(int j=0;j<this.dim;j++)
					if(i==j) S.el[i][j]=((1+kms)*B.el[i]*B.el[j]-kms*Bn2)/Bn2;
					else
						S.el[i][j]=((1+kms)*B.el[i]*B.el[j])/Bn2;
			}
			else{
				
				//== plane stress 
				double rr=(v-.5)/(1-v);
				
				for(int i=0;i<this.dim;i++)
					for(int j=0;j<this.dim;j++)
						if(i==j) S.el[i][j]=(1.5*B.el[i]*B.el[j]+rr*Bn2)/Bn2;
						else
							S.el[i][j]=1.5*B.el[i]*B.el[j]/Bn2;
				
			}

			double a;

			
				a=-model.lamB[nLam].getLam(Bn);	

				T=S.times(G*a);


			
			}


			return T;
		}
		
		
		public Mat Tth(Model model,int ie)
		{

			double dT=model.element[ie].getDeltaT();
				double ct=model.region[model.element[ie].getRegion()].getThermalCoef();
				
				double E=model.element[ie].getYng().el[0];
				double v=model.element[ie].getPois().el[0];

				double ss=-dT*ct*E/(1-v);
				if(dim==3) ss=-dT*ct*E/((1+v)*(1-2*v));
		
				 Mat Tth=new Mat(dim,dim);
				 for(int i=0;i<Tth.nCol;i++)
					 Tth.el[i][i]=ss;
				
				 if(dim==3) Tth.el[2][2]=0;		


			return Tth;
		}
		
		

		public void setStressForce(Model model){


			for(int i=1;i<=this.numberOfNodes;i++)
				if(model.node[i].isDeformable())
				model.node[i].Fms=new Vect(model.dim);
				
			for(int i=1;i<=this.numberOfElements;i++){
				if(model.element[i].getStress()==null) continue; 
				int[] vertNumb=model.element[i].getVertNumb();
				
		
				Vect[] nodalForce=this.nodalStressForce(model,i);
				for(int j=0;j<this.nElVert;j++){
					int nn=vertNumb[j];		
					if(model.node[nn].Fms==null)
				model.node[nn].Fms=nodalForce[j].deepCopy();
					else
						model.node[nn].Fms=model.node[nn].Fms.add(nodalForce[j]);
				}
			}
			
				
			}
		
		
		
		private Vect[] nodalStressForce(Model model,int i){


			double[][] PW;
			PW=this.PW;

			Node[] vertexNode=model.elementNodes(i);
			Vect[] F=new Vect[this.nElVert];
			for(int j=0;j<this.nElVert;j++)
				F[j]=new Vect(this.dim);
			

			Vect[] gradN=new Vect[this.nElVert];

			Mat T=model.element[i].getStressTensor().times(1e6);


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
		



	}
