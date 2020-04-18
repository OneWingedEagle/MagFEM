package fem;

import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import javax.swing.JOptionPane;

import math.Mat;
import math.Vect;
import math.util;

public class Calculator {
	double[][] PW;
	double[][] PWNL;
	double[][] PW3ang;
	double[][] PWtetra;
	double[][] PWpyr;
	boolean centerNu;
	public int dim,elCode,nElVert,nElEdge,numberOfElements,numberOfNodes,numberOfRegions;
	public int struc2D=0; // 0:plane stress, 1: plain strain

	public Calculator(){	}

	public Calculator(Model model)
	{
		this.nElVert=model.nElVert;
		this.nElEdge=model.nElEdge;
		this.numberOfElements=model.numberOfElements;
		this.numberOfNodes=model.numberOfNodes;
		this.numberOfRegions=model.numberOfRegions;
		this.dim=model.dim;
		this.elCode=model.elCode;
		this.PW=gaussInteg(2);
		this.PWNL=gaussInteg(3);
		this.PW3ang=gaussInteg3(4);
		this.PWtetra=gaussIntegTetra(5);
		this.PWpyr=gaussIntegPyramid(8);
		
		this.centerNu=false;
	}


	public Vect getElementB(Model model,int i, Vect[] rotNe){

		Edge[] elEdge=model.elementEdges(i);
		Vect B=new Vect(model.dim);
		for(int j=0;j<model.nElEdge;j++)		{
			B=B.add(rotNe[j].times(elEdge[j].A));
		}

		return B;

	}

	public double getElementAQ(Model model,int i){
		Edge[] elEdge=model.elementEdges(i);
		double[] Ae=new double[4];
		for(int j=0;j<4;j++)
			Ae[j]=elEdge[j].A;
		double A=0;
		Vect zero=new Vect(2);

		double[] Ne=NeQuad(zero);

		for(int j=0;j<model.nElEdge;j++)	{		
			A= A+Ne[j]*Ae[j];
		}
		return  A;	

	}


	public Mat nodalMat(Vect nu,Vect gN1, Vect gN2){
		Vect gNnuN2=nu.times(gN2);
		double p=gN1.dot(gNnuN2);
		Mat M=new Mat(this.dim,this.dim);
		for(int i=0;i<this.dim;i++)
			for(int j=0;j<this.dim;j++)
				if(i==j)
					M.el[i][j]=p-gN1.el[j]*gNnuN2.el[j];
				else
					M.el[i][j]=-gN1.el[j]*gNnuN2.el[j];
		return M;

	}

	public double[][] He(Model model,int nBH,int nLam,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){

		if(this.elCode==0 && !model.axiSym) return He3ang(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==0&& model.axiSym) return He3angAxi(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==1&& !model.axiSym) return HeQuad(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==1&& model.axiSym) return HeQuadAxi(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==3) return HePrism(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==5) return HePyramid(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		else 	if(this.elCode==2) {
			 return HeTetra(model,nBH,ie,nonLinear,eddy,hasJ,hasM);
		//	String msg="Tetrahedral edge element not implemented!";
			//JOptionPane.showMessageDialog(null, msg," ", JOptionPane. ERROR_MESSAGE);
		}
		//else if(this.elCode==2) return He3ang2nd(model,nBH,ie,nonLinear,eddy,hasJ,hasM);

		Vect nu=new Vect(this.dim),nuVar=new Vect(this.dim);
		Vect B=new Vect();
		Vect  M=new Vect();		


		if(hasM){ M=model.element[ie].getM();}
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();

		int n;

		if(!nonLinear){
			nu=model.element[ie].getNu();
			n=this.PW[0].length; 
		}
		else{
			n=this.PWNL[0].length; 
			if(this.centerNu){
				double nux,nuv=0;
				B=model.element[ie].getB();
				double Bn2=B.norm2();
				double Bn=sqrt(Bn2);

				nux=model.BH[nBH].getNu(Bn);
				nu=new Vect(nux,nux,nux);
				if(Bn>0){
					nuv=model.BH[nBH].getdHdB(Bn);	
					nuVar.el[0]=(nuv-nux)/Bn2;
					nuVar.el[1]=nuVar.el[0];
					nuVar.el[2]=nuVar.el[0];

				}
			}
		}

		Node[] vertexNode=model.elementNodes(ie);
		
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		Vect[] Cj=new Vect[model.nElEdge];
		for(int i=0;i<model.nElEdge;i++)
			Cj[i]=new Vect(3);

		double detJac,ws=1,wsJ=0,term1,term2,term3;

		Mat jac;

		

		Vect[] rotNe=new Vect[model.nElEdge];


		Vect localCo=new Vect(this.dim);

		Vect[] Ne=new Vect[model.nElEdge];

		Vect sigma=new Vect();

		if(eddy)
			sigma=model.element[ie].getSigma();

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){


					if(nonLinear){
						localCo.el[0]=this.PWNL[0][p];
						localCo.el[1]=this.PWNL[0][q];
						localCo.el[2]=this.PWNL[0][r];

						if(n!=2)
							ws=this.PWNL[1][p]*this.PWNL[1][q]*this.PWNL[1][r];
						else
							ws=1;
						B=this.getBAt(model,ie,localCo);
						double Bn2=B.norm2();
						double Bn=sqrt(Bn2);

						if(!this.centerNu){
							double nux,nuv=0;


							nux=model.BH[nBH].getNu(Bn);
							nu=new Vect(nux,nux,nux);
							if(Bn>0){
								nuv=model.BH[nBH].getdHdB(Bn);	
								nuVar.el[0]=(nuv-nux)/Bn2;
								nuVar.el[1]=nuVar.el[0];
								nuVar.el[2]=nuVar.el[0];

							}
						}

					}
					else
					{
						localCo.el[0]=this.PW[0][p];
						localCo.el[1]=this.PW[0][q];
						localCo.el[2]=this.PW[0][r];

						if(n!=2)
							ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r];
						else
							ws=1;


					}

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());
		

					wsJ=ws*detJac;


					rotNe=rotNe(jac,localCo,edgeDir);
					
	
					
					if(hasJ || eddy)
						Ne=Ne(jac,localCo,edgeDir);
					
					for(int i=0;i<model.nElEdge;i++){
						
						
				if(hasJ){

							Cj[i]=Cj[i].add(Ne[i].times(wsJ));
						}

						if(hasM){
							C[i]+=wsJ*rotNe[i].dot(M);

						}

						for(int j=0;j<=i;j++){

							term1=wsJ*rotNe[i].times(nu).dot(rotNe[j]);
							H1[i][j]+=term1;
						
			
							if(nonLinear){
								term2=wsJ*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
								H2[i][j]+=term2;


							}

							if(eddy){
								term3=wsJ*Ne[i].dot(Ne[j].times(sigma));
								H3[i][j]+=term3;

							}


						}
					}
				}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy){
			lowSym(H3);
		}
		
		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj=Cj;


/*		if(ie==1 ||model.element[ie].getRegion()==-1){
		//	util.pr("elem "+ie);
		//	util.hshow(model.element[ie].getEdgeNumb());
		//	util.show(H1);
		//	new Mat(H1).diagonal().show();
		}
*/
		return H1;


	}


	
	public Vect elemVectPhiCoil(Model model,int ie, double[] loss,double conductivity){
		
		if(model.elCode==1) return elemVectPhiCoilQuad(model,ie,loss,conductivity);
		if(model.elCode==2) return elemVectPhiCoilTet(model,ie,loss,conductivity);
		
		Vect elemV=new Vect(model.nElEdge);
		
		double heat=0;
		int n;

		n=this.PW[0].length; 

		Node[] vertexNode=model.elementNodes(ie);
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();
		
		double detJac,ws=1,wsJ=0;

		Mat jac;

		Vect localCo=new Vect(this.dim);

		Vect[] Ne=new Vect[model.nElEdge];

		Mat G;
		Vect Je=new Vect(this.dim);
		Vect[] localGradsN=new Vect[model.nElVert];;
		Vect[] gradN=new Vect[model.nElVert];
		Vect[] gradPhi=new Vect[model.nElVert];;
	

		double[] nodePhi=new double[this.nElVert];
		for(int j=0;j<model.nElVert;j++){
			nodePhi[j]=vertexNode[j].getPhi();
		}

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){

					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r];
					else
						ws=1;

					jac=jacobian(vertexNode,localCo);

					G=jac.inv3();

					localGradsN=localGradN(localCo);

					detJac=abs(jac.determinant());

					wsJ=ws*detJac;

					Ne=Ne(jac,localCo,edgeDir);
					
					Je.zero();
					for(int j=0;j<model.nElVert;j++){
						gradN[j]=G.mul(localGradsN[j]);
						gradPhi[j]= gradN[j].times(nodePhi[j]);
						Je=Je.add(gradPhi[j]);
					}
	
					for(int i=0;i<model.nElEdge;i++){

												
						elemV.el[i]+=Ne[i].dot(Je)*wsJ;

						}
					
					heat+=Je.dot(Je)*wsJ;
				}

		elemV=elemV.times(-1);
		
		loss[0]=heat;

		return elemV;

	}


	public Vect elemVectPhiCoilTet(Model model,int ie, double[] loss,double conductivity){

		Vect elemV=new Vect(model.nElEdge);
		
		Vect[] Ne=new Vect[model.nElEdge];
		
		double heat=0;
		
		Node[] vertexNode=model.elementNodes(ie);

		boolean[] edgeDir=model.element[ie].getEdgeReverse();
		
		double[] nodePhi=new double[this.nElVert];
		for(int j=0;j<model.nElVert;j++){
			nodePhi[j]=vertexNode[j].getPhi();
		}

		Mat jac=this.jacobianTetra(vertexNode, new Vect(dim));

		Vect[] gradN=this.gradNTet(jac);
		Vect[] gradPhi=new Vect[this.nElVert];
		double detJac=abs(jac.determinant());
		double ws0=detJac/6;
		
		Vect Je=new Vect(3);
		for(int j=0;j<model.nElVert;j++){

			gradPhi[j]= gradN[j].times(nodePhi[j]);
			Je=Je.add(gradPhi[j]);
		}

		heat=Je.dot(Je)*ws0;

		
		Vect lc=new Vect(4);

		double ws;
		
		int n=this.PWtetra.length; 					
		for(int p=0;p<n;p++)
		{

			lc.el[0]=this.PWtetra[p][0];
			lc.el[1]=this.PWtetra[p][1];
			lc.el[2]=this.PWtetra[p][2];
			lc.el[3]=this.PWtetra[p][3];
			
			Ne=NeTetra(jac,lc,edgeDir);

			ws=this.PWtetra[p][4]*detJac;
	
			for(int i=0;i<model.nElEdge;i++){
									
			elemV.el[i]+=Ne[i].dot(Je)*ws;

			}
		}
		
		
		elemV=elemV.times(-1);
		
		loss[0]=heat/conductivity;

		return elemV;
	
	}
	

	public Vect elemVectPhiCoilQuad(Model model,int ie, double[] loss,double conductivity){
		
		Vect elemV=new Vect(model.nElEdge);
		
		double heat=0;

		int n;

		n=this.PW[0].length; 

		Node[] vertexNode=model.elementNodes(ie);
		Edge[] elemEdges=model.elementEdges(ie);


		double detJac,ws=1,wsJ=0;

		Mat jac;

		Vect localCo=new Vect(this.dim);

		double[] Ne=new double[model.nElEdge];

		Mat G;
		double Jez=0;



		double[] nodePhi=new double[this.nElVert];
		for(int j=0;j<model.nElVert;j++){
			nodePhi[j]=vertexNode[j].getPhi();
		}
		
		double jezphi=0;
		for(int j=0;j<model.nElVert;j++){
			jezphi+=nodePhi[j]/model.height;
		}

		
		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	{
			A[j]=elemEdges[j].getA();
		}
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){
			

					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q];
					else
						ws=1;

					jac=jacobian(vertexNode,localCo);
					
					detJac=abs(jac.determinant());

					wsJ=ws*detJac;

					Ne=NeQuad(localCo);
				
					Jez=0;
					for(int j=0;j<model.nElEdge;j++){
						Jez+=Ne[j]*A[j];
					
					}
					
					Jez+=jezphi/(n*n);

			

					for(int i=0;i<model.nElEdge;i++){

												
						elemV.el[i]+=Ne[i]*Jez*wsJ;

						}
					
					heat+=Jez*Jez*wsJ;

				}
				
		
		loss[0]=heat/conductivity;
		
		return elemV;


	}
	
	public Vect elemVectCLN(Model model,int ie, double[] loss){
		
		if(this.elCode==1) return elemVectCLNQuad(model, ie,loss);
		if(this.elCode==2) return elemVectCLNTetra(model, ie,loss);

		Vect elemV=new Vect(model.nElEdge);
		
		double heat=0;

		int n;

		n=this.PW[0].length; 

		Node[] vertexNode=model.elementNodes(ie);
		Edge[] elemEdges=model.elementEdges(ie);
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();
		
		double sigma=model.element[ie].getSigma().el[0];

		double detJac,ws=1,wsJ=0;

		Mat jac;

		Vect localCo=new Vect(this.dim);

		Vect[] Ne=new Vect[model.nElEdge];

		Mat G;
		Vect Ee=new Vect(this.dim);
		Vect[] localGradsN=new Vect[model.nElVert];;
		Vect[] gradN=new Vect[model.nElVert];
		Vect[] gradPhi=new Vect[model.nElVert];;
	

		double[] nodePhi=new double[this.nElVert];
		for(int j=0;j<model.nElVert;j++){
			nodePhi[j]=vertexNode[j].getPhi();
		}
		
		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	{
			A[j]=elemEdges[j].getA();
		}
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){

					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r];
					else
						ws=1;

					jac=jacobian(vertexNode,localCo);

					G=jac.inv3();

					localGradsN=localGradN(localCo);

					detJac=abs(jac.determinant());

					wsJ=ws*detJac;

					Ne=Ne(jac,localCo,edgeDir);

					for(int j=0;j<model.nElVert;j++){
						gradN[j]=G.mul(localGradsN[j]);
						gradPhi[j]= gradN[j].times(nodePhi[j]);
					}

					
					Ee.zero();
					for(int j=0;j<model.nElEdge;j++){
						Ee=Ee.add(Ne[j].times(A[j]));
					
					}

					for(int j=0;j<model.nElVert;j++){

						Ee=Ee.add(gradPhi[j]);
					}

					for(int i=0;i<model.nElEdge;i++){

												
						elemV.el[i]+=Ne[i].dot(Ee)*wsJ;

						}
					
					heat+=Ee.dot(Ee)*wsJ;
	
				}
				
		elemV=elemV.times(sigma);

		
		loss[0]=heat*sigma;

		return elemV;


	}

	public Vect elemVectCLNTetra(Model model,int ie, double[] loss){
		Vect elemV=new Vect(model.nElEdge);
			
			Vect[] Ne=new Vect[model.nElEdge];
			
			double sigma=model.element[ie].getSigma().el[0];

			double heat=0;
			
			Node[] vertexNode=model.elementNodes(ie);

			boolean[] edgeDir=model.element[ie].getEdgeReverse();
			
			double[] nodePhi=new double[this.nElVert];
			for(int j=0;j<model.nElVert;j++){
				nodePhi[j]=vertexNode[j].getPhi();
			}

			Mat jac=this.jacobianTetra(vertexNode, new Vect(dim));

			Vect[] gradN=this.gradNTet(jac);
			Vect[] gradPhi=new Vect[this.nElVert];
			double detJac=abs(jac.determinant());
			double ws0=detJac/6;
			
			Vect Je=new Vect(3);
			for(int j=0;j<model.nElVert;j++){

				gradPhi[j]= gradN[j].times(nodePhi[j]);
				Je=Je.add(gradPhi[j]);
			}

			heat=Je.dot(Je)*ws0;

			
			Vect lc=new Vect(4);

			double ws;
			
			int n=this.PWtetra.length; 					
			for(int p=0;p<n;p++)
			{

				lc.el[0]=this.PWtetra[p][0];
				lc.el[1]=this.PWtetra[p][1];
				lc.el[2]=this.PWtetra[p][2];
				lc.el[3]=this.PWtetra[p][3];
				
				Ne=NeTetra(jac,lc,edgeDir);

				ws=this.PWtetra[p][4]*detJac;
		
				for(int i=0;i<model.nElEdge;i++){
										
				elemV.el[i]+=Ne[i].dot(Je)*ws;

				}
			}
			
			
			elemV=elemV.times(sigma);

			
			loss[0]=heat*sigma;
			return elemV;
		
		}
	

	public Vect elemVectCLNQuad(Model model,int ie, double[] loss){
		
		Vect elemV=new Vect(model.nElEdge);
		
		double heat=0;

		int n;

		n=this.PW[0].length; 

		Node[] vertexNode=model.elementNodes(ie);
		Edge[] elemEdges=model.elementEdges(ie);
		Vect[] localGradsN=new Vect[model.nElVert];;

		double sigma=model.element[ie].getSigma().el[0];

		double detJac,ws=1,wsJ=0;

		Mat jac;

		Vect localCo=new Vect(this.dim);

		double[] Ne=new double[model.nElEdge];

		Mat G;
		double Jez=0;



		double[] nodePhi=new double[this.nElVert];
		for(int j=0;j<model.nElVert;j++){
			nodePhi[j]=vertexNode[j].getPhi();
		}
		
		double jezphi=0;
		for(int j=0;j<model.nElVert;j++){
			jezphi+=nodePhi[j]/model.height;
		}

		
		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	{
			A[j]=elemEdges[j].getA();
		}
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){
			

					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q];
					else
						ws=1;

					jac=jacobian(vertexNode,localCo);
					
					detJac=abs(jac.determinant());

					wsJ=ws*detJac;

					Ne=NeQuad(localCo);
				
					Jez=0;
					for(int j=0;j<model.nElEdge;j++){
						Jez+=Ne[j]*A[j];
					
					}
					
					Jez+=jezphi/(n*n);

			

					for(int i=0;i<model.nElEdge;i++){

												
						elemV.el[i]+=Ne[i]*Jez*wsJ;

						}
					
					heat+=Jez*Jez*wsJ;

				}
				
		elemV=elemV.times(sigma);

		
		loss[0]=heat*sigma;
		return elemV;


	}

	public double[][] HeQuad(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){


		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect M=new Vect();		

		boolean[] edgeDir=model.element[ie].getEdgeReverse();



		if(hasM){ M=model.element[ie].getM();}


		if(!nonLinear)
			nu=model.element[ie].getNu();
		else if(this.centerNu){
			B=model.element[ie].getB();
			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=model.BH[nBH].getNu(Bn);
			nu=new Vect(nux,nux);
			if(Bn>0){
				double nuv=model.BH[nBH].getdHdB(Bn);	
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}

		}

		double sigma=0;
		if(eddy){
			sigma=model.element[ie].getSigma().el[2];			


		}

		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		double[] Cj=new double[model.nElEdge];

		double detJac,ws=1,wsJ=0,term1,term2,term3;

		Mat jac;

		Vect[] rotNe=new Vect[model.nElEdge];
		double [] Ne=new double[model.nElEdge];

		int n=0;
		Vect localCo=new Vect(2);
		if(nonLinear)
			n=this.PWNL[0].length; 
		else
			n=this.PW[0].length; 

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){

				if(nonLinear){
					localCo.el[0]=this.PWNL[0][p];
					localCo.el[1]=this.PWNL[0][q];

					if(n!=2)
						ws=this.PWNL[1][p]*this.PWNL[1][q];
					else
						ws=1;
					B=this.getBAt(model,ie,localCo);
					double Bn2=B.norm2();
					double Bn=sqrt(Bn2);

					if(!this.centerNu){

						double nux=model.BH[nBH].getNu(Bn);
						nu=new Vect(nux,nux);
						if(Bn>0){
							double nuv=model.BH[nBH].getdHdB(Bn);	
							nuVar.el[0]=(nuv-nux)/Bn2;
							nuVar.el[1]=nuVar.el[0];


						}
					}


				}
				else
				{
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q];
					else
						ws=1;

				}

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				wsJ=ws*detJac;

				rotNe=rotNe(jac,localCo,edgeDir);
				if(hasJ || eddy)
					Ne=NeQuad(localCo);

				for(int i=0;i<model.nElEdge;i++){

					if(hasJ){
						Cj[i]+=wsJ*Ne[i];
					}

					if(hasM)
						C[i]+=wsJ*rotNe[i].dot(M);

					for(int j=0;j<=i;j++){

						term1=wsJ*rotNe[i].times(nu).dot(rotNe[j]);
						H1[i][j]+=term1;


						if(nonLinear){
							term2=wsJ*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
							H2[i][j]+=term2;

						}

						if(eddy){

							term3=wsJ*Ne[i]*Ne[j]*sigma;
							H3[i][j]+=term3;

						}
					}
				}
			}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);


		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;

		return H1;

	}

	public double[][] HeQuadAxi(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){


		boolean centerNu=true;

		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect  M=new Vect();		


		if(hasM){ M=model.element[ie].getM();}


		if(!nonLinear)
			nu=model.element[ie].getNu();
		else if(centerNu){
			B=model.element[ie].getB();
			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=model.BH[nBH].getNu(Bn);
			nu=new Vect(nux,nux);
			if(Bn>0){
				double nuv=model.BH[nBH].getdHdB(Bn);	
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}

		}

		double sigma=0;
		if(eddy){
			sigma=model.element[ie].getSigma().el[2];			


		}


		double rr=0;
		boolean rcent=true;

		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		double[] Cj=new double[model.nElEdge];

		double[] rrj=new double[model.nElEdge];
		for(int i=0;i<model.nElEdge;i++){
			double r1=vertexNode[i].getCoord(0)+1e-5;
			rrj[i]=1.0/r1;
		}


		double detJac,ws=1,wsJ=0,term1,term2,term3;

		Mat jac;

		Vect[] rotNe=new Vect[model.nElEdge];
		double [] Ne=new double[model.nElEdge];

		int n=0;
		Vect localCo=new Vect(2);
		if(nonLinear)
			n=this.PWNL[0].length; 
		else
			n=this.PW[0].length; 

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){

				if(nonLinear){
					localCo.el[0]=this.PWNL[0][p];
					localCo.el[1]=this.PWNL[0][q];

					if(n!=2)
						ws=this.PWNL[1][p]*this.PWNL[1][q];
					else
						ws=1;
					B=this.getBAt(model,ie,localCo);

					double Bn2=B.norm2();

					double Bn=sqrt(Bn2);

					if(!centerNu){

						double nux=model.BH[nBH].getNu(Bn);
						nu=new Vect(nux,nux);
						if(Bn>0){
							double nuv=model.BH[nBH].getdHdB(Bn);	
							nuVar.el[0]=(nuv-nux)/Bn2;
							nuVar.el[1]=nuVar.el[0];


						}
					}


				}
				else
				{
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q];
					else
						ws=1;

				}

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				wsJ=ws*detJac;

				rotNe=gradN(jac,localCo);

				Ne=NQuad(localCo);


				rr=0;
				if(!rcent){

					for(int i=0;i<model.nElEdge;i++)
						rr+=Ne[i]*rrj[i];
				}
				else{

					for(int i=0;i<model.nElEdge;i++)
						rr+=.25*rrj[i];
				}


				for(int i=0;i<model.nElEdge;i++){

					if(hasJ){
						Cj[i]+=wsJ*Ne[i];
					}

					if(hasM)
						C[i]+=wsJ*rotNe[i].dot(M);




					for(int j=0;j<4;j++){

						term1=wsJ*rotNe[i].times(nu).dot(rotNe[j]);
						H1[i][j]+=term1;


						if(nonLinear){
							term2=wsJ*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
							H2[i][j]+=term2;
						}

						if(eddy){
							term3=wsJ*Ne[i]*Ne[j]*sigma;
							H3[i][j]+=term3;

						}
					}
				}
			}


		for(int i=0;i<model.nElEdge;i++)
			for(int j=0;j<model.nElEdge;j++){
				H1[i][j]*=rr;
				if(nonLinear){
					H2[i][j]*=rr;
				}
				if(eddy)
					H3[i][j]*=rr;

			}


		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;


		return H1;

	}

	public double[][] HeTetra(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){
		Vect nu=new Vect(3),
		nuVar=new Vect(3);
		Vect M=new Vect();		
		Vect B=model.element[ie].getB();
		Vect lc=new Vect().ones(4).times(.25);

		boolean[] edgeDir=model.element[ie].getEdgeReverse();


		if(hasM){ M=model.element[ie].getM();}

		Vect sigma=null;
		if(eddy)
			sigma=model.element[ie].getSigma();

		if(nonLinear){

			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=0,nuv;

			nux=model.BH[nBH].getNu(Bn);

			nu=new Vect(nux,nux,nux);
			if(Bn>0){
				nuv=model.BH[nBH].getdHdB(Bn);
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
				nuVar.el[2]=nuVar.el[0];
			}
		}

		else
			nu=model.element[ie].getNu();


		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		
		Vect[] Cj=new Vect[model.nElEdge];
		for(int i=0;i<model.nElEdge;i++)
			Cj[i]=new Vect(3);

	//	double vol=model.getElementVolumeTetra(ie);
		double detJac,ws=1;
		Mat jac;
		Node[] vertexNode=model.elementNodes(ie);

		Vect[] Ne=new Vect[model.nElEdge];
		jac=jacobianTetra(vertexNode,lc);
		detJac=abs(jac.determinant());
		
		double vol=detJac/6;

		Vect[] rotNe=rotNeTetra(jac,lc,edgeDir);
		


		if(hasM){
			for(int i=0;i<model.nElEdge;i++)
				C[i]=rotNe[i].dot(M)*vol;
		}
		if(hasJ){
			Ne=NeTetra(jac,lc,edgeDir);
		
			for(int i=0;i<model.nElEdge;i++)
			Cj[i]=Ne[i].times(vol);
		}

		for(int i=0;i<model.nElEdge;i++)
			for(int j=0;j<=i;j++)	{
	
				H1[i][j]=vol*rotNe[i].times(nu).dot(rotNe[j]);

				if(nonLinear)
					H2[i][j]=vol*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));

				if(eddy){
					int n=this.PWtetra.length; 					
					for(int p=0;p<n;p++)
					{

						lc.el[0]=this.PWtetra[p][0];
						lc.el[1]=this.PWtetra[p][1];
						lc.el[2]=this.PWtetra[p][2];
						lc.el[3]=this.PWtetra[p][3];
		
						Ne=NeTetra(jac,lc,edgeDir);
						ws=this.PWtetra[p][4]*detJac;
				
						H3[i][j]+=ws*Ne[i].dot(Ne[j].times(sigma));
								
				
					}
				}


			}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj=Cj;

		return H1;
	}
	
	public double[][] HePrism(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){

		//return new double[1][1];

		Vect nu=new Vect(this.dim),nuVar=new Vect(this.dim);
		Vect B=new Vect();
		Vect  M=new Vect();		

		if(hasM){ M=model.element[ie].getM();}
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();


		int n3ang=this.PW3ang.length; 		

		int nGauss;

		if(!nonLinear){
			nu=model.element[ie].getNu();
			nGauss=this.PW[0].length; 
		}
		else{
			nGauss=this.PWNL[0].length; 
			if(this.centerNu){
				double nux,nuv=0;
				B=model.element[ie].getB();
				double Bn2=B.norm2();
				double Bn=sqrt(Bn2);

				nux=model.BH[nBH].getNu(Bn);
				nu=new Vect(nux,nux,nux);
				if(Bn>0){
					nuv=model.BH[nBH].getdHdB(Bn);	
					nuVar.el[0]=(nuv-nux)/Bn2;
					nuVar.el[1]=nuVar.el[0];
					nuVar.el[2]=nuVar.el[0];

				}
			}
		}

		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		Vect[] Cj=new Vect[model.nElEdge];
		for(int i=0;i<model.nElEdge;i++)
			Cj[i]=new Vect(3);

		double detJac,ws=1,wsJ=0,term1,term2,term3;

		Mat jac;


		Vect[] rotNe=new Vect[model.nElEdge];


		Vect localCo=new Vect(this.dim);

		Vect[] Ne=new Vect[model.nElEdge];




		///


		Vect sigma=new Vect();

		if(eddy)
			sigma=model.element[ie].getSigma();

		for(int p=0;p<n3ang;p++)
			for(int q=0;q<nGauss;q++)
			{
				if(nonLinear){
					localCo.el[0]=this.PW3ang[p][0];
					localCo.el[1]=this.PW3ang[p][1];
					localCo.el[2]=this.PWNL[0][q];

					if(nGauss!=2)
						ws=this.PWNL[1][q]*this.PW3ang[p][2];
					else
						ws=this.PW3ang[p][2];

					B=this.getBAt(model,ie,localCo);
					double Bn2=B.norm2();
					double Bn=sqrt(Bn2);

					if(!this.centerNu){
						double nux,nuv=0;


						nux=model.BH[nBH].getNu(Bn);
						nu=new Vect(nux,nux,nux);

						if(Bn>0){
							nuv=model.BH[nBH].getdHdB(Bn);	
							nuVar.el[0]=(nuv-nux)/Bn2;
							nuVar.el[1]=nuVar.el[0];
							nuVar.el[2]=nuVar.el[0];

						}
					}

				}
				else
				{
					localCo.el[0]=this.PW3ang[p][0];
					localCo.el[1]=this.PW3ang[p][1];
					localCo.el[2]=this.PW[0][q];

					if(nGauss!=2)
						ws=this.PW[1][p]*this.PW3ang[p][2];
					else
						ws=this.PW3ang[p][2];


				}

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				wsJ=ws*detJac;


				rotNe=rotNePrism(jac,localCo,edgeDir);
				if(hasJ || eddy)
					Ne=NePrism(jac,localCo,edgeDir);

				for(int i=0;i<model.nElEdge;i++){

					if(hasJ){
						Cj[i]=Cj[i].add(Ne[i].times(wsJ));
					}

					if(hasM){
						C[i]+=wsJ*rotNe[i].dot(M);

					}

					for(int j=0;j<=i;j++){

						term1=wsJ*rotNe[i].times(nu).dot(rotNe[j]);
						H1[i][j]+=term1;

						if(nonLinear){
							term2=wsJ*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
							H2[i][j]+=term2;


						}

						if(eddy){
							term3=wsJ*Ne[i].dot(Ne[j].times(sigma));
							H3[i][j]+=term3;

						}


					}
				}
			}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy){
			lowSym(H3);
		}

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj=Cj;


		return H1;


	}

	
	public double[][] HePyramid(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ,boolean hasM){


		double vol=0;
		Vect nu=new Vect(this.dim),nuVar=new Vect(this.dim);
		Vect B=new Vect();
		Vect  M=new Vect();		

		if(hasM){ M=model.element[ie].getM();}
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();



		double[][] PWpyramid=this.PWpyr;
		int nGauss=PWpyramid.length;

		if(!nonLinear){
			nu=model.element[ie].getNu();
			nGauss=PWpyramid.length; 
		}
		else{
			nGauss=PWpyramid.length; 
			if(this.centerNu){
				double nux,nuv=0;
				B=model.element[ie].getB();
				double Bn2=B.norm2();
				double Bn=sqrt(Bn2);

				nux=model.BH[nBH].getNu(Bn);
				nu=new Vect(nux,nux,nux);
				if(Bn>0){
					nuv=model.BH[nBH].getdHdB(Bn);	
					nuVar.el[0]=(nuv-nux)/Bn2;
					nuVar.el[1]=nuVar.el[0];
					nuVar.el[2]=nuVar.el[0];

				}
			}
		}

		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		Vect[] Cj=new Vect[model.nElEdge];
		for(int i=0;i<model.nElEdge;i++)
			Cj[i]=new Vect(3);

		double detJac,ws=1,wsJ=0,term1,term2,term3;

		Mat jac;


		Vect[] rotNe=new Vect[model.nElEdge];


		Vect localCo=new Vect(this.dim);

		Vect[] Ne=new Vect[model.nElEdge];



		Vect sigma=new Vect();

		if(eddy)
			sigma=model.element[ie].getSigma();

		for(int p=0;p<nGauss;p++){

				if(nonLinear){
					localCo.el[0]=PWpyramid[p][0];
					localCo.el[1]=PWpyramid[p][1];
					localCo.el[2]=PWpyramid[p][2];
					
					ws=PWpyramid[p][3];

	
					B=this.getBAt(model,ie,localCo);
					double Bn2=B.norm2();
					double Bn=sqrt(Bn2);

					if(!this.centerNu){
						double nux,nuv=0;


						nux=model.BH[nBH].getNu(Bn);
						nu=new Vect(nux,nux,nux);

						if(Bn>0){
							nuv=model.BH[nBH].getdHdB(Bn);	
							nuVar.el[0]=(nuv-nux)/Bn2;
							nuVar.el[1]=nuVar.el[0];
							nuVar.el[2]=nuVar.el[0];

						}
					}

				}
				else
				{
					localCo.el[0]=PWpyramid[p][0];
					localCo.el[1]=PWpyramid[p][1];
					localCo.el[2]=PWpyramid[p][2];
					ws=PWpyramid[p][3];

				}

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				wsJ=ws*detJac;


				rotNe=rotNePyramid(jac,localCo,edgeDir);
				//if(hasJ || eddy)
				//	Ne=NePyramid(jac,localCo,edgeDir);

				for(int i=0;i<model.nElEdge;i++){

					if(hasJ){
						Cj[i]=Cj[i].add(Ne[i].times(wsJ));
					}

					if(hasM){
						C[i]+=wsJ*rotNe[i].dot(M);

					}

					for(int j=0;j<=i;j++){

						term1=wsJ*rotNe[i].times(nu).dot(rotNe[j]);
						H1[i][j]+=term1;

						if(nonLinear){
							term2=wsJ*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
							H2[i][j]+=term2;


						}

						if(eddy){
							term3=wsJ*Ne[i].dot(Ne[j].times(sigma));
							H3[i][j]+=term3;

						}


					}
				}
		}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy){
			lowSym(H3);
		}

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj=Cj;

//util.pr("volume of el "+ie+" = "+vol);
		return H1;


	}


	public double[][] He3ang(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ, boolean hasM){

		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect M=new Vect();		
		B=model.element[ie].getB();
		Vect lc=new Vect(2);



		if(hasM){ M=model.element[ie].getM();}

		double sigma=0;
		if(eddy)
			sigma=model.element[ie].getSigma().el[2];

		if(nonLinear){

			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=0,nuv;

			nux=model.BH[nBH].getNu(Bn);

			nu=new Vect(nux,nux);
			if(Bn>0){
				nuv=model.BH[nBH].getdHdB(Bn);
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}
		}

		else
			nu=model.element[ie].getNu();


		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[3];
		double[] Cj=new double[3];

		double S=el3angArea(model,ie);
		double detJac,ws=1;
		Mat jac;
		Node[] vertexNode=model.elementNodes(ie);
		double[] localNe=new double[3];
		Vect[] rotNe=rotNe3ang(model,ie);


		if(hasJ){
			for(int i=0;i<3;i++)
				Cj[i]=S/3;

		}

		if(hasM){
			for(int i=0;i<3;i++)
				C[i]=S*rotNe[i].dot(M);
		}

		for(int i=0;i<this.nElEdge;i++)
			for(int j=0;j<=i;j++)	{
				H1[i][j]=S*rotNe[i].times(nu).dot(rotNe[j]);

				if(nonLinear)
					H2[i][j]=S*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));

				if(eddy){
					int n=this.PW3ang.length; 					
					for(int p=0;p<n;p++)
					{

						lc.el[0]=this.PW3ang[p][0];
						lc.el[1]=this.PW3ang[p][1];
						localNe[0]=lc.el[0];
						localNe[1]=lc.el[1];
						localNe[2]=1-lc.el[0]-lc.el[1];
						jac=jacobian3ang(vertexNode,lc);

						detJac=abs(jac.determinant());
						ws=this.PW3ang[p][2]*detJac;
						H3[i][j]+=ws*localNe[i]*localNe[j]*sigma;	
					}
				}


			}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;


		return H1;
	}

	public double[][] He3angAxi(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ, boolean hasM){

		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect M=new Vect();		
		B=model.element[ie].getB();
		Vect lc=new Vect(2);



		if(hasM){ M=model.element[ie].getM();}

		double sigma=0;
		if(eddy)
			sigma=model.element[ie].getSigma().el[2];

		if(nonLinear){

			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=0,nuv;

			nux=model.BH[nBH].getNu(Bn);

			nu=new Vect(nux,nux);
			if(Bn>0){
				nuv=model.BH[nBH].getdHdB(Bn);
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}
		}

		else
			nu=model.element[ie].getNu();



		Node[] vertexNode=model.elementNodes(ie);
		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];
		double[] Cj=new double[model.nElEdge];



		double rr=0;

		for(int i=0;i<model.nElEdge;i++){
			double r1=vertexNode[i].getCoord(0);
			rr+=0.33333333/r1;
		}



		double S=el3angArea(model,ie);
		double detJac,ws=1;
		Mat jac;
		double[] localNe=new double[3];
		Vect[] rotNe=rotNe3ang(model,ie);


		if(hasJ){
			for(int i=0;i<3;i++)
				Cj[i]=S/3;

		}

		if(hasM){
			for(int i=0;i<3;i++)
				C[i]=S*rotNe[i].dot(M);
		}

		for(int i=0;i<this.nElEdge;i++)
			for(int j=0;j<3;j++)	{
				H1[i][j]=S*rotNe[i].times(nu).dot(rotNe[j]);

				if(nonLinear)
					H2[i][j]=S*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));

				if(eddy){
					int n=this.PW3ang.length; 					
					for(int p=0;p<n;p++)
					{

						lc.el[0]=this.PW3ang[p][0];
						lc.el[1]=this.PW3ang[p][1];
						localNe[0]=lc.el[0];
						localNe[1]=lc.el[1];
						localNe[2]=1-lc.el[0]-lc.el[1];
						jac=jacobian3ang(vertexNode,lc);

						detJac=abs(jac.determinant());
						ws=this.PW3ang[p][2]*detJac;
						H3[i][j]+=ws*localNe[i]*localNe[j]*sigma;	
					}
				}


			}

		for(int i=0;i<model.nElEdge;i++)
			for(int j=0;j<model.nElEdge;j++){
				H1[i][j]*=rr;
				if(nonLinear)
					H2[i][j]*=rr;
				if(eddy)
					H3[i][j]*=rr;

			}


		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;


		return H1;
	}



	public double[][] He3ang1st(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ, boolean hasM){

		Vect nu=new Vect(2),nuVar=new Vect(2);
		Vect B=new Vect();
		Vect M=new Vect();		
		B=model.element[ie].getB();
		Vect lc=new Vect(2);

		if(hasM) M=model.element[ie].getM();

		double sigma=0;
		if(eddy)
			sigma=model.element[ie].getSigma().el[2];

		if(!nonLinear)
			nu=model.element[ie].getNu();
		else if(this.centerNu){
			B=model.element[ie].getB();
			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=model.BH[nBH].getNu(Bn);
			nu=new Vect(nux,nux);
			if(Bn>0){
				double nuv=model.BH[nBH].getdHdB(Bn);	
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
			}

		}

		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[3];
		double[] Cj=new double[3];

		double detJac,ws=1;
		Mat jac;
		Node[] vertexNode=model.elementNodes(ie);
		double[] localNe;
		Vect[] rotNe;


		int n=this.PW3ang.length; 					
		for(int p=0;p<n;p++){

			lc.el[0]=this.PW3ang[p][0];
			lc.el[1]=this.PW3ang[p][1];

			localNe=Ne3ang1st(lc);


			jac=jacobian3ang(vertexNode,lc);

			rotNe=rotNe3ang1st(jac,ie);

			detJac=abs(.5*jac.determinant());

			ws=this.PW3ang[p][2]*detJac;

			if(nonLinear){

				double Bn2=B.norm2();
				double Bn=sqrt(Bn2);

				if(!this.centerNu){

					double nux=model.BH[nBH].getNu(Bn);
					nu=new Vect(nux,nux);
					if(Bn>0){
						double nuv=model.BH[nBH].getdHdB(Bn);	
						nuVar.el[0]=(nuv-nux)/Bn2;
						nuVar.el[1]=nuVar.el[0];


					}
				}

			}


			if(hasJ){
				for(int i=0;i<3;i++)
					Cj[i]+=ws*localNe[i];
			}

			if(hasM){
				for(int i=0;i<3;i++)
					C[i]+=ws*rotNe[i].dot(M);
			}

			for(int i=0;i<this.nElEdge;i++)
				for(int j=0;j<=i;j++)	{
					H1[i][j]+=ws*rotNe[i].times(nu).dot(rotNe[j]);
					if(nonLinear){
						H2[i][j]+=ws*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
					}



					if(eddy){
						H3[i][j]+=ws*localNe[i]*localNe[j]*sigma;	
					}

				}
		}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;
		model.Cj2d=Cj;


		return H1;
	}

		/*public double[][] HeTetra(Model model,int nBH,int ie, boolean nonLinear,boolean eddy,boolean hasJ, boolean hasM){

		Vect nu=new Vect(3),nuVar=new Vect(3);
		Vect B=new Vect();
		Vect J=new Vect(), M=new Vect();		
		B=model.element[ie].getB();
		Vect lc=new Vect(3);

		if( hasJ)	{ J=model.element[ie].getJ().times(model.region[model.element[ie].getRegion()].getFactJ());}
		if(hasM) M=model.element[ie].getM();

		Vect sigma=new Vect(3);
		if(eddy)
			sigma=model.element[ie].getSigma().times(sigma);

		if(!nonLinear)
			nu=model.element[ie].getNu();
		else if(this.centerNu){
			B=model.element[ie].getB();
			double Bn2=B.norm2();
			double Bn=sqrt(Bn2);
			double nux=model.BH[nBH].getNu(Bn);
			nu=new Vect(nux,nux);
			if(Bn>0){
				double nuv=model.BH[nBH].getdHdB(Bn);	
				nuVar.el[0]=(nuv-nux)/Bn2;
				nuVar.el[1]=nuVar.el[0];
				nuVar.el[2]=nuVar.el[0];
			}

		}

		double[][] H1=new double[model.nElEdge][model.nElEdge];
		double[][] H2=new double[model.nElEdge][model.nElEdge];
		double[][] H3=new double[model.nElEdge][model.nElEdge];
		double[] C=new double[model.nElEdge];

		double detJac,ws=1;
		Mat jac;
		Node[] vertexNode=model.elementNodes(ie);
		Vect[] localNe;
		Vect[] rotNe;


		int n=this.PWtetra.length; 					
		for(int p=0;p<n;p++){

			lc.el[0]=this.PWtetra[p][0];
			lc.el[1]=this.PWtetra[p][1];
			lc.el[2]=this.PWtetra[p][2];

			localNe=NeTetra(lc);

			jac=jacobianTetra(vertexNode,lc);

			rotNe=rotNeTetra(jac,ie);

			detJac=abs(.5*jac.determinant());

			ws=this.PWtetra[p][3]*detJac;

			if(nonLinear){

				double Bn2=B.norm2();
				double Bn=sqrt(Bn2);

				if(!this.centerNu){

					double nux=model.BH[nBH].getNu(Bn);
					nu=new Vect(nux,nux);
					if(Bn>0){
						double nuv=model.BH[nBH].getdHdB(Bn);	
						nuVar.el[0]=(nuv-nux)/Bn2;
						nuVar.el[1]=nuVar.el[0];
						nuVar.el[2]=nuVar.el[0];

					}
				}

			}


			if(hasJ){
				for(int i=0;i<model.nElEdge;i++)
					C[i]+=localNe[i].dot(J)*ws;
			}

			if(hasM){
				for(int i=0;i<model.nElEdge;i++)
					C[i]+=ws*rotNe[i].dot(M);
			}

			for(int i=0;i<this.nElEdge;i++)
				for(int j=0;j<=i;j++)	{
					H1[i][j]+=ws*rotNe[i].times(nu).dot(rotNe[j]);
					if(nonLinear){
						H2[i][j]+=ws*rotNe[i].dot(B)*rotNe[j].dot(B.times(nuVar));
					}



					if(eddy){
						H3[i][j]+=ws*localNe[i].dot(localNe[j].times(sigma));	
					}

				}
		}

		lowSym(H1);
		if(nonLinear)
			lowSym(H2);
		if(eddy)
			lowSym(H3);

		model.H2=H2;
		model.H3=H3;
		model.C=C;


		return H1;
	}
	 */

	public Vect[] gradN(Mat jac,Vect localCo){

		if(this.elCode==3) return gradNPrism(jac,localCo);

		Mat invJac=jac.inv();
		Vect[] gradN=new Vect[this.nElVert];
		Vect[] localGradN=localGradN(localCo);
		for(int i=0;i<this.nElVert;i++){
			gradN[i]=invJac.mul(localGradN[i]);
		}

		return gradN;
	}

	public Vect[] gradNPrism(Mat jac,Vect localCo){


		Mat invJac=jac.inv();
		Vect[] gradN=new Vect[this.nElVert];
		Vect[] localGradN=localGradNPrism(localCo);
		for(int i=0;i<this.nElVert;i++){
			gradN[i]=invJac.mul(localGradN[i]);
		}

		return gradN;
	}

	

	public  Vect[] rotNe(Mat jac, Vect localCo, boolean[] edgeReverse){

		if(this.elCode==1) return rotNeQuad(jac,localCo);
		if(this.elCode==2) return rotNeTetra(jac,localCo,edgeReverse);
		if(this.elCode==3) return rotNePrism(jac,localCo,edgeReverse);
		if(this.elCode==5) return rotNePyramid(jac,localCo,edgeReverse);
		Vect[] rotNe=new Vect[this.nElEdge];

		double a=localCo.el[0];
		double b=localCo.el[1];
		double c=localCo.el[2];

		Mat invJac=jac.inv3();

		Vect[] grad=new Vect[3];

		for(int j=0;j<3;j++)
			grad[j]=invJac.getColVect(j);



		rotNe[0]= grad[1].times(-(1-c)).add(grad[2].times(-(1-b))).times(0.125).cross(grad[0]); 
		rotNe[1]= grad[1].times((1-c)).add(grad[2].times(-(1+b))).times(0.125).cross(grad[0]);
		rotNe[2]= grad[1].times(-(1+c)).add(grad[2].times(+(1-b))).times(0.125).cross(grad[0]);
		rotNe[3]= grad[1].times((1+c)).add(grad[2].times(+(1+b))).times(0.125).cross(grad[0]);
		rotNe[4]= grad[0].times(-(1-c)).add(grad[2].times(-(1-a))).times(0.125).cross(grad[1]);
		rotNe[5]= grad[0].times((1-c)).add(grad[2].times(-(1+a))).times(0.125).cross(grad[1]);
		rotNe[6]= grad[0].times(-(1+c)).add(grad[2].times(+(1-a))).times(0.125).cross(grad[1]);
		rotNe[7]= grad[0].times((1+c)).add(grad[2].times(+(1+a))).times(0.125).cross(grad[1]);
		rotNe[8]= grad[0].times(-(1-b)).add(grad[1].times(-(1-a))).times(0.125).cross(grad[2]);
		rotNe[9]= grad[0].times((1-b)).add(grad[1].times(-(1+a))).times(0.125).cross(grad[2]);
		rotNe[10]= grad[0].times(-(1+b)).add(grad[1].times(+(1-a))).times(0.125).cross(grad[2]);
		rotNe[11]= grad[0].times((1+b)).add(grad[1].times(+(1+a))).times(0.125).cross(grad[2]);

		for(int k=0;k<rotNe.length;k++)
			if(edgeReverse[k])
				rotNe[k].timesVoid(-1);


		return rotNe;
	}


	public double[] NQuad(Vect localCo){
		double a,b;
		a=localCo.el[0];b=localCo.el[1];;
		double[] N=new double[this.nElVert];
		double c=.25;
		N[0]=(1+a)*(1+b)*c;
		N[1]=(1-a)*(1+b)*c;
		N[2]=(1-a)*(1-b)*c;
		N[3]=(1+a)*(1-b)*c;

		return N;
	}
	public Vect[] rotNeQuad(Mat jac, Vect localCo){


		double a=localCo.el[0];
		double b=localCo.el[1];

		Vect[] rotNe=new Vect[this.nElVert];

		Mat invJac=jac.inv();

		Vect[] grad=new Vect[2];
		grad[0]=invJac.getColVect(0);
		grad[1]=invJac.getColVect(1);

		Vect v;
		v=grad[0].times(1+b).add(grad[1].times(1+a)).times(.25);
		rotNe[0]= new Vect(v.el[1],-v.el[0]);
		v=grad[0].times(-(1+b)).add(grad[1].times(1-a)).times(.25);
		rotNe[1]= new Vect(v.el[1],-v.el[0]);
		v=grad[0].times(b-1).add(grad[1].times(a-1)).times(.25);
		rotNe[2]= new Vect(v.el[1],-v.el[0]);
		v=grad[0].times(1-b).add(grad[1].times(-(1+a))).times(.25);
		rotNe[3]= new Vect(v.el[1],-v.el[0]);

		return rotNe;
	}




	public Vect[] rotNePrism(Mat jac,Vect localCo,boolean[] edgeReverse){

		double c=localCo.el[2];


		Mat invJac=jac.inv3();

		Vect[] grad=new Vect[3];

		for(int j=0;j<3;j++)
			grad[j]=invJac.getColVect(j);

		Vect[] rotNe=new Vect[9];


		rotNe[0]=grad[0].cross(grad[1]).times(2*c);
		rotNe[1]=rotNe[0];
		rotNe[2]=rotNe[0];

		rotNe[3]=grad[0].cross(grad[1]).times(2*(1-c));
		rotNe[4]=rotNe[3];
		rotNe[5]=rotNe[3];

		rotNe[6]=grad[0].cross(grad[2]);
		rotNe[7]=grad[1].cross(grad[2]);
		rotNe[8]=grad[0].add(grad[1]).times(-1).cross(grad[2]);

		for(int k=0;k<rotNe.length;k++)
			if(edgeReverse[k])
				rotNe[k].timesVoid(-1);


		return rotNe;
	}
	
	public Vect[] rotNePyramid(Mat jac,Vect localCo,boolean[] edgeReverse){

		double u=localCo.el[0];
		double v=localCo.el[1];
		double w=localCo.el[2];

		Mat invJac=jac.inv3();

		Vect[] grad=new Vect[3];

		for(int j=0;j<3;j++)
			grad[j]=invJac.getColVect(j);

		Vect[] rotNe=new Vect[8];

		Vect cross12=grad[1].cross(grad[2]);
		Vect cross20=grad[2].cross(grad[0]);
		Vect cross01=grad[0].cross(grad[1]);

		rotNe[0] = cross12.times(u/(4*(w - 1))).add(cross20.times(- v/(4*(w - 1)) - 0.5)).add(cross01.times(.25));
		rotNe[1] = cross12.times(u/(4*(w - 1)) + 0.5).add(cross20.times(-v/(4*(w - 1)))).add(cross01.times(-.25));	
		rotNe[2] = cross12.times(-0.5*((u+w-1)/(w-1))).add(cross20.times(0.5*(v+w-1)/(w-1)));
		rotNe[3] = cross12.times(0.5 - u/(4*(w - 1))).add(cross20.times(v/(4*(w - 1)))).add(cross01.times(.25));
		rotNe[4] = cross12.times((u/(4*(w - 1)))).add(cross20.times(-0.5*((v+w-1)/(w-1))));	
		rotNe[5] = cross12.times(u/(4*(w - 1))).add(cross20.times(0.5 - v/(4*(w - 1)))).add(cross01.times(.25));
		rotNe[6] = cross12.times(-0.5*((u-w+1)/(w-1))).add(cross20.times(+0.5*((v-w+1)/(w-1))));
		rotNe[7] = cross12.times(0.5*((u+w-1)/(w-1))).add(cross20.times(-0.5*((v-w+1)/(w-1))));
		

		for(int k=0;k<rotNe.length;k++)
			if(edgeReverse[k])
				rotNe[k].timesVoid(-1);


		return rotNe;
	}
	
	public Vect[] rotNeTetra(Mat jac,Vect localCo,boolean[] edgeReverse){



		Vect[] rotNe=new Vect[6];

		Vect[] gradN=gradNTet(jac);

		rotNe[0]=gradN[0].cross(gradN[1]).times(2);
		rotNe[1]=gradN[1].cross(gradN[2]).times(2);
		rotNe[2]=gradN[2].cross(gradN[0]).times(2);
		rotNe[3]=gradN[0].cross(gradN[3]).times(2);
		rotNe[4]=gradN[1].cross(gradN[3]).times(2);
		rotNe[5]=gradN[2].cross(gradN[3]).times(2);


		for(int k=0;k<rotNe.length;k++)
			if(edgeReverse[k])
				rotNe[k].timesVoid(-1);


		return rotNe;
	}
	

	public double[] NeQuad(Vect localCo){

		double a=localCo.el[0];
		double b=localCo.el[1];

		double[] Ne=new double[this.nElVert];

		Ne[0]= (1+a)*(1+b)*0.25; 
		Ne[1]= (1-a)*(1+b)*0.25; 
		Ne[2]= (1-a)*(1-b)*0.25; 
		Ne[3]= (1+a)*(1-b)*0.25; 


		return Ne;
	}


	Vect[] localGradN(Vect localCo){
		if(this.elCode==0) return localGradN3ang();
		else if(this.elCode==1) return localGradNQuad(localCo);
		else if(this.elCode==3) return localGradNPrism(localCo);
		else if(this.elCode==5) return localGradNPyramid(localCo);
		double a=localCo.el[0];
		double b=localCo.el[1];
		double c=localCo.el[2];

		Vect[] gradN=new Vect[this.nElVert];

		gradN[0]=new Vect((1+b)*(1+c),(1+a)*(1+c),(1+a)*(1+b));
		gradN[1]=new Vect(-(1+b)*(1+c),(1-a)*(1+c),(1-a)*(1+b));
		gradN[2]=new Vect(-(1-b)*(1+c),-(1-a)*(1+c),(1-a)*(1-b));
		gradN[3]=new Vect((1-b)*(1+c),-(1+a)*(1+c),(1+a)*(1-b));
		gradN[4]=new Vect((1+b)*(1-c),(1+a)*(1-c),-(1+a)*(1+b));
		gradN[5]=new Vect(-(1+b)*(1-c),(1-a)*(1-c),-(1-a)*(1+b));
		gradN[6]=new Vect(-(1-b)*(1-c),-(1-a)*(1-c),-(1-a)*(1-b));
		gradN[7]=new Vect((1-b)*(1-c),-(1+a)*(1-c),-(1+a)*(1-b));


		for(int i=0;i<this.nElVert;i++)
			gradN[i].timesVoid(0.125);
		return gradN;
	}

	public Vect[] localGradN3ang(){


		Vect[] gradN=new Vect[3];

		gradN[0]=new Vect(1,0,0);
		gradN[1]=new Vect(0,1,0);
		gradN[2]=new Vect(0,0,1);

		return gradN;
	}

	public Vect[] localGradNtet(){


		Vect[] gradN=new Vect[4];

		gradN[0]=new Vect(1,0,0);
		gradN[1]=new Vect(0,1,0);
		gradN[2]=new Vect(0,0,1);
		gradN[3]=new Vect(-1,-1,-1);

		return gradN;
	}



	public Vect[] localGradNPrism(Vect localCo){

		double a=localCo.el[0]/2;
		double b=localCo.el[1]/2;
		double c=.5-a-b;

		double h1=(1-localCo.el[2])/2;
		double h2=(1+localCo.el[2])/2;
		Vect[] gradN=new Vect[6];

		gradN[0]=new Vect(h2,0,a);
		gradN[1]=new Vect(0,h2,b);
		gradN[2]=new Vect(-h2,-h2,c);
		gradN[3]=new Vect(h1,0,-a);
		gradN[4]=new Vect(0,h1,-b);
		gradN[5]=new Vect(-h1,-h1,-c);

		return gradN;
	}
	
	public Vect[] localGradNPyramid(Vect localCo){

		double u=localCo.el[0];
		double v=localCo.el[1];
		double w=localCo.el[2];

	
		Vect[] gradN=new Vect[6];
		
		double val3=u*v/((1-w)*(1-w));

		gradN[0]=new Vect(-(1-v)+v*w/(1-w),-(1-u)+u*w/(1-w),-1+val3).times(.25);
		gradN[1]=new Vect((1-v)+v*w/(1-w),-(1+u)+u*w/(1-w),-1-val3).times(.25);
		gradN[2]=new Vect((1+v)+v*w/(1-w),(1+u)+u*w/(1-w),-1+val3).times(.25);
		gradN[3]=new Vect(-(1+v)+v*w/(1-w),(1-u)+u*w/(1-w),-1-val3).times(.25);
		gradN[4]=new Vect(0,0,1);

		return gradN;
	}

	public Vect[] localGradNQuad(Vect localCo){
		double a=localCo.el[0];
		double b=localCo.el[1];

		Vect[] gradN=new Vect[this.nElVert];

		gradN[0]=new Vect((1+b),(1+a));
		gradN[1]=new Vect(-(1+b),(1-a));
		gradN[2]=new Vect(-(1-b),-(1-a));
		gradN[3]=new Vect((1-b),-(1+a));


		for(int i=0;i<this.nElVert;i++)
			gradN[i].timesVoid(0.25);
		return gradN;
	}

	Vect[] gradN3ang(Model model,int ie){
		Node[] vertexNode=model.elementNodes(ie);

		Vect v1=vertexNode[0].getCoord();
		Vect v2=vertexNode[1].getCoord();
		Vect v3=vertexNode[2].getCoord();

		double rS=.5/el3angArea(model,ie);

		Vect[] gradN=new Vect[3];

		gradN[0]=new Vect(v2.el[1]-v3.el[1],v3.el[0]-v2.el[0]).times(rS);
		gradN[1]=new Vect(v3.el[1]-v1.el[1],v1.el[0]-v3.el[0]).times(rS);
		gradN[2]=new Vect(v1.el[1]-v2.el[1],v2.el[0]-v1.el[0]).times(rS);

		return gradN;
	}

	Vect[] gradNTet(Mat jac){
		
		Mat grads=jac.inv().transp();
		
		Vect[] gradN=new Vect[4];

		for(int i=0;i<3;i++)
			gradN[i]=grads.getColVect(i);
		
		gradN[3]=gradN[0].add(gradN[1]).add(gradN[2]).times(-1);
				
		return gradN;
	}
	
	Vect[] gradNTetOlder(Mat jac){
	//	return gradNTetSimple(jac);

		Vect[] locgN=this.localGradNtet();
		Mat grads=new Mat(3,4) ;
		for(int i=0;i<4;i++)
			grads.setCol(locgN[i], i);

		grads=jac.inv().transp().mul(grads);// correct
		
		Vect[] gradN=new Vect[4];

		for(int i=0;i<4;i++)
			gradN[i]=grads.getColVect(i);
		return gradN;
	}

	public double getElementArea(Model model,int i){
		if(this.elCode==0) return el3angArea(model,i);
		else 	if(this.elCode==1) return elQuadArea(model,i);
		else throw new NullPointerException(" Element is not 2D. ");

	}

	public double el3angArea(Model model,int i){
		Node[] vertexNode=model.elementNodes(i);

		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=abs(v1.cross(v2).norm())/2;
		return S;
	}

	public double elPrismVol(Model model,int i){
		Node[] vertexNode=model.elementNodes(i);
		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=v1.cross(v2).norm()/2;
		double h=abs(vertexNode[0].getCoord(2)-vertexNode[3].getCoord(2));
		return S*h;
	}

	public double elPrismBase(Model model,int i){
		Node[] vertexNode=model.elementNodes(i);

		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[2].getCoord().sub(vertexNode[0].getCoord());
		double S=v1.cross(v2).norm()/2;

		return S;
	}

	public double elQuadArea(Model model,int i){
		Node[] vertexNode=model.elementNodes(i);

		Vect v1=vertexNode[1].getCoord().sub(vertexNode[0].getCoord());
		Vect v2=vertexNode[3].getCoord().sub(vertexNode[0].getCoord());
		Vect v3=vertexNode[1].getCoord().sub(vertexNode[2].getCoord());
		Vect v4=vertexNode[3].getCoord().sub(vertexNode[2].getCoord());
		double S=(v1.cross(v2).norm()+v4.cross(v3).norm())/2;
		return S;
	}

	double[] N(Vect localCo){
		if(elCode==0) return N3ang1st(localCo);
		else if(elCode==1) return NQuad(localCo);
		double a,b,c;
		a=localCo.el[0];b=localCo.el[1];c=localCo.el[2];
		double[] N=new double[this.nElVert];
		N[0]=(1+a)*(1+b)*(1+c)*0.125;
		N[1]=(1-a)*(1+b)*(1+c)*0.125;
		N[2]=(1-a)*(1-b)*(1+c)*0.125;
		N[3]=(1+a)*(1-b)*(1+c)*0.125;
		N[4]=(1+a)*(1+b)*(1-c)*0.125;
		N[5]=(1-a)*(1+b)*(1-c)*0.125;
		N[6]=(1-a)*(1-b)*(1-c)*0.125;
		N[7]=(1+a)*(1-b)*(1-c)*0.125;

		return N;
	}

	public  Vect[] Ne(Mat jac,Vect localCo, boolean[] edgeReverse){
		if(elCode==2) return NeTetra(jac,localCo,  edgeReverse);
	//	if(elCode==5) return NePyramid(jac,localCo,  edgeReverse);
		
		double a=localCo.el[0];
		double b=localCo.el[1];
		double c=localCo.el[2];
		Vect[] Ne=new Vect[this.nElEdge];
		Mat invJac=jac.inv3();

		Vect[] grad=new Vect[3];

		for(int j=0;j<3;j++)
			grad[j]=invJac.getColVect(j);

		Ne[0]= grad[0].times((1-b)*(1-c)*0.125); 
		Ne[1]= grad[0].times((1+b)*(1-c)*0.125); 
		Ne[2]= grad[0].times((1-b)*(1+c)*0.125); 
		Ne[3]= grad[0].times((1+b)*(1+c)*0.125); 
		Ne[4]= grad[1].times((1-a)*(1-c)*0.125); 
		Ne[5]= grad[1].times((1+a)*(1-c)*0.125); 
		Ne[6]= grad[1].times((1-a)*(1+c)*0.125); 
		Ne[7]= grad[1].times((1+a)*(1+c)*0.125); 
		Ne[8]= grad[2].times((1-a)*(1-b)*0.125); 
		Ne[9]= grad[2].times((1+a)*(1-b)*0.125); 
		Ne[10]= grad[2].times((1-a)*(1+b)*0.125); 
		Ne[11]= grad[2].times((1+a)*(1+b)*0.125); 
		
		for(int k=0;k<Ne.length;k++)
			if(edgeReverse[k])
				Ne[k].timesVoid(-1);


		return Ne;

	}
	public  Vect[] NeTetra(Mat jac,Vect localCo, boolean[] edgeReverse)
	{

		
		double a=localCo.el[0];
		double b=localCo.el[1];
		double c=localCo.el[2];
		double d=localCo.el[3];
		
		Vect[] Ne=new Vect[this.nElEdge];


		Vect[] gradN=gradNTet(jac);
		double[] N=new double [4];
		N[0]=a;
		N[1]=b;
		N[2]=c;
		N[3]=d;

		
		 Ne[0]= gradN[1].times(N[0]).sub(gradN[0].times(N[1])) ;
		 Ne[1]= gradN[2].times(N[1]).sub(gradN[1].times(N[2])) ;
		 Ne[2]= gradN[0].times(N[2]).sub(gradN[2].times(N[0])) ;
		 
		 Ne[3]= gradN[3].times(N[0]).sub(gradN[0].times(N[3])) ;
		 Ne[4]= gradN[3].times(N[1]).sub(gradN[1].times(N[3])) ;
		 Ne[5]= gradN[3].times(N[2]).sub(gradN[2].times(N[3])) ;

         
 		for(int k=0;k<Ne.length;k++)
			if(edgeReverse[k])
				Ne[k].timesVoid(-1);

	         
	  return Ne;

	}
	

	public  Vect[] NePrism(Mat jac,Vect localCo, boolean[] edgeReverse)
	{

		double a=localCo.el[0];
		double b=localCo.el[1];
		double c=localCo.el[2];

		Vect[] Ne=new Vect[this.nElEdge];

		Mat invJac=jac.inv3();

		Vect[] grad=new Vect[3];

		for(int j=0;j<3;j++)
			grad[j]=invJac.getColVect(j);


		Ne[0]= grad[1].times(a).sub( grad[0].times(b)).times(c);
		Ne[1]= Ne[0].sub( grad[1]).times(c);
		Ne[2]= Ne[0].add( grad[0]).times(c);


		Ne[3]= grad[1].times(a).sub( grad[0].times(b)).times(1-c);
		Ne[4]= Ne[3].sub( grad[1]).times(1-c);
		Ne[5]= Ne[3].add( grad[0]).times(1-c);

		Ne[6]= grad[2].times(a);
		Ne[7]= grad[2].times(b);
		Ne[8]= grad[2].times(1-a-b);
		
		for(int k=0;k<Ne.length;k++)
			if(edgeReverse[k])
				Ne[k].timesVoid(-1);

		return Ne;

	}

	double[] localNPrims(Vect lc){

		double[] N=new double[6];


		double c1=(1+lc.el[2])/2;
		double c2=(1-lc.el[2])/2;

		N[0]=lc.el[0]*c1;
		N[1]=lc.el[1]*c1;
		N[2]=(1-lc.el[0]-lc.el[1])*c1;
		N[3]=lc.el[0]*c2;
		N[4]=lc.el[1]*c2;
		N[5]=(1-lc.el[0]-lc.el[1])*c2;

		return N;


	}

	double[] NPrims(Node[] vertexNode, Vect globalCo){

		Vect P=new Vect(globalCo.el[0],globalCo.el[1]);
		Vect v3D1=vertexNode[3].getCoord();
		Vect v3D2=vertexNode[4].getCoord();
		Vect v3D3=vertexNode[5].getCoord();
		Vect v3D4=vertexNode[0].getCoord();

		double H=v3D1.el[2]-v3D4.el[2];
		double h=globalCo.el[2]-v3D4.el[2];

		Vect v1=new Vect(v3D1.el[0],v3D1.el[1]).sub(P);
		Vect v2=new Vect(v3D2.el[0],v3D2.el[1]).sub(P);
		Vect v3=new Vect(v3D3.el[0],v3D3.el[1]).sub(P);

		double c1=h/H;
		double c2=(1-h/H);
		double[] N=new double[6];

		double[] s=new double[3];
		s[0]=v2.cross(v3).norm();
		s[1]=v3.cross(v1).norm();
		s[2]=v1.cross(v2).norm();

		for(int i=0;i<6;i++)
			if(i<3)
				N[i]=s[i]*c1/2;
			else
				N[i]=s[i-3]*c2/2;

		return N;


	}

	double[] localNPyramid(Vect lc){
		

		
		double u=lc.el[0];
		double v=lc.el[1];
		double w=lc.el[2];
		
		double r=0;
		   if(w != 1.) r = u*v*w / (1.-w) ;
	

		
		double[] N=new double[5];
	
	    
	     N[0]=1/4*((1-u)*(1-v)-w+r);
	     N[1]=1/4*((1+u)*(1-v)-w-r);
	     N[2]=1/4*((1+u)*(1+v)-w+r);
	     N[3]=1/4*((1-u)*(1+v)-w-r);
	     N[4]=w;
	     
	     return N;
	}


	//============

	public double[] N3ang(Node[] vertexNode, Vect globalCo){

		Vect v1=vertexNode[0].getCoord().sub(globalCo);
		Vect v2=vertexNode[1].getCoord().sub(globalCo);
		Vect v3=vertexNode[2].getCoord().sub(globalCo);

		double[] N=new double[3];

		N[0]=v2.cross(v3).norm();
		N[1]=v3.cross(v1).norm();
		N[2]=v1.cross(v2).norm();

		double S=N[0]+N[1]+N[2];

		for(int i=0;i<3;i++)
			N[i]/=S;

		return N;


	}
	
	public double[] N3ang(Model model,int ie, Vect globalCo){
		Node[] vertexNode=model.elementNodes(ie);

		return N3ang(vertexNode,globalCo);

	}

	public double[] N3ang(Model model,int ie){
		Node[] vertexNode=model.elementNodes(ie);
		Vect globalCo=vertexNode[0].getCoord().add(vertexNode[1].getCoord()).add(vertexNode[2].getCoord()).times(1.0/3);

		return N3ang(vertexNode,globalCo);


	}

	public double[] N3ang(Node[] vertexNode){

		Vect globalCo=vertexNode[0].getCoord().add(vertexNode[1].getCoord()).add(vertexNode[2].getCoord()).times(1.0/3);
		return N3ang(vertexNode,globalCo);


	}


	public double[] Ne3ang(Node[] vertexNode){
		return N3ang(vertexNode);

	}

	double[] Ne3ang(Node[] vertexNode, Vect globalCo){

		return N3ang(vertexNode, globalCo);
	}

	double[] N3ang1st(Vect localCo){
		double[] N=new double[3];
		N[0]=localCo.el[0];
		N[1]=localCo.el[1];
		N[2]=1-N[0]-N[1];
		return N;
	}

	public double [] localNTetra(Vect localCo){
		double[] N=new double[3];
		N[0]=localCo.el[0];
		N[1]=localCo.el[1];
		N[2]=localCo.el[2];
		N[3]=1-N[0]-N[1]-N[2];
		return N;
	}


	public double[] Ne3ang1st(Vect localCo){
		return N3ang1st(localCo);
	}


	double[] Ne3ang(Model model,int ie, Vect globalCo){
		return Ne3ang(model.elementNodes(ie), globalCo);
	}


	public Vect[] rotNe3ang(Node[] vertexNode){

		Vect v1=vertexNode[0].getCoord();
		Vect v2=vertexNode[1].getCoord();
		Vect v3=vertexNode[2].getCoord();

		double rS=1.0/v2.sub(v1).cross(v3.sub(v1)).norm();

		Vect[] rotNe=new Vect[3];

		rotNe[0]=new Vect(v3.el[0]-v2.el[0], v3.el[1]-v2.el[1]).times(rS);
		rotNe[1]=new Vect(v1.el[0]-v3.el[0], v1.el[1]-v3.el[1]).times(rS);
		rotNe[2]=new Vect(v2.el[0]-v1.el[0],v2.el[1]-v1.el[1]).times(rS);

		return rotNe;
	}

	private Vect[] rotNe3ang1st(Mat jac,int ie){


		Vect[] rotNe=new Vect[this.nElEdge];


		double Jr=1.0/(jac.determinant());
		Mat P=new Mat(3,2);
		P.el[0][0]=jac.el[2][1]-jac.el[2][2];
		P.el[0][1]=jac.el[1][2]-jac.el[1][1];
		P.el[1][0]=jac.el[2][2]-jac.el[2][0];
		P.el[1][1]=jac.el[1][0]-jac.el[1][2];
		P.el[2][0]=jac.el[2][0]-jac.el[2][1];
		P.el[2][1]=jac.el[1][1]-jac.el[1][0];

		P=P.times(Jr);

		Vect[] grad=new Vect[3];

		for(int i=0;i<3;i++){
			grad[i]=P.rowVect(i);

			rotNe[i]=new Vect(grad[i].el[1], -grad[i].el[0]);

		}

		return rotNe;


	}

	private Vect[] rotNe3ang2nd(Mat jac, Vect lc){

		double[] tu=lc.el;

		Vect[] rotNe=new Vect[6];


		double Jr=1.0/(jac.determinant());
		Mat P=new Mat(3,2);
		P.el[0][0]=jac.el[2][1]-jac.el[2][2];
		P.el[0][1]=jac.el[1][2]-jac.el[1][1];
		P.el[1][0]=jac.el[2][2]-jac.el[2][0];
		P.el[1][1]=jac.el[1][0]-jac.el[1][2];
		P.el[2][0]=jac.el[2][0]-jac.el[2][1];
		P.el[2][1]=jac.el[1][1]-jac.el[1][0];

		P=P.times(Jr);


		Vect[] grad=new Vect[3];

		for(int i=0;i<3;i++){
			grad[i]=P.rowVect(i).times(4*tu[i]-1);
		}

		grad[3]=P.rowVect(0).times(4*tu[1]).add(P.rowVect(1).times(4*tu[0]));
		grad[4]=P.rowVect(1).times(4*tu[2]).add(P.rowVect(2).times(4*tu[1]));
		grad[5]=P.rowVect(0).times(4*tu[2]).add(P.rowVect(2).times(4*tu[0]));



		for(int i=0;i<3;i++){
			rotNe[i]=new Vect(grad[i].el[1], -grad[i].el[0]);

		}


		return rotNe;


	}


	private Vect[] rotNe3ang(Model model,int ie){
		Node[] vertexNode=model.elementNodes(ie);
		return rotNe3ang(vertexNode);

	}

	private double[] Ne3ang(Model model,int ie){
		Node[] vertexNode=model.elementNodes(ie);
		return Ne3ang(vertexNode);

	}

	public double[] N3ang2nd(Vect localCo){

		double[] N=new double[6];
		double[] tu=localCo.el;

		for(int i=0;i<3;i++)
			N[i]=tu[i]*(2*tu[i]-1);

		N[3]=4*tu[0]*tu[1];
		N[4]=4*tu[1]*tu[2];
		N[5]=4*tu[0]*tu[2];

		return N;
	}

	public double[] Ne3ang2nd(Vect localCo){
		return N3ang2nd(localCo);
	}

	public Vect[] rotNe3ang2nd(Node[] vertexNode,Vect globalCo){

		Vect v1=vertexNode[0].getCoord().sub(globalCo);
		Vect v2=vertexNode[1].getCoord().sub(globalCo);
		Vect v3=vertexNode[2].getCoord().sub(globalCo);

		double[] N=new double[6];

		double[] su=new double [3];
		su[0]=v2.cross(v3).norm();
		su[1]=v3.cross(v1).norm();
		su[2]=v1.cross(v2).norm();

		double S=su[0]+su[1]+su[2];
		double rS=1.0/S;
		double[] tu=new double [3];
		for(int i=0;i<3;i++)
			tu[i]=su[i]*rS;

		Vect[] rotNe=new Vect[6];

		Vect[] rn=new Vect[3];
		rn[0]=new Vect(v3.el[0]-v2.el[0], v3.el[1]-v2.el[1]).times(rS);
		rn[1]=new Vect(v1.el[0]-v3.el[0], v1.el[1]-v3.el[1]).times(rS);
		rn[2]=new Vect(v2.el[0]-v1.el[0],v2.el[1]-v1.el[1]).times(rS);
		rotNe[0]=rn[0].times(4*tu[0]-1);
		rotNe[1]=rn[1].times(4*tu[1]-1);
		rotNe[2]=rn[2].times(4*tu[2]-1);
		rotNe[3]=rn[0].times(4*tu[1]).add(rn[1].times(4*tu[0]));
		rotNe[4]=rn[1].times(4*tu[2]).add(rn[2].times(4*tu[1]));
		rotNe[5]=rn[2].times(4*tu[0]).add(rn[0].times(4*tu[2]));


		return rotNe;
	}


	public Mat jacobian(Node[] vertexNode,Vect localCo){
		if(this.elCode==0) return jacobian3ang(vertexNode,localCo);
		if(this.elCode==2) return jacobianTetra(vertexNode,localCo);
		if(this.elCode==3) return jacobianPrism(vertexNode,localCo);
		if(this.elCode==5) return jacobianPyramid(vertexNode,localCo);


		Mat J=new Mat(this.dim,this.dim);
		Vect[] gN=new Vect[this.nElVert];

		gN=localGradN(localCo);

		for(int i=0;i<this.nElVert;i++) {
			for(int j=0;j<this.dim;j++)
				for(int k=0;k<this.dim;k++)
					J.el[j][k]+=gN[i].el[j]* vertexNode[i].getCoord(k);
		}

		return J;
	}



	public Mat jacobian3ang(Node[] vertexNode,Vect localCo){

		Mat J=new Mat(3,3);
		Vect[] gN=new Vect[3];

		gN=localGradN3ang();

		for(int i=0;i<3;i++) 
			J.el[0][i]=1;

		for(int i=0;i<3;i++) {
			for(int j=1;j<3;j++)
				for(int k=0;k<3;k++)
					J.el[j][k]+=gN[i].el[k]* vertexNode[i].getCoord(j-1);
		}


		return J;
	}

	public Mat jacobianTetra(Node[] vertexNode,Vect localCo){

		Mat J=new Mat(3,3);

	

		J.el[0][0]=vertexNode[0].getCoord(0)-vertexNode[3].getCoord(0);
		J.el[1][0]=vertexNode[0].getCoord(1)-vertexNode[3].getCoord(1);
		J.el[2][0]=vertexNode[0].getCoord(2)-vertexNode[3].getCoord(2);
		J.el[0][1]=vertexNode[1].getCoord(0)-vertexNode[3].getCoord(0);
		J.el[1][1]=vertexNode[1].getCoord(1)-vertexNode[3].getCoord(1);
		J.el[2][1]=vertexNode[1].getCoord(2)-vertexNode[3].getCoord(2);		
		J.el[0][2]=vertexNode[2].getCoord(0)-vertexNode[3].getCoord(0);
		J.el[1][2]=vertexNode[2].getCoord(1)-vertexNode[3].getCoord(1);
		J.el[2][2]=vertexNode[2].getCoord(2)-vertexNode[3].getCoord(2);



		return J;
	}

	public Mat jacobianPrism(Node[] vertexNode,Vect localCo){

		Mat J=new Mat(3,3);

		Vect[] gN=new Vect[this.nElVert];

		gN=localGradNPrism(localCo);

		for(int i=0;i<this.nElVert;i++) {
			for(int j=0;j<this.dim;j++)
				for(int k=0;k<this.dim;k++)
					J.el[j][k]+=gN[i].el[j]* vertexNode[i].getCoord(k);
		}

		return J;
	}

	
	public Mat jacobianPyramid(Node[] vertexNode,Vect localCo){

		Mat J=new Mat(3,3);

		Vect[] gN=new Vect[this.nElVert];
	

		gN=localGradNPyramid(localCo);

		for(int i=0;i<this.nElVert;i++) {
			for(int j=0;j<this.dim;j++)
				for(int k=0;k<this.dim;k++)
					J.el[j][k]+=gN[i].el[j]* vertexNode[i].getCoord(k);
		}

		return J;
	}

	public Mat Qe(Model model,int ie){

	//	if(model.elCode==0) return new Mat(Qe3ang(model,ie));
		if(model.elCode==2) return  QeTetra(model,ie);


		
		Vect sigma=model.element[ie].getSigma();
			
		Node[] vertexNode=model.elementNodes(ie);
		Mat M=new Mat(this.nElVert,this.nElVert);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN=new Vect[this.nElVert];
		Vect localCo=new Vect(3);

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;


					gradN=gradN(jac,localCo);


					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<=i;j++)	
							M.el[i][j]+=ws*gradN[i].dot(gradN[j].times(sigma));
				}	


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<i;j++)	
			{
				M.el[j][i]=M.el[i][j];
			}
		
		return M;

	}
	
	public Mat QeTetra(Model model,int ie){

		
		Vect sigma=model.element[ie].getSigma();
			

		Node[] vertexNode=model.elementNodes(ie);


		
		Mat M=new Mat(this.nElVert,this.nElVert);



		Mat jac=this.jacobianTetra(vertexNode, new Vect(dim));
		Vect[] gradN=this.gradNTet(jac);

		double ws=abs(jac.determinant())/6;


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<=i;j++)	
				M.el[i][j]=ws*gradN[i].dot(gradN[j].times(sigma));



		for(int i=0;i<this.nElVert;i++)
		for(int j=0;j<i;j++)	
		{
			M.el[j][i]=M.el[i][j];
		}
			
		return M;

	}


	public Mat elemPhiMat(Model model,int ie){

		if(model.elCode==1) return elemPhiMatQuad(model, ie);
		if(model.elCode==2) return elemPhiMatTet(model, ie);
			
		Node[] vertexNode=model.elementNodes(ie);
		Mat M=new Mat(this.nElVert,this.nElVert);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN=new Vect[this.nElVert];
		Vect localCo=new Vect(3);

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;


					gradN=gradN(jac,localCo);


					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<=i;j++)	
							M.el[i][j]+=ws*gradN[i].dot(gradN[j]);
				}	


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<=i;j++)	
			{
				M.el[j][i]=M.el[i][j];
			}

		return M;

	}
	
	public Mat elemPhiMatTet(Model model,int ie){

		Node[] vertexNode=model.elementNodes(ie);


	
		Mat M=new Mat(this.nElVert,this.nElVert);



		Mat jac=this.jacobianTetra(vertexNode, new Vect(dim));
		Vect[] gradN=this.gradNTet(jac);

		double ws=abs(jac.determinant())/6;




		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){
	
					M.el[i][j]=ws*gradN[i].dot(gradN[j]);
			}
		



		return M;
	}
	
	public Mat elemPhiMatQuad(Model model,int ie){

			
		Node[] vertexNode=model.elementNodes(ie);
		Mat M=new Mat(this.nElVert,this.nElVert);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN=new Vect[this.nElVert];
		Vect localCo=new Vect(model.dim);

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				{
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];


					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*detJac;
					else
						ws=detJac;


					gradN=gradN(jac,localCo);


					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<=i;j++)	
							M.el[i][j]+=ws*gradN[i].dot(gradN[j]);
				}	


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<=i;j++)	
			{
				M.el[j][i]=M.el[i][j];
			}

		return M;

	}
	

	public Vect elemPhiVect(Model model,int ie){

		if(model.elCode==1) return elemPhiVectQuad(model, ie);
		
		Node[] vertexNode=model.elementNodes(ie);
		Edge[] elemEdges=model.elementEdges(ie);
		Vect elemVec=new Vect(this.nElVert);
		
		boolean[] edgeDir=model.element[ie].getEdgeReverse();


		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN=new Vect[this.nElVert];
		
		Vect localCo=new Vect(3);
		Vect[] Ne=new Vect[this.nElEdge];
		
		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	
			A[j]=elemEdges[j].getA();
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];

					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;


					gradN=gradN(jac,localCo);

					Ne=Ne(jac,localCo,edgeDir);
					
					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<this.nElEdge;j++)	
							elemVec.el[i]+=ws*gradN[i].dot(Ne[j].times(A[j]));
				}	

		return elemVec;

	}
	
	public Vect elemPhiVectQuad(Model model,int ie){

		
		Node[] vertexNode=model.elementNodes(ie);
		Edge[] elemEdges=model.elementEdges(ie);
		Vect elemVec=new Vect(this.nElVert);

		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect localCo=new Vect(model.dim);
		double[] N=new double[this.nElEdge];
		
		double[] A=new double[this.nElEdge];
		for(int j=0;j<this.nElEdge;j++)	
			A[j]=elemEdges[j].getA();
		
		double gradN=1./model.height;
		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];


					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*detJac;
					else
						ws=detJac;


					N=N(localCo);
					

					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<this.nElEdge;j++)	
									elemVec.el[i]+=ws*gradN*N[j]*A[j];
				}	
		

		return elemVec;

	}

/*	public double[][] Qe3ang(Model model, int ie){
		double rdt=1.0/model.dt;

		double[][] Qe=new double[this.nElVert][this.nElEdge];

		double S=el3angArea(model,ie);
		double sigmaZ=model.element[ie].getSigma().el[2]*rdt;


		Vect[] gradN=gradN3ang(model,ie);


		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElEdge;j++)	{
				Qe[i][j]=S*gradN[i].dot(gradN[j].times(sigmaZ));

			}



		return Qe;


	}
*/

	public double[][] Pe(Model model, int ie){
		
		if(model.elCode==2) return PeTetra(model,ie);

		Vect sigma=model.element[ie].getSigma();

		Node[] vertexNode=model.elementNodes(ie);
		double[][] M=new double[this.nElVert][this.nElEdge];

		boolean[] edgeDir=model.element[ie].getEdgeReverse();

		double detJac,ws=1;

		Mat jac;

		Vect[] Ne=new Vect[this.nElEdge];
		Vect[] gradN=new Vect[this.nElVert];

		Vect localCo=new Vect(3);
		int n=this.PW[0].length; 
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){
					localCo.el[0]=this.PW[0][p];
					localCo.el[1]=this.PW[0][q];
					localCo.el[2]=this.PW[0][r];


					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					if(n!=2)
						ws=this.PW[1][p]*this.PW[1][q]*this.PW[1][r]*detJac;
					else
						ws=detJac;

					Ne=Ne(jac,localCo,edgeDir);


					gradN=gradN(jac,localCo);

					for(int i=0;i<this.nElVert;i++)
						for(int j=0;j<this.nElEdge;j++)	
							M[i][j]+=ws*gradN[i].dot(Ne[j].times(sigma));

				}	

		return M;
	}
	
	public double[][] PeTetra(Model model, int ie){
		
		Vect sigma=model.element[ie].getSigma();

		Node[] vertexNode=model.elementNodes(ie);
		double[][] M=new double[this.nElVert][this.nElEdge];

		boolean[] edgeDir=model.element[ie].getEdgeReverse();

		double ws=1;

		Vect lc=new Vect().ones(4).times(.25);
		

		Vect[] Ne=new Vect[this.nElEdge];


		Mat jac=jacobianTetra(vertexNode,lc);
		
		double detJ=abs(jac.determinant());

		Vect[] gradN=this.gradNTet(jac);
		
					int n=this.PWtetra.length; 	

					for(int p=0;p<n;p++)
					{

						lc.el[0]=this.PWtetra[p][0];
						lc.el[1]=this.PWtetra[p][1];
						lc.el[2]=this.PWtetra[p][2];
						lc.el[3]=this.PWtetra[p][3];
						
						Ne=NeTetra(jac,lc,edgeDir);
	
						ws=this.PWtetra[p][4]*detJ;;
						//if(ie==56){
					//	for(int i=0;i<this.nElVert;i++)
					//	vol+=Ne[i].norm();
						//util.pr("vol  "+vol);
					//	}

						for(int i=0;i<this.nElVert;i++)
							for(int j=0;j<this.nElEdge;j++)	
								M[i][j]+=ws*gradN[i].dot(Ne[j].times(sigma));

					}					

			//if(ie==56)
			//	util.pr("vol  "+vol);
			//	util.pr(ie+":    "+new Mat(M).norm());
			
			return M;
			
	}

	public Vect gradPhi(Node[] vertexNode,double[] nodePhi,Vect localCo){

		Vect gradPhi=new Vect(3);
		Vect[] gradN=new Vect[this.nElVert];
		Mat G=jacobian(vertexNode,localCo).inv3();
		Vect[] localGradN=localGradN(localCo);
		for(int i=0;i<this.nElVert;i++){
			gradN[i]=G.mul(localGradN[i]);
			gradPhi= gradPhi.add(gradN[i].times(nodePhi[i]));
		}

		return  gradPhi;	

	}
	public Vect gradPhiTet(Node[] vertexNode,double[] nodePhi){

		Vect zero=new Vect(3);
		Vect gradPhi=new Vect(3);
		Mat jac=jacobianTetra(vertexNode,zero);
		Vect[] gradN=this.gradNTet(jac);
		for(int i=0;i<this.nElVert;i++){
			gradPhi= gradPhi.add(gradN[i].times(nodePhi[i]));
		}

		return  gradPhi;	

	}



	public Mat Se(Model model,int ie){

		if(this.elCode==0) return Se3ang(model,ie);
		else if(this.elCode==1) return SeQuad(model,ie);
		else return null;
	}

	public Mat Se(Model model,int ie,Mat Me){

		if(this.elCode==0) { fillMe3angScalar(model,ie,Me); return Se3ang(model,ie);}

		if(this.elCode==1) { fillMeQuadScalar(model,ie,Me); return SeQuad(model,ie);}
		else return null;
	}



	public Mat SeQuad(Model model,int ie){


		Node[] vertexNode=model.elementNodes(ie);
		Mat Ke=new Mat(this.nElVert,this.nElVert);


		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;


		Vect kt=model.element[ie].getSigma();

		Vect[] gradN;
		Vect localCo=new Vect(dim);
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;


				gradN=gradN(jac,localCo);


				for(int i=0;i<this.nElVert;i++){
					for(int j=0;j<this.nElVert;j++){
						double BtkB=gradN[i].dot(kt.times(gradN[j]))*ws;
						Ke.el[i][j]+=BtkB;

					}
				}

			}	


		return Ke;
	}	



	public Mat gradNi_gradNjQ(Model model,int ie){


		Node[] vertexNode=model.elementNodes(ie);
		Mat Ke=new Mat(this.nElVert,this.nElVert);


		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;


		Vect[] gradN;
		Vect localCo=new Vect(dim);
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];

				jac=jacobian(vertexNode,localCo);

				detJac=abs(jac.determinant());

				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;


				gradN=gradN(jac,localCo);


				for(int i=0;i<this.nElVert;i++){
					for(int j=0;j<this.nElVert;j++){
						double BtkB=gradN[i].dot(gradN[j])*ws;
						Ke.el[i][j]+=BtkB;

					}
				}

			}	


		return Ke;
	}

	public Mat NiNjQ(Model model,int ie){

		Node[] vertexNode=model.elementNodes(ie);
		Mat Me=new Mat(this.nElVert,this.nElVert);

		double detJac,ws=1;
		Mat jac;
		Mat M1=new Mat(this.nElVert,this.nElVert);
		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW[0].length;	
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{

				lc.el[0]=this.PW[0][p];
				lc.el[1]=this.PW[0][q];

				localNe=NQuad(lc);

				for(int j=0;j<this.nElVert;j++){
					Nv.el[0][j]=localNe[j];
				}

				jac=jacobian(vertexNode,lc);
				detJac=abs(jac.determinant());

				ws=this.PW[1][p]*this.PW[1][q]*detJac;

				M1=M1.add(Nv.transp().mul(Nv).times(ws));



			}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Me.el[i][j]=M1.el[i][j];

		return Me;
	}



	public Mat Se3ang(Model model,int ie){

		Mat Ke=new Mat(this.nElVert,this.nElVert);	

		Vect kt=model.element[ie].getSigma();
		double S=el3angArea(model,ie);


		Vect[] gradN=gradN3ang(model,ie);

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){

				Ke.el[i][j]=gradN[i].dot(kt.times(gradN[j]))*S;

			}


		return Ke;
	}


	private void fillMe3angScalar(Model model,int ie,Mat Me){

		Node[] vertexNode=model.elementNodes(ie);

		double detJac,ws=1;
		Mat jac;
		double gw=model.gw;
		double he=model.getElementT(ie);
		double lam=gw*Math.exp(-.05*he);
		Mat M1=new Mat(this.nElVert,this.nElVert);

		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW3ang.length;	
		for(int p=0;p<n;p++)
		{
			lc.el[0]=this.PW3ang[p][0];
			lc.el[1]=this.PW3ang[p][1];
			localNe[0]=lc.el[0];
			localNe[1]=lc.el[1];
			localNe[2]=1-lc.el[0]-lc.el[1];


			for(int j=0;j<this.nElVert;j++)
				Nv.el[0][j]=localNe[j];

			jac=jacobian3ang(vertexNode,lc);

			detJac=abs(jac.determinant());
			ws=this.PW3ang[p][2]*detJac;
			M1=M1.add(Nv.transp().mul(Nv).times(ws));
		}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++){

				Me.el[i][j]=lam*M1.el[i][j];

			}

	}

	private void fillMeQuadScalar(Model model,int ie,Mat Me){

		Node[] vertexNode=model.elementNodes(ie);


		double detJac,ws=1;
		Mat jac;
		double gw=model.gw;
		double he=model.getElementT(ie);
		double lam=gw*Math.exp(-.05*he);
		Mat M1=new Mat(this.nElVert,this.nElVert);
		double[] localNe=new double[this.nElVert];
		Mat Nv=new Mat(1,this.nElVert);
		Vect lc=new Vect(dim);		
		int n=this.PW[0].length;	
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{

				lc.el[0]=this.PW[0][p];
				lc.el[1]=this.PW[0][q];

				localNe=NQuad(lc);

				for(int j=0;j<this.nElVert;j++){
					Nv.el[0][j]=localNe[j];
				}

				jac=jacobian(vertexNode,lc);
				detJac=abs(jac.determinant());

				ws=this.PW[1][p]*this.PW[1][q]*detJac;

				M1=M1.add(Nv.transp().mul(Nv).times(ws));



			}

		for(int i=0;i<this.nElVert;i++)
			for(int j=0;j<this.nElVert;j++)
				Me.el[i][j]=lam*M1.el[i][j];



	}

	public Vect getElementCenter(Model model,int ie){


		Vect center=new Vect(this.dim);
		int[] vertNumb=model.element[ie].getVertNumb();
		for(int j=0;j<this.nElVert;j++)
			center=center.add(model.node[vertNumb[j]].getCoord());


		return center.times(1.0/this.nElVert);
	}



	public Vect globalCo(Node[] vertex,double[]  localN){
		Vect v=new Vect(this.dim);
		for(int i=0;i<this.nElVert;i++) 
			v=v.add( vertex[i].getCoord().times(localN[i]));
		return v;
	}

	public  Vect localCo(Model model,int[] m,Vect globalCo){

		if(model.elCode==1) return localCoQuad(model,m[0],globalCo);
		Vect lc=new Vect(3);
		Vect dlc=new Vect(3);
		Vect gc=new Vect(3);
		double[] localN;

		Node[] vertex=model.elementNodes(m[0]);
		double resNorm=1;
		double error=1e-6;

		if(m[1]==0) lc.el[0]=-1;
		if(m[2]==0) lc.el[0]=1;
		if(m[3]==0) lc.el[1]=-1;
		if(m[4]==0) lc.el[1]=1;
		if(m[5]==0) lc.el[2]=-1;
		if(m[6]==0) lc.el[2]=1;

		boolean[] known=new boolean[3];
		for(int j=0;j<3;j++)
			if(lc.el[j]!=0) known[j]=true;

		for(int i=0;(i<1 && resNorm>error);i++){
			localN=N(lc);
			gc=globalCo(vertex,localN);

			Vect res=globalCo.sub(gc);

			resNorm=res.norm();
			for(int j=0;j<3;j++)
				if(!known[j])
					lc.el[j]+=dlc.el[j];
		}

		return lc;
	}

	private Vect localCoQuad(Model model,int ie,Vect globalCo){
		Vect lc=new Vect(2);
		Vect dlc=new Vect(2);
		Vect gc=new Vect(2);
		double[] localN;

		Node[] vertex=model.elementNodes(ie);
		double resNorm=1;
		double error=1e-12*model.scaleFactor;

		Vect[] v=new Vect[4];


		int[] vertNumb=model.element[ie].getVertNumb();
		for(int j=0;j<4;j++)
			v[j]=model.node[vertNumb[j]].getCoord();


		if(v[0].sub(globalCo).cross(v[1].sub(globalCo)).norm()<error) lc.el[1]=1;
		if(v[1].sub(globalCo).cross(v[2].sub(globalCo)).norm()<error) lc.el[0]=-1;
		if(v[2].sub(globalCo).cross(v[3].sub(globalCo)).norm()<error) lc.el[1]=-1;
		if(v[3].sub(globalCo).cross(v[0].sub(globalCo)).norm()<error) lc.el[0]=1;

		boolean[] known=new boolean[2];
		for(int j=0;j<2;j++)
			if(lc.el[j]!=0) known[j]=true;

		for(int i=0;(i<1 && resNorm>error);i++){
			localN=NQuad(lc);
			gc=globalCo(vertex,localN);

			Vect res=globalCo.sub(gc);

			resNorm=res.norm();
			for(int j=0;j<2;j++)
				if(!known[j])
					lc.el[j]+=dlc.el[j];
		}

		return lc;
	}

	public double[] getApAnAt(Model model,Vect globalCo){
		if(dim==3){
			throw new IllegalArgumentException(" for 3D not set up yet.");
		}

		if(model.elCode==0){
			throw new IllegalArgumentException(" for triangular element not set up yet.");
		}

		double na[]=new double[2];
		na[0]=-1e10;
		na[1]=-1e10;

		int[] m=model.getContainingElement(globalCo);
		if(m[0]<=0) return na;


		Vect	lc=localCoQuad(model,m[0],globalCo);


		return getApAnAt(model, m[0], lc);

	}

	public Vect getBAt(Model model,Vect globalCo){
		if(this.elCode==3) return new Vect(dim);

		Vect na=new Vect(this.dim);
		na.el[0]=1e10;
		na.el[1]=-1e10;
		int[] m=model.getContainingElement(globalCo);
		if(m[0]<=0) return na;

		if(this.elCode==0) return model.element[m[0]].getB();
		Vect lc;
		if(this.dim==3){
			lc=localCo(model,m,globalCo);
		}
		else{
			lc=localCoQuad(model,m[0],globalCo);

		}
		return getBAt(model, m[0], lc);

	}

	public double[] getApAnAt(Model model,int i,Vect lc){

		Edge[] elEdges=model.elementEdges(i);

		double[] Ne=NeQuad(lc);


		double[] ApAn=new double[3];

		for(int j=0;j<this.nElEdge;j++)	{	
			ApAn[0]= ApAn[0]+Ne[j]*elEdges[j].Ap;
			ApAn[1]= ApAn[1]+Ne[j]*elEdges[j].A;
		}


		return 	ApAn;
	}


	public Vect getBAt(Model model,int i,Vect lc){
		if(model.fluxLoaded || model.elCode==0) return model.element[i].getB();
		else if(this.elCode==3) model.element[i].getB();

		boolean[] edgeDir=model.element[i].getEdgeReverse();

		Vect B=new Vect(this.dim);
		Mat jac=new Mat(this.dim,this.dim);
		Node[] vertex=model.elementNodes(i);
		Edge[] elEdges=model.elementEdges(i);
		Vect[] rotNe=new Vect[this.nElEdge];
		if(this.dim==3){
			jac=jacobian(vertex,lc);
			rotNe=rotNe(jac,lc,edgeDir);
		}
		else if(this.elCode==1){
			jac=jacobian(vertex,lc);

			rotNe=rotNeQuad(jac,lc);

		}


		for(int j=0;j<this.nElEdge;j++)	{	
			B= B.add(rotNe[j].times(elEdges[j].A));
		}

		if(model.axiSym){
			double rr=0;
			int L=vertex.length;
			for(int j=0;j<L;j++){
				double r1=vertex[j].getCoord(0);
				rr+=1.0/L/r1;
			}
			B=B.times(rr);
			//B.show();
		}

		return B;
	}

	public double[][] Te(Model model,int ie){



		boolean[] edgeDir=model.element[ie].getEdgeReverse();


		int n=this.PWNL[0].length; 


		Node[] vertexNode=model.elementNodes(ie);
		double[][] H=new double[model.nElEdge][model.nElEdge];


		double detJac,ws=1,wsJ=0,term1;

		Mat jac;


		Vect[] rotNe=new Vect[model.nElEdge];
		Vect localCo=new Vect(model.dim);


		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
				for(int r=0;r<n;r++){

					localCo.el[0]=this.PWNL[0][p];
					localCo.el[1]=this.PWNL[0][q];
					localCo.el[2]=this.PWNL[0][r];


					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					wsJ=ws*detJac;


					rotNe=rotNe(jac,localCo,edgeDir);


					for(int i=0;i<model.nElEdge;i++){

						for(int j=0;j<=i;j++){

							term1=wsJ*rotNe[i].dot(rotNe[j]);
							H[i][j]+=term1;


						}
					}
				}

		lowSym(H);


		return H;


	}


	public Mat SePCQuad(Model model,int ie,Mat Te,Mat Pe,double kw){



		Node[] vertexNode=model.elementNodes(ie);
		Mat Se=new Mat(this.nElVert,this.nElVert);


		double detJac,ws=1;
		int n=this.PW[0].length;	
		Mat jac;

		Vect[] gradN;
		double[] Nshp;
		Vect lc=new Vect(dim);

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++)
			{
				lc.el[0]=this.PW[0][p];
				lc.el[1]=this.PW[0][q];

				jac=jacobian(vertexNode,lc);

				detJac=abs(jac.determinant());

				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q]*detJac;
				else
					ws=detJac;


				gradN=gradN(jac,lc);
				Nshp=NQuad(lc);

				for(int i=0;i<this.nElVert;i++){
					for(int j=0;j<this.nElVert;j++){
						double BtkB=gradN[i].dot(gradN[j])*ws;
						Se.el[i][j]+=BtkB;
						Te.el[i][j]+=Nshp[i]*Nshp[j]*ws;

					}
				}

			}	


		return Se;
	}




	private void lowSym(double[][] H) {
		for(int i=0;i<H.length;i++)
			for(int j=i+1;j<H[0].length;j++)
				H[i][j]= H[j][i];
	}
	
	
	public double obtainElementLoss(Model model,int ie){
	
		if(model.elCode==0) return 0;
		if(model.elCode==2) return obtainElementLossTet(model,ie);
		
		double loss=0;
		
		int n;

		n=this.PW[0].length; 

		Node[] vertexNode=model.elementNodes(ie);
		
		
		double detJac,ws=1,wsJ=0;

		Mat jac;

		Vect localCo=new Vect(this.dim);


		Vect Je=new Vect(this.dim);

		Vect sigmaInv=model.element[ie].getSigma().inv();

		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];
				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q];
				else
					ws=1;
				
				if(model.dim==3)
				for(int r=0;r<n;r++){

					localCo.el[2]=this.PW[0][r];

					if(n!=2)
						ws*=this.PW[1][r];
					
					jac=jacobian(vertexNode,localCo);

					detJac=abs(jac.determinant());

					wsJ=ws*detJac;
					
					Je=model.getElementJe(ie, localCo);
							
					loss+=Je.dot(sigmaInv.times(Je))*wsJ;
				}
			}
	
		return loss;
	}
	
	public double obtainElementEnergy(Model model,int ie){
		if(model.elCode==0) return 0;
		if(model.elCode==2) obtainElementEnergyTet(model,ie);
		
		double energy=0;
		
		int n;

		n=this.PW[0].length; 

		Node[] vertexNode=model.elementNodes(ie);
		
		
		double detJac,ws=1,wsJ=0;

		Mat jac;

		Vect localCo=new Vect(this.dim);


		Vect B=null;

		Vect nu=model.element[ie].getNu();

		
		for(int p=0;p<n;p++)
			for(int q=0;q<n;q++){
				localCo.el[0]=this.PW[0][p];
				localCo.el[1]=this.PW[0][q];
				if(n!=2)
					ws=this.PW[1][p]*this.PW[1][q];
				else
					ws=1;
				if(model.dim==3)
				for(int r=0;r<n;r++){

					localCo.el[2]=this.PW[0][r];

					if(n!=2)
						ws*=this.PW[1][r];
			
					
					jac=jacobian(vertexNode,localCo);


					detJac=abs(jac.determinant());

					wsJ=ws*detJac;
					B=model.getElementB(ie, localCo);

							
					energy+=B.dot(nu.times(B))*wsJ;
				}
	}
	
			energy/=2;
		
		return energy;
	}
	
	public double obtainElementEnergyTet(Model model,int ie){


		double energy=0;
		
		Node[] vertexNode=model.elementNodes(ie);
		
		Vect lc=new Vect(4);

		double wsJ;
		
		Vect B=null;

		Vect nu=model.element[ie].getNu();
		

		Mat jac=this.jacobianTetra(vertexNode, new Vect(dim));

		double detJac=abs(jac.determinant());
		
		int n=this.PWtetra.length; 					
		for(int p=0;p<n;p++)
		{

			lc.el[0]=this.PWtetra[p][0];
			lc.el[1]=this.PWtetra[p][1];
			lc.el[2]=this.PWtetra[p][2];
			lc.el[3]=this.PWtetra[p][3];
			wsJ=this.PWtetra[p][4]*detJac;
			
			B=model.getElementB(ie, lc);

			
			energy+=B.dot(nu.times(B))*wsJ;

			

		}
		
		energy/=2;
		
		return energy;
	
	}

	
	public double obtainElementLossTet(Model model,int ie){


		double loss=0;
		
		Node[] vertexNode=model.elementNodes(ie);
		
		Vect lc=new Vect(4);

		double wsJ;
		
		Vect Je=new Vect(3);
		
		Vect sigmaInv=model.element[ie].getSigma().inv();
		

		Mat jac=this.jacobianTetra(vertexNode, new Vect(dim));

		double detJac=abs(jac.determinant());
		
		int n=this.PWtetra.length; 					
		for(int p=0;p<n;p++)
		{

			lc.el[0]=this.PWtetra[p][0];
			lc.el[1]=this.PWtetra[p][1];
			lc.el[2]=this.PWtetra[p][2];
			lc.el[3]=this.PWtetra[p][3];
			
			Je=model.getElementJe(ie, lc);


			wsJ=this.PWtetra[p][4]*detJac;
			
			loss+=Je.dot(sigmaInv.times(Je))*wsJ;

		}
		


		return loss;
	
	}
	

	public double[][] gaussInteg(int n){
		double[] pg;
		double[] wg;
		if(n==1) {pg=new double[1]; wg=new double[1]; wg[0]=2;}
		else if(n==2) {double a=1.0/sqrt(3.0); pg=new double[2]; pg[0]=-a;pg[1]=a; wg=new double[2]; wg[0]=1;wg[1]=1;}
		else if(n==3) {double a=sqrt(3.0/5.0); pg=new double[3]; pg[0]=-a;pg[1]=0; pg[2]=a;wg=new double[3]; wg[0]=5.0/9;wg[1]=8.0/9;wg[2]=5.0/9;}
		else if(n==4) {double a1=sqrt((3.0-2*sqrt(6.0/5))/7.0); double a2=sqrt((3.0+2*sqrt(6.0/5))/7.0);
		pg=new double[4]; pg[0]=-a2;pg[1]=-a1; pg[2]=a1; pg[3]=a2;wg=new double[4]; double w1=(18-sqrt(30))/36; double w2=(18+sqrt(30))/36;
		wg[0]=w1;wg[1]=w2;wg[2]=w2;wg[3]=w1;}
		else if(n==5) {double a1=sqrt(5.0-2*sqrt(10.0/7))/3.0; double a2=sqrt(5.0+2*sqrt(10.0/7))/3.0;
		pg=new double[5]; pg[0]=-a2;pg[1]=-a1; pg[2]=0;pg[3]=a1; pg[4]=a2;wg=new double[5]; double w1=(322-13*sqrt(70))/900; double w2=(322+13*sqrt(70))/900;
		wg[0]=w1;wg[1]=w2; wg[2]=128.0/225; wg[3]=w2;wg[4]=w1;}
		else {  throw new NullPointerException("Gauss integration for degree "+ n+" is not defined in this software");};		

		double[][] WP=new double[2][wg.length];
		WP[0]=pg;
		WP[1]=wg;

		return WP;
	}

	public double[][] gaussInteg3(int n){

		double[][] PW=new double[n][3];

		if(n==1) {
			PW[0][0]=1.0/3;
			PW[0][1]=1.0/3;
			PW[0][2]=.5;
		}
		else if(n==3) { 
			PW[0][0]=.5; 
			PW[0][1]=.5;
			PW[1][0]=.5; 
			PW[1][1]=0; 
			PW[2][0]=0; 
			PW[2][1]=0.5; 
			double a=1.0/6;
			for(int i=0;i<3;i++)
				PW[i][2]=a;
		}

		else if(n==4) { 
			PW[0][0]=1./3; 
			PW[0][1]=1./3;
			PW[1][0]=.2; 
			PW[1][1]=.6; 
			PW[2][0]=.2; 
			PW[2][1]=0.2; 
			PW[3][0]=.6; 
			PW[3][1]=0.2;

			PW[0][2]=-27./96;

			double a=25./96;
			for(int i=1;i<4;i++)
				PW[i][2]=a;
		}

		else if(n==7) { 
			double c1=1.0/3;
			double c2=1.0/15;
			double c3=1.0/40;

			PW[0][0]=c1; 
			PW[0][1]=c1;
			PW[1][0]=0.5; 
			PW[1][1]=0; 
			PW[2][0]=0.5; 
			PW[2][1]=0.5; 
			PW[3][0]=0; 
			PW[3][1]=0.5;
			PW[4][0]=1; 
			PW[4][1]=0; 
			PW[5][0]=0; 
			PW[5][1]=1; 
			PW[6][0]=0; 
			PW[6][1]=0;

			PW[0][2]=9*c3;


			for(int i=1;i<4;i++)
				PW[i][2]=c2;
			for(int i=4;i<7;i++)
				PW[i][2]=c3;
		}
		else {  throw new NullPointerException("Gauss integration for degree "+ n+" is not defined in this software");};		


		return PW;
	}

	public double[][] gaussIntegPyramid(int n){
	
		double[][] PW=new double[n][4];

		if(n==4){
		int ix=0;
		
		PW[ix][0]=0.750000;
		PW[ix][1]=0;
		PW[ix][2]=1./2*(-0.500000000000000+1);
		PW[ix][3]=0.699453551912568;
		
		ix++;
		PW[ix][0]=-0.266666666666667;
		PW[ix][1]=0.700991361493770;
		PW[ix][2]=1/2*(-0.500000000000000+1);
		PW[ix][3]=0.699453551912568;
		
		
		ix++;
		PW[ix][0]=-0.266666666666667;
		PW[ix][1]=-0.386753854617252;
		PW[ix][2]=1/2*(+0.131034482758621+1);
		PW[ix][3]=0.560442489670798;
		
		ix++;
		PW[ix][0]=-0.266666666666667;
		PW[ix][1]=-0.386753854617252;
		PW[ix][2]=1/2*(-1.000000000000000+1);
		PW[ix][3]=0.699453551912568;
		
	}
		else if(n==8){
			

	 double upyr8[] = {0.2631840555694285,-0.2631840555694285,
				  0.2631840555694285,-0.2631840555694285,
				  0.5066163033492386,-0.5066163033492386,
				  0.5066163033492386,-0.5066163033492386};
	 double vpyr8[] = {0.2631840555694285,0.2631840555694285,
				  -0.2631840555694285,-0.2631840555694285,
				  0.5066163033492386,0.5066163033492386,
				  -0.5066163033492386,-0.5066163033492386};
	 double wpyr8[] = {0.544151844011225,0.544151844011225,
				  0.544151844011225,0.544151844011225,
				  0.122514822655441,0.122514822655441,
				  0.122514822655441,0.122514822655441};
	 
	 double ppyr8[] = {0.100785882079825,0.100785882079825,
				  0.100785882079825,0.100785882079825,
				  0.232547451253508,0.232547451253508,
				  0.232547451253508,0.232547451253508};
	 
	 for(int i=0;i<n;i++){
		 PW[i][0]=upyr8[i];
		 PW[i][1]=vpyr8[i];
		 PW[i][2]=wpyr8[i];
		 PW[i][3]=ppyr8[i];
		 }
	}


		return PW;
	}
	public double[][] gaussIntegTetra(int n){


		double[][] PW=new double[n][5];

		if(n==1) {
			double a=1.0/4;
			PW[0][0]=a;
			PW[0][1]=a;
			PW[0][2]=a;
			PW[0][3]=a;
			PW[0][4]=1.0/6;
		}
		else if(n==4) { 
			double a=(5+3*sqrt(5))/20;
			double b=(5-sqrt(5))/20;
			double w=1.0/24;

			int ix=0;
			PW[ix][0]=a; 
			PW[ix][1]=b;
			PW[ix][2]=b;
			PW[ix][3]=b;

			ix++;
			PW[ix][0]=b; 
			PW[ix][1]=a;
			PW[ix][2]=b;
			PW[ix][3]=b;
			
			ix++;
			PW[ix][0]=b; 
			PW[ix][1]=b;
			PW[ix][2]=a;
			PW[ix][3]=b;
			
			ix++;
			PW[ix][0]=b; 
			PW[ix][1]=b;
			PW[ix][2]=b;
			PW[ix][3]=a;

			for(int k=0;k<4;k++)
				PW[k][4]=w;

		}

		else if(n==5) { 

			double a=.25;
			double b=1.0/6;
			double w=3.0/40;

			int ix=0;
			
	
			PW[ix][0]=a; 
			PW[ix][1]=a;
			PW[ix][2]=a;
			PW[ix][3]=a;
			
			ix++;
			PW[ix][0]=.5; 
			PW[ix][1]=b;
			PW[ix][2]=b;
			PW[ix][3]=b;
			
			ix++;
			PW[ix][0]=b; 
			PW[ix][1]=.5;
			PW[ix][2]=b;
			PW[ix][3]=b;
			
			
			ix++;
			PW[ix][0]=b; 
			PW[ix][1]=b;
			PW[ix][2]=.5;
			PW[ix][3]=b;
			
			
			ix++;
			PW[ix][0]=b; 
			PW[ix][1]=b;
			PW[ix][2]=b;
			PW[ix][3]=.5;
			
			for(int k=0;k<5;k++)
				if(k==0)
					PW[k][4]=-4.0/30;
				else
					PW[k][4]=w;

		}

		else {  throw new NullPointerException("Gauss integration for degree "+ n+" is not defined in this software");};		

		return PW;
	}

}
