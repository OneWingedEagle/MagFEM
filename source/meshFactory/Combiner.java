package meshFactory;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import io.Writer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import fem.BoundarySet;
import fem.Calculator;
import fem.Element;
import fem.Force;
import main.Main;
import fem.Model;
import math.Mat;
import math.Vect;
import math.util;



public class Combiner {


	BoundarySet bs=new BoundarySet();

	public static void mainxxxx(String[] args){


	}


	public void combineOneNode(Model motor,String rotf, String statf,double ang,double step){
		combineTwoNode(motor,rotf,statf,ang,step,true);

		for(int ir=1;ir<=motor.numberOfRegions;ir++){
			if(!motor.region[ir].rotor) continue;

			for( int i=motor.region[ir].getFirstEl();i<=motor.region[ir].getLastEl();i++)
			{
				int[] vertNumb=motor.element[i].getVertNumb();

				for(int j=0;j<motor.nElVert;j++){
					if(motor.node[vertNumb[j]].common){
						motor.element[i].setVertNumb(j,motor.node[vertNumb[j]].getMap());
					}
				}
			}

		}		
		
	}

	public void combineTwoNode(Model motor,String rotf, String statf,double ang){
		combineTwoNode(motor,rotf,statf,ang,1,false);
	}

	public void combineTwoNode(Model motor,String rotf, String statf)
	{
		combineTwoNode(motor,rotf,statf,0);
	}

	public void combineTwoNode(Model motor,String rotf, String statf,double ang,double step,boolean forWriteMesh){}


	public int[] getBorderNodeSorted(Model model,int nb){

		return bs.getBorderNodeSorted(model,nb);

	}



	public void setSliceBounds(Model model)
	{
		bs.setSliceBounds(model);
	}


	// set indices for rotation coordinates change
	
/*	public void setRotorPosition(Model motor, double newAng)
	{
		setRotorPosition(motor,newAng,motor.meshAngStep);
	}
	
	public void setRotorPosition(Model motor, double newAng, double angStep0){
		if(abs(newAng)>=360) newAng=newAng%360;
		if(newAng<0) newAng=360+newAng;
		
			
		int L=motor.commonNodes[0].length;


		boolean full=motor.fullMotor;
		int cpb=motor.cpb;
		int cpb2=cpb*cpb;
		 
	if(full){setRotorPositionFull(motor,newAng,angStep0); return;}
		
		double dangNew=newAng-motor.rotAng;
		
		if(dangNew==0) return;

		int nSteps=(int)Math.round(newAng/angStep0);

		double diff=newAng-nSteps*angStep0;

		
		double angComm=-diff+motor.rotSubAng;
		


		motor.rotSubAng=diff;
		

		motor.rotAng=newAng;
		double dangx=dangNew;
		
		int[] commonNodeIndSorted=new int[motor.numberOfNodes+1];
		int Ld=L-1;

		for(int i=0;i<L;i++)
			commonNodeIndSorted[motor.commonNodes[0][i]]=i+1;
		
		Mat R=util.rotMat2D(dangx*Math.PI/180);
		Mat Rcomm=R.mul(util.rotMat2D(angComm*Math.PI/180));
		
		
	//	double rx=1.0001*motor.node[908].getR();

		
		boolean[] nodeRotated=new boolean[motor.numberOfNodes+1];

		for(int ir=1;ir<=motor.numberOfRegions;ir++)
			if(motor.region[ir].rotor)
				for(int i=motor.region[ir].getFirstEl();i<=motor.region[ir].getLastEl();i++){

					Vect Mr=R.mul(motor.element[i].getM());
					if(motor.element[i].hasM())
					motor.element[i].setM(Mr);

					int[] vertNumb=motor.element[i].getVertNumb();
					for(int j=0;j<motor.nElVert;j++){

						int nn=vertNumb[j];

						if(nodeRotated[nn]) continue;
						nodeRotated[nn]=true;
						if(!motor.node[nn].common){
							
							if(ir==7  && motor.node[nn].getR()>rx )
							motor.node[nn].setCoord(Rcomm.mul(motor.node[nn].getCoord()));	
						
							else if(ir==7  && motor.node[nn].getR()>.999*rx ){
								 
								double a=(motor.node[nn].getR()-motor.node[816].getR())/.0003;
								Mat Rx=util.rotMat2D((dangx+a*angComm)*Math.PI/180);
								motor.node[nn].setCoord(Rx.mul(motor.node[nn].getCoord()));	

							}
							else 
								motor.node[nn].setCoord(R.mul(motor.node[nn].getCoord()));	
							
								
					}

						else{
							int ind=commonNodeIndSorted[nn];
							int k=nSteps+ind-1;
						
							if(k<L){
								int nx=motor.commonNodes[1][k];
								motor.node[nn].setMap(nx);	
								motor.node[nn].setCoord(motor.node[nx].getCoord());

							
								if(k==Ld){
									motor.node[nn].setPBC(motor.cpb);
									if(motor.node[nn].hasPBC())
										motor.node[nn].setMap(motor.commonNodes[1][0]);
							
								}
							}
							else
							{
								
								
								int kx=k/Ld;
								int p=(k)%Ld;

								if(kx%2==1){
									motor.node[nn].setPBC(cpb);
								}
								else{
									motor.node[nn].setPBC(cpb2);

								}								
			
								motor.node[nn].setCoord(Rcomm.mul(motor.node[nn].getCoord()));	
								
								motor.node[nn].common=false;
								motor.node[motor.node[nn].getMap()].common=false;

								motor.node[nn].setMap(0);

								if(motor.hasPBC){
									int nmp=motor.commonNodes[1][p];
									int nmp2=motor.node[nmp].getMap();
									if(nmp2==0)
									motor.node[nn].setMap(nmp);
									else
									motor.node[nn].setMap(motor.node[nmp].getMap());	
									
								}
						}
							
						
							}
					}
				}
		
		


		
	}
	

	
	// set indices for rotation but no coordinate changes
	
	public void setRotorIndex(Model motor, int steps){
		
		motor.rotIndex=true;

		
		boolean full=motor.fullMotor;
		int cpb=motor.cpb;
		int cpb2=cpb*cpb;

	if(full){setRotorIndexFull(motor,steps); return;}
		
		int nstep=steps;
		
		int L=motor.commonNodes[0].length;
		

	
		int[] commonNodeIndSorted=new int[motor.numberOfNodes+1];
		
		int Ld=L-1;

		for(int i=0;i<L;i++)
			commonNodeIndSorted[motor.commonNodes[0][i]]=i+1;


		boolean[] nodeRotated=new boolean[motor.numberOfNodes+1];

		for(int ir=1;ir<=motor.numberOfRegions;ir++)
			if(motor.region[ir].rotor)
				for(int i=motor.region[ir].getFirstEl();i<=motor.region[ir].getLastEl();i++){


					int[] vertNumb=motor.element[i].getVertNumb();
					for(int j=0;j<motor.nElVert;j++){

						int nn=vertNumb[j];

						if(!motor.node[nn].common) continue;
						
						if(nodeRotated[nn]) continue;
						
							int ind=commonNodeIndSorted[nn];
							int k=nstep+ind-1;
							
							if(k<L){
								int nx=motor.commonNodes[1][k];
								motor.node[nn].setMap(nx);	
								if(k==Ld){
									motor.node[nn].setPBC(motor.cpb);
									if(motor.node[nn].hasPBC())
										motor.node[nn].setMap(motor.commonNodes[1][0]);
							
								}
							}
							else
							{
																
								int kx=k/Ld;
								int p=(k)%Ld;
								
								
								if(kx%2==1){
									motor.node[nn].setPBC(cpb);
								}
								else{
									motor.node[nn].setPBC(cpb2);

								}
								
								motor.node[nn].common=false;
								motor.node[motor.node[nn].getMap()].common=false;

								motor.node[nn].setMap(0);

								if(motor.hasPBC){
									int nmp=motor.commonNodes[1][p];
									int nmp2=motor.node[nmp].getMap();
									if(nmp2==0)
									motor.node[nn].setMap(nmp);
									else
									motor.node[nn].setMap(motor.node[nmp].getMap());	
									
								}
						}
							
						
					}
				}
		
		
	}
	
	public void setRotorPositionFull(Model motor, double newAng,double angStep0){


		double dangNew=newAng-motor.rotAng;
		
		if(dangNew==0) return;

		int nSteps=(int)Math.round(newAng/angStep0);

		double diff=newAng-nSteps*angStep0;

		motor.rotSubAng=diff;

		motor.rotAng=newAng;
		double dangx=dangNew;

		int L=motor.commonNodes[0].length;
		int[] commonNodeIndSorted=new int[motor.numberOfNodes+1];

		for(int i=0;i<L;i++)
			commonNodeIndSorted[motor.commonNodes[0][i]]=i+1;
		Mat R=util.rotMat2D(dangx*Math.PI/180);

		int Ld=L;
	
		boolean[] nodeRotated=new boolean[motor.numberOfNodes+1];

		for(int ir=1;ir<=motor.numberOfRegions;ir++)
			if(motor.region[ir].rotor)
				for(int i=motor.region[ir].getFirstEl();i<=motor.region[ir].getLastEl();i++){

					Vect Mr=R.mul(motor.element[i].getM());
					if(motor.element[i].hasM())
					motor.element[i].setM(Mr);

					int[] vertNumb=motor.element[i].getVertNumb();
					for(int j=0;j<motor.nElVert;j++){

						int nn=vertNumb[j];

						if(nodeRotated[nn]) continue;

						if( motor.hasTwoNodeNumb){
							motor.node[nn].setCoord(R.mul(motor.node[nn].getCoord()));	
							nodeRotated[nn]=true;
							if(motor.node[nn].common)
							{
								int ind=commonNodeIndSorted[nn];
								int k=nSteps+ind-1;
								int p=(k)%Ld;						
								motor.node[nn].setMap(motor.commonNodes[1][p]);
								

							}
						}

					}
				}
	}
		
	
	
	public void setGearsPositionFull(Model motor, double newAng,double angStep0){


		double dangNew=newAng-motor.rotAng;
		
		if(dangNew==0) return;

		int nSteps=(int)Math.round(newAng/angStep0);

		double diff=newAng-nSteps*angStep0;

		motor.rotSubAng=diff;

		motor.rotAng=newAng;
		double dangx=dangNew;

		int L=motor.commonNodes[0].length;
		int[] commonNodeIndSorted=new int[motor.numberOfNodes+1];

		for(int i=0;i<L;i++)
			commonNodeIndSorted[motor.commonNodes[0][i]]=i+1;
		
		Mat R=util.rotMat2D(dangx*Math.PI/180);

		int Ld=L;
	
		boolean[] nodeRotated=new boolean[motor.numberOfNodes+1];

		for(int ir=1;ir<=motor.numberOfRegions;ir++)
			if(motor.region[ir].rotor)
				for(int i=motor.region[ir].getFirstEl();i<=motor.region[ir].getLastEl();i++){

					Vect Mr=R.mul(motor.element[i].getM());
					if(motor.element[i].hasM())
					motor.element[i].setM(Mr);

					int[] vertNumb=motor.element[i].getVertNumb();
					for(int j=0;j<motor.nElVert;j++){

						int nn=vertNumb[j];

						if(nodeRotated[nn]) continue;

						if( motor.hasTwoNodeNumb){
							motor.node[nn].setCoord(R.mul(motor.node[nn].getCoord()));	
							nodeRotated[nn]=true;
							if(motor.node[nn].common)
							{
								int ind=commonNodeIndSorted[nn];
								int k=nSteps+ind-1;
								int p=(k)%Ld;						
								motor.node[nn].setMap(motor.commonNodes[1][p]);
								

							}
						}

					}
				}
	}
		

	
	public void setRotorIndexFull(Model motor, int steps){
		

		int L=motor.commonNodes[0].length;
		int[] commonNodeIndSorted=new int[motor.numberOfNodes+1];

		for(int i=0;i<L;i++)
			commonNodeIndSorted[motor.commonNodes[0][i]]=i+1;

		int Ld=L;
	
		boolean[] nodeRotated=new boolean[motor.numberOfNodes+1];
		for(int ir=1;ir<=motor.numberOfRegions;ir++)
			if(motor.region[ir].rotor)
				for(int i=motor.region[ir].getFirstEl();i<=motor.region[ir].getLastEl();i++){

					
					int[] vertNumb=motor.element[i].getVertNumb();
					for(int j=0;j<motor.nElVert;j++){

						int nn=vertNumb[j];

						if(nodeRotated[nn]) continue;

						if(!motor.node[nn].common) continue;
				
						
							nodeRotated[nn]=true;
					
								int ind=commonNodeIndSorted[nn];
								int k=steps+ind-1;
								int p=(k)%Ld;						
								motor.node[nn].setMap(motor.commonNodes[1][p]);
						
						
					}
				}
	}

	public void combineFull(Model motor,String rotf, String statf, double angx){
	
		angx=angx-(int)(angx/360)*360;
		Model rotor=new Model(rotf);

		if(angx!=0)
			rotor.rotate(angx);
		
		Model stator=new Model(statf);


		double rmax=0;

		int nRotNodes=rotor.numberOfNodes;

		boolean[] rotNode=new boolean[nRotNodes+1];
		for(int i=1;i<=rotor.numberOfElements;i++)
			for(int j=0;j<rotor.nElVert;j++)
				rotNode[rotor.element[i].getVertNumb(j)]=true;


		for(int i=1;i<=nRotNodes;i++)
		{
			if(!rotNode[i]) continue;
			Vect z=rotor.node[i].getCoord();
			double r=new Vect(z.el[0],z.el[1]).norm();

			if(r>rmax) rmax=r;

		}

		int ix=0;
		boolean[] commonNode=new boolean[nRotNodes+1];
		for(int i=1;i<=nRotNodes;i++){

			if(!rotNode[i]) continue;

			if(rotor.node[i].getCoord().norm()>rmax*.99999) 
			{
				ix++;
				commonNode[i]=true;
			}

		}

		int L=ix;
		int ang=(int)(Math.round(angx)+1e-6);

		int[] commonNodeNumb=new int[L+1];
		int[] commonNodeInd=new int[nRotNodes+1];
		Vect angs=new Vect(L+1);
		angs.el[0]=-1000;
		ix=0;
		for(int i=1;i<=nRotNodes;i++)
			if(commonNode[i]) {
				ix++;
				commonNodeInd[i]=ix;
				commonNodeNumb[ix]=i;
				angs.el[ix]=util.getAng(rotor.node[i].getCoord());
			}

		int[] indSorted=angs.bubble();
angs.show();

		int[] BorderNode=new int[L+1];
		for(int i=1;i<=L;i++)
			BorderNode[i]= commonNodeNumb[indSorted[i]];


		int[] commonNodeIndSorted=new int[nRotNodes+1];

		for(int i=1;i<=L;i++)
			commonNodeIndSorted[BorderNode[i]]=i;

		Mat R=util.rotMat2D(ang*Math.PI/180);
		boolean[] nodeCounted=new boolean[nRotNodes+1];

		for(int i=1;i<=rotor.numberOfElements;i++){
			int[] vertNumb=rotor.element[i].getVertNumb();
			for(int j=0;j<rotor.nElVert;j++){

				if(nodeCounted[vertNumb[j]]) continue;

				int ind=commonNodeIndSorted[vertNumb[j]];
				if(ind==0){
					rotor.node[vertNumb[j]].setCoord(R.mul(rotor.node[vertNumb[j]].getCoord()));	
					nodeCounted[vertNumb[j]]=true;
				}

				else {	
					int k=ang+ind;
					if(k<=L)
						rotor.element[i].setVertNumb(j,BorderNode[k]);
					else{
						rotor.element[i].setVertNumb(j,BorderNode[k-L]);
					}

				}


			}
		}




		int nRegions1=rotor.numberOfRegions;
		int nRegions2=stator.numberOfRegions;
		int nElements1=rotor.numberOfElements;
		int nElements2=stator.numberOfElements;
		int nNodes1=rotor.numberOfNodes;
		int nNodes2=stator.numberOfNodes;



		int nRegions=nRegions1+nRegions2;
		int nElements=nElements1+nElements2;

		motor.alloc(nRegions,nElements,nNodes2,rotor.elType);

		for(int i=1;i<=nRegions1;i++){
			motor.region[i]=rotor.region[i];
			motor.region[i].setMaterial("region"+i+"R");
		}

		for(int i=1;i<=nRegions2;i++)	{

			motor.region[i+nRegions1].setMaterial("region"+i+"S");
			motor.region[i+nRegions1].setFirstEl(stator.region[i].getFirstEl()+nElements1);
			motor.region[i+nRegions1].setLastEl(stator.region[i].getLastEl()+nElements1);
		}


		for(int i=1;i<=nElements1;i++)	{
			motor.element[i]=rotor.element[i];
			motor.element[i].setVertNumb(rotor.element[i].getVertNumb());
		}

		ix=0;
		for(int i=1;i<=nElements2;i++)	{
			motor.element[i+nElements1]=stator.element[i];


		}

		for(int i=1;i<=nNodes1;i++)	{

			if(rotor.node[i].getCoord().norm()==0 && stator.node[i].getCoord().norm()>0)
				motor.node[i]=stator.node[i];
			else{
				motor.node[i].setCoord(rotor.node[i].getCoord());	
			}
		}



		motor.scaleFactor=rotor.scaleFactor;
		motor.rotAng=ang;
		motor.setElType(rotor.elType);
		motor.setFemCalc();
		motor.setForceCalc();



	}


	



*/

	public void setNodeOnMotorBound(Model model){

		bs.setNodeOnMotorBound(model);

	}

	public int[][] getBoundaryNodesMap(Model model,int nb1,int nb2){

		return  bs.mapBorderNodes(model,nb1,nb2);

	}


}
