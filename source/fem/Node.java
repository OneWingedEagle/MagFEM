package fem;
import math.Mat;
import math.Vect;
import math.util;

public class Node {

	public Vect F,Fms,Fstr,u;
	private Vect coord;
	private double phi,nodalMass;
	public double stress,T;
	private boolean[] uKnown;
	public boolean[] onBound;
	private boolean deformable,phiKnown,phiVar,hasF,hasFms,Tknown;
	public boolean common,sPBC,aPBC,inUse,rotor;
	private byte dim;
	private int map;
	public int id;

	public Node(int id1,int dim)
	{
		this.id=id1;
	this.dim=(byte)dim;
	coord=new Vect(dim);
	uKnown=new boolean[dim];
	if(dim==2) onBound=new boolean[4];
	}

	public void setPhi(double phi){
		this.phi=phi;
	}
		
	public double getPhi(){
		return this.phi;
		
	}
	
	public void setNodalMass(double mass){
		nodalMass=mass;		
	}
	
	public double getNodalMass(){
		return this.nodalMass;
		
	}
	
	public Vect getNodalVect(int k){

		if(k==-1)
		{
		if(u==null) return null;
		else
		return this.u.deepCopy();
		}
		else if(k==1)
		{
			
			if(F==null) return null;
			else{
			return this.F.deepCopy();
			}
		
		}
		else if(k==2) {
			if(Fms==null) return null;
			else
			return this.Fms.deepCopy();
		}
		else if(k==3){
			if(F==null && Fms==null) return null;
			else
			if(F==null) return Fms;
			else
			if(Fms==null) return F;	
			else{

				return this.F.add(this.Fms);
	}
		}
		else
			return null;
		
	}


	
	public void setDeformable(boolean b){
		this.deformable=b;
		if(b){
			//F=new Vect(dim);
			u=new Vect(dim);
		//	Fms=new Vect(dim);
		}
		else
		{
			F=null;
			u=null;
			Fms=null;
		}
	}

	public boolean isDeformable(){
		return this.deformable;
	}
	
	public void setPhiVar(boolean b){
		this.phiVar=b;
	}

	public boolean isPhiVar(){
		return this.phiVar;
	}

	public void setPhiKnown(boolean b){
		this.phiKnown=b;
	}

	public boolean isPhiKnown(){
		return this.phiKnown;
	}

	public void setHasF(boolean b){
		 this.hasF=b;
	}
	public boolean hasF(){
		return this.hasF;
	}
	
	public void setHasFms(boolean b){
		 this.hasFms=b;
	}
	public boolean hasFms(){
		return this.hasFms;
	}
	
	public void setU_is_known(boolean b){
		for(int k=0;k<dim;k++)
		uKnown[k]=b;
	}
	
	
	public void setU_is_known(int i,boolean b){
		uKnown[i]=b;
	}
	
	public boolean is_U_known(){

		for(byte k=0;k<dim;k++)
			if(!uKnown[k]) return false;
		
		return true;
	}
	
	public boolean is_T_known(){

		return this.Tknown;
	}
	
	public void setKnownT(double T){

		this.Tknown=true;
		this.T=T;
	}
	
	public boolean is_A_known(){

		for(byte k=0;k<dim;k++)
			if(!uKnown[k]) return false;
		
		return true;
	}

	public boolean has_U_known(){
		for(byte k=0;k<dim;k++)
			if(uKnown[k]) return true;
		return false;
	}
	
	public boolean is_U_known(int i){
		return uKnown[i];
	}
	
	
	public void setKnownU(Vect u){
		this.u=u.deepCopy();
		setU_is_known(true);
	}
	

	
	public void setKnownU(int k,double a){
		this.u.el[k]=a;
		setU_is_known(k,true);
	}
	
	public void setKnownU(double a, double b, double c){
		setKnownU(new Vect(a,b,c));
		
	}
	
	
	
	public void setU(Vect u){
		this.u=u.deepCopy();
	}

	
	public void setU(int k,double a){
		this.u.el[k]=a;

	}


	public Vect getU(){
		return	this.u.deepCopy();

		}

	
	public double getU(int k){
		return	this.u.el[k];

		}
	
	public void setCoord(Vect coord){
			this.coord=coord.deepCopy();

		}
	
	public Vect getCoord(){
		return	this.coord.deepCopy();

		}
	
	public double getR(){
		return	this.coord.norm();

		}
	
	public void setCoord(int i, double u){
		this.coord.el[i]=u;

	}
	public double getCoord(int i){
		return this.coord.el[i];

	}
	
	
	public void setF(Vect F){
		this.F=F.deepCopy();
		
	}
	
	public void setF(int k, double Fu){
		if(F==null) F=new Vect(dim);
		this.F.el[k]=Fu;
		
	}
	
	public Vect getF(){
		return F.deepCopy();
		
	}
	public void setFms(Vect F){
		this.Fms=F.deepCopy();
		
	}

	public void setPBC(int nPBC){
		this.sPBC=(nPBC==1);
		this.aPBC=(nPBC==-1);
		
	}
	
	public boolean hasPBC(){
		return (this.sPBC || this.aPBC);
		
	}

	public int getnPBC(){
		if(this.aPBC) return -1;
		 return 1;

		}
	
	public Mat getRPBC(){
		if(this.aPBC) return util.rotMat2D(-Math.PI/2);
		Mat I=new Mat();
		I.eye(2);
		 return I;
		}
	
	public void setMap(int map){
		this.map=map;
		
	}
	
	public int getMap(){
		return this.map;
		}
	
	public double dist(Node n1){
		return Math.sqrt(dist2(n1));
		
	}
	public double dist2(Node n1){
		return this.getCoord().sub(n1.getCoord()).norm2();
		
	}
}
