package fem;

import math.util;

public class Edge {
	public byte direction;
	public double  length,A,Ap,T;
	public Node[] node=new Node[2];
	public boolean edgeKnown,edgeKnownT,hasJ,sPBC,aPBC,common;
	public int map;
	
	
	public Edge(Node n1,Node n2)
	{

		node[0]=n1;
		node[1]=n2;
	if(n1.id<n2.id){
		node[0]=n1;
		node[1]=n2;
	}
	else{
		node[0]=n2;
		node[1]=n1;
	}
	}

/*	public Edge(int n1,int n2)
	{

	if(n1<n2){
	endNodeNumber[0]=n1;
	endNodeNumber[1]=n2;
	}
	else{
		endNodeNumber[0]=n2;
		endNodeNumber[1]=n1;
	}

	}
	*/
	public void setKnownA(double A){
		edgeKnown=true;
		this.A=A;

	}
	
	public void setSolvedAL(double A){
		this.A=A;
		
	}
	
	public void saveAp(){
		this.Ap=this.A;
	}

	public void setLength(double length){
				
		this.length=length;
	}
	
	public void setDirection(byte i){
		
		this.direction=i;
	}
	
	public double getDiffA(){

		return A-Ap;
	}

	public void setA(double A) {

		this.A=A;
		
	}
	
	public double getA() {

		return this.A;
		
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
	
}
