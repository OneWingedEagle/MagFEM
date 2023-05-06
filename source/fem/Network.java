package fem;


import java.io.BufferedReader;

import java.io.IOException;

import fem.Network.ElemType;
import io.Loader;
import main.Main;

import math.Mat;

import math.Vect;

import math.util;


public class Network {
	
	public enum ElemType {
	    R, L, C, VPS,
	    CPS, FEM, XX
	}
	
	public class Elem {
		public int id;
		public int index;
		public Elem(){}
		public Elem(ElemType tp){
			type=tp;
		}
		
		 public Node[] nodes=new Node[2];
		 private ElemAndSign[] dependent_elements;
		 public ElemType type;
		 public double R;
		 public double V,I;
		 public double L;
		 public int fem_id,time_id;
		 public int fem_index;

		 public boolean tree;
		 private boolean passed=false;
		public  int unknown_seq_no;

	}
	
	public class Node {
		public Node(){}
		public Node(int id1){ id=id1;}
		public Elem [] connectingElems;
		public int id;
		public Node parent;
		private boolean passed=false;
		private int treeIndex=0;
		private ElemAndSign linkFromParent;
		private ElemAndSign [] linksFromRoot;
	}


	public class ElemAndSign {
	
		public  ElemAndSign(){}
		public  ElemAndSign(Elem elem1, double sign1){elem=elem1;sign=sign1;}
		public double sign;
		public Elem elem;
	}
	public int numElements,numNodes;
	public Elem[] elems;
	public Node[] nodes;
	public Elem[] indep_elems;
	Node nullNode=null;
	
	 public int no_unknown_currents;
	public  int[] unknownCurrentAddress;
	 
	 public Mat PRPt;
	 public Mat tiesetMat;
	
	public Network(){}

	public void print(Model model, Main main){}


	public void read(Loader loader, BufferedReader br) throws IOException{
		String[] sp;
		String line;
		Elem[] elems1=new Elem[1000];
		Node[] nodes1=new Node[1000];
		int ie=0;
		int in=0;
		line=" ";
		
		while(line!=null && line.length()>0){

			line=loader.getNextDataLine(br,"// TYPE ID  NODE1 NODE2 VALUE ");

			line=util.dropLeadingSpaces(line);
			
			if(line.equals("END")) break;
			
			sp=line.split(loader.regex2);

				elems1[ie]=new Elem();
				elems1[ie].id=Integer.parseInt(sp[1]);
				int n1=	Integer.parseInt(sp[2]);
				int nodeIndex=-1;
				for(int j=0;j<in;j++){
					if( nodes1[j].id==n1) {
						nodeIndex=j;
						elems1[ie].nodes[0]=nodes1[j];

					break;
					}
				}
				if(nodeIndex==-1){
					

					nodes1[in]=new Node(n1);
					elems1[ie].nodes[0]=nodes1[in];
					in++;


				}
				
				
				int n2=	Integer.parseInt(sp[3]);
				 nodeIndex=-1;
				for(int j=0;j<in;j++){
					if( nodes1[j].id==n2) {
						nodeIndex=j;
						elems1[ie].nodes[1]=nodes1[j];
					break;
					}
				}
				if(nodeIndex==-1){
					nodes1[in]=new Node(n2);
					elems1[ie].nodes[1]=nodes1[in];
					in++;
				}
				if(sp[0].equals("VPS")){
					System.out.println("!!!!!!!!!!! VPS not available! Job stops.!!!!!!!!!!");
					
					try {
						wait(10000*10000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					elems1[ie].type=ElemType.VPS;
					elems1[ie].time_id=Integer.parseInt(sp[4]);
				}
				else if(sp[0].equals("CPS")){
					elems1[ie].type=ElemType.CPS;
					elems1[ie].time_id=Integer.parseInt(sp[4]);
				}
				else if(sp[0].equals("R")){
					elems1[ie].type=ElemType.R;
					elems1[ie].R=Double.parseDouble(sp[4]);
				}
				else if(sp[0].equals("L")){
					elems1[ie].type=ElemType.L;
				}	
				else if(sp[0].equals("FEM")){
					elems1[ie].type=ElemType.FEM;
					elems1[ie].fem_id=Integer.parseInt(sp[4]);
				}else{
					elems1[ie].type=ElemType.XX;
				}
				
				elems1[ie].index = ie;
				ie++;		
		}
		
		elems=new Elem[ie];
		for(int i=0;i<ie;i++){
			elems[i]=elems1[i];
		}
		
		nodes=new Node[in];
		for(int i=0;i<in;i++)
			nodes[i]=nodes1[i];
		
		numElements=ie;
		numNodes=in;

	
	}
	
	public void constructNetwork(){
		
		AssignConnectingNetworkElementsToNodes();

		//bool invalid = CheckForFreeElementEnds();
		//if (invalid) EXIT(1);

		SetInitialCotree(true);

		SetParentToNodes();

		//invalid = CheckNetworkValidity();
		//if (invalid) EXIT(1);

		TraceNodesToRoot();

		SetFinalCotree();

		SetElementDependency();

		SetTiesetMatrix();
		util.pr("\n Tieset Matrix:\n" );

		tiesetMat.show("%1.0f");
		
	
		SetPRPt();

		//util.pr(this.no_unknown_currents);

//		SetFieldSources();
//
		//SetVariableType();

		//OutputNetwork(ffile_print_out);	
	}


void AssignConnectingNetworkElementsToNodes()
{
	Node node;
	Elem elem;
	
	

	for (int k = 0; k<numNodes; k++)
	{
		node = nodes[k];
		Elem [] connecting=new Elem[100];
		int ix=0;
		for (int j = 0; j<numElements;j++)
		{
			elem = elems[j];
			if (node.id == elem.nodes[0].id || node.id == elem.nodes[1].id )
			{
				connecting[ix++]=elem;
			}
		}
		
		node.connectingElems=new Elem[ix];
		for(int i=0;i<ix;i++)
			node.connectingElems[i]=connecting[i];
	}
}

void SetInitialCotree(boolean forced_cotree) {

	Elem elem;
	
	for (int k = 0; k<numElements;k++){
		elem = elems[k];
		elem.tree = true;

	if (forced_cotree) {
			if (elem.type == ElemType.CPS)
				elem.tree = false;


			if (elem.type == ElemType.VPS)
				elem.tree = false;

		}
	}

}



public void SetParentToNodes()
{

	Node node;
	Elem elem;
	

	nullNode = new Node();
	nullNode.id = -1;
	int lastElemIndex = numElements - 1;
	Node rootNode=null;
	
	int numNetworkSeparateParts = 0;

	while (true) {
		rootNode = null;
		for (int k = lastElemIndex; k >= 0; k--) {
			elem = elems[k];
			if (elem.nodes[0].parent==null) {
				rootNode = elem.nodes[0];
				break;
			}
			else if (elem.nodes[1].parent==null) {
				rootNode = elem.nodes[1];
				break;
			}
		}
		if (rootNode==null) break;

		rootNode.passed = true;
		rootNode.parent = nullNode;
		rootNode.treeIndex = numNetworkSeparateParts;

		SetParentToNodes(rootNode);

		numNetworkSeparateParts++;
	}

}

public void SetParentToNodes(Node node) {

	Elem elem;
	
	double sign;
	
	for (int k = 0; k<node.connectingElems.length; k++)
	{
		elem = node.connectingElems[k];
	
		if (elem.tree == false) continue;

		if (node.id == elem.nodes[0].id) {
		
			if (elem.nodes[1].passed==false && elem.nodes[1].connectingElems!=null) {

				elem.nodes[1].parent = node;
				elem.nodes[1].treeIndex = node.treeIndex;
				sign = 1;
				
				elem.nodes[1].linkFromParent = new ElemAndSign(elem, sign);
				elem.passed = true;
				elem.nodes[1].passed = true;
				
				SetParentToNodes(elem.nodes[1]);
				

			}

		}
		else {
			
			if (elem.nodes[0].passed==false && elem.nodes[0].connectingElems!=null) {

				elem.nodes[0].parent = node;
				elem.nodes[0].treeIndex = node.treeIndex;
				sign = -1;
				
				elem.nodes[0].linkFromParent = new ElemAndSign(elem, sign);
				elem.passed = true;
				elem.nodes[0].passed = true;

				SetParentToNodes(elem.nodes[0]);
			}
				

		
	}

}
}



void TraceNodesToRoot() {


	Node node,node1;
	ElemAndSign connect;

	for (int k = 0; k < numNodes; k++)
	{
		node = nodes[k];
		if (node.linksFromRoot!=null) continue;

		ElemAndSign[] links= new ElemAndSign[100];
		int ix=0;
		node1 = node;

	while (node1.parent.id!= nullNode.id) {

			connect = node1.linkFromParent;

			links[ix++]=connect;

			if (node1.parent.linksFromRoot!=null) {

				int count = node1.parent.linksFromRoot.length;
				for (int j = 0; j < count; j++)
					links[ix++]=node1.parent.linksFromRoot[j];


				break;
			}
			else {

				node1 = node1.parent;
			}
			
		}
	
	//util.pr(">>>>>>>>>>>>.           ----     >>>"+ix);
	
	node.linksFromRoot = new ElemAndSign[ix];
	for (int j = 0; j < ix; j++)
		node.linksFromRoot[j]=links[j];
		
	}
	
}
	public  void solveStatic(Model model){

		double time=model.getCurrentTime();
		
		Vect indepCurrents=new Vect(indep_elems.length);
		for (int j = 0; j<indep_elems.length; ++j)
		{
			if(indep_elems[j].type==ElemType.CPS){
				int time_id=indep_elems[j].time_id;

				double value=0;
				if(model.timeFunctions!=null && time_id>0)
				  value=model.timeFunctions[time_id].getValue(time);
				indepCurrents.el[j]=value;
				elems[j].I=value;
			}
		}
		
		Vect allCurrents=tiesetMat.transp().mul(indepCurrents);

		for (int j = 0; j<numElements; ++j)
		{
			//if(network.elems[j].type!=ElemType.CPS) {
				elems[j].I=allCurrents.el[j];
				if(elems[j].type==ElemType.FEM) {
					int femIndex=elems[j].fem_index;
					model.phiCoils[femIndex].current=elems[j].I;
				}
				util.pr("element: "+elems[j].id+"   current: "+elems[j].I);
			//}
		}
		
	}



void SetFinalCotree() {

	Elem elem;
	
	for (int k = 0; k<numElements; k++)
	{
		elem = elems[k];

		if (elem.tree == false) continue;

		if (!elem.passed)
			elem.tree = false;
	}

}

void SetElementDependency() {


	Elem elem;
	ElemAndSign connect;
	double sign;
	int ix = 0;
	Elem[] independent_elements=new Elem[1000];
	for (int k = 0; k<numElements; k++)
	{
		elem = elems[k];


		if (elem.tree == false)
		{
			independent_elements[ix++]=elem;
			sign = 1.;
			AddDependentElementToElement(elem);

		}
	}
	
	util.pr("Dependent Network Elements:" );
	
	indep_elems=new Elem[ix];
	for (int k = 0; k<ix; k++){
		indep_elems[k]=independent_elements[k];
		
		util.pr(k+" :  "+indep_elems[k].id);

	}

	for (int k = 0; k < numElements; k++)
	{
		elem = elems[k];

		elem.unknown_seq_no = -1;
	}



	int i=0;

	for (int k = 0; k < indep_elems.length; k++)
	{
		elem = indep_elems[k];


		if (elem.type == ElemType.CPS) continue;


		elem.unknown_seq_no = i++;
	}

	no_unknown_currents = i;

	if(no_unknown_currents>0){

	unknownCurrentAddress = new int[no_unknown_currents];

	int jx = 0;
	for (int k = 0; k < indep_elems.length; k++)
	{
		elem = indep_elems[k];
		if (elem.type == ElemType.CPS)
			
		unknownCurrentAddress[jx++] = elem.index;
	}

	}
}


void AddDependentElementToElement(Elem elem){



	for (int k = 0; k<numElements; k++)
	{
		elems[k].passed = false;
	}

	TraceLoop(elem);


}


void TraceLoop(Elem elem) {

	ElemAndSign[] trace1 = elem.nodes[0].linksFromRoot;
	ElemAndSign[] trace2 = elem.nodes[1].linksFromRoot;

	//util.pr(trace1.length);
	//util.pr(trace2.length);

	int numCommonLinks = 0;
	int count1 = trace1.length;
	int count2 = trace2.length;
	ElemAndSign lastConnect1, lastConnect2;
	while (numCommonLinks < Math.min(count1, count2)) {
		lastConnect1 = trace1[count1 - 1 - numCommonLinks];
		lastConnect2 = trace2[count2 - 1 - numCommonLinks];

		if (lastConnect1.elem.id == lastConnect2.elem.id) {
			numCommonLinks++;
		}
		else {
			break;
		}

	}

	ElemAndSign connect;


	ElemAndSign[] loopLinks=new ElemAndSign[1000];
	int ix=0;

	double elemSign = 1.0;

	connect = new ElemAndSign(elem, elemSign);
	


	loopLinks[ix++]=connect;
	for (int k = 0; k<count1 - numCommonLinks; k++)
	{
		connect = trace1[k];

		loopLinks[ix++]=new ElemAndSign(connect.elem, connect.sign*elemSign);

	}


	for (int k = 0; k<count2 - numCommonLinks; k++)
	{
		connect = trace2[k];

		loopLinks[ix++]=new ElemAndSign(connect.elem, -connect.sign*elemSign);
	}
	

	elem.dependent_elements=new ElemAndSign[ix];
	int jx=0;

	Vect ids = new Vect(ix);
	for (int k = 0; k < ix; k++)
		ids.el[k]=loopLinks[k].elem.id;

	int [] indices = ids.bubble();

	for (int k = 0; k < ix; k++)
	{
		int ind = indices[k];
		elem.dependent_elements[k]=loopLinks[ind];
	}

}

void SetTiesetMatrix() {

	int numRows =indep_elems.length;
	int dimension = numElements;
	
	tiesetMat=new Mat(numRows,dimension);

	Elem row_elem;
	Elem elem;;


	for (int i = 0; i<numRows; ++i) {

		row_elem = this.indep_elems[i];

		for (int j = 0; j<row_elem.dependent_elements.length; ++j)
		{
			int col=row_elem.dependent_elements[j].elem.index;
	
			double sign=row_elem.dependent_elements[j].sign;
			tiesetMat.el[i][col]=sign;

		}

	}


}


private void SetPRPt(){
	
Mat R=new Mat(numElements,numElements);

for (int k = 0; k < numElements; k++)
{
	if(elems[k].type==ElemType.R) R.el[k][k]=elems[k].R;
}


PRPt =tiesetMat.mul(R.mul(tiesetMat.transp()));
//PRPt.show("%2.2e");

}

}

